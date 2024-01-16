import msprime
import argparse
import pandas as pd

## PARSE ARGUMENTS ##

parser = argparse.ArgumentParser(description='Generate vcf files from msprime simulation')
parser.add_argument('-s','--nstart', help='ID of starting replicate', required=True, type=int)
parser.add_argument('-e','--nend', help='ID of ending replicate', required=True, type=int)
parser.add_argument('-t', '--type', help="model ID", required=True, type=int)
parser.add_argument('-l','--nsites', help='number of sites', required=True, type=int)
args = parser.parse_args()

nstart = args.nstart
nend = args.nend
type = args.type
nsites = args.nsites

param_tab = pd.read_table(f'../models/M{type}.model', sep = ' ', names = ['param', 'value'])
param_dct = {param_tab.param[i]:float(param_tab.value[i]) for i in range(len(param_tab))}

## DEFINE PARAMETERS ##

mu = param_dct['mu']
r = param_dct['r']
g = param_dct['g']

### For models 1 and 2

# Effective population size
N_A = param_dct['N_A']
N_B = param_dct['N_B']
N_C = param_dct['N_C']
N_D = param_dct['N_D']
N_E = param_dct['N_E']
N_AB = param_dct['N_AB']
N_ABC = param_dct['N_ABC']
N_ABCD = param_dct['N_ABCD']
N_ABCDE = param_dct['N_ABCDE']
# Time, in generations, from present to...
#    A/B split for A
t_A = param_dct['t_A']/g
#    A/B split for B
t_B = param_dct['t_B']/g
#    AB/C split for C
t_C = param_dct['t_C']/g
#    ABC/D split for D
t_D = param_dct['t_D']/g
#    ABCD/E split for E
t_E = param_dct['t_E']/g
# Time, in generations, between...
#    A/B and AB/C splits
t_AB = param_dct['t_AB']/g # t_C-t_B
#    AB/C and ABC/D splits
t_ABC = param_dct['t_ABC']/g # t_D-t_C
#    ABC/D and ABCD/E splits
t_ABCD = param_dct['t_ABCD']/g # t_E-t_D
# This is the maximum time in generations from the ABCD/E split to present
t_1 = max([t_A+t_AB+t_ABC+t_ABCD, t_B+t_AB+t_ABC+t_ABCD, t_C+t_ABC+t_ABCD, t_D+t_ABCD, t_E])

### For admixture model only

if 'm' in param_dct:
    # Effective population size
    N_left = param_dct['N_left']
    # Time, in generations, between...
    #    left/B split and ABC/C split
    t_split = param_dct['t_split']/g
    # Admixture proportion
    m = param_dct['m']

## DEFINE MSPRIME MODEL ##
    
model = msprime.Demography()
model.add_population(name="A", initial_size=N_A, default_sampling_time=t_1-(t_A+t_AB+t_ABC+t_ABCD))
model.add_population(name="B", initial_size=N_B, default_sampling_time=t_1-(t_B+t_AB+t_ABC+t_ABCD))
if 'm' in param_dct:
    model.add_population(name="left", initial_size=N_left, initially_active=False)
model.add_population(name="C", initial_size=N_C, default_sampling_time=t_1-(t_C+t_ABC+t_ABCD))
model.add_population(name="D", initial_size=N_D, default_sampling_time=t_1-(t_D+t_ABCD))
model.add_population(name="E", initial_size=N_E, default_sampling_time=t_1-t_E)
model.add_population(name="AB", initial_size=N_AB)
model.add_population(name="ABC", initial_size=N_ABC)
model.add_population(name="ABCD", initial_size=N_ABCD)
model.add_population(name="ABCDE", initial_size=N_ABCDE)
if 'm' in param_dct:
    model.add_admixture(time = t_1-(t_AB+t_ABC+t_ABCD), derived="A", ancestral=["left", "B"], proportions=(m, 1-m))
    model.add_population_split(time=t_1-(t_split+t_ABC+t_ABCD), derived=["left", "B"], ancestral="AB")
else:
    model.add_population_split(time=t_1-(t_AB+t_ABC+t_ABCD), derived=["A", "B"], ancestral="AB")
model.add_population_split(time=t_1-(t_ABC+t_ABCD), derived=["AB", "C"], ancestral="ABC")
model.add_population_split(time=t_1-t_ABCD, derived=["ABC", "D"], ancestral="ABCD")
model.add_population_split(time=t_1, derived=["ABCD", "E"], ancestral="ABCDE")

# samples = [
#         msprime.SampleSet(1, "A", ploidy=1),
#         msprime.SampleSet(1, "B", ploidy=1),
#         msprime.SampleSet(1, "C", ploidy=1),
#         msprime.SampleSet(1, "D", ploidy=1),
#         msprime.SampleSet(1, "E", ploidy=1)
#     ]

replicates = msprime.sim_ancestry(
    samples = {"A": 1, "B": 1, "C": 1, "D": 1, "E" : 1},
    demography = model,
    recombination_rate = r,
    sequence_length = nsites,
    ploidy = 1,
    num_replicates=int(nend-nstart+1)
    )

for idx, ts in enumerate(replicates):
    mts = msprime.sim_mutations(ts, rate = mu)
    rep = list(range(nstart, nend+1))[idx]
    with open(f"../results/01_vcf_files/model_{type}.{rep}.vcf", "w") as vcf:
        vcf.writelines(mts.as_vcf())