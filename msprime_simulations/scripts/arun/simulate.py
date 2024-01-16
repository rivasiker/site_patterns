import msprime
import numpy as np
import pandas
import scipy.stats as ss
from scipy.stats import uniform
from tqdm.auto import tqdm
import sys
import argparse
import pickle
import pprint
sys.path.insert(1, '../scripts')
from utils import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate plots and data to see the effects of ')
    parser.add_argument('-s','--sim_file', help='two column space delimited file', required=True)
    parser.add_argument('-m', '--model', help="model ID (will load the file in ../models/<model_id>")
    parser.add_argument('-o','--out', help='output VCF', required=True)
    args = parser.parse_args()


    sim_params = get_sim_params(args.sim_file)
    print("Simulating with parameters:")
    pprint.pprint(sim_params)
    if sim_params['RECOMBINATION_RATE_MAP'] == True:
        rec_param = msprime.RateMap(
        position=sim_params['RECOMBINATION_RATE_MAP_POSITIONS'],
        rate=sim_params['RECOMBINATION_RATE_MAP_RATES'])
    else:
        rec_param = sim_params['RECOMBINATION_RATE']
    
    demography = create_demography("../models/"+args.model)
    samples = [
        msprime.SampleSet(1, "A", ploidy=1),
        msprime.SampleSet(1, "B", ploidy=1),
        msprime.SampleSet(1, "C", ploidy=1),
        msprime.SampleSet(1, "D", ploidy=1),
        msprime.SampleSet(1, "E", ploidy=1)
    ]
    replicates = msprime.sim_ancestry(
        samples,
        demography=demography,
        recombination_rate=rec_param,
        sequence_length=distances[-1]+1,
        ploidy=1,
        num_replicates=int(sim_params['NUM_REPLICATES']),
        gene_conversion_rate=sim_params['GENE_CONVERSION_RATE'],
        gene_conversion_tract_length=sim_params['GENE_CONVERSION_LENGTH']
    )

    for idx,ts in enumerate(tqdm(replicates)):
        mts = msprime.sim_mutations(ts, rate=sim_params['MU'])
        mts.write_vcf(args.out+"."+str(idx)+".vcf")
        
