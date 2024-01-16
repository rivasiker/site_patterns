import msprime
import pandas
import pickle

def create_demography(file):
    # read a file specifying the demography
    GEN_TIME=20
    with open(file) as f:
        demography = msprime.Demography()
        for line in f:
            l = line.split()
            if l[0] == "pop":
                demography.add_population(name=l[1], initial_size=int(l[2]))
            elif l[0] == "bottleneck":
                demography.add_instantaneous_bottleneck(time=int(l[1]) / GEN_TIME, strength=int(l[2]), population=l[3])
            elif l[0] == "split":
                derived = list(l[2].split("/"))
                demography.add_population_split(time=int(l[1]) / GEN_TIME, derived=derived, ancestral=l[3])
            elif l[0] == 'admixture':
                derived = l[2]
                ancestral = list(l[3].split("/"))
                proportions = [float(props) for props in list(l[4].split("/"))]
                demography.add_admixture(time=int(l[1]) / GEN_TIME, derived=derived, ancestral=ancestral, proportions=proportions)
    return(demography)

def get_sim_params(file):
    out = {}
    with open(file) as f:
        for line in f:
            l = line.rstrip().split()
            if l[1] == "None":
                out[l[0]] = None
                continue
            if l[1] == "True":
                out[l[0]] = True
                continue
            if l[1] == "False":
                out[l[0]] = False
                continue
            if l[0] == "RECOMBINATION_RATE_MAP_POSITIONS":
                out[l[0]] = [int(i) for i in l[1].split(",")]
                continue
            if l[0] == "RECOMBINATION_RATE_MAP_RATES":
                out[l[0]] = [float(i) for i in l[1].split(",")]
                continue
            try:
                out[l[0]] = float(l[1])
            except ValueError:
                out[l[0]] = l[1]
    return(out)

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    with open(filename, 'rb') as inp:
        loaded = pickle.load(inp)
    return(loaded)

def get_tree_type(tree):
    if tree.parent(0) == tree.parent(1):
        type="HC"
    elif tree.parent(0) == tree.parent(2):
        type="HG"
    elif tree.parent(1) == tree.parent(2):
        type="CG"
    else:
        type="other" # should not happen. draw tree to debug
        print(tree.draw_text())
    return(type)

def compute_counts(tree, ps, pd):
    tree_type = get_tree_type(tree)
    # exit early if we get an unexpected toplogy. output the correctly formatted df

    if tree_type == 'other':
        names = ['H', 'C', 'G', 'O', 'M', 'HC', 'HG', 'CG', 'HCG', 'HO', 'CO', 'GO', 'HCO', 'HGO', 'CGO']
        patterns = [np.nan * len(names)]
        d = pandas.DataFrame({'names':names, 'pattern_counts': patterns})
        d['prop'] = d['pattern_counts'] /  d['pattern_counts'].sum()
        return(d)

    tH = tree.branch_length(0)
    tC = tree.branch_length(1)
    tG = tree.branch_length(2)
    tO = tree.branch_length(3)
    tM = tree.branch_length(4) + tree.branch_length(tree.parent(3)) # tM reflects the length of the unrooted tree
    # so it is the length of the M branch plus the parent of HCGO (which we access as the parent of O)
    tHCG = tree.branch_length(tree.parent(tree.parent(0)))
    if tree_type == 'HC':   
        tHC = tree.branch_length(tree.parent(0))
        tHG = 0
        tCG = 0
    elif tree_type == 'HG':
        tHC = 0
        tHG = tree.branch_length(tree.parent(0))
        tCG = 0
    elif tree_type == 'CG':
        tHC = 0
        tHG = 0
        tCG = tree.branch_length(tree.parent(1))

    H   = ps*tH   + pd*tHC*tC + pd*tHCG *tCG + pd*tHG*tG 
    C   = ps*tC   + pd*tHC*tH + pd*tHCG *tHG + pd*tCG*tG 
    G   = ps*tG   + pd*tHCG*tHC + pd*tHG*tG + pd*tCG*tC 
    O   = ps*tO   + pd*tM*tHCG
    M   = ps*tM   + pd*tHCG *tO
    HC  = ps*tHC  + pd*tH*tC  + pd*tHCG*tG
    HG  = ps*tHG  + pd*tH*tG  + pd*tHCG*tC
    CG  = ps*tCG  + pd*tC*tG  + pd*tHCG*tH
    HCG = ps*tHCG + pd*tHC*tG + pd*tHG*tC + pd*tCG*tH + pd*tO*tM
    HO  = pd*tH*tO + pd*tCG*tM
    CO  = pd*tC*tO + pd*tHG*tM
    GO  = pd*tG*tO + pd*tHC*tM
    HCO = pd*tHC*tO + pd*tG*tM
    HGO = pd*tHG*tO + pd*tC*tM
    CGO = pd*tCG*tO + pd*tH*tM

    patterns = [H, C, G, O, M, HC, HG, CG, HCG, HO, CO, GO, HCO, HGO, CGO]
    names = ['H', 'C', 'G', 'O', 'M', 'HC', 'HG', 'CG', 'HCG', 'HO', 'CO', 'GO', 'HCO', 'HGO', 'CGO']
    d = pandas.DataFrame({'names':names, 'pattern_counts': patterns})
    d['prop'] = d['pattern_counts'] /  d['pattern_counts'].sum()
    return(d)
