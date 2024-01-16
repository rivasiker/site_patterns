import io
import pandas as pd
import numpy as np
import argparse


def read_vcf(path, human, chimp, goril, orang, macaq):
    """
    This function reads a (possibly gzipped) vcf file and parses it into a 
    pandas data frame. 

    Variables:
        path (str): the path of the gzipped vcf file
        human, chimp, goril, orang, macaq (str): the name of the 
            column corresponding to the human, chimp, gorilla, 
            orangutan and macaque sites in the vcf file, 
            respectively. 
    """
    with open(path, 'r') as f:
        if path.split('.')[-1] == 'gz':
            lines = [l.decode("utf-8") for l in f if not l.startswith(b'##')]
        else:
            lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str,
               'FORMAT':str,  human:str, chimp:str, 
               goril:str, orang:str, macaq:str,},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})



parser = argparse.ArgumentParser(description='Generate dist files from vcf')
parser.add_argument('-t', '--type', help="model ID", required=True, type=int)
parser.add_argument('-s','--nstart', help='ID of starting replicate', required=True, type=int)
parser.add_argument('-e','--nend', help='ID of ending replicate', required=True, type=int)
parser.add_argument('-hu', '--human', help="human name in vcf", required=True, type=str)
parser.add_argument('-ch', '--chimp', help="human name in vcf", required=True, type=str)
parser.add_argument('-go', '--goril', help="human name in vcf", required=True, type=str)
parser.add_argument('-or', '--orang', help="human name in vcf", required=True, type=str)
parser.add_argument('-ma', '--macaq', help="human name in vcf", required=True, type=str)
args = parser.parse_args()

type = args.type
nstart = args.nstart
nend = args.nend
human = args.human
chimp = args.chimp
goril = args.goril
orang = args.orang
macaq = args.macaq

for i in range(nstart, nend+1):
    tab = read_vcf(f'../results/01_vcf_files/model_{type}.{i}.vcf', human, chimp, goril, orang, macaq)
    # Create grouping variable
    tab['group'] = tab[human] +tab[chimp] +tab[goril] +tab[orang] +tab[macaq]
    tab = tab[['CHROM', 'POS', 'group']]
    # Remove missing data or >2-allelic sites
    tab = tab[['.' not in i for i in tab['group']]]
    tab = tab[['2' not in i for i in tab['group']]]
    tab = tab[['3' not in i for i in tab['group']]]
    tab = tab[['4' not in i for i in tab['group']]]
    # Force genotype in outgroup to be 0
    lst_variants = []
    gr = list(tab['group'])
    for geno in gr:
        if geno[4] != '0':
            dct = {geno[4]:'0', '0':geno[4]}
            geno = [dct[i] if i in list(dct.keys()) else i for i in geno]
        lst_variants.append(''.join(geno))
    tab['group'] = lst_variants
    # Convert group to pattern
    spstr = 'HCGOM'
    tab['group'] = [''.join([spstr[j] for j in range(len(i)) if i[j]=='1']) for i in tab['group']]
    # Create dictionary {chrom:{pos:pattern}}
    dct = dict(zip(tab['POS'], tab['group']))
    # Find pattern at pos+d distance
    for d in sorted(list(set([int(i) for i in np.logspace(0, 5, num=30).round()]))):
        print(d)
        tab['group_%s' % d] = [dct.get(i+d) for i in tab['POS']]
    tab.to_csv(f'../results/02_dist_files/dist_model_{type}.{i}.csv', index = False)