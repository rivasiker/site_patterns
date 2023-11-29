import io
import pandas as pd
import gzip
import sys
import numpy as np


def read_vcf(path, human, chimp, goril, orang, macac):
    with gzip.open(path, 'r') as f:
        lines = [l.decode("utf-8") for l in f if not l.startswith(b'##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str,
               'FORMAT':str,  human:str, chimp:str, 
               goril:str, orang:str, macac:str,},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

i = sys.argv[1]
human = sys.argv[2]
chimp = sys.argv[3]
goril = sys.argv[4]
orang = sys.argv[5]
macac = sys.argv[6]

tab = read_vcf('../results/vcf_files/%s.vcf.gz' % i, human, chimp, goril, orang, macac)
# Create grouping variable
tab['group'] = tab[human] +tab[chimp] +tab[goril] +tab[orang] +tab[macac]
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

tab.to_csv('../results/dist_files/dist_%s.csv' % i, index = False)