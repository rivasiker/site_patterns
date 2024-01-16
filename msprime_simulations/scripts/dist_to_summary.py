import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Generate dist files from vcf')
parser.add_argument('-t', '--type', help="model ID", required=True, type=int)
parser.add_argument('-s','--nstart', help='ID of starting replicate', required=True, type=int)
parser.add_argument('-e','--nend', help='ID of ending replicate', required=True, type=int)
args = parser.parse_args()

type = args.type
nstart = args.nstart
nend = args.nend


for i in range(nstart, nend+1):
    # Load data
    dat = pd.read_csv(f'../results/02_dist_files/dist_model_{type}.{i}.csv')
    # Pivot longer
    dat2 = dat.melt(['CHROM', 'POS', 'group'])
    # Remove null values
    dat2 = dat2[-dat2['value'].isnull()]
    # Count sites per combination of distance, focal site and target site
    dat3 = dat2.groupby(['variable', 'group', 'value']).size().reset_index()
    # Save distance as numeric variable
    dat3['variable'] = [int(i.replace('group_', '')) for i in dat3['variable']]
    # Rename columns
    dat3.columns = ['distance', 'from', 'to', 'n']
    # Save file
    dat3.to_csv(f'../results/03_summary_files/summary_model_{type}.{i}.csv', index = False)
