import pandas as pd
import sys

i = sys.argv[1]

# Load data
dat = pd.read_csv('../results/dist_files/dist_%s.csv' % i)
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
dat3.to_csv('../results/summary_files/summary_%s.csv' % i, index = False)
