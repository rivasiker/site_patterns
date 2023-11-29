import pandas as pd
import sys

i = sys.argv[1]

dat = pd.read_csv('../results/dist_files/dist_%s.csv' % i)
dat2 = dat.melt(['CHROM', 'POS', 'group'])
dat2 = dat2[-dat2['value'].isnull()]
dat3 = dat2.groupby(['variable', 'group', 'value']).size().reset_index()
dat3['variable'] = [int(i.replace('group_', '')) for i in dat3['variable']]
dat3.columns = ['distance', 'from', 'to', 'n']
dat3.to_csv('../results/summary_files/summary_%s.csv' % i, index = False)
