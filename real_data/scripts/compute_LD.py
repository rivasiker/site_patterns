import argparse
from cyvcf2 import VCF
import pandas as pd
import numpy as np
from tqdm import tqdm
import datetime

def log_progress(message):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'[{timestamp}] {message}')
    
genostring2pattern = {
    '00001':'M',
    '00010':'O',
    '00011':'HCG',
    '00100':'G',
    '01000':'C',
    '01111':'H',
    '00111':'HC',
    '00101':'HCO',
    '01100':'CG',
    '01011':'HG',
    '01001':'HGO',
    '01110':'CGO',
    '00110':'GO',
    '01010':'CO',
    '01101':'HO',
}

def parse_pattern(l):
    s = ''.join(map(str,l))
    try:
        return(genostring2pattern[s])
    except:
        return("NA")

def get_TMRCA(df2, window_start, length):
    X = window_start
    Y = window_start+length
    selected_rows = df2[(df2['pos'] > X) & (df2['pos'] < Y)]
    count = selected_rows['pattern'].value_counts()[0]
    type = selected_rows['pattern'].value_counts().index.tolist()[0]
    return((count/length), type)

def get_data(chromosome):
    genomat = []
    pos = []
    log_progress("Reading VCF for chr"+str(chromosome))
    for variant in VCF('../vcf/chr'+str(chromosome)+'.vcf.gz'): 
        curr_geno = [str(i[0]) for i in variant.genotypes]
        x = parse_pattern(curr_geno)
        if x == "NA":
            continue
        pos.append(int(variant.start))
        genomat.append(x)
    df = pd.DataFrame({'pos':pos, 'pattern':genomat})
    keep = ['HC', 'CG', 'HG']
    df2 = df[df['pattern'].isin(keep)]
    return(df2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute LD curve for a specified number of points. Run from the scripts directory. Note: does not check for missing data. Make sure to filter out windows with lots of missingness beforehand.")
    parser.add_argument("-o", "--output", help="Output filename")
    parser.add_argument("-w", "--window_size", default=1000, type=int, help="Size of window used to estimate TMRCA at a locus.")
    #parser.add_argument("-s", "--step_size", default=10000, type=int, help="How far apart each starting locus is.")
    parser.add_argument("-c", "--chromosome", help="Chromosome to analyze")

    args = parser.parse_args()
    
    data = get_data(args.chromosome)
    window_size = args.window_size
    #step_size = args.step_size
    #distance = args.distance
    distances = [1, 2, 3, 4, 5, 10, 50, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 50000, 75000, 100000] # 19
    step_sizes = [10000, 10000, 10000, 10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000] # smaller step sizes for bigger distances
    log_progress("Computing pi and estimating TMRCA")
    output_l = []
    for jdx,step_size in enumerate(step_sizes):
        starts = np.arange(data['pos'].iloc[0], data['pos'].iloc[-1], step_size)
        d = distances[jdx]
        for idx,start in enumerate(tqdm(starts)):
            try:
                pi1, type1 = (get_TMRCA(data, start, window_size))
                pi2, type2 = (get_TMRCA(data, start+d, window_size))
                output_l.append([d, pi1, type1, pi2, type2])
            except IndexError:
                continue
    paired_df = pd.DataFrame(output_l, columns=['distance', 'pi1', 'site1', 'pi2', 'site2'])
    correlation_df = paired_df.groupby(['distance', 'site1', 'site2'])[['pi1', 'pi2']].corr()
    correlation_df.to_csv(args.output+'.csv')
