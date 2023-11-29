from gwf import Workflow
from os.path import exists

gwf = Workflow()

chr_lst = [f'chr{i}' for i in list(range(1, 23))+['X']]
human = 'hg38'
chimp = 'panTro4'
goril = 'gorGor3'
orang = 'ponAbe2'
macac = 'rheMac3'

for chrom in chr_lst:
    # if not exists('../results/filtered_maf_files/{}.filtered.maf.gz'.format(chrom)):
    if False:
        gwf.target('maf2filteredmaf_{}'.format(chrom),
               inputs=['../data/{}.maf.gz'.format(chrom)],
               outputs=['../results/filtered_maf_files/{}.filtered.maf.gz'.format(chrom)],
               cores=1,
               memory='1g',
               queue='short',
               walltime= '00:30:00') << """
        maffilter param=options_filt.txt CHR={}
        """.format(chrom)
    # if not exists('../results/summary_files/summary_{}.csv'.format(chrom)):
    # if not exists('../results/vcf_files/{}.vcf.gz'.format(chrom)):
    if True:
        gwf.target('filteredmaf2vfc_{}'.format(chrom),
               inputs=['../data/{}.maf.gz'.format(chrom)],
               outputs=['../results/vcf_files/{}.vcf.gz'.format(chrom)],
               cores=1,
               memory='1g',
               queue='short',
               walltime= '00:30:00') << """
        maffilter param=options_vcf.txt CHR={}
        """.format(chrom)
    # if not exists('../results/dist_files/dist_{}.csv'.format(chrom)):
    # if True:
        gwf.target(f'vcf2dist_{chrom}',
                   inputs=['../results/vcf_files/{}.vcf.gz'.format(chrom)],
                   outputs=['../results/dist_files/dist_{}.csv'.format(chrom)],
                   cores=1,
                   memory='8g',
                   queue='short',
                   walltime= '1:00:00') << f"""
            python calc_dist_pat_2.py {chrom} {human} {chimp} {goril} {orang} {macac}
        """
    # if not exists('../results/summary_files/summary_{}.csv'.format(chrom)):
    # if True:
        gwf.target(f'dist2summary_{chrom}',
                   inputs=['../results/dist_files/dist_{}.csv'.format(chrom)],
                   outputs=['../results/summary_files/summary_{}.csv'.format(chrom)],
                   cores=1,
                   memory='36g',
                   queue='short',
                   walltime= '0:30:00') << f"""
            python calc_summary_pat_2.py {chrom}
            """

    
