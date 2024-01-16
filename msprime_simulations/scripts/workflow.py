from gwf import Workflow, AnonymousTarget

gwf = Workflow()

nstart = 8
nend = 20
nsites = 10_000_000

def msprime_to_vcf(type, nstart, nend, nsites):
    """Simulate msprime and create vcf file"""
    inputs = []
    outputs = [f'../results/01_vcf_files/model_{type}.{i}.vcf' for i in range(nstart, nend+1)]
    options = {'cores': 1, 'memory': '10g', 'walltime': "01:00:00"}
    spec = f"""
        python msprime_to_vcf.py -t {type} -s {nstart} -e {nend} -l {nsites}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_to_dist(type, nstart, nend, human, chimp, goril, orang, macaq):
    """Convert vcf to dist file"""
    inputs = [f'../results/01_vcf_files/model_{type}.{i}.vcf' for i in range(nstart, nend+1)]
    outputs = [f'../results/02_dist_files/dist_model_{type}.{i}.csv' for i in range(nstart, nend+1)]
    options = {'cores': 1, 'memory': '10g', 'walltime': "01:00:00"}
    spec = f"""
        python vcf_to_dist.py -t {type} -s {nstart} -e {nend} -hu {human} -ch {chimp} -go {goril} -or {orang} -ma {macaq}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def dist_to_summary(type, nstart, nend):
    """Convert dist file to summary file"""
    inputs = [f'../results/02_dist_files/dist_model_{type}.{i}.csv' for i in range(nstart, nend+1)]
    outputs = [f'../results/03_summary_files/summary_model_{type}.{i}.csv' for i in range(nstart, nend+1)]
    options = {'cores': 1, 'memory': '10g', 'walltime': "01:00:00"}
    spec = f"""
        python dist_to_summary.py -t {type} -s {nstart} -e {nend}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

for t in range(1, 7):
    gwf.target_from_template(
        f"msprime_to_vcf_{t}", 
        msprime_to_vcf(
            type = t,
            nstart = nstart,
            nend = nend,
            nsites = nsites
            )
        )
    
for t in range(1, 7):
    gwf.target_from_template(
        f"vcf_to_dist_{t}", 
        vcf_to_dist(
            type = t,
            nstart = nstart,
            nend = nend,
            human = 'tsk_0',
            chimp = 'tsk_1',
            goril = 'tsk_2',
            orang = 'tsk_3',
            macaq = 'tsk_4'
            )
        )
    
for t in range(1, 7):
    gwf.target_from_template(
        f"dist_to_summary_{t}", 
        dist_to_summary(
            type = t,
            nstart = nstart,
            nend = nend
            )
        )
    

