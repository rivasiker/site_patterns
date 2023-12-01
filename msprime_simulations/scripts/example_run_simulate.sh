#!/bin/bash
#SBATCH --account=durvasul_1174
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sim-%x.%j.out
#SBATCH --error=logs/sim-%x.%j.err
#SBATCH --array=1-10

module purge
eval "$(conda shell.bash hook)"

conda activate msprime-env
I_LINE=$SLURM_ARRAY_TASK_ID
pair=$(head -n $I_LINE  list_to_run2.txt | tail -n 1)
demographic_model_id=`echo $pair | cut -f 1 -d ' '`
param_model_id=`echo $pair | cut -f 2 -d ' '`

python scripts/simulate.py -s sim_params/${param_model_id}.txt -m ${demographic_model_id} -o outputs/${demographic_model_id}-${param_model_id}
