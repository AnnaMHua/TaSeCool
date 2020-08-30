#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --account=jeffreyteo
#SBATCH -o myRprog.out
#SBATCH --partition=standard
#SBATCH --mail-type=END


module load anaconda/2019.10-py3.7
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source venv/bin/activate
python TB_findWPs.py 0