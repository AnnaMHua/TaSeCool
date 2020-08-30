#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --account=jeffreyteo
#SBATCH -o myRprog.out
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --array=0-1


module load anaconda/2019.10-py3.7

source venv/bin/activate
python TB_rivannafindWPs.py ${SLURM_ARRAY_TASK_ID}