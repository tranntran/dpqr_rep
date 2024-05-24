#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=26
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=48:00:00
#SBATCH --array=1-20
#SBATCH --partition=open
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[your_psuid]@psu.edu
#SBATCH --error=[your_directory]/error/kng_%j.err
#SBATCH --output=[your_directory]/log/kng_%A_%a.out

cd [your_directory]

module load r/4.3.2

Rscript src/simulation_kng_heteroscedastic.R $SLURM_ARRAY_TASK_ID > out/simulation_kng_heteroscedastic_$SLURM_ARRAY_TASK_ID.Rout 
