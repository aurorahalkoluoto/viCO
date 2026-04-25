#!/bin/bash -l
#SBATCH --job-name=tradeSeq           # Job name
#SBATCH --account=project_2007686         # Billing project, has to be defined!
#SBATCH --time=08:00:00             # Max. duration of the job
#SBATCH --mem=80G            # Memory saved
#SBATCH --partition=small           # Job queue (partition)
#SBATCH --output=log/tradeSeq_GlutaNeuro.out
#SBATCH --error=log/tradeSeq_GlutaNeuro.err

module purge
module load tykky

export PATH="/scratch/tykky_env/renv/bin:$PATH"

Rscript pseudotime_tradeSeq_run_GlutaNeuro.R

