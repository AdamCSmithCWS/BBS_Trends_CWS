#!/bin/bash
#SBATCH --job-name=bbs
#SBATCH --account=eccc_cws
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2G
#
# -------------------------
# Logs
# -----------------------
#SBATCH --output=/home/acs001/BBS_Trends_CWS/logs/%x_%A_%a.out
#SBATCH --error=/home/acs001/BBS_Trends_CWS/logs/%x_%A_%a.err


# echo "Node: $(hostname)"
# echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
# echo "Array ID: $SLURM_ARRAY_TASK_ID"

# -------------------------
# Parameters
# -------------------------

aou=${aou:-5450} # default Baird's Sparrow


# -------------------------
# Call R processing script run with container
# -------------------------
apptainer exec /home/acs001/r-bbsbayes2.sif Rscript /home/acs001/BBS_Trends_CWS/gpsc_model_fit.R \
        "$aou"

