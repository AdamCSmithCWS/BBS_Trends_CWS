#!/bin/bash
#SBATCH --job-name=bbs_gamye
#SBATCH --account=eccc_cws
#SBATCH --array=1-3                 # one task per species (1-indexed line numbers in aou_list.txt)
#SBATCH --cpus-per-task=4            # matches chains=4 / parallel_chains=4 in the R script
#SBATCH --mem-per-cpu=2G                  # adjust based on largest species dataset
#SBATCH --time=01:00:00              # adjust based on expected runtime per species
#SBATCH --partition=standard
#SBATCH --output=/home/acs001/BBS_Trends_CWS/logs/bbs_%A_%a.out
#SBATCH --error=/home/acs001/BBS_Trends_CWS/logs/bbs_%A_%a.err

# --- setup ---
set -euo pipefail

# make sure log dir exists (harmless if it already does)
mkdir -p /home/acs001/BBS_Trends_CWS/logs

# load apptainer (adjust module name/version to match your cluster's module system)
# module load apptainer

# container image and path to the R script it will run
CONTAINER="/home/acs001/BBS_Trends_CWS/r-bbsbayes2.sif"
RSCRIPT_PATH="/home/acs001/BBS_Trends_CWS/gpsc_model_fit.R"

# path to the file listing all aou codes, one per line
AOU_LIST="/home/acs001/BBS_Trends_CWS/aou_list.txt"

# pull the aou code corresponding to this array task
aou=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${AOU_LIST}")

if [ -z "${aou}" ]; then
  echo "ERROR: no aou code found for array task ${SLURM_ARRAY_TASK_ID} in ${AOU_LIST}"
  exit 1
fi

export aou

echo "Running species aou=${aou} on array task ${SLURM_ARRAY_TASK_ID} with ${SLURM_CPUS_PER_TASK} cpus"

# --- run the R script inside the apptainer container ---
# --bind mounts /home/acs001 (and thus BBS_Trends_CWS) into the container at the same path
# --env aou passes the SLURM-exported variable through to the container's environment
apptainer exec \
  --env aou="${aou}" \
  "${CONTAINER}" \
  Rscript "${RSCRIPT_PATH}"
