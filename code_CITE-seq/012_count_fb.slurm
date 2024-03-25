#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH --mem=120gb
#SBATCH --tmp=50gb
#SBATCH --error=%x.e%j
#SBATCH --output=%x.o%j
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
#SBATCH --partition=agsmall


#######################################################################
# Script preliminaries
#######################################################################

# Exit script immediately upon error
set -o errexit -o errtrace -o pipefail -o functrace


function trap_my_error {
    >&2 echo "ERROR: \"${BASH_COMMAND}\" command failed with exit code $? [$(date)]"
}

function trap_my_exit {
    echo "[$(date)] Script exit."
    # Print env variables
    declare -p
    # Print slurm job details
    if [ -n "${SLURM_JOB_ID+x}" ]; then
        scontrol show job "${SLURM_JOB_ID}"
        sstat -j "${SLURM_JOB_ID}" --format=JobID,MaxRSS,MaxVMSize,NTasks,MaxDiskWrite,MaxDiskRead
    fi
}
# Execute these functions after any error (i.e. nonzero exit code) or 
# when exiting the script (i.e. with zero or nonzero exit code).
trap trap_my_error ERR
trap trap_my_exit EXIT

# If not a slurm job, set THREADS to 1
THREADS=$([ -n "${SLURM_CPUS_PER_TASK+x}" ] && echo "${SLURM_CPUS_PER_TASK}" || echo 1)
export THREADS

echo "[$(date)] Script start."

#######################################################################
# Script
#######################################################################

# ---------------------------------------------------------------------
# Project variables
# ---------------------------------------------------------------------

proj="cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix="012_"
out="${prefix}count_fb"
group="farrarm"
proj_dir="/home/${group}/shared/riss/knut0297/${proj}"
out_dir="${proj_dir}/code_out/${out}"

mkdir -p $out_dir
cd $out_dir




module load /home/lmnp/knut0297/software/modulesfiles/R/4.2.0
module load /home/lmnp/knut0297/software/modulesfiles/python3/3.7.10
module load /home/lmnp/knut0297/software/modulesfiles/kallisto/0.48.0
module load /home/lmnp/knut0297/software/modulesfiles/bustools/0.39.3
module load /home/lmnp/knut0297/software/modulesfiles/kite/aafc8c8
# Activate python virtual environment
. /home/lmnp/knut0297/software/python_venvs/kallisto_bustools/bin/activate


# Copy the 10X Genomics barcodes whitelist (provided by the software)
# Copy the version 2 barcodes used in this project
module load /home/lmnp/knut0297/software/modulesfiles/cellranger/6.1.2
cp "/home/lmnp/knut0297/software/modules/cellranger/6.1.2/cellranger-6.1.2/lib/python/cellranger/barcodes/737K-august-2016.txt" "."

Rscript --vanilla ${proj_dir}/code/${out}.R
