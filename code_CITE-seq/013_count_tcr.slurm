#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=96gb
#SBATCH --tmp=16gb
#SBATCH --error=%x.e%j
#SBATCH --output=%x.o%j
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
#SBATCH --partition=agsmall,ag2tb


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
prefix="013_"
out="${prefix}count_tcr"
group="farrarm"
proj_dir="/home/${group}/shared/riss/knut0297/${proj}"
out_dir="${proj_dir}/code_out/${out}"

mkdir -p $out_dir
cd $out_dir

module load /home/lmnp/knut0297/software/modulesfiles/cellranger/7.1.0


# ---------------------------------------------------------------------
# Get input data
# ---------------------------------------------------------------------

read1="$(
awk '
BEGIN{FS="\t"}
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{ if (NR > 1 && $(f["seq_library"]) == "TCR") { print $(f["fastq_r1"]) }}
' $proj_dir/code_out/010_samples/010_samples.txt
)"

read2="$(
awk '
BEGIN{FS="\t"}
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{ if (NR > 1 && $(f["seq_library"]) == "TCR") { print $(f["fastq_r2"]) }}
' $proj_dir/code_out/010_samples/010_samples.txt
)"

echo $read1
echo $read2

fastq_dir=$(dirname ${read1})


# Find the sample prefix
input="$(basename $read1)"
if [[ $input =~ ^(.+)_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+\.fastq\.gz$ ]]; then
    sample_base="${BASH_REMATCH[1]}"
    echo $sample_base
else
    echo "No match"
fi



# ---------------------------------------------------------------------
# Cellranger vdj (TCR dataset)
# ---------------------------------------------------------------------


cellranger vdj \
--id "Farrar_Project_069" \
--fastqs "${fastq_dir}" \
--reference "$CELLRANGER_REFS/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0" \
--sample "${sample_base}" \
--chain "TR" \
--localcores 48 \
--localmem 96 \
--localvmem 96 \
--disable-ui
