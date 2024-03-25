#!/usr/bin/env Rscript

rm(list = ls(all.names = TRUE))

set.seed(20)

options(
    rlang_backtrace_on_error = "full",
    error = rlang::entrace,
    menu.graphics = FALSE,
    repos = c("CRAN" = "https://mirror.las.iastate.edu/CRAN"),
    mc.cores = as.integer(system("[ ! -z ${THREADS+x} ] && echo ${THREADS} || echo 1", intern = TRUE))
)

curr_threads <- getOption("mc.cores")
if (curr_threads < 1 | is.null(curr_threads)) {stop("Error: The bash THREADS variable is less than 1 or null.")}
print(paste0("number of threads: ", curr_threads))



#######################################################################
# Load R packages
#######################################################################


library(toddr)
library(tidyverse)
library(glue)




#######################################################################
# Script parameters
#######################################################################

proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "999_"
out <- glue("{prefix}permissions")
group <- "farrarm"
proj_dir <- glue("/home/{group}/shared/riss/knut0297/{proj}")
out_dir <- glue("{proj_dir}/code_out/{out}")


if (!dir.exists(glue("{out_dir}"))) {
    dir.create(glue("{out_dir}"), recursive = TRUE)
}
setwd(glue("{out_dir}"))




#######################################################################
# Analysis
#######################################################################

# ---------------------------------------------------------------------
# Find paths
# ---------------------------------------------------------------------

# Small number of paths and files
paths_small <- Sys.glob(glue("{proj_dir}/*")) %>% 
    str_subset("code_out", negate = TRUE)

for (i in seq_along(paths_small)) {
    system(glue("find {paths_small[i]} -type f -print0 | xargs -0 -I[] chmod ug+r,g-w,o-rwx []"))
    system(glue("find {paths_small[i]} -type d -print0 | xargs -0 -I[] chmod ug+rxs,g-w,o-rwx []"))
}

# Large number of paths and files
paths <- Sys.glob(glue("{proj_dir}/code_out/*"))
names(paths) <- basename(paths)


# ---------------------------------------------------------------------
# Create SLURM
# ---------------------------------------------------------------------



for (j in seq_along(paths)) {
    curr_path <- paths[j]
    curr_path_name <- names(paths)[j]

    if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

    # Use the following for SLURM and R filenames
    slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}permissions_{curr_path_name}.slurm")


# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
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



CURR_PATH="<<curr_path>>"

# ---------------------------------------------------------------------
# Open permissions for group
# ---------------------------------------------------------------------


find $CURR_PATH -type f -print0 | xargs -0 chmod ug+r,g-w,o-rwx
find $CURR_PATH -type d -print0 | xargs -0 chmod ug+rxs,g-w,o-rwx




', .open = "<<", .close = ">>")

long_string_slurm_vect <- str_split(long_string_slurm, "\\n")[[1]]

con <- file(slurm_file, "w")
for (m  in seq_along(long_string_slurm_vect)) {
	writeLines(str_trim(long_string_slurm_vect[m], side = "right"), con = con)
}
close(con)




    # Launch the SLURM file
    setwd(glue("{proj_dir}/code/{out}.slurm_jobs"))
    system(glue("sbatch {prefix}permissions_{curr_path_name}.slurm"))
}






#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
