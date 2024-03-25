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
prefix <- "051_"
out <- glue("{prefix}cluster_de")
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
# Get samples
# ---------------------------------------------------------------------


gemset_name <- read_tsv(glue("{proj_dir}/code_out/031_gemsets/031_gemsets_readme.txt")) %>%
    pull(gemset_name)



cell_cycle_norm_method <- c("cc_none", "cc_s_and_g2m", "cc_s_minus_g2m")
gene_norm_method <- c("gene_none")
norm_method <- tidyr::crossing(cell_cycle_norm_method, gene_norm_method) %>% 
    tidyr::unite(col = "unite", everything()) %>%
    pull()


# resolutions <- c(0.3, 0.5, 0.7, 0.9)
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)


# ---------------------------------------------------------------------
# Create SLURM and R job files
# ---------------------------------------------------------------------



for (i in seq_along(gemset_name)) {
    for (j in seq_along(norm_method)) {
        for (k in seq_along(resolutions)) {
            curr_gemset_name <- gemset_name[i]
            curr_norm_method <- norm_method[j]
            curr_resolution <- resolutions[k]
            curr_res_name <- glue("SCT_snn_res_{curr_resolution}")
            
            if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

            # Use the following for SLURM and R filenames
            slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm")
            r_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.R")



# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=6:00:00
#SBATCH --mem=100gb
#SBATCH --tmp=60gb
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
    # if [ -n "${SLURM_JOB_ID+x}" ]; then
    #     scontrol show job "${SLURM_JOB_ID}"
    #     sstat -j "${SLURM_JOB_ID}" --format=JobID,MaxRSS,MaxVMSize,NTasks,MaxDiskWrite,MaxDiskRead
    # fi
}
# Execute these functions after any error (i.e. nonzero exit code) or 
# when exiting the script (i.e. with zero or nonzero exit code).
trap trap_my_error ERR
trap trap_my_exit EXIT


# Prepend MODULPATH with personal software modules
module use /home/lmnp/knut0297/software/modulesfiles

# If not a slurm job, set THREADS to 1
THREADS=$([ -n "${SLURM_CPUS_PER_TASK+x}" ] && echo "${SLURM_CPUS_PER_TASK}" || echo 1)
export THREADS

echo "[$(date)] Script start."

#######################################################################
# Script
#######################################################################



PROJ_DIR="<<proj_dir>>"


module load /home/lmnp/knut0297/software/modulesfiles/R/4.2.0


Rscript --vanilla ${PROJ_DIR}/code/<<out>>.slurm_jobs/<<basename(r_file)>>



', .open = "<<", .close = ">>")

long_string_slurm_vect <- str_split(long_string_slurm, "\\n")[[1]]

con <- file(slurm_file, "w")
for (m  in seq_along(long_string_slurm_vect)) {
	writeLines(str_trim(long_string_slurm_vect[m], side = "right"), con = con)
}
close(con)






# ---------------------------------------------------------------------
# R file
# ---------------------------------------------------------------------
long_string_r <- glue('#!/usr/bin/env Rscript

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
library(openxlsx)
library(cowplot)
library(Seurat)
library(SingleR)
library(celldex)
library(patchwork)
library(plotly)
library(gridExtra)
library(NMF)
library(RColorBrewer)
library(parallel)
library(viridis)
library(multtest)
library(metap)
library(ggrepel)





#######################################################################
# Script parameters
#######################################################################


proj <- "<<proj>>"
prefix <- "<<prefix>>"
out <- "<<out>>"
group <- "<<group>>"
proj_dir <- "<<proj_dir>>"
out_dir <- "<<out_dir>>"


if (!dir.exists(glue("{out_dir}"))) {
    dir.create(glue("{out_dir}"), recursive = TRUE)
}
setwd(glue("{out_dir}"))


#######################################################################
# Analysis
#######################################################################



# ---------------------------------------------------------------------
# Variable names and data
# ---------------------------------------------------------------------


# Get sample names info
curr_gemset_name <- "<<curr_gemset_name>>"
curr_norm_method  <- "<<curr_norm_method>>"
curr_resolution <- <<curr_resolution>>
curr_res_name <- "<<curr_res_name>>"



if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


# ---------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------


curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/041_pca_umap/{curr_gemset_name}/{curr_norm_method}/041_seurat_object.rds"))

# Make sure the variable is set to orig.ident
curr_seurat_object@meta.data$orig.ident <- curr_seurat_object@meta.data$sample_name

# ---------------------------------------------------------------------
# Seurat SNN based clustering
# ---------------------------------------------------------------------




DefaultAssay(curr_seurat_object) <- "SCT"
curr_seurat_object <- FindNeighbors(curr_seurat_object, reduction = "pca", dims = 1:30, verbose = FALSE)
curr_seurat_object <- FindClusters(curr_seurat_object, verbose = FALSE, resolution = curr_resolution)





# ---------------------------------------------------------------------
# Merge clusters based on DE gene criteria
# ---------------------------------------------------------------------


setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"))


DefaultAssay(curr_seurat_object) <- "SCT"

# Setting num_genes_diff_between_clusters_threshold = 0 will not attempt any cluster merging based on DE analysis.
source(glue("{proj_dir}/code/002_r_functions.R"))
cluster_merge_list <- tk_citeseq_cluster_merge2(seurat_object = curr_seurat_object, seurat_object_name = curr_seurat_object@project.name, lfc_threshold = 0.25, num_genes_diff_between_clusters_threshold = 0)




# ---------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------


setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


saveRDS(cluster_merge_list, glue("{prefix}cluster_merge_list.rds"))
saveRDS(curr_seurat_object, file = glue("{prefix}seurat_object.rds"))
saveRDS(curr_seurat_object@meta.data, glue("{prefix}cell_metadata.rds"))
write_tsv(as_tibble(curr_seurat_object@meta.data, rownames = "cell_id"), glue("{prefix}cell_metadata.txt"))






#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{prefix}"))



', .open = "<<", .close = ">>")
long_string_r_vect <- str_split(long_string_r, "\\n")[[1]]

con <- file(r_file, "w")
for (m  in seq_along(long_string_r_vect)) {
	writeLines(str_trim(long_string_r_vect[m], side = "right"), con = con)
}
close(con)


        # Launch the SLURM file
        setwd(glue("{proj_dir}/code/{out}.slurm_jobs"))
        system(glue("sbatch {prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm"))
        }
    }
}





#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
