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
prefix <- "041_"
out <- glue("{prefix}pca_umap")
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





# ---------------------------------------------------------------------
# Create SLURM and R job files
# ---------------------------------------------------------------------

for (i in seq_along(gemset_name)) {
    for (j in seq_along(norm_method)) {
        curr_gemset_name <- gemset_name[i]
        curr_norm_method <- norm_method[j]

        if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

        # Use the following for SLURM and R filenames
        slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}.slurm")
        r_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}.R")


# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00
#SBATCH --mem=72gb
#SBATCH --tmp=40gb
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
library(Seurat)





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



if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"))


# ---------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------


curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/031_gemsets/{curr_gemset_name}/{curr_norm_method}/031_seurat_object.rds"))





# ---------------------------------------------------------------------
# Generate PCA, TSNE, and UMAP
# ---------------------------------------------------------------------



# All of these methods use the "scale.data" slot of the respective assay for input data


# SCT assay
# These are now standard steps in the Seurat workflow for visualization and clustering
curr_seurat_object <- RunPCA(curr_seurat_object, assay = "SCT", verbose = FALSE, reduction.name = "pca")
#curr_seurat_object <- RunTSNE(curr_seurat_object, assay = "SCT", reduction = "pca", dims = 1:30, reduction.name = "tsne")
curr_seurat_object <- RunUMAP(curr_seurat_object, assay = "SCT", reduction = "pca", dims = 1:30, verbose = FALSE, reduction.name = "umap")
# might get lots of warnings from the umap-learn, but it seems to work!




# RNA assay
# These are now standard steps in the Seurat workflow for visualization and clustering
curr_seurat_object <- RunPCA(curr_seurat_object, assay = "RNA", verbose = FALSE, reduction.name = "pca_rna", reduction.key = "rnaPC_")
#curr_seurat_object <- RunTSNE(curr_seurat_object, assay = "RNA", reduction = "pca_rna", dims = 1:30, reduction.name = "tsne_rna", reduction.key = "rnaTSNE_")
curr_seurat_object <- RunUMAP(curr_seurat_object, assay = "RNA", reduction = "pca_rna", dims = 1:30, verbose = FALSE, reduction.name = "umap_rna", reduction.key = "rnaUMAP_")
# might get lots of warnings from the umap-learn, but it seems to work!





# ---------------------------------------------------------------------
# SCT assay plot
# ---------------------------------------------------------------------


data <- as.data.frame(curr_seurat_object@reductions$pca@cell.embeddings)
pdf(glue("{prefix}sct_pca.pdf"), width = 10, height = 10)
ggplot(data, aes(x = .panel_x, y = .panel_y)) +
	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
	ggforce::geom_autodensity() +
	geom_density2d() +
	ggforce::facet_matrix(vars(glue("PC_{seq(1:4)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
	labs(title = "PCA (derived from SCT scaled data)")
dev.off()


# data <- as.data.frame(curr_seurat_object@reductions$tsne@cell.embeddings)
# pdf(glue("{prefix}sct_tsne.pdf"), width = 6, height = 6)
# ggplot(data, aes(x = .panel_x, y = .panel_y)) +
# 	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
# 	ggforce::geom_autodensity() +
# 	geom_density2d() +
# 	ggforce::facet_matrix(vars(glue("tSNE_{seq(1:2)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
# 	labs(title = "TSNE (derived from SCT scaled data)")
# dev.off()


data <- as.data.frame(curr_seurat_object@reductions$umap@cell.embeddings)
pdf(glue("{prefix}sct_umap.pdf"), width = 6, height = 6)
ggplot(data, aes(x = .panel_x, y = .panel_y)) +
	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
	ggforce::geom_autodensity() +
	geom_density2d() +
	ggforce::facet_matrix(vars(glue("UMAP_{seq(1:2)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
	labs(title = "UMAP (derived from SCT scaled data)")
dev.off()







# ---------------------------------------------------------------------
# RNA assay plot
# ---------------------------------------------------------------------


data <- as.data.frame(curr_seurat_object@reductions$pca_rna@cell.embeddings)
pdf(glue("{prefix}rna_pca.pdf"), width = 10, height = 10)
ggplot(data, aes(x = .panel_x, y = .panel_y)) +
	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
	ggforce::geom_autodensity() +
	geom_density2d() +
	ggforce::facet_matrix(vars(glue("rnaPC_{seq(1:4)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
	labs(title = "PCA (derived from RNA scaled data)")
dev.off()


# data <- as.data.frame(curr_seurat_object@reductions$tsne_rna@cell.embeddings)
# pdf(glue("{prefix}rna_tsne.pdf"), width = 6, height = 6)
# ggplot(data, aes(x = .panel_x, y = .panel_y)) +
# 	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
# 	ggforce::geom_autodensity() +
# 	geom_density2d() +
# 	ggforce::facet_matrix(vars(glue("rnaTSNE_{seq(1:2)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
# 	labs(title = "TSNE (derived from RNA scaled data)")
# dev.off()


data <- as.data.frame(curr_seurat_object@reductions$umap_rna@cell.embeddings)
pdf(glue("{prefix}rna_umap.pdf"), width = 6, height = 6)
ggplot(data, aes(x = .panel_x, y = .panel_y)) +
	geom_point(alpha = 0.2, shape = 16, size = 0.5) +
	ggforce::geom_autodensity() +
	geom_density2d() +
	ggforce::facet_matrix(vars(glue("rnaUMAP_{seq(1:2)}")), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
	labs(title = "UMAP (derived from RNA scaled data)")
dev.off()







# ---------------------------------------------------------------------
# Save Seurat object
# ---------------------------------------------------------------------

saveRDS(curr_seurat_object, file = glue("{prefix}seurat_object.rds"))




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
        system(glue("sbatch {prefix}{curr_gemset_name}_{curr_norm_method}.slurm"))
    }
}






#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
