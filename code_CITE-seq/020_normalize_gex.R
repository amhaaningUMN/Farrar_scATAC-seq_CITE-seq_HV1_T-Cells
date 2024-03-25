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
library(cowplot)
library(Seurat)
library(msigdbr)





#######################################################################
# Script parameters
#######################################################################


proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "020_"
out <- glue("{prefix}normalize_gex")
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
# Get sample sheet
# ---------------------------------------------------------------------


samples <- readRDS(glue("{proj_dir}/code_out/010_samples/010_samples.rds"))

cell_cycle_norm_method <- c("cc_none", "cc_s_and_g2m", "cc_s_minus_g2m")
gene_norm_method <- c("gene_none")
norm_method <- tidyr::crossing(cell_cycle_norm_method, gene_norm_method) %>%
    tidyr::unite(col = "unite", everything()) %>%
    pull()




# ---------------------------------------------------------------------
# Get cell cycle genes
# ---------------------------------------------------------------------


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
# The cc.genes vector is from the Seurat package
s_phase <- tibble(human_gene_symbol = Seurat::cc.genes$s.genes,
    marker_class = "s_phase")
g2m_phase <- tibble(human_gene_symbol = Seurat::cc.genes$g2m.genes,
    marker_class = "g2m_phase")
cell_cycle_human <- bind_rows(s_phase, g2m_phase)



# Use the msigdbr R package because it includes annotated orthologs
human_mouse_homologs <- msigdbr::msigdbr(species = "mouse")

human_mouse_homologs_s_phase <- human_mouse_homologs %>%
    dplyr::filter(human_gene_symbol %in% s_phase$human_gene_symbol) %>%
    distinct(gene_symbol, human_gene_symbol) %>%
    dplyr::rename(mouse_gene_symbol = "gene_symbol")

human_mouse_homologs_g2m_phase <- human_mouse_homologs %>%
    dplyr::filter(human_gene_symbol %in% g2m_phase$human_gene_symbol) %>%
    distinct(gene_symbol, human_gene_symbol) %>%
    dplyr::rename(mouse_gene_symbol = "gene_symbol")




write_tsv(human_mouse_homologs_s_phase, glue("{prefix}human_mouse_homologs_s_phase.txt"))
write_tsv(human_mouse_homologs_g2m_phase, glue("{prefix}human_mouse_homologs_g2m_phase.txt"))





#
# # ---------------------------------------------------------------------
# # Get list of TCR genes
# # ---------------------------------------------------------------------
#
#
# # Remove TCR genes from the variable gene list
#
# tcr_features <- read_tsv(glue("{proj_dir}/code_out/011_count_gex/kallisto_bustools/t2g.txt"), col_names = FALSE) %>%
#     dplyr::rename(ensembl_id_txpt = "X1", ensembl_id = "X2", symbol = "X3") %>%
#     dplyr::distinct() %>%
#     dplyr::filter(str_detect(symbol, "^Tr[ab][vdj]")) %>%
#     dplyr::filter(!str_detect(symbol, "^Tradd$"))
#
# write_tsv(tcr_features, glue("{prefix}tcr_features.txt"))
#
#

# ---------------------------------------------------------------------
# Create SLURM and R job files
# ---------------------------------------------------------------------



for (j in seq_along(norm_method)) {
    curr_norm_method <- norm_method[j]

    if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

    # Use the following for SLURM and R filenames
    slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_norm_method}.slurm")
    r_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_norm_method}.R")


# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=8:00:00
#SBATCH --mem=128gb
#SBATCH --tmp=60gb
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
library(gridExtra)
library(cowplot)
library(ggforce)
library(future)

# This should respect the cgroups for the slurm job
future::availableCores()
future::plan("multicore")
# Set max size to 24 GB
options(future.globals.maxSize = 24000 * 1024^2)


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
# Variable names
# ---------------------------------------------------------------------


# Get sample names info
curr_norm_method  <- "<<curr_norm_method>>"




curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/019_prepare_datasets_for_gex/019_seurat_object.rds"))



# ---------------------------------------------------------------------
# Run standard log normalization on the RNA assay
# ---------------------------------------------------------------------


if (!dir.exists(glue("{out_dir}/{curr_norm_method}"))) {dir.create(glue("{out_dir}/{curr_norm_method}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_norm_method}"))


# No restriction on genes potentially included in variable gene list
curr_seurat_object <- NormalizeData(curr_seurat_object, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
curr_seurat_object <- FindVariableFeatures(curr_seurat_object, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = TRUE)

# # Remove TCR genes from the variable gene list
# tcr_features <- read_tsv(glue("{proj_dir}/code_out/020_normalize_gex/020_tcr_features.txt"))
# tcr_features_logical <- curr_seurat_object@assays$RNA@var.features %in% tcr_features$symbol
# curr_seurat_object@assays$RNA@var.features <- curr_seurat_object@assays$RNA@var.features[!tcr_features_logical]


curr_seurat_object <- ScaleData(curr_seurat_object, assay = "RNA", features = rownames(curr_seurat_object), vars.to.regress = c("nCount_RNA", "percent_mito"))




# ---------------------------------------------------------------------
# Run SCT transformation (first time)
# ---------------------------------------------------------------------

# How to use SCT with cell cycle regression
# https://github.com/satijalab/seurat/issues/1679#issuecomment-557781838


# Run sctransform
# Results are saved in a new assay (named SCT by default) with counts being (corrected)
# counts, data being log1p(counts), scale.data being pearson
# residuals; sctransform::vst intermediate results are saved in misc
# slot of new assay.

# Will probably get warnings: "In theta.ml(y = y, mu = fit$fitted) : iteration limit reached"
# This is not a problem:
# https://github.com/satijalab/seurat/issues/1378

# UNRELIABLE VALUE: Future (future_lapply-2) unexpectedly generated random numbers without specifying argument future.seed...
# This is not a problem:
# https://github.com/satijalab/seurat/issues/3622


# No restriction on genes potentially included in variable gene list
# residual.features = NULL allows use of all genes as potential variable genes
curr_seurat_object <- SCTransform(curr_seurat_object, assay = "RNA", vars.to.regress = "percent_mito", verbose = FALSE, return.only.var.genes = FALSE, residual.features = NULL)








# ---------------------------------------------------------------------
# Get cell cycle scores
# ---------------------------------------------------------------------


# Get cell cycle gene list
cell_cycle_s_genes_homologs <- read_tsv(glue("{out_dir}/{prefix}human_mouse_homologs_s_phase.txt"))
cell_cycle_g2m_genes_homologs <- read_tsv(glue("{out_dir}/{prefix}human_mouse_homologs_g2m_phase.txt"))



# We assign each cell a score, based on its expression of G2/M and S phase markers. These
# marker sets should be anticorrelated in their expression levels, and cells expressing neither are
# likely not cycling and in G1 phase.


# CellCycleScoring needs to use scaled data (which has already been done above)
# Find cell cycle scores
curr_seurat_object <- CellCycleScoring(curr_seurat_object, assay = "SCT", s.features = cell_cycle_s_genes_homologs$mouse_gene_symbol, g2m.features = cell_cycle_g2m_genes_homologs$mouse_gene_symbol, set.ident = FALSE)

# Find the difference between S-Phase and G2/M Phase.
curr_seurat_object$CC.Difference <- curr_seurat_object$S.Score - curr_seurat_object$G2M.Score








# ---------------------------------------------------------------------
# Run SCT transformation (second time)
# ---------------------------------------------------------------------

# Regressing cell cycle scores AND/OR remove genes from the variable list




# Determine which variables to include in regression
cc_sct_vars <- NULL
cc_rna_vars <- NULL
if (str_detect(curr_norm_method, "cc_s_and_g2m")) {
    cc_sct_vars <- c("percent_mito", "S.Score", "G2M.Score")
    cc_rna_vars <- c("nCount_RNA", "percent_mito", "S.Score", "G2M.Score")
} else if (str_detect(curr_norm_method, "cc_s_minus_g2m")) {
    cc_sct_vars <- c("percent_mito", "CC.Difference")
    cc_rna_vars <- c("nCount_RNA", "percent_mito", "CC.Difference")
}

gene_vars <- rownames(curr_seurat_object)
if (str_detect(curr_norm_method, "gene_tcr")) {
    # Do not allow TCR genes to be included in the possible set of variable features list.
    tcr_features_logical <- rownames(curr_seurat_object) %in% tcr_features$symbol
    gene_vars <- rownames(curr_seurat_object)[!tcr_features_logical]
}


# Re-normalize data if necessary
if (str_detect(curr_norm_method, "cc_none") & str_detect(curr_norm_method, "gene_none")) {
    # Do nothing
    print("SCTransform normalization was NOT rerun using cell cycle variables -- or TCR-removed variable genes.")
    print("RNA assay ScaleData was NOT rerun using cell cycle variables -- or TCR-removed variable genes.")
} else {
    # Overwrite previous SCT assay
    print(glue("SCTransform normalization was rerun using variables: {cc_sct_vars}"))
    curr_seurat_object <- SCTransform(curr_seurat_object, assay = "RNA", new.assay.name = "SCT", vars.to.regress = cc_sct_vars, verbose = FALSE, return.only.var.genes = FALSE, residual.features = gene_vars)
    print(glue("RNA assay ScaleData was rerun using variables: {cc_sct_vars}"))
    curr_seurat_object <- ScaleData(curr_seurat_object, assay = "RNA", features = gene_vars, vars.to.regress = cc_rna_vars)
}






# ---------------------------------------------------------------------
# Difference in gene tally between SCT and RNA
# ---------------------------------------------------------------------


# NOTE:
# The SCT assay might have fewer genes than RNA@counts.
# https://github.com/ChristophH/sctransform/issues/27
# The vst function has a parameter min_cells set to 5 by default. This means that genes that
# are detected in fewer than 5 cells are not considered during normalization and are not part of the output.

print("Genes not expressed in more than 5 cells are dropped from SCT assay:")
print(paste0("Number of features in RNA assay: ", dim(curr_seurat_object@assays$RNA@counts)[1]))
print(paste0("Number of features in SCT assay: ", dim(curr_seurat_object@assays$SCT@counts)[1]))





# ---------------------------------------------------------------------
# Update seurat object
# ---------------------------------------------------------------------

# Set the default assay
DefaultAssay(curr_seurat_object) <- "SCT"

# Update the project name
new_proj_name <- as.character(glue("Farrar069_{curr_norm_method}"))
curr_seurat_object@project.name <- new_proj_name
curr_seurat_object@meta.data$orig.ident <- new_proj_name
curr_seurat_object@meta.data$orig.ident <- factor(curr_seurat_object@meta.data$orig.ident)







#######################################################################
# Save session info
#######################################################################

# ---------------------------------------------------------------------
# Save Seurat object only
# ---------------------------------------------------------------------

saveRDS(curr_seurat_object, file = glue("{prefix}seurat_object.rds"))
saveRDS(curr_seurat_object@meta.data, glue("{prefix}cell_metadata.rds"))
write_tsv(as_tibble(curr_seurat_object@meta.data), glue("{prefix}cell_metadata.txt"))


# ---------------------------------------------------------------------
# Save all R objects
# ---------------------------------------------------------------------

# save.image(file = glue("{prefix}r_objects.RData"))
write_tsv(data.frame(RData_objects = ls()), glue("{prefix}r_objects.txt"))



# ---------------------------------------------------------------------
# Write out session info
# ---------------------------------------------------------------------


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
    system(glue("sbatch {prefix}{curr_norm_method}.slurm"))
}






#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
