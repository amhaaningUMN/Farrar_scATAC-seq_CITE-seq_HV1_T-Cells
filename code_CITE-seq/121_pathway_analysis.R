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
if (curr_threads < 1 | is.null(curr_threads)) {
    stop("Error: The bash THREADS variable is less than 1 or null.")
}
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
prefix <- "121_"
out <- glue("{prefix}pathway_analysis")
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


gemset_name <- read_tsv(glue("{proj_dir}/code_out/031_gemsets/031_gemsets_readme.txt")) %>%
    pull(gemset_name)


cell_cycle_norm_method <- c("cc_none", "cc_s_and_g2m", "cc_s_minus_g2m")
gene_norm_method <- c("gene_none")
norm_method <- tidyr::crossing(cell_cycle_norm_method, gene_norm_method) %>%
    tidyr::unite(col = "unite", everything()) %>%
    pull()



resolutions <- c(0.3, 0.5, 0.7, 0.9)



# ---------------------------------------------------------------------
# Provide list of interesting datasets to investigate
# ---------------------------------------------------------------------

# After discussion with the group, Brian Fife has decided to persue the
# following datasets for pathways analysis:

dataset_list <- list(
    # list(gemset_name = "adtCD4neg", norm_method = "cc_none_gene_none", resolutions = "0.5"),
    list(gemset_name = "adtCD4pos", norm_method = "cc_none_gene_none", resolutions = "0.5")
    # list(gemset_name = "all", norm_method = "cc_none_gene_none", resolutions = "0.5")
    )



# ---------------------------------------------------------------------
# Create SLURM and R job files
# ---------------------------------------------------------------------


for (i in seq_along(dataset_list)) {
    curr_gemset_name <- dataset_list[[i]]$gemset_name
    curr_norm_method <- dataset_list[[i]]$norm_method
    curr_resolution <- dataset_list[[i]]$resolutions
    curr_res_name <- glue("SCT_snn_res_{curr_resolution}")


    cluster_merge_list <- readRDS(glue("{proj_dir}/code_out/051_cluster_de/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/051_cluster_merge_list.rds"))

    # ---------------------------------------------------------------------
    # Find individual comparisons to test
    # ---------------------------------------------------------------------


    # Select only some comparisons to test. By default, test all possible comparisons. Thus,
    # make sure the "my_fav_comparisons" object does not exist. If you want to only test a few
    # comparisons, create a "my_favorite_comparisons" vector.

    # my_fav_comparisons <- c("cluster_0_vs_6", "cluster_1_vs_6")
    my_fav_comparisons <- c(str_subset(names(cluster_merge_list[[1]]$de_data_list), "^cluster"), "1.untreated_leuk_vs_2.nilotinib_pdl1", "1.nilotinib_vs_2.nilotinib_pdl1")


    if ("de_data_list" %in% names(cluster_merge_list[[1]])) {
        # This is a multi-sample comparison (i.e. de_across_samples)
        if (exists("my_fav_comparisons")) {
            # Never run the "conserved across samples"
            my_fav_comparisons <- my_fav_comparisons[!stringr::str_detect(my_fav_comparisons, fixed("_conserved_across_samples"))]
            de_comparisons_to_test <- names(cluster_merge_list[[1]]$de_data_list) %in% my_fav_comparisons
        } else {
            # Never run the "conserved across samples"
            de_comparisons_to_test <- !stringr::str_detect(names(cluster_merge_list[[1]]$de_data_list), fixed("_conserved_across_samples"))
        }
        all_de_data <- cluster_merge_list[[1]]$de_data_list
    } else {
        # This is a single-sample comparison (i.e. de_across_samples not avail)
        if (exists("my_fav_comparisons")) {
            de_comparisons_to_test <- names(cluster_merge_list[[1]]$merge_data_list) %in% my_fav_comparisons
        } else {
            de_comparisons_to_test <- !logical(length = length(names(cluster_merge_list[[1]]$merge_data_list)))
        }
        all_de_data <- cluster_merge_list[[1]]$merge_data_list
    }

    for (k in which(de_comparisons_to_test)) {
        curr_name <- all_de_data[[k]]$name
        curr_de <- all_de_data[[k]]$de
        
        # Write out the necessary de_data files in the appropriate dir
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"))) {{dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"), recursive = TRUE)}}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"))
        
        write_tsv(curr_de, glue("{curr_name}_de.txt"))


        #######################################################################
        # Create job files for each individual comparison
        #######################################################################

        if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

        # Use the following for SLURM and R filenames
        slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}_{curr_name}.slurm")
        r_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}_{curr_name}.R")



# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=32gb
#SBATCH --tmp=16gb
#SBATCH --error=%x.e%j
#SBATCH --output=%x.o%j
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
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
library(viridis)
library(pathwayPCA)
library(cowplot)
library(Seurat)
library(NMF)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(openxlsx)
library(parallel)
library(gtools)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)



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



# Contains clusterProfiler custom functions
source(glue("{proj_dir}/code/001_r_functions.R"))

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
curr_name <- "<<curr_name>>"


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"))) {{dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"), recursive = TRUE)}}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}"))

# ---------------------------------------------------------------------
# Get genesets
# ---------------------------------------------------------------------


mm_msigdb_long <- read_tsv(glue("{proj_dir}/code_out/026_clean_gene_sets_and_msigdb/026_mm_msigdb_long.tsv")) %>%
    # Factor the database type, run analyses by database type
    purrr::modify_at(vars("database"), factor)

gene_info <- read_tsv(glue("{proj_dir}/code_out/015_tidy_count_tables/015_gex_filtered/features.tsv.gz"), col_names = FALSE, col_types = cols(.default = "c")) %>%
    dplyr::rename(ensembl_id = "X1")
# symbol, feature_type

# columns(org.Mm.eg.db)

gene_info$entrez_id <- mapIds(org.Mm.eg.db,
					keys = str_remove(gene_info$ensembl_id, "\\\\..*$"),
					keytype = "ENSEMBL",
					column = "ENTREZID",
					multiVals = "first")

gene_info$symbol <- mapIds(org.Mm.eg.db,
					keys = str_remove(gene_info$ensembl_id, "\\\\..*$"),
					keytype = "ENSEMBL",
					column = "SYMBOL",
					multiVals = "first")
		    
		    
		    
curr_de <- read_tsv(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_name}/{curr_name}_de.txt"))






# ---------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------




print(curr_name)


# Cluster profiler parameters
de <- curr_de %>%
    # Add ensembl and entrez ids
    dplyr::left_join(dplyr::select(.data = gene_info, ensembl_id, entrez_id, symbol), by = "symbol")

log2fc_cutoff <- 0.25
padj_cutoff <- 0.1
logFC_col <- 3
padj_col <- 6
ensembl_id_col <- 7
entrez_id_col <- 8
geneset_to_gene_db_long <- mm_msigdb_long
database_types <- c("custom", "c2_cgp", "h", "c7_immunesigdb", "c8")
# database_types <- c("c2_cgp", "c2_cp_biocarta", "c2_cp_kegg", "c2_cp_pid", "c2_cp_reactome", "c3_tft", "c7", "h")
# database_types <- levels(mm_msigdb_long$database)
plot_title <- glue("{curr_gemset_name}, {curr_norm_method}\n{curr_res_name}, {curr_name}")
run_ora <- FALSE
run_go <- FALSE
run_kegg <- FALSE
run_gsea <- TRUE


source(glue("{proj_dir}/code/001_r_functions.R"))
cluster_profiler <- tk_cluster_profiler(de = de, log2fc_cutoff = log2fc_cutoff, padj_cutoff = padj_cutoff, logFC_col = logFC_col, padj_col = padj_col, ensembl_id_col = ensembl_id_col, entrez_id_col = entrez_id_col, geneset_to_gene_db_long = geneset_to_gene_db_long, database_types = database_types, plot_title = plot_title, run_ora = run_ora, run_go = run_go, run_kegg = run_kegg, run_gsea = run_gsea)
names(cluster_profiler) <- curr_name



#######################################################################
# Save session info
#######################################################################

# ---------------------------------------------------------------------
# Save Cluster Profiler results object only
# ---------------------------------------------------------------------


saveRDS(cluster_profiler, file = glue("{prefix}cluster_profiler.rds"))




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
        system(glue("sbatch {prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}_{curr_name}.slurm"))
    }
}





#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
