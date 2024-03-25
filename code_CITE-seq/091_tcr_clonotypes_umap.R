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



# Set R working directory
proj_dir <- "/home/fifeb/shared/ris/knut0297/cd4_ins_pd1_scrnaseq_20230501"



# Set up output directory and filename prefix
out_dir_name <- "091_tcr_clonotypes_umap"
out_dir <- glue("{proj_dir}/code_out/{out_dir_name}")
if (!dir.exists(glue("{out_dir}"))) {dir.create(glue("{out_dir}"), recursive = TRUE)}
setwd(glue("{out_dir}"))

prefix <- "091_"




#######################################################################
# Analysis
#######################################################################



gemset_name <- read_tsv(glue("{proj_dir}/code_out/031_gemsets/031_gemsets_readme.txt")) %>%
    pull(gemset_name)



cell_cycle_norm_method <- c("cc_none", "cc_s_and_g2m", "cc_s_minus_g2m")
gene_norm_method <- c("gene_none", "gene_tcr")
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
    list(gemset_name = "lymall", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "lymctrl", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "lympd1", norm_method = "cc_none_gene_none", resolutions = "0.7"),
    list(gemset_name = "lympdl1", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "panall", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "panctrl", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "panpd1", norm_method = "cc_none_gene_none", resolutions = "0.9"),
    list(gemset_name = "panpdl1", norm_method = "cc_none_gene_none", resolutions = "0.5")
    )


#######################################################################
# Create SLURM and R job files
#######################################################################


for (h in seq_along(dataset_list)) {
    curr_gemset_name <- dataset_list[[h]]$gemset_name
    curr_norm_method <- dataset_list[[h]]$norm_method
    curr_resolution <- dataset_list[[h]]$resolutions
    curr_res_name <- glue("SCT_snn_res_{curr_resolution}")

    if (!dir.exists(glue("{proj_dir}/code/{out_dir_name}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out_dir_name}.slurm_jobs"), recursive = TRUE)}

    # Use the following for SLURM and R filenames
    slurm_file <- glue("{proj_dir}/code/{out_dir_name}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm")
    r_file <- glue("{proj_dir}/code/{out_dir_name}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.R")



# ---------------------------------------------------------------------
# Create SLURM job file
# ---------------------------------------------------------------------
long_string_slurm <- glue('#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem=32gb
#SBATCH --tmp=16gb
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



PROJ_DIR="<<proj_dir>>"


module load /home/lmnp/knut0297/software/modulesfiles/R/4.2.0


Rscript --vanilla ${PROJ_DIR}/code/<<out_dir_name>>.slurm_jobs/<<basename(r_file)>>



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
library(forcats)
library(glue)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(viridis)
library(scales)
library(openxlsx)
library(Seurat)




#######################################################################
# Script parameters
#######################################################################


# Set R working directory
proj_dir <- "<<proj_dir>>"


# Set up output directory and filename prefix
out_dir_name <- "<<out_dir_name>>"
out_dir <- glue("{proj_dir}/code_out/{out_dir_name}")
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
setwd(out_dir)

prefix <- "<<prefix>>"



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

curr_res_name_final <- glue("{curr_res_name}_de_merge_final")
col_cluster <- paste0(curr_res_name, "_de_merge_final")
col_sample <- "sample_name"


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}"))



# ---------------------------------------------------------------------
# Get input data
# ---------------------------------------------------------------------


hto <- readRDS(glue("{proj_dir}/code_out/016_demultiplex_hto/hto/match_gex_filtered/016_gmm_demux.rds")) %>%
    dplyr::select(barcode, hto_class1, hto_class2, hto_class3, gem_type, sample_name)

tet <- readRDS(glue("{proj_dir}/code_out/017_normalize_tet/017_tet_dsb_normalized.rds")) %>%
    # Transpose the matrix
    t(.) %>%
    tibble::as_tibble(rownames = "barcode")

meta <- readRDS(glue("{proj_dir}/code_out/081_tetramer_and_cluster/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/081_cell_metadata.rds")) %>%
    dplyr::select(barcode, Cluster, starts_with("ord_tet"), starts_with("pos_tet"), tet_features_combined) %>%
    magrittr::set_rownames(.$barcode)

contigs2 <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_contigs.rds")) %>%
    dplyr::left_join(hto, by = "barcode") %>%
    dplyr::mutate(across(c(hto_class1, hto_class2, hto_class3, gem_type, sample_name), as.character)) %>%
    dplyr::mutate(across(c(hto_class1, hto_class2, hto_class3, gem_type, sample_name), ~ replace_na(.x, "negative"))) %>%
    dplyr::left_join(meta, by = "barcode") %>%
    dplyr::left_join(tet, by = "barcode")


clonotypes2 <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_clonotypes.rds"))





# ---------------------------------------------------------------------
# Restrict contigs and clonotypes to current gemset
# ---------------------------------------------------------------------


contigs3 <- contigs2 %>%
    dplyr::filter(barcode %in% rownames(meta))


clonotypes3 <- clonotypes2 %>%
    dplyr::filter(clonotype_id %in% contigs3$raw_clonotype_id)

# ---------------------------------------------------------------------
# Morista
# ---------------------------------------------------------------------

vdj_groups <- contigs3 %>%
    drop_na(Cluster) %>%
    dplyr::filter(Cluster != "all") %>%
    group_by(Cluster)

            contig_list1 <- as.list(group_split(vdj_groups))
            contig_list1_names <- group_keys(vdj_groups) %>%
                pull(Cluster) %>%
                as.character()

            combined <- combineTCR(contig_list1, 
                samples = rep("clust", times = length(contig_list1_names)),
                ID = contig_list1_names,
                cells = "T-AB")

            table <- clonalOverlap(combined, cloneCall = "aa", method = "morisita", exportTable = TRUE) %>%
                tidyr::pivot_longer(cols = -names) %>%
                as_tibble() %>%
                dplyr::rename(group1 = "names", group2 = "name") %>%
                dplyr::mutate(group1 = str_remove(group1, "^_")) %>%
                dplyr::mutate(group2 = str_remove(group2, "^_")) %>%
                dplyr::filter(!is.na(value))

            p <- table %>%
                dplyr::mutate(group1 = forcats::fct_relevel(group1, gtools::mixedsort)) %>%
                dplyr::mutate(group2 = forcats::fct_relevel(group2, gtools::mixedsort)) %>%
                ggplot(aes(x = group1, y = group2, fill = value)) +
                geom_tile() +
                geom_text(aes(label = round(value, 3)), color = "white", size = 4) +
                scale_fill_gradientn(name = "Morisita Index", colors = c("gray95", rev(viridis::plasma(100)))) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                labs(title = "morisita by cluster, CDR3_aa")
            pdf(glue("{prefix}clonal_overlap_morisita_cluster.pdf"), width = 8, height = 8)
            print(p)
            dev.off()




# ---------------------------------------------------------------------
# Unique aa clonotypes
# ---------------------------------------------------------------------

# Get a list of only the unique CDR3aa clonotypes and seqs
clonotype_id_cdr3aa_unique <- clonotypes3 %>%
    dplyr::select(clonotype_id_cdr3aa, cdr3s_aa) %>%
    distinct() %>%
    arrange(order(gtools::mixedorder(clonotype_id_cdr3aa)))


# Get a list of only the unique CDR3nt clonotypes and seqs
clonotype_id_cdr3nt_unique <- clonotypes3 %>%
    dplyr::select(clonotype_id_cdr3nt, cdr3s_nt) %>%
    distinct() %>%
    arrange(order(gtools::mixedorder(clonotype_id_cdr3nt)))


saveRDS(clonotype_id_cdr3aa_unique, glue("{prefix}clonotype_id_cdr3aa_unique.rds"))
write_tsv(clonotype_id_cdr3aa_unique, glue("{prefix}clonotype_id_cdr3aa_unique.txt"))




# ---------------------------------------------------------------------
# Find top clonotypes by clonotype_id_cdr3aa, split by sample
# ---------------------------------------------------------------------


# Find the top number of clonotypes by sample
gem_tally_by_clonotype_id_cdr3aa_sample_name_top <- contigs3 %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(hto_class2, clonotype_id_cdr3aa) %>%
    tally() %>%
    arrange(desc(n)) %>%
    slice_head(n = 5) %>%
    dplyr::rename(n_gems_with_clonotype = "n") %>%
    dplyr::rename(sample_name = "hto_class2")


saveRDS(gem_tally_by_clonotype_id_cdr3aa_sample_name_top, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.rds"))

write_tsv(gem_tally_by_clonotype_id_cdr3aa_sample_name_top, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_clonotype_id_cdr3aa_sample_name_top), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)



gem_tally_by_v_gene_sample_name_top <- contigs3 %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(hto_class2, v_gene) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    arrange(desc(n)) %>%
    slice_head(n = 10) %>%
    dplyr::rename(n_gems_with_clonotype = "n") %>%
    dplyr::rename(freq_within_sample_name = "freq") %>%
    dplyr::rename(sample_name = "hto_class2")

saveRDS(gem_tally_by_v_gene_sample_name_top, glue("{prefix}gem_tally_by_v_gene_sample_name_top.rds"))

write_tsv(gem_tally_by_v_gene_sample_name_top, glue("{prefix}gem_tally_by_v_gene_sample_name_top.txt"))
out_filename <- glue("{prefix}gem_tally_by_v_gene_sample_name_top.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_v_gene_sample_name_top), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)


# ---------------------------------------------------------------------
# Find top 20 clonotypes expressed in tetramer specific cells
# ---------------------------------------------------------------------

# Get only single TRA and single TRB clonotypes. Within each set of tetramer
# positive cells, tally the number of GEMs per clonotype. Find the top 21 clonotypes.

tet_groups_names <- meta %>%
    dplyr::select(starts_with("pos_tet")) %>%
    colnames()

tet_groups_list <- list()
for (i in seq_along(tet_groups_names)) {
    curr_tet_group <- tet_groups_names[i]

    tet_groups_list[[i]] <- contigs3 %>%
        distinct(barcode, .keep_all = TRUE) %>%
        dplyr::mutate(count_tra = str_count(cdr3s_aa, "TRA:")) %>%
        dplyr::mutate(count_trb = str_count(cdr3s_aa, "TRB:")) %>%
        dplyr::mutate(count_tra_trb = factor(glue("tra_{count_tra}_trb_{count_trb}"))) %>%
        dplyr::filter(count_tra_trb == "tra_1_trb_1") %>%
        dplyr::filter(!!rlang::sym(curr_tet_group) == "pos") %>%
        group_by(clonotype_id_cdr3aa) %>%
        tally() %>%
        arrange(desc(n)) %>%
        slice_head(n = 20) %>%
        dplyr::rename(n_gems_with_clonotype = "n") %>%
        dplyr::mutate(tet_group = curr_tet_group)

    names(tet_groups_list)[i] <- curr_tet_group
}
tet_groups <- bind_rows(tet_groups_list)





saveRDS(tet_groups, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_tetramer_pos_top20.rds"))
write_tsv(tet_groups, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_tetramer_pos_top20.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_tetramer_pos_top20.xlsx")
wb <- openxlsx::write.xlsx(list(tally = tet_groups), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)




# ---------------------------------------------------------------------
# Clonotypes by cluster
# ---------------------------------------------------------------------


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))

# ---------------------------------------------------------------------
# Public clonotypes in TET positive cells only
# ---------------------------------------------------------------------

# Determine a list of clonotypes that appear
# to be public (shared across mice). Plot the normalized expression of the TET tetramers
# in cells expressing public vs non-public clonotpes.

# Find the public clonotypes that are tet pos 2.5HP.
sample_name_public_pos_tet2.5HP <- contigs3 %>%
    dplyr::filter(pos_tet2.5HP == "pos") %>%
    group_by(clonotype_id_cdr3aa) %>%
    count(sample_name) %>%
    tidyr::pivot_wider(clonotype_id_cdr3aa, names_from = "sample_name", values_from = "n") %>%
    tidyr::drop_na()



tet_by_clonotype_id_raw <- contigs2 %>%
    dplyr::mutate(public_tet2.5HP = case_when(
        clonotype_id_cdr3aa %in% sample_name_public_pos_tet2.5HP$clonotype_id_cdr3aa ~ "yes",
        TRUE ~ "no"))






p <- tet_by_clonotype_id_raw %>%
    dplyr::select(barcode, Cluster, pos_tet2.5HP, public_tet2.5HP) %>%
    ggplot(aes(x = Cluster, fill = factor(public_tet2.5HP))) +
    geom_bar(position = position_dodge2(preserve = "single")) +
    scale_fill_discrete(name = "tet2.5HP public\\nclonotypes") +
    theme_minimal() +
    labs(title = "Number of GEMs (cells) in each cluster,\\nsplit by clonotype membership (tet2.5HP public/private)",
        y = "Number of GEMs (i.e. cells)")
pdf(glue("{prefix}gem_tally_by_cluster_and_tet2.5HP_public.pdf"))
print(p)
dev.off()

p <- tet_by_clonotype_id_raw %>%
    dplyr::select(barcode, Cluster, tet2.5HP, public_tet2.5HP) %>%
    group_by(Cluster, public_tet2.5HP) %>%
    ggplot(aes(x = Cluster, y = tet2.5HP, fill = factor(public_tet2.5HP))) +
    geom_boxplot() +
    scale_fill_discrete(name = "tet2.5HP public\\nclonotypes") +
    theme_minimal() +
    labs(title = "Expression of tet2.5HP in each cluster,\\nsplit by clonotype membership (tet2.5HP public/private)",
        y = "Normalized (dsb method) expression\\nof tet2.5HP tet")
pdf(glue("{prefix}tet2.5HP_tet_by_cluster_and_tet2.5HP_public.pdf"))
print(p)
dev.off()


# ---------------------------------------------------------------------
# Plot gem tallies by cluster and sample
# ---------------------------------------------------------------------


curr_file <- glue("{proj_dir}/code_out/051_cluster_de/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/final_cluster_labels.txt")
if (file.exists(curr_file)) {
    gem_tally_sample_name_cluster <- contigs3 %>%
        distinct(barcode, .keep_all = TRUE) %>%
        dplyr::mutate(sample_name = forcats::fct_drop(sample_name, "all")) %>%
        group_by(clonotype_id_cdr3aa, sample_name, Cluster, .drop = FALSE) %>%
        tally() %>%
        dplyr::mutate(sample_name_cluster = glue("sample_name_{sample_name}_cluster_{Cluster}")) %>%
        ungroup() %>%
        # Remove redundant columns
        dplyr::select(-sample_name, -Cluster) %>%
        tidyr::pivot_wider(names_from = sample_name_cluster, values_from = n) %>%
        dplyr::relocate(clonotype_id_cdr3aa, .before = 1) %>%
        dplyr::left_join(clonotype_id_cdr3aa_unique, by = "clonotype_id_cdr3aa") %>%
        arrange(order(gtools::mixedorder(clonotype_id_cdr3aa))) %>%
        dplyr::mutate(count_tra = str_count(cdr3s_aa, "TRA:")) %>%
        dplyr::mutate(count_trb = str_count(cdr3s_aa, "TRB:")) %>%
        dplyr::mutate(count_tra_trb = factor(glue("tra_{count_tra}_trb_{count_trb}"))) %>%
        dplyr::mutate(clonotype_id_cdr3aa_number = as.integer(str_remove(clonotype_id_cdr3aa, fixed("clonotype_cdr3aa")))) %>%
        dplyr::relocate(clonotype_id_cdr3aa_number, .after = clonotype_id_cdr3aa)



    # Write this file out
    write_tsv(gem_tally_sample_name_cluster, glue("{prefix}gem_tally_sample_name_cluster.txt"))
    out_filename <- glue("{prefix}gem_tally_sample_name_cluster.xlsx")
    wb <- openxlsx::write.xlsx(list(tally = gem_tally_sample_name_cluster), file = out_filename, rowNames = FALSE)
    openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
    openxlsx::setColWidths(wb, sheet = 1, cols = 2, widths = 6)
    openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)


    p <- gem_tally_sample_name_cluster %>%
        dplyr::select(-matches("NA$", perl = TRUE)) %>%
        tidyr::pivot_longer(cols = starts_with("sample_name_")) %>%
        dplyr::filter(value > 1 & !is.na(value)) %>%
        dplyr::mutate(clonotype_id_cdr3aa = fct_rev(fct_relevel(clonotype_id_cdr3aa, gtools::mixedsort))) %>%
        dplyr::mutate(name = fct_relevel(name, gtools::mixedsort)) %>%
        ggplot(aes(x = name, y = clonotype_id_cdr3aa, fill = value)) +
        scale_x_discrete(position = "top") +
        geom_tile() +
        scale_fill_gradientn(name = "Number of\\nGEMs (cells)", colors = rev(c(viridis::plasma(100)))) +
        theme_minimal() +
        # hjust = 0 makes the labels align at the top, next to the plot
        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, size = 3),
            axis.text.y = element_text(size = 0.5),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
        labs(title = glue("{curr_gemset_name}\\n{curr_norm_method}\\n{curr_res_name}"),
            subtitle = "Only clonotypes, sample_names, and clusters \\nwith >1 GEM (cell) are shown",
            x = "Expt. sample_name and Cluster",
            y = "Amino acid based clonotype\\n(clonotype_id_cdr3aa)")

    pdf(glue("{prefix}gem_tally_sample_name_cluster.pdf"), width = 6, height = 11)
    print(p)
    dev.off()

} else {
    readr::write_lines("final_cluster_labels.txt file is missing.", glue("{prefix}gem_tally_arm_cluster_MISSING.txt"))
}










# ---------------------------------------------------------------------
# UMAP highlight by clonotypes
# ---------------------------------------------------------------------


file1 <- glue("{proj_dir}/code_out/051_cluster_de/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/051_seurat_object.rds")
file2 <- glue("{proj_dir}/code_out/081_tetramer_and_cluster/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/081_cell_metadata.rds")

if (file.exists(file1) & file.exists(file2)) {
    col_cluster <- paste0(curr_res_name, "_de_merge_final")
    col_sample <- "sample_name"

    curr_seurat_object <- readRDS(file1)
    curr_seurat_object@meta.data <- readRDS(file2) %>%
        # Duplicate final cluster ids with simple colname
        dplyr::mutate(Cluster = !!sym(col_cluster)) %>%
        dplyr::mutate(sample_name = fct_relevel(sample_name, gtools::mixedsort))

    # Add the clonotypes when they match a GEX barcode
    umap_data <- cbind(curr_seurat_object@reductions$umap@cell.embeddings, curr_seurat_object@meta.data) %>%
        as_tibble() %>%
        left_join(dplyr::select(contigs3, barcode, clonotype_id_cdr3aa), by = "barcode")
    rm(curr_seurat_object)
    toddr::clean_mem()

    # Duplicate the umap df, adding facet variable and data
    umap_data2 <- list()
    for (m in seq_len(nrow(gem_tally_by_clonotype_id_cdr3aa_sample_name_top))) {
        curr_sample_name <- as.character(gem_tally_by_clonotype_id_cdr3aa_sample_name_top$sample_name[m])
        curr_clonotype <- as.character(gem_tally_by_clonotype_id_cdr3aa_sample_name_top$clonotype_id_cdr3aa[m])
        curr_facet_group <- glue("{curr_sample_name}\\n{curr_clonotype}")

        umap_data2[[m]] <- umap_data %>%
            dplyr::mutate(facet_group = curr_facet_group) %>%
            dplyr::mutate(facet_value = case_when(
                sample_name == curr_sample_name & clonotype_id_cdr3aa == curr_clonotype ~ "yes",
                TRUE ~ NA_character_))
    }

    umap_data3 <- bind_rows(umap_data2) %>%
        dplyr::mutate(facet_group = fct_relevel(as.character(facet_group), gtools::mixedsort)) %>%
        dplyr::mutate(facet_value = fct_relevel(as.character(facet_value), gtools::mixedsort))

    n_plot_rows <- length(levels(umap_data3$sample_name))


    centers <- umap_data %>%
        dplyr::select(Cluster, starts_with("UMAP_")) %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarize_all(median)

    p <- umap_data3 %>%
        # Order by facet_value
        # arrange always sorts NAs to the end, regardless of desc()
        dplyr::arrange(!is.na(facet_value), facet_value) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = facet_value)) +
            geom_point(size = 0.5) +
            scale_color_manual(name = "Clonotype Positive", na.value = "grey95", values = c("firebrick1")) +
            labs(title = paste0("Top 5 most frequent clonotypes per sample_name\\n", curr_sample_name, "\\n", curr_res_name), subtitle = "Cluster numbers inside plot") +
            geom_text(data = centers, mapping = aes(label = Cluster), size = 2, color = "grey40") +
            facet_wrap(vars(facet_group), nrow = n_plot_rows) +
            theme_cowplot() +
            theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                legend.position = "none",
                strip.background = element_blank(),
                strip.text = element_text(size = 6))

    pdf(glue("{prefix}umap_by_top_clonotypes.pdf"), width = 5 * 2, height = n_plot_rows * 2)
    print(p)
    dev.off()
} else {
    readr::write_lines("seurat_object is missing.", glue("{prefix}umap_by_top_clonotypes_MISSING.txt"))
}




# ---------------------------------------------------------------------
# Top 20 clonotypes
# ---------------------------------------------------------------------





if (file.exists(file1) & file.exists(file2)) {
    col_cluster <- paste0(curr_res_name, "_de_merge_final")
    col_sample <- "sample_name"

    curr_seurat_object <- readRDS(file1)
    curr_seurat_object@meta.data <- readRDS(file2) %>%
        # Duplicate final cluster ids with simple colname
        dplyr::mutate(Cluster = !!sym(col_cluster)) %>%
        dplyr::mutate(sample_name = fct_relevel(sample_name, gtools::mixedsort))

    # Add the clonotypes when they match a GEX barcode
    umap_data <- cbind(curr_seurat_object@reductions$umap@cell.embeddings, curr_seurat_object@meta.data) %>%
        as_tibble() %>%
        left_join(dplyr::select(contigs3, barcode, clonotype_id_cdr3aa), by = "barcode")
    rm(curr_seurat_object)
    toddr::clean_mem()


    top20_groups <- unique(tet_groups$tet_group)


    for (m in seq_along(top20_groups)) {
        top_clonotypes <- tet_groups %>%
            dplyr::filter(tet_group == top20_groups[m])

        curr_clonotypes <- top_clonotypes$clonotype_id_cdr3aa
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{top20_groups[m]}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{top20_groups[m]}"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{top20_groups[m]}"))
        for (n in seq_along(curr_clonotypes)) {
            curr_clonotype <- as.character(curr_clonotypes[n])


            umap_data2 <- umap_data %>%
                dplyr::mutate(facet_value = case_when(
                    clonotype_id_cdr3aa == curr_clonotype ~ "yes",
                    TRUE ~ NA_character_)) %>%
                dplyr::mutate(facet_value = fct_relevel(as.character(facet_value), gtools::mixedsort))

            centers <- umap_data %>%
                dplyr::select(Cluster, starts_with("UMAP_")) %>%
                dplyr::group_by(Cluster) %>%
                dplyr::summarize_all(median)

            p <- umap_data2 %>%
                # Order by facet_value
                # arrange always sorts NAs to the end, regardless of desc()
                dplyr::arrange(!is.na(facet_value), facet_value) %>%
                ggplot(aes(x = UMAP_1, y = UMAP_2, color = facet_value)) +
                    geom_point(size = 0.5) +
                    scale_color_manual(name = "Clonotype Positive", na.value = "grey95", values = c("firebrick1")) +
                    labs(title = glue("{top20_groups[m]}\\n{curr_clonotype}\\nMost frequent clonotypes per tetramer pos cells\\n{curr_res_name}"), subtitle = "Cluster numbers inside plot") +
                    geom_text(data = centers, mapping = aes(label = Cluster), size = 2, color = "grey40") +
                    theme_cowplot() +
                    theme(aspect.ratio = 1,
                        plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5),
                        legend.position = "none",
                        strip.background = element_blank(),
                        strip.text = element_text(size = 6))

            pdf(glue("{prefix}umap_by_top_clonotypes_{top20_groups[m]}_{curr_clonotype}.pdf"), width = 6, height = 6)
            print(p)
            dev.off()
        }
    }
} else {
    readr::write_lines("seurat_object is missing.", glue("{prefix}umap_by_top_clonotypes_MISSING.txt"))
}

















#######################################################################
# Save session info
#######################################################################



toddr::write_session_info(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{prefix}"))



', .open = "<<", .close = ">>")
long_string_r_vect <- str_split(long_string_r, "\\n")[[1]]

con <- file(r_file, "w")
for (m  in seq_along(long_string_r_vect)) {
	writeLines(str_trim(long_string_r_vect[m], side = "right"), con = con)
}
close(con)


    # Launch the SLURM file
    setwd(glue("{proj_dir}/code/{out_dir_name}.slurm_jobs"))
    system(glue("sbatch {prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm"))
}




#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))










