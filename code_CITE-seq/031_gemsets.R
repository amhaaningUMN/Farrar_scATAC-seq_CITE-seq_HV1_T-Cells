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
library(Seurat)
library(gtools)
library(openxlsx)



#######################################################################
# Script parameters
#######################################################################


proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "031_"
out <- glue("{prefix}gemsets")
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
# GEM (cell-barcode) set group #1
# ---------------------------------------------------------------------

# Create subsets of the original data that we call GEM sets (or "gs"). It's possible that
# we'll want to make more (further refined) splits of the data later that might be based
# on things we learn from these GEM sets. Future GEM sets will have a different suffix
# (e.g. gs2).


# For the "gs1" GEM set group, split the data into the sequencing libraries.

# Create readme describing the GEM sets

readme <- tibble::tribble(
  ~gemset_group, ~gemset_name, ~seq_library, ~gemset_description,
"gs1", "all", "GEX", "GEMs/cells from the GEX sequencing library, which includes naive HV1, untreated leukemia, nilotinib, nilotinib+anti-PDL1, and nilotinib+anti-ILR treated samples.",
"gs2", "adtCD4pos", "GEX", "GEMs/cells from the GEX sequencing library that were classified as adtCD4 positive.",
"gs3", "adtCD4neg", "GEX", "GEMs/cells from the GEX sequencing library that were classified as adtCD4 negative."
)


write_tsv(readme, glue("{prefix}gemsets_readme.txt"))




# ---------------------------------------------------------------------
# Get sample sheet
# ---------------------------------------------------------------------


samples <- readRDS(glue("{proj_dir}/code_out/010_samples/010_samples.rds"))




cell_cycle_norm_method <- c("cc_none", "cc_s_and_g2m", "cc_s_minus_g2m")
gene_norm_method <- c("gene_none")
norm_method <- tidyr::crossing(cell_cycle_norm_method, gene_norm_method) %>%
    tidyr::unite(col = "unite", everything()) %>%
    pull()




for (i in seq_along(norm_method)) {
    for (j in seq_len(nrow(readme))) {
        curr_norm_method <- norm_method[i]
        curr_gemset_name <- readme$gemset_name[j]

        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}"))

        # Get normalized data
        curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/020_normalize_gex/{curr_norm_method}/020_seurat_object.rds"))
        meta_data <- glue("{proj_dir}/code_out/081_adt_and_cluster/all/{curr_norm_method}/SCT_snn_res_0.3/081_cell_metadata.rds")

        # Subset data
        if (curr_gemset_name == "all") {
            curr_seurat_object_subset <- curr_seurat_object
        } else if (curr_gemset_name == "adtCD4pos") {
            keep_cells <- readRDS(meta_data) %>%
                dplyr::select(barcode, pos_adtCD4) %>%
                dplyr::filter(pos_adtCD4 == "pos")

            curr_seurat_object_subset <- curr_seurat_object[, curr_seurat_object@meta.data$barcode %in% keep_cells$barcode]
        } else if (curr_gemset_name == "adtCD4neg") {
            keep_cells <- readRDS(meta_data) %>%
                dplyr::select(barcode, pos_adtCD4) %>%
                dplyr::filter(pos_adtCD4 == "neg")

            curr_seurat_object_subset <- curr_seurat_object[, curr_seurat_object@meta.data$barcode %in% keep_cells$barcode]
        } else {
            stop("gemset name is wrong.")
        }

        # Export data
        saveRDS(curr_seurat_object_subset, file = glue("{prefix}seurat_object.rds"))
        saveRDS(curr_seurat_object_subset@meta.data, file = glue("{prefix}cell_metadata.rds"))
        write_tsv(as_tibble(curr_seurat_object_subset@meta.data), glue("{prefix}cell_metadata.txt"))
    }
}




#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
