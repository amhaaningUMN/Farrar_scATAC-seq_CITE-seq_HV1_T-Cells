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
prefix <- "019_"
out <- glue("{prefix}prepare_datasets_for_gex")
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


# ---------------------------------------------------------------------
# Prepare datasets for GEX analysis with Seurat
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Process GEX feature data
# ---------------------------------------------------------------------



# Fix the duplicate symbol issue.
# Are there any original symbols that contain: "--1", etc? No.
# Seurat cannot handle underscores in feature names and many genes already have hyphens.
# Therefore, use "-dup" as our duplicate separator

txpt_genes_symbol <- read_tsv(glue("{proj_dir}/code_out/011_count_gex/kallisto_bustools/t2g.txt"), col_names = FALSE) %>%
    dplyr::rename(ensembl_id_txpt = "X1", ensembl_id = "X2", symbol_with_dups = "X3")

txpt_genes_symbol_distinct <- txpt_genes_symbol %>%
    dplyr::select(-ensembl_id_txpt) %>%
    dplyr::distinct()

features_gex1 <- read_tsv(glue("{proj_dir}/code_out/015_tidy_count_tables/015_gex_filtered/features.tsv.gz"), col_names = FALSE) %>%
    dplyr::rename(ensembl_id = "X1") %>%
    dplyr::left_join(txpt_genes_symbol_distinct, by = "ensembl_id") %>%
    dplyr::mutate(symbol = base::make.unique(symbol_with_dups, sep = "-dup"))



# There appear to be duplicated symbol names. Search Ensembl to find the dups and pick the
# best gene id to represent the locus. Searched for info in GRCh38.p13 (release 104).
features_gex1_symbol_dups <- features_gex1 %>%
    dplyr::group_by(symbol_with_dups) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n > 1) %>%
    dplyr::arrange(desc(n), symbol_with_dups)

write_tsv(features_gex1_symbol_dups, glue("{prefix}features_gex_symbol_dups.txt"))



# ---------------------------------------------------------------------
# Process GEX raw counts
# ---------------------------------------------------------------------

gex_counts1 <- readRDS(glue("{proj_dir}/code_out/017_normalize_adt/017_gex_raw_matrix_cell_gems.rds"))


# This is a numeric vector.
# The length equals the number of genes in dataset (nrow).
# The values indicate the number of cells that express this particular gene (with a RNA count > zero).
number_of_cells_that_exp_this_gene <- Matrix::rowSums(gex_counts1 > 0)

# This is a numeric vector.
# The length equals the number of cells in dataset (ncol).
# The values indicate the number of genes that are expressed this particular cell (with a RNA count > zero).
number_of_genes_that_are_exp_in_this_cell <- Matrix::colSums(gex_counts1 > 0)


# Subset the data
min_cells_exp_the_gene <- 3
min_genes_exp_in_the_cell <- 300

# Require a cell to express at least 300 genes.
# Require a gene to be expressed in at least 3 cells.
gex_counts2 <- gex_counts1[which(number_of_cells_that_exp_this_gene >= min_cells_exp_the_gene), which(number_of_genes_that_are_exp_in_this_cell >= min_genes_exp_in_the_cell)]


# Filter the GEX features to match the GEX table
features_gex2 <- features_gex1 %>%
    dplyr::filter(ensembl_id %in% rownames(gex_counts2)) %>%
    dplyr::arrange(factor(ensembl_id, levels = rownames(gex_counts2)))


# Convert the GEX rownames from ensembl id to unique symbols
stopifnot(identical(features_gex2$ensembl_id, rownames(gex_counts2)))
rownames(gex_counts2) <- features_gex2$symbol



# ---------------------------------------------------------------------
# Import the ADT data
# ---------------------------------------------------------------------


# Raw, dsb-normalized, and clr-normalized data
adt_counts1 <- readRDS(glue("{proj_dir}/code_out/017_normalize_adt/017_adt_raw_matrix_cell_gems.rds"))
adt_dsb_norm1 <- readRDS(glue("{proj_dir}/code_out/017_normalize_adt/017_adt_dsb_normalized.rds"))
adt_clr_norm1 <- readRDS(glue("{proj_dir}/code_out/017_normalize_adt/017_adt_clr_transformed.rds"))

# ---------------------------------------------------------------------
# Process HTO metadata
# ---------------------------------------------------------------------

# Keep only GEM barcodes that represent "single sample GEMs" (no multiplets) from the
# HTO demultiplexing process.
hto_meta1 <- readRDS(glue("{proj_dir}/code_out/016_demultiplex_hto/hto/match_gex_filtered/016_gmm_demux.rds")) %>%
    dplyr::filter(gem_type == "same_sample_gem") %>%
    dplyr::select(barcode, gem_type, sample_name, sample_name_long, treatment, starts_with("hto_counts_clr_")) %>%
    dplyr::mutate(across(c(gem_type, sample_name, sample_name_long, treatment), as.character))


# Keep only GEMs with GEX data that intersect with HTO gems (i.e. also have a sample label)
keep_barcodes <- intersect(colnames(gex_counts2), hto_meta1$barcode)

hto_counts1 <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_hto_raw.rds"))


hto_clr_norm1 <- hto_meta1 %>%
    dplyr::select(starts_with("hto_counts_clr")) %>%
    dplyr::rename_all(~ stringr::str_remove(.,"^hto_counts_clr_")) %>%
    as.matrix() %>%
    magrittr::set_rownames(hto_meta1$barcode) %>%
    t()




# ---------------------------------------------------------------------
# # Filter all the datasets to all contain the same GEMs
# ---------------------------------------------------------------------


hto_meta2 <- hto_meta1 %>%
    dplyr::select(-starts_with("hto_counts_clr")) %>%
    dplyr::filter(barcode %in% keep_barcodes) %>%
    dplyr::arrange(factor(barcode, levels = keep_barcodes))

hto_counts2 <- hto_counts1$matrix[, keep_barcodes]
hto_clr_norm2 <- hto_clr_norm1[, keep_barcodes]



gex_counts3 <- gex_counts2[, keep_barcodes]



adt_counts2 <- adt_counts1[, keep_barcodes]
adt_dsb_norm2 <- adt_dsb_norm1[, keep_barcodes]
adt_clr_norm2 <- adt_clr_norm1[, keep_barcodes]





# ---------------------------------------------------------------------
# Verify GEM barcodes
# ---------------------------------------------------------------------

stopifnot(identical(hto_meta2$barcode, colnames(gex_counts3)))
stopifnot(identical(hto_meta2$barcode, colnames(hto_counts2)))
stopifnot(identical(hto_meta2$barcode, colnames(hto_clr_norm2)))
stopifnot(identical(hto_meta2$barcode, colnames(adt_counts2)))
stopifnot(identical(hto_meta2$barcode, colnames(adt_dsb_norm2)))
stopifnot(identical(hto_meta2$barcode, colnames(adt_clr_norm2)))

rownames(gex_counts3) %>% head()
rownames(hto_counts2) %>% head()
rownames(hto_clr_norm2) %>% head()
rownames(adt_counts2) %>% head()
rownames(adt_dsb_norm2) %>% head()
rownames(adt_clr_norm2) %>% head()



# ---------------------------------------------------------------------
# Create Seurat object with GEX counts
# ---------------------------------------------------------------------

curr_seurat_object <- CreateSeuratObject(counts = gex_counts3, assay = "RNA", project = "Farrar_Project_069")


# Add the HTO dataset
curr_seurat_object[["HTO"]] <- CreateAssayObject(counts = hto_counts2)
curr_seurat_object[["HTO"]] <- CreateAssayObject(data = hto_clr_norm2)

# Add the ADT dataset
curr_seurat_object[["ADT"]] <- CreateAssayObject(counts = adt_counts2)
curr_seurat_object[["ADT"]] <- CreateAssayObject(data = adt_dsb_norm2)


# Add RNA features metadata to seurat object
curr_seurat_object@assays$RNA@meta.features$ensembl_id <- features_gex2$ensembl_id
curr_seurat_object@assays$RNA@meta.features$symbol <- features_gex2$symbol


curr_seurat_object@meta.data$barcode <- colnames(gex_counts3)

curr_seurat_object@meta.data <- curr_seurat_object@meta.data %>%
    dplyr::relocate(barcode) %>%
    dplyr::left_join(hto_meta2, by = "barcode") %>%
    magrittr::set_rownames(.$barcode)

stopifnot(identical(rownames(curr_seurat_object@meta.data), colnames(curr_seurat_object)))

# ---------------------------------------------------------------------
# Add mito percent per cell
# ---------------------------------------------------------------------


mito_genes <- features_gex2 %>%
    dplyr::filter(str_detect(symbol, "^mt-")) %>%
    pull(symbol)


percent_mito <- Matrix::colSums(curr_seurat_object@assays$RNA@counts[mito_genes, ])/Matrix::colSums(curr_seurat_object@assays$RNA@counts)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
curr_seurat_object <- AddMetaData(object = curr_seurat_object, metadata = percent_mito, col.name = "percent_mito")





# ---------------------------------------------------------------------
# Export data
# ---------------------------------------------------------------------


saveRDS(curr_seurat_object, file = glue("{prefix}seurat_object.rds"))
saveRDS(curr_seurat_object@meta.data, file = glue("{prefix}cell_metadata.rds"))
write_tsv(as_tibble(curr_seurat_object@meta.data), glue("{prefix}cell_metadata.txt"))





#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
