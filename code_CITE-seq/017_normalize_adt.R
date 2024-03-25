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
library(Matrix)
library(dsb)



#######################################################################
# Script parameters
#######################################################################


proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "017_"
out <- glue("{prefix}normalize_adt")
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


#######################################################################
# Normalize ADT data
#######################################################################


# Use the dsb R package for ADT data normalization

# Plan:
# get raw counts
# subset ADT and rawGEX to only include only GEM barcodes that intersect
# generate some QC stats (mito, total RNA umis, total ADT umis)




# This section was run interactively first to manually define the parameters below
thresholds <- list(adt_raw_lib_size_lower = 0.2, adt_raw_lib_size_upper = 1.5, gex_raw_lib_size_upper = 2.0)


# Seurat CLR function
# https://github.com/satijalab/seurat/blob/a1294c4d363780548dbf9cc4a4abb3a6078a6d64/R/preprocessing.R#L2484-L2486
clr_function <- function(x) {
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}



# Get data
adt_raw <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_adt_raw.rds"))
gex_raw <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_gex_raw.rds"))
gex_filtered <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_gex_filtered.rds"))


# Keep only barcodes in ADT dataset that are also found in "raw" GEX data
# The raw GEX barcode list is larger
adt_barcodes_found_in_gex_raw <- intersect(colnames(adt_raw$matrix), colnames(gex_raw$matrix))
# Make sure order is same for both, by specifying the order
gex_raw_matrix <- gex_raw$matrix[, adt_barcodes_found_in_gex_raw]
adt_raw_matrix <- adt_raw$matrix[, adt_barcodes_found_in_gex_raw]


# ---------------------------------------------------------------------
# Droplet quality control metadata
# ---------------------------------------------------------------------

# Create metadata of droplet QC stats used in standard scRNAseq processing

# Find the mitochondrial genes and ensembl ids
gex_features_mt <- gex_raw$features %>%
    dplyr::filter(str_detect(symbol, "^mt-"))

meta_data <- tibble(barcode = adt_barcodes_found_in_gex_raw) %>%
    dplyr::mutate(gex_raw_lib_size = log10(Matrix::colSums(gex_raw_matrix))) %>%
    dplyr::mutate(adt_raw_lib_size = log10(Matrix::colSums(adt_raw_matrix))) %>%
    dplyr::mutate(gex_raw_n_genes_expr = Matrix::colSums(gex_raw_matrix > 0)) %>%
    dplyr::mutate(prop_of_umi_from_mito = Matrix::colSums(gex_raw_matrix[gex_features_mt$unique_id, ]) / Matrix::colSums(gex_raw_matrix)) %>%
    # Mark GEMs if they were classified to contain a cell by kallisto
    dplyr::mutate(gem_type = if_else(adt_barcodes_found_in_gex_raw %in% colnames(gex_filtered$matrix), "gex_cell", "empty_gem")) %>%
    # Keep only barcodes with some evidence of capture in the experiment
    dplyr::filter(gex_raw_lib_size > 0 & adt_raw_lib_size > 0)




# ---------------------------------------------------------------------
# Isolate GEMs with only background ADT levels (i.e. empty gems)
# ---------------------------------------------------------------------


# This section was run interactively first to manually define the parameters below for each seq library

adt_lower_bound <- thresholds$adt_raw_lib_size_lower
adt_upper_bound <- thresholds$adt_raw_lib_size_upper
gex_upper_bound <- thresholds$gex_raw_lib_size_upper

p <- meta_data %>%
    ggplot(aes(x = adt_raw_lib_size, y = gex_raw_lib_size)) +
    facet_grid(cols = vars(gem_type)) +
    geom_bin2d(bins = 25) +
    scale_fill_gradientn(name = "Number of GEMs\nin bin", colors = c("gray95", viridis::viridis(100))) +
    theme_bw() +
    geom_hline(yintercept = gex_upper_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_lower_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_upper_bound, linetype = "dashed")

pdf(glue("{prefix}gem_qc_metrics.pdf"))
print(p)
dev.off()


p <- meta_data %>%
    ggplot(aes(x = adt_raw_lib_size, y = gex_raw_lib_size, color = gex_raw_n_genes_expr)) +
    geom_point() +
    facet_grid(cols = vars(gem_type)) +
    scale_color_gradientn(name = "Number of Genes\nExpressed", colors = c("gray95", viridis::viridis(100))) +
    theme_bw() +
    geom_hline(yintercept = gex_upper_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_lower_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_upper_bound, linetype = "dashed")

pdf(glue("{prefix}gem_qc_metrics_number_of_genes_expr.pdf"))
print(p)
dev.off()



p <- meta_data %>%
    ggplot(aes(x = adt_raw_lib_size, y = gex_raw_lib_size, color = prop_of_umi_from_mito)) +
    geom_point() +
    facet_grid(cols = vars(gem_type)) +
    scale_color_gradientn(name = "Proportion of UMIs\nfrom mito genes", colors = c("gray95", viridis::viridis(100))) +
    theme_bw() +
    geom_hline(yintercept = gex_upper_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_lower_bound, linetype = "dashed") +
    geom_vline(xintercept = adt_upper_bound, linetype = "dashed")

pdf(glue("{prefix}gem_qc_metrics_mito_prop.pdf"))
print(p)
dev.off()


background_gems <- meta_data %>%
    dplyr::filter(gem_type == "empty_gem") %>%
    dplyr::filter(adt_raw_lib_size > adt_lower_bound & adt_raw_lib_size < adt_upper_bound) %>%
    dplyr::filter(gex_raw_lib_size < gex_upper_bound) %>%
    pull(barcode)

adt_raw_matrix_background_gems <- adt_raw_matrix[, background_gems] %>%
    as.matrix()



# ---------------------------------------------------------------------
# Isolate GEX-cell-containing GEMs (i.e. non-empty gems)
# ---------------------------------------------------------------------



cell_gems <- meta_data %>%
    dplyr::filter(gem_type == "gex_cell") %>%
    # Calculate library QC stats
    # Remove GEMs that fall outside 4 median absolute deviations from the median library size
    dplyr::mutate(gex_lower = median(gex_raw_lib_size) - (4 * mad(gex_raw_lib_size))) %>%
    dplyr::mutate(gex_upper = median(gex_raw_lib_size) + (4 * mad(gex_raw_lib_size))) %>%
    dplyr::mutate(adt_lower = median(adt_raw_lib_size) - (4 * mad(adt_raw_lib_size))) %>%
    dplyr::mutate(adt_upper = median(adt_raw_lib_size) + (4 * mad(adt_raw_lib_size))) %>%
    dplyr::filter(gex_raw_lib_size > gex_lower) %>%
    dplyr::filter(gex_raw_lib_size < gex_upper) %>%
    dplyr::filter(adt_raw_lib_size > adt_lower) %>%
    dplyr::filter(adt_raw_lib_size < adt_upper) %>%
    dplyr::filter(prop_of_umi_from_mito < 0.2) %>%
    pull(barcode)

adt_raw_matrix_cell_gems <- adt_raw_matrix[, cell_gems] %>%
    as.matrix()

gex_raw_matrix_cell_gems <- gex_raw_matrix[, cell_gems] %>%
    as.matrix()


# ---------------------------------------------------------------------
# CLR normalization
# ---------------------------------------------------------------------

adt_raw_matrix_cell_gems_clr1 <- apply(X = adt_raw_matrix_cell_gems,
    MARGIN = 1,
    FUN = clr_function)
adt_raw_matrix_cell_gems_clr <- Matrix::t(adt_raw_matrix_cell_gems_clr1)


# ---------------------------------------------------------------------
# dsb normalization
# ---------------------------------------------------------------------

adt_raw_matrix_cell_gems_dsb <- dsb::DSBNormalizeProtein(cell_protein_matrix = adt_raw_matrix_cell_gems,
    empty_drop_matrix = adt_raw_matrix_background_gems,
    denoise.counts = TRUE,
    use.isotype.control = FALSE)



# ---------------------------------------------------------------------
# Histogram
# ---------------------------------------------------------------------


adt_raw_matrix_cell_gems_clr_long <- adt_raw_matrix_cell_gems_clr %>%
    tibble::as_tibble(rownames = "unique_id") %>%
    tidyr::pivot_longer(!unique_id) %>%
    dplyr::mutate(norm_method = "clr")



adt_raw_matrix_cell_gems_dsb_long <- adt_raw_matrix_cell_gems_dsb %>%
    tibble::as_tibble(rownames = "unique_id") %>%
    tidyr::pivot_longer(!unique_id) %>%
    dplyr::mutate(norm_method = "dsb")

adt_norm_long <- bind_rows(adt_raw_matrix_cell_gems_clr_long, adt_raw_matrix_cell_gems_dsb_long)

p <- adt_norm_long %>%
    ggplot(aes(x = value)) +
    geom_histogram() +
    facet_grid(rows = vars(unique_id), cols = vars(norm_method), scales = "free") +
    labs(title = glue("Distribution of normalized or transformed UMI counts split by ADT feature"),
        x = "Normalized or transformed expression (30 bins)",
        y = "Number of GEMs with expression value")

pdf(glue("{prefix}dsb_clr_adt_histogram.pdf"), width = 8, height = 8)
print(p)
dev.off()



# ---------------------------------------------------------------------
# Export data
# ---------------------------------------------------------------------

setwd(glue("{out_dir}"))
saveRDS(adt_raw_matrix_cell_gems_clr, glue("{prefix}adt_clr_transformed.rds"))
if (!dir.exists(glue("{out_dir}/{prefix}adt_clr_transformed"))) {dir.create(glue("{out_dir}/{prefix}adt_clr_transformed"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}adt_clr_transformed"))
write_tsv(tibble(rownames(adt_raw_matrix_cell_gems_clr)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(adt_raw_matrix_cell_gems_clr)), glue("barcodes.tsv.gz"), col_names = FALSE)
sparse_mtx <- Matrix::Matrix(adt_raw_matrix_cell_gems_clr, sparse = TRUE)
Matrix::writeMM(sparse_mtx, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))

setwd(glue("{out_dir}"))
saveRDS(adt_raw_matrix_cell_gems_dsb, glue("{prefix}adt_dsb_normalized.rds"))
if (!dir.exists(glue("{out_dir}/{prefix}adt_dsb_normalized"))) {dir.create(glue("{out_dir}/{prefix}adt_dsb_normalized"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}adt_dsb_normalized"))
write_tsv(tibble(rownames(adt_raw_matrix_cell_gems_dsb)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(adt_raw_matrix_cell_gems_dsb)), glue("barcodes.tsv.gz"), col_names = FALSE)
sparse_mtx <- Matrix::Matrix(adt_raw_matrix_cell_gems_dsb, sparse = TRUE)
Matrix::writeMM(sparse_mtx, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))

setwd(glue("{out_dir}"))
saveRDS(cell_gems, glue("{prefix}cell_gems_meta_data.rds"))
saveRDS(adt_raw_matrix_cell_gems, glue("{prefix}adt_raw_matrix_cell_gems.rds"))
saveRDS(gex_raw_matrix_cell_gems, glue("{prefix}gex_raw_matrix_cell_gems.rds"))


















#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))










