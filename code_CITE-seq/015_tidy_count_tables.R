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


#######################################################################
# Script parameters
#######################################################################



proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "015_"
out <- glue("{prefix}tidy_count_tables")
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

# Purpose:
# Export nicely formatted count tables.



# ---------------------------------------------------------------------
# Get sample sheet
# ---------------------------------------------------------------------


samples <- readRDS(glue("{proj_dir}/code_out/010_samples/010_samples.rds"))


# ---------------------------------------------------------------------
# Subset the HTO and ADT data by features of each sequencing library
# ---------------------------------------------------------------------

# Read in raw data
fb_raw_allfeatures <- readRDS(glue("{proj_dir}/code_out/012_count_fb/explore_counts/012_fb_raw_allfeatures.rds"))


hto_features <- samples %>%
    dplyr::filter(application == "hashtag") %>%
    dplyr::pull(sample_name)

adt_features <- samples %>%
    dplyr::filter(application == "adt") %>%
    dplyr::pull(sample_name)


hto_raw <- list()
adt_raw <- list()



# HTO
curr_features <- fb_raw_allfeatures$features %>%
    dplyr::filter(unique_id %in% hto_features) %>%
    dplyr::arrange(factor(unique_id, levels = hto_features))
curr_barcodes <- fb_raw_allfeatures$barcodes
curr_matrix <- fb_raw_allfeatures$matrix[hto_features, ]
# Remove any GEM barcodes with zero counts for all features
keep_barcodes <- which(colSums(curr_matrix) != 0)

hto_raw$features <- curr_features
hto_raw$barcodes <- curr_barcodes[keep_barcodes, ]
hto_raw$matrix <- curr_matrix[, keep_barcodes]

# ADT
curr_features <- fb_raw_allfeatures$features %>%
    dplyr::filter(unique_id %in% adt_features) %>%
    dplyr::arrange(factor(unique_id, levels = adt_features))
curr_barcodes <- fb_raw_allfeatures$barcodes
curr_matrix <- fb_raw_allfeatures$matrix[adt_features, ]
# Remove any GEM barcodes with zero counts for all features
keep_barcodes <- which(colSums(curr_matrix) != 0)

adt_raw$features <- curr_features
adt_raw$barcodes <- curr_barcodes[keep_barcodes, ]
adt_raw$matrix <- curr_matrix[, keep_barcodes]



dim(fb_raw_allfeatures$matrix)
dim(hto_raw$matrix)
dim(adt_raw$matrix)






#######################################################################
# Import the GEX kallisto data
#######################################################################


gex_raw <- list()
gex_filtered <- list()


txpt_genes_symbol <- read_tsv(glue("{proj_dir}/code_out/011_count_gex/kallisto_bustools/t2g.txt"), col_names = FALSE) %>%
   magrittr::set_colnames(c("ensembl_id_txpt", "ensembl_id", "symbol"))

genes_symbol <- txpt_genes_symbol %>%
    dplyr::distinct(ensembl_id, symbol)


# Import raw GEX matrix
kallisto_gex_raw_dir <- glue("{proj_dir}/code_out/011_count_gex/kallisto_bustools/counts_unfiltered")

curr_features <- read_tsv(glue("{kallisto_gex_raw_dir}/cells_x_genes.genes.txt"), col_names = FALSE) %>%
    dplyr::rename(unique_id = "X1") %>%
    dplyr::left_join(genes_symbol, by = c("unique_id" = "ensembl_id")) %>%
    dplyr::mutate(library_type = "Gene Expression")

curr_barcodes <- read_tsv(glue("{kallisto_gex_raw_dir}/cells_x_genes.barcodes.txt"), col_names = FALSE) %>%
    dplyr::rename(barcode = "X1")

curr_matrix <- Matrix::readMM(glue("{kallisto_gex_raw_dir}/cells_x_genes.mtx")) %>%
    as(object = ., Class = "CsparseMatrix") %>%
    magrittr::set_colnames(curr_features$unique_id) %>%
    magrittr::set_rownames(curr_barcodes$barcode) %>%
    t()

# Remove any GEM barcodes with zero counts for all features
keep_barcodes <- which(colSums(curr_matrix) != 0)

gex_raw$features <- curr_features
gex_raw$barcodes <- curr_barcodes[keep_barcodes, ]
gex_raw$matrix <- curr_matrix[, keep_barcodes]




# Import filtered GEX matrix
# NOTE: The filtered GEX matrix is simply a subset of the raw matrix
gex_filtered$matrix <- readRDS(glue("{proj_dir}/code_out/011_count_gex/kallisto_bustools/011_res_mat3.rds"))
gex_filtered$features <- rownames(gex_filtered$matrix)
gex_filtered$barcodes <- colnames(gex_filtered$matrix)




dim(gex_raw$matrix)
dim(gex_filtered$matrix)



#######################################################################
# Export tables
#######################################################################



# GEX (all barcodes: raw)
saveRDS(gex_raw, glue("{prefix}gex_raw.rds"))
curr_matrix <- gex_raw$matrix
if (!dir.exists(glue("{out_dir}/{prefix}gex_raw"))) {dir.create(glue("{out_dir}/{prefix}gex_raw"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}gex_raw"))
write_tsv(tibble(rownames(curr_matrix)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(curr_matrix)), glue("barcodes.tsv.gz"), col_names = FALSE)
Matrix::writeMM(curr_matrix, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))



# GEX (filtered barcodes: correspond to cell-containing-GEMs)
setwd(glue("{out_dir}"))
saveRDS(gex_filtered, glue("{prefix}gex_filtered.rds"))
curr_matrix <- gex_filtered$matrix
if (!dir.exists(glue("{out_dir}/{prefix}gex_filtered"))) {dir.create(glue("{out_dir}/{prefix}gex_filtered"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}gex_filtered"))
write_tsv(tibble(rownames(curr_matrix)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(curr_matrix)), glue("barcodes.tsv.gz"), col_names = FALSE)
Matrix::writeMM(curr_matrix, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))



# HTO
setwd(glue("{out_dir}"))
saveRDS(hto_raw, glue("{prefix}hto_raw.rds"))
curr_matrix <- hto_raw$matrix
if (!dir.exists(glue("{out_dir}/{prefix}hto_raw"))) {dir.create(glue("{out_dir}/{prefix}hto_raw"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}hto_raw"))
write_tsv(tibble(rownames(curr_matrix)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(curr_matrix)), glue("barcodes.tsv.gz"), col_names = FALSE)
Matrix::writeMM(curr_matrix, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))



# ADT
setwd(glue("{out_dir}"))
saveRDS(adt_raw, glue("{prefix}adt_raw.rds"))
curr_matrix <- adt_raw$matrix
if (!dir.exists(glue("{out_dir}/{prefix}adt_raw"))) {dir.create(glue("{out_dir}/{prefix}adt_raw"), recursive = TRUE)}
setwd(glue("{out_dir}/{prefix}adt_raw"))
write_tsv(tibble(rownames(curr_matrix)), glue("features.tsv.gz"), col_names = FALSE)
write_tsv(tibble(colnames(curr_matrix)), glue("barcodes.tsv.gz"), col_names = FALSE)
Matrix::writeMM(curr_matrix, file = glue("matrix.mtx"))
system(glue("gzip matrix.mtx"))








# ---------------------------------------------------------------------
# Compare barcode overlap
# ---------------------------------------------------------------------

setwd(out_dir)



# Function
make_all_combinations <- function(set){
      unlist(lapply(seq_along(set), function(size){
            apply(combn(set, size), 2, paste0, collapse="-")
      }))
}



barcode_tbl_gex_raw <- tibble(barcode = colnames(gex_raw$matrix),
    dataset = "gex_raw")

barcode_tbl_gex_filtered <- tibble(barcode = colnames(gex_filtered$matrix),
    dataset = "gex_filtered")

barcode_tbl_hto <- tibble(barcode = colnames(hto_raw$matrix),
    dataset = "hto_raw")

barcode_tbl_adt <- tibble(barcode = colnames(adt_raw$matrix),
    dataset = "adt_raw")

barcode_tbl_tcr <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_contigs.rds")) %>%
    distinct(barcode) %>%
    dplyr::mutate(dataset = "tcr")



# Combine tables
barcode_tbl <- barcode_tbl_gex_raw %>%
    bind_rows(barcode_tbl_gex_filtered) %>%
    bind_rows(barcode_tbl_hto) %>%
    bind_rows(barcode_tbl_adt) %>%
    bind_rows(barcode_tbl_tcr) %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarize(dataset_list = list(dataset)) %>%
    dplyr::mutate(comparisons = lapply(dataset_list, make_all_combinations)) %>%
    tidyr::unnest(comparisons)

# Find order, largest set to smallest
barcode_tbl_order <- barcode_tbl %>%
    group_by(comparisons) %>%
    tally() %>%
    arrange(desc(n)) %>%
    pull(comparisons)


p <- barcode_tbl %>%
    ggplot(aes(x = comparisons)) +
    geom_bar(fill = "dodgerblue") +
    # Place the number 5% higher than the top of the bar
    geom_text(stat = "count", size = 2, aes(label = after_stat(count), y = after_stat(count) + (0.05 * after_stat(count)))) +
    ggupset::axis_combmatrix(sep = "-") +
    scale_y_continuous(breaks = scales::breaks_extended(n = 10), labels = scales::label_comma()) +
    scale_x_discrete(limits = barcode_tbl_order) +
    labs(title = glue("Number of GEMs that intersect between datasets"),
        x = "Datasets\n(linked dots indicate set intersections of datasets)",
        y = "Number of GEMs") +
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
        combmatrix.panel.line.size = 1,
        combmatrix.label.make_space = TRUE)

pdf(glue("{prefix}gem_overlap_by_dataset.pdf"), width = 10, height = 7)
print(p)
dev.off()



#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
