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
library(BUSpaRse)
library(Seurat)
library(DropletUtils)
library(Matrix)
library(ggpointdensity)
library(scico)
library(scales)




#######################################################################
# Script parameters
#######################################################################



proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "011_"
out <- glue("{prefix}count_gex")
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
# Get fastq paths
# ---------------------------------------------------------------------


fastqs <- samples %>%
    dplyr::filter(seq_library == "GEX")



#######################################################################
# Use kallisto and bustools (kb) to count GEX
#######################################################################

# Follow this tutorial:
# https://www.kallistobus.tools/tutorials/kb_getting_started/r/kb_intro_2_r/


if (!dir.exists(glue("{out_dir}/kallisto_bustools"))) {dir.create(glue("{out_dir}/kallisto_bustools"), recursive = TRUE)}
setwd(glue("{out_dir}/kallisto_bustools"))


# Download mouse index
system_kb_index <- toddr::robust_system(glue("kb ref -d mouse -i index.idx -g t2g.txt"))
toddr::robust_system_check(system_kb_index)


# https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-
# The 10X Genomics kit used, was:
# "Single Cell 5' R2-only", "SC5P-R2",
# "whitelist": "737K-august-2016"
# Version 2 chemistry
# # 737k-august-2016.txt:	Single Cell 3' v2, Single Cell 5' v1 and v2, Single Cell 5' HT v2
# I1 means sample index, R1 means barcode and UMI, and R2 means the piece of cDNA -- only R1 and R2 are needed for analysis below.
# Count reads per gene
system_kb_count <- toddr::robust_system(glue("kb count -i index.idx -g t2g.txt -x 10xv2 -t {curr_threads} -o . {fastqs$fastq_r1} {fastqs$fastq_r2}"))
toddr::robust_system_check(system_kb_count)



#######################################################################
# Basic QC
#######################################################################


# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
    dir <- normalizePath(dir, mustWork = TRUE)
    m <- readMM(paste0(dir, "/", name, ".mtx"))
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # The matrix read has cells in rows
    ge <- ".genes.txt"
    genes <- readLines(file(paste0(dir, "/", name, ge)))
    barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
    colnames(m) <- barcodes
    rownames(m) <- genes
    return(m)
}

res_mat <- read_count_output("counts_unfiltered", name = "cells_x_genes")
saveRDS(res_mat, glue("{prefix}res_mat.rds"))
dim(res_mat)
tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)
bc_rank <- DropletUtils::barcodeRanks(res_mat, lower = 10)



# ---------------------------------------------------------------------
# Knee plot
# ---------------------------------------------------------------------

#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(rank, total)) +
    geom_line() +
    geom_vline(aes(xintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_hline(aes(yintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "UMI Counts", x = "Cell Barcodes", title = glue("Estimated Number of cells: {round(annot$rank_cutoff, 0)}"))
  return(p)
}



options(repr.plot.width = 9, repr.plot.height = 6)
pdf(glue("{prefix}knee_plot.pdf"))
print(knee_plot(bc_rank))
dev.off()

res_mat2 <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
res_mat3 <- res_mat2[Matrix::rowSums(res_mat2) > 0,]
print(dim(res_mat3))

saveRDS(res_mat2, glue("{prefix}res_mat2.rds"))
saveRDS(res_mat3, glue("{prefix}res_mat3.rds"))


# seu <- CreateSeuratObject(res_mat3, min.cells = 3, min.features = 200)







#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
