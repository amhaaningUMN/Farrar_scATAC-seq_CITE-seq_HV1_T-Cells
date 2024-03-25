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
prefix <- "071_"
out <- glue("{prefix}plot_and_tally")
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
# Create SLURM and R job files
# ---------------------------------------------------------------------



for (i in seq_along(gemset_name)) {
    for (j in seq_along(norm_method)) {
        for (k in seq_along(resolutions)) {
            curr_gemset_name <- gemset_name[i]
            curr_norm_method <- norm_method[j]
            curr_resolution <- resolutions[k]
            curr_res_name <- glue("SCT_snn_res_{curr_resolution}")

            if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {
                dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)
            }

            # Use the following for SLURM and R filenames
            slurm_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm")
            r_file <- glue("{proj_dir}/code/{out}.slurm_jobs/{prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.R")



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


Rscript --vanilla ${PROJ_DIR}/code/<<out>>.slurm_jobs/<<basename(r_file)>>



', .open = "<<", .close = ">>")

            long_string_slurm_vect <- str_split(long_string_slurm, "\\n")[[1]]

            con <- file(slurm_file, "w")
            for (m in seq_along(long_string_slurm_vect)) {
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
library(cowplot)
library(openxlsx)
library(viridis)
library(ggrepel)
library(plotly)



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

curr_res_name_final <- glue("{curr_res_name}_de_merge_final")
col_cluster <- paste0(curr_res_name, "_de_merge_final")
col_sample <- "sample_name"


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


# ---------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------


curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/041_pca_umap/{curr_gemset_name}/{curr_norm_method}/041_seurat_object.rds"))

curr_seurat_object@meta.data <- readRDS(glue("{proj_dir}/code_out/051_cluster_de/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/cell_metadata.rds")) %>%
    dplyr::mutate(Cluster = !!sym(col_cluster)) %>%
    dplyr::mutate(Cluster = fct_relevel(Cluster, gtools::mixedsort)) %>%
    dplyr::mutate(sample_name = fct_relevel(sample_name, gtools::mixedsort)) %>%
    dplyr::mutate(sample_name_long = fct_relevel(sample_name_long, gtools::mixedsort))


DefaultAssay(curr_seurat_object) <- "SCT"
Idents(curr_seurat_object) <- col_cluster


if (length(levels(curr_seurat_object@meta.data$sample_name)) == 1) {
    single_sample <- TRUE
} else {
    single_sample <- FALSE
}


# ---------------------------------------------------------------------
# Seurat style heatmap
# ---------------------------------------------------------------------


# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- Seurat::FindAllMarkers(curr_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 10, order_by = avg_log2FC)

pdf(glue("{prefix}top_10_genes_per_cluster_PurYel.pdf"), width = 15, height = 10)
print(
DoHeatmap(curr_seurat_object, assay = "SCT", features = top10$gene, raster = FALSE) +
    scale_color_discrete(name = "Cluster") +
    theme(legend.position = "bottom")
)
dev.off()


pdf(glue("{prefix}top_10_genes_per_cluster_RedBlu.pdf"), width = 15, height = 10)
print(
DoHeatmap(curr_seurat_object, assay = "SCT", features = top10$gene, raster = FALSE) +
    scale_color_discrete(name = "Cluster") +
    theme(legend.position = "bottom") +
    scale_fill_gradientn(colors = CustomPalette(low = "dodgerblue", high = "firebrick1", mid = "white", k = 50))
)
dev.off()




# ---------------------------------------------------------------------
# 3D UMAP
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


curr_seurat_object3D <- RunUMAP(curr_seurat_object, assay = DefaultAssay(curr_seurat_object), reduction = "pca", dims = 1:30, verbose = FALSE, reduction.name = "umap", n.components = 3L)

umap3d <- cbind(curr_seurat_object3D@reductions$umap@cell.embeddings, curr_seurat_object3D@meta.data) %>%
    as_tibble() %>%
    # Duplicate final cluster ids with simple colname
    dplyr::mutate(Cluster = !!sym(col_cluster))
# Consider writing out data for using with Embedding Projector website
# write_tsv(umap3d, glue("{prefix}umap_3D_data.tsv"))
rm(curr_seurat_object3D)
toddr::clean_mem()
gg_color_hue <- function(n, alpha = 1) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
n <- length(levels(umap3d$Cluster))
my_colors <- gg_color_hue(n, alpha = 1)
p <- plot_ly(data = umap3d, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    type = "scatter3d", mode = "markers", marker = list(size = 2),
    color = ~Cluster, colors = my_colors) %>%
    layout(title = glue("{col_cluster}"))
htmlwidgets::saveWidget(as_widget(p, height = 500, width = 1000), glue("{prefix}umap_3d.html"))






# ---------------------------------------------------------------------
# UMAP highlighted by clusters
# ---------------------------------------------------------------------



umap_data <- cbind(curr_seurat_object@reductions$umap@cell.embeddings, curr_seurat_object@meta.data) %>%
    as_tibble()



gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(umap_data$Cluster)))


centers <- umap_data %>%
    dplyr::select(Cluster, starts_with("UMAP_")) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize_all(median)

umap_cluster_ggplot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.5) +
    labs(color = "Cluster", title = paste0(curr_gemset_name, "\\n", curr_res_name)) +
    geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
    geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
    theme_cowplot() +
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))


# Generate a plot for each cluster, where the points are highlighted for each cluster
umap_clusters_highlight <- umap_cluster_ggplot +
    gghighlight::gghighlight(keep_scales = TRUE, use_direct_label = FALSE, unhighlighted_params = list(color = "gray95")) +
    facet_wrap(~ Cluster) +
    # Remove the facet labels
    theme(strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          strip.background = element_blank())

pdf(glue("{prefix}umap_highlight_by_cluster.pdf"), width = 8, height = 8)
print(umap_clusters_highlight)
dev.off()




# Individual plots
if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_cluster"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_cluster"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_cluster"))

umap_clusters_highlight_ind_plots <- list()
for (i in seq_len(length(levels(umap_data$Cluster)))) {
    curr_cluster <- levels(umap_data$Cluster)[i]
    curr_color <- cluster_colors[i]
    unhighlight_color <- "gray95"

    umap_data2 <- umap_data %>%
        dplyr::mutate(Cluster2 = case_when(
            Cluster != curr_cluster ~ "",
            TRUE ~ as.character(Cluster))) %>%
        dplyr::mutate(Cluster2 = factor(Cluster2, levels = c(curr_cluster, ""))) %>%
        dplyr::select(-Cluster) %>%
        dplyr::rename(Cluster = "Cluster2") %>%
        dplyr::arrange(desc(Cluster))

    centers2 <- centers %>%
        dplyr::filter(Cluster == curr_cluster)

    umap_clusters_highlight_ind_plots[[i]] <- ggplot(umap_data2, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = c(curr_color, unhighlight_color)) +
        labs(color = "Cluster", title = paste0(curr_gemset_name, "\\n", curr_res_name, "\\n", "Cluster: ", curr_cluster)) +
        geom_point(data = centers2, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
        geom_text(data = centers2, mapping = aes(label = Cluster), size = 4, color = "black") +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5),
            # Ensure the plots stay square
            aspect.ratio = 1,
            legend.position = "none")

    pdf(glue("{prefix}umap_highlight_by_cluster{curr_cluster}.pdf"), width = 6, height = 6)
    print(umap_clusters_highlight_ind_plots[[i]])
    dev.off()
}










# ---------------------------------------------------------------------
# sample_name, cluster
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))

if (single_sample) {
    curr_width <- 6
    curr_height <- 5
} else {
    curr_width <- 22
    curr_height <- 5
}

p <- umap_data %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
        geom_point(size = 0.5) +
        labs(title = paste0(curr_gemset_name, "\\n", curr_res_name), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_grid(cols = vars(sample_name)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank(),
            legend.position = "none")
pdf(glue("{prefix}umap_by_sample_name_cluster.pdf"), width = curr_width, height = curr_height)
print(p)
dev.off()





# ---------------------------------------------------------------------
# cell cycle
# ---------------------------------------------------------------------

if (single_sample) {
    curr_width <- 6
    curr_height <- 5
} else {
    curr_width <- 22
    curr_height <- 5
}



p <- umap_data %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = Phase)) +
        geom_point(size = 0.5) +
        labs(title = paste0(curr_gemset_name, "\\n", curr_res_name), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_grid(cols = vars(sample_name)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_by_sample_name_cell_cycle_phase.pdf"), width = curr_width, height = curr_height)
print(p)
dev.off()







d1 <- umap_data


# Duplicate the df, adding facet variable and data
s_phase <- d1 %>%
    dplyr::mutate(cell_cycle_facet = "S.Score") %>%
    dplyr::mutate(cell_cycle_facet_score = S.Score)

g2m_phase <- d1 %>%
    dplyr::mutate(cell_cycle_facet = "G2M.Score") %>%
    dplyr::mutate(cell_cycle_facet_score = G2M.Score)

d2 <- bind_rows(s_phase, g2m_phase)


p <- d2 %>%
    # Order by high cell cycle score levels
    dplyr::arrange(cell_cycle_facet_score) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell_cycle_facet_score)) +
        geom_point(size = 0.5) +
        scale_color_viridis(name = "Score") +
        labs(title = paste0(curr_gemset_name, "\\n", curr_res_name), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "grey80") +
        facet_grid(cols = vars(cell_cycle_facet)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

pdf(glue("{prefix}umap_by_cell_cycle_scores.pdf"), width = 11, height = 6)
print(p)
dev.off()



if (single_sample) {
    curr_width <- 6
    curr_height <- 9
} else {
    curr_width <- 22
    curr_height <- 9
}

p <- d2 %>%
    # Order by high cell cycle score levels
    dplyr::arrange(cell_cycle_facet_score) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell_cycle_facet_score)) +
        geom_point(size = 0.5) +
        scale_color_viridis(name = "Score") +
        labs(title = paste0(curr_gemset_name, "\\n", curr_res_name), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "grey80") +
        facet_grid(cols = vars(sample_name), rows = vars(cell_cycle_facet)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

pdf(glue("{prefix}umap_by_cell_cycle_scores_sample_name.pdf"), width = curr_width, height = curr_height)
print(p)
dev.off()








# ---------------------------------------------------------------------
# UMAP by UMI COUNTS
# ---------------------------------------------------------------------

# Use a NULL slot because you are not plotting expression data, but instead a column from meta.data
# nCount_RNA equals the colSums() of the Assay = RNA, slot = counts
p1 <- FeaturePlot(curr_seurat_object, features = "nCount_RNA", order = TRUE, cols = viridis(20))
# nCount_SCT what will you just have to equals the colSums() of the Assay = SCT, slot = counts
p2 <- FeaturePlot(curr_seurat_object, features = "nCount_SCT", order = TRUE, cols = viridis(20))
pdf(glue("{prefix}umap_umi_count_lognorm_vs_sct.pdf"), width = 11, height = 6)
print(
cowplot::plot_grid(plotlist = list(p1, p2), align = "hv", ncol = 2)
)
dev.off()



#######################################################################
# Cell tallies
#######################################################################

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/cell_counts"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/cell_counts"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/cell_counts"))



# ---------------------------------------------------------------------
# Bar plot
# ---------------------------------------------------------------------


out_filename_prefix <- glue("{prefix}cells_per_sample_name_and_cluster")
cell_numbers_tbl1 <- umap_data %>%
    dplyr::mutate(cluster = factor(!!rlang::sym("Cluster"))) %>%
    dplyr::group_by(sample_name, cluster, .drop = FALSE) %>%
    tally(name = "cells") %>%
    dplyr::mutate(prop_cells_per_cluster = round(cells / sum(cells), 3))

write_tsv(cell_numbers_tbl1, glue("{out_filename_prefix}.txt"))
wb <- openxlsx::write.xlsx(list(cell_numbers = cell_numbers_tbl1), file = glue("{out_filename_prefix}.xlsx"))
setColWidths(wb, sheet = 1, cols = 3, widths = 12)
setColWidths(wb, sheet = 1, cols = 4, widths = 12)
saveWorkbook(wb, glue("{out_filename_prefix}.xlsx"), overwrite = TRUE)





tally_data <- umap_data %>%
    group_by(sample_name, Cluster, .drop = FALSE) %>%
    tally() %>%
    ungroup() %>%
    group_by(sample_name) %>%
    dplyr::mutate(sample_name_sum = sum(n)) %>%
    dplyr::mutate(sample_name_prop = n / sample_name_sum)

p <- tally_data %>%
    ggplot(aes(x = Cluster, y = sample_name_prop)) +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
        facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 6)) +
        coord_flip() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank(),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5),
            panel.border = element_blank()) +
            labs(title = "Percent of cells for each sample_name", subtitle = "(Sum of all bars within a sample_name, equals 100%)",
                x = "Cluster", y = "Percent of cells")
pdf(glue("{out_filename_prefix}_prop.pdf"), width = 16, height = 5)
print(p)
dev.off()


p <- tally_data %>%
    ggplot(aes(x = Cluster, y = n)) +
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
        facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
        coord_flip() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank(),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5),
            panel.border = element_blank()) +
            labs(title = "Number of cells at each time point",
                x = "Cluster", y = "Number of cells")
pdf(glue("{out_filename_prefix}_counts.pdf"), width = 16, height = 5)
print(p)
dev.off()






#######################################################################
# UMAP highlight
#######################################################################

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))



# ---------------------------------------------------------------------
# UMAP, highlight by sample_name
# ---------------------------------------------------------------------

n <- length(levels(umap_data$sample_name))
my_colors <- gg_color_hue(n, alpha = 1)

p <- umap_data %>%
    dplyr::arrange(sample_name) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, group = Cluster)) +
    geom_point(aes(color = sample_name), size = 0.5, alpha = 0.5, stroke = 0) +
    labs(color = "sample_name", title = paste0(curr_gemset_name, "\\n", curr_res_name)) +
    theme_cowplot() +
    gghighlight::gghighlight(keep_scales = TRUE, use_direct_label = FALSE, unhighlighted_params = list(color = "gray95")) +
    scale_color_manual(name = "sample_name", values = my_colors) +
    facet_wrap(~ Cluster) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    # Remove the facet labels
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "right")

pdf(glue("{prefix}umap_highlight_by_sample_name_cluster.pdf"), width = 8, height = 10)
print(p)
dev.off()



# Individual plots
if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_sample_name_cluster"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_sample_name_cluster"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_highlight_by_sample_name_cluster"))

plot_list <- list()
for (i in seq_len(length(levels(umap_data$Cluster)))) {
    curr_cluster <- levels(umap_data$Cluster)[i]

    plot_list[[i]] <- umap_data %>%
        dplyr::mutate(sample_name2 = case_when(
            Cluster != curr_cluster ~ NA_character_,
            TRUE ~ as.character(sample_name))) %>%
        dplyr::mutate(sample_name2 = addNA(fct_relevel(sample_name2, gtools::mixedsort))) %>%
        dplyr::arrange(desc(sample_name2)) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = sample_name2), size = 0.5) +
            scale_color_discrete(na.value = "gray95") +
            labs(title = paste0(curr_gemset_name, "\\n", curr_res_name, "\nCluster: ", curr_cluster)) +
            theme_cowplot() +
            facet_wrap(vars(sample_name), nrow = 1) +
            theme(plot.title = element_text(hjust = 0.5),
                aspect.ratio = 1,
                legend.position = "none")

    pdf(glue("{prefix}umap_highlight_by_sample_name_cluster_{curr_cluster}.pdf"), width = 16, height = 3)
    print(plot_list[[i]])
    dev.off()
}









#######################################################################
# Save session info
#######################################################################



# ---------------------------------------------------------------------
# Write out session info
# ---------------------------------------------------------------------


toddr::write_session_info(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{prefix}"))



', .open = "<<", .close = ">>")
            long_string_r_vect <- str_split(long_string_r, "\\n")[[1]]

            con <- file(r_file, "w")
            for (m in seq_along(long_string_r_vect)) {
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
