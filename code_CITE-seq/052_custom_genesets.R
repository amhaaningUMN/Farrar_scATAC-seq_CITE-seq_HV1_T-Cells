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
prefix <- "052_"
out <- glue("{prefix}custom_genesets")
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
# Gemset
# ---------------------------------------------------------------------


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
# Which gene set to use for these plots
# ---------------------------------------------------------------------

# This should be the name (first column) of the geneset row, in the gmt file 
# This name gets incorpterated into "mm_msigdb_long"
custom_genesets_title <- "Farrar_custom_genesets"




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



# Get a list of genes that are interesting and will want to be included in plots
mm_msigdb_long <- read_tsv(glue("{proj_dir}/code_out/026_clean_gene_sets_and_msigdb/026_mm_msigdb_long.tsv"))




# ---------------------------------------------------------------------
# UMAP of cluster ids
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


plot_data <- as.data.frame(curr_seurat_object@reductions$umap@cell.embeddings)
plot_data$active.ident <- as.factor(curr_seurat_object@active.ident)

gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(plot_data$active.ident)))



centers <- plot_data %>%
    dplyr::group_by(active.ident) %>%
    dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

umap_cluster_ggplot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
    geom_point(size = 0.5) +
    labs(color = "Cluster", title = paste0(curr_gemset_name, "\\n", curr_res_name)) +
    geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
    geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
    theme_cowplot() +
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))



# ---------------------------------------------------------------------
# UMAP feature plot
# ---------------------------------------------------------------------



genes <- mm_msigdb_long %>%
    dplyr::filter(ont == custom_genesets_title) %>%
    # Add NA columns because I dont have this info
    dplyr::mutate(dataset = NA_character_) %>%
    dplyr::mutate(cluster_resolution_name = NA_character_) %>%
    dplyr::mutate(cluster_id = NA_character_) %>%
    # Keep only genes that were specified for this datasets/resolutions/clusters, or NA (i.e. any dataset)
    dplyr::filter(is.na(dataset) | dataset == curr_gemset_name) %>%
    dplyr::filter(is.na(cluster_resolution_name) | cluster_resolution_name == curr_res_name) %>%
    group_by(dataset, cluster_resolution_name, cluster_id) %>%
    dplyr::mutate(group_idx = cur_group_id() - 1) %>%
    dplyr::arrange(group_idx, symbol)



tk_get_subgroups <- function(group_idx, n = 15) {
    data <- tibble(group_idx = group_idx) %>%
        group_by(group_idx) %>%
        summarize(colors = list(group_idx))
    all_out <- list()
    for (i in seq_len(nrow(data))) {
        mylist <- split(data$colors[[i]], ceiling(seq_along(data$colors[[i]])/n))
        if (length(mylist) > 1) {
            out <- list()
            for (j in seq_along(names(mylist))) {
                out[[j]] <- paste0(mylist[[j]], "_", names(mylist)[j])
            }
            all_out[[i]] <- unname(unlist(out))
        } else {
            all_out[[i]] <- unname(unlist(mylist))
        }
    }
    unname(unlist(all_out))
}
genes$group_idx_sub <- tk_get_subgroups(genes$group_idx)

genes2 <- genes %>%
    dplyr::mutate(group_name = case_when(
        !is.na(cluster_id) ~ paste0("cluster_", group_idx_sub),
        is.na(cluster_id) ~ paste0("set_", str_remove(group_idx_sub, "\\\\d+_")),
        TRUE ~ "misc"))


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))

write_tsv(genes2, glue("{prefix}genes_by_cluster_with_labels.txt"))
openxlsx::write.xlsx(as.data.frame(genes2), file = glue("{prefix}genes_by_cluster_with_labels.xlsx"))




for (i in seq_along(unique(genes2$group_idx_sub))) {
    curr_group <- unique(genes2$group_idx_sub)[i]
    curr_genes <- genes2 %>%
        dplyr::filter(group_idx_sub == curr_group) %>%
        dplyr::pull(symbol)

    curr_group_name <- genes2 %>%
        dplyr::filter(group_idx_sub == curr_group) %>%
        dplyr::pull(group_name) %>%
        unique(.)

    # UMAP
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_umap"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_umap"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_umap"))
    tk_gene_plots(curr_seurat_object, curr_genes, paste0(curr_group_name), all_together = TRUE, number_of_ind_plots = Inf, umap_ggplot_object = umap_cluster_ggplot, type = "feature_plot", facet_by = NULL)


    # Violin
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_violin"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_violin"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_violin"))
    tk_gene_plots(curr_seurat_object, curr_genes, paste0(curr_group_name), all_together = TRUE, number_of_ind_plots = Inf, umap_ggplot_object = umap_cluster_ggplot, type = "violin_plot", facet_by = NULL)


    # Box
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_box"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_box"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_box"))
    tk_gene_plots(curr_seurat_object, curr_genes, paste0(curr_group_name), all_together = TRUE, number_of_ind_plots = Inf, umap_ggplot_object = umap_cluster_ggplot, type = "box_plot", facet_by = NULL)


    # Ridge
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_ridge"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_ridge"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_ridge"))
    tk_gene_plots(curr_seurat_object, curr_genes, paste0(curr_group_name), all_together = TRUE, number_of_ind_plots = Inf, umap_ggplot_object = umap_cluster_ggplot, type = "ridge_plot", facet_by = NULL)


    # Dot
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_dot"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_dot"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_dot"))
    tk_gene_plots(curr_seurat_object, curr_genes, paste0(curr_group_name), all_together = TRUE, number_of_ind_plots = Inf, umap_ggplot_object = umap_cluster_ggplot, type = "dot_plot", facet_by = NULL)



    # Heatmap
    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"))
    pdf(glue("{curr_group_name}_PurYel.pdf"), width = 15, height = 7)
    print(
    DoHeatmap(curr_seurat_object, assay = "SCT", features = curr_genes, raster = FALSE) + scale_color_discrete(name = "Cluster") + theme(legend.position = "bottom")
    )
    dev.off()

    pdf(glue("{curr_group_name}_RedBlu.pdf"), width = 15, height = 7)
    print(
    DoHeatmap(curr_seurat_object, assay = "SCT", features = curr_genes, raster = FALSE) + scale_color_discrete(name = "Cluster") + scale_fill_gradientn(colors = CustomPalette(low = "dodgerblue", high = "firebrick1", mid = "white", k = 50)) + theme(legend.position = "bottom")
    )
    dev.off()

}


# Heatmap (all genes from all sets)
if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/gene_expr_heatmap"))


curr_genes <- genes2 %>%
    dplyr::pull(symbol)

pdf(glue("all_custom_genes_PurYel.pdf"), width = 15, height = 10)
print(
DoHeatmap(curr_seurat_object, assay = "SCT", features = curr_genes, raster = FALSE) + scale_color_discrete(name = "Cluster") + theme(legend.position = "bottom")
)
dev.off()

pdf(glue("all_custom_genes_RedBlu.pdf"), width = 15, height = 10)
print(
DoHeatmap(curr_seurat_object, assay = "SCT", features = curr_genes, raster = FALSE) + scale_color_discrete(name = "Cluster") + scale_fill_gradientn(colors = CustomPalette(low = "dodgerblue", high = "firebrick1", mid = "white", k = 50)) + theme(legend.position = "bottom")
)
dev.off()



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
