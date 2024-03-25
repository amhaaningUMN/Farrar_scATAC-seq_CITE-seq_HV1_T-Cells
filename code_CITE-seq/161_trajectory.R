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
prefix <- "161_"
out <- glue("{prefix}trajectory")
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
    list(gemset_name = "adtCD4pos", norm_method = "cc_none_gene_none", resolutions = "0.5")
    )



# ---------------------------------------------------------------------
# Create SLURM and R job files
# ---------------------------------------------------------------------


for (i in seq_along(dataset_list)) {
    curr_gemset_name <- dataset_list[[i]]$gemset_name
    curr_norm_method <- dataset_list[[i]]$norm_method
    curr_resolution <- dataset_list[[i]]$resolutions
    curr_res_name <- glue("SCT_snn_res_{curr_resolution}")


    #######################################################################
    # Create job files for each individual comparison
    #######################################################################

    if (!dir.exists(glue("{proj_dir}/code/{out}.slurm_jobs"))) {dir.create(glue("{proj_dir}/code/{out}.slurm_jobs"), recursive = TRUE)}

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
library(openxlsx)
library(Seurat)
library(SeuratWrappers)
library(leidenbase)
library(GenomicRanges)
library(slingshot)
library(mclust)
library(RColorBrewer)
library(tidyverse)
library(scales)
library(tradeSeq)
library(scater)
library(NMF)
library(BiocParallel)
library(cowplot)





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


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))



curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/041_pca_umap/{curr_gemset_name}/{curr_norm_method}/041_seurat_object.rds"))

curr_seurat_object@meta.data <- readRDS(glue("{proj_dir}/code_out/051_cluster_de/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/051_cell_metadata.rds")) %>%
    dplyr::mutate(Cluster = seurat_clusters) %>%
    dplyr::mutate(Cluster = fct_relevel(Cluster, gtools::mixedsort)) %>%
    dplyr::mutate(sample_name = fct_relevel(sample_name, gtools::mixedsort)) %>%
    dplyr::mutate(sample_name_long = fct_relevel(sample_name_long, gtools::mixedsort))


DefaultAssay(curr_seurat_object) <- "SCT"
Idents(curr_seurat_object) <- "sample_name"


# ---------------------------------------------------------------------
# Determine barcodes that represent the root of differentiation
# ---------------------------------------------------------------------

# Get a vector of barcodes that should be the root of tree
naive_root_cluster <- 4
naive_root_cell_barcodes <- read_tsv(glue("{proj_dir}/code_out/051_cluster_de/adtCD4pos/cc_none_gene_none/SCT_snn_res_0.5/final_cluster_labels.txt")) %>%
    dplyr::select(cell_id, clust_de_merge_final) %>%
    dplyr::rename(barcode = "cell_id", Cluster = "clust_de_merge_final") %>%
    dplyr::filter(Cluster == naive_root_cluster) %>%
    dplyr::pull(barcode)

# Which cluster (from this dataset) contains the majority of these naive cells (originally determined from cluster 4, SCT_snn_res_0.5, adtCD4pos)?
num_barcode_overlap_with_naive <- vector()
for (i in seq_along(levels(curr_seurat_object@meta.data$Cluster))) {
    curr_cluster_barcodes <- curr_seurat_object@meta.data %>%
        dplyr::filter(Cluster == levels(curr_seurat_object@meta.data$Cluster)[i]) %>%
        dplyr::pull(barcode)
    
    num_barcode_overlap_with_naive[i] <- length(intersect(naive_root_cell_barcodes, curr_cluster_barcodes))
}
# This is the corresponding cluster in this dataset
root_cluster <- levels(curr_seurat_object@meta.data$Cluster)[which.max(num_barcode_overlap_with_naive)]
root_cell_barcodes <- curr_seurat_object@meta.data %>%
    dplyr::filter(Cluster == root_cluster) %>%
    dplyr::pull(barcode)



write_lines(glue("The root_cluster used in this dataset was: {root_cluster}. This cluster is most similar to the niave cluster (cluster 5) from the SCT_snn_res_0.8, TcellExpr, dataset."), glue("root_cluster_{root_cluster}.txt"))
 



       
#######################################################################
# SingleCellExperiment object and Slingshot trajectory
#######################################################################

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/slingshot"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/slingshot"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/slingshot"))


DefaultAssay(curr_seurat_object) <- "SCT"

sce1 <- Seurat::as.SingleCellExperiment(curr_seurat_object)
identical(assays(sce1)$counts, curr_seurat_object@assays$SCT@counts)
# [1] TRUE
identical(assays(sce1)$logcounts, curr_seurat_object@assays$SCT@data)
# [1] TRUE



sce1_na_clusters <- sce1
# Set certain cells to NA clusters (i.e. clusters we do not want to allow the tool to view.) 
x <- colData(sce1_na_clusters) %>% 
    as.data.frame() %>%
    dplyr::mutate(Cluster_null = factor(case_when(
        as.character(Cluster) %in% as.character(c(9)) ~ NA_character_,
        TRUE ~ as.character(Cluster))))
colData(sce1_na_clusters)$Cluster_null <- x$Cluster_null




# Get a trajectory
# https://bustools.github.io/BUS_notebooks_R/slingshot.html#slingshot
sce2 <- slingshot(sce1, 
    clusterLabels = "Cluster", 
    reducedDim = "UMAP",
    start.clus = root_cluster,
    stretch = 0)



sce2_na_clusters <- slingshot(sce1_na_clusters,
    clusterLabels = "Cluster_null",
    reducedDim = "UMAP",
    start.clus = root_cluster,
    stretch = 0)

    
class(sce2)
# class(sce2_na_clusters)
# [1] "SingleCellExperiment"
# attr(,"package")
# [1] "SingleCellExperiment"



# ---------------------------------------------------------------------
# UMAP sce2
# ---------------------------------------------------------------------

# reducedDim(sce2, "UMAP")
umap_data1 <- sce2@int_colData@listData$reducedDims$UMAP
# The DFrame contains a non-standard column that cannot be coerced
umap_data2 <- subset(colData(sce2), select = -slingshot)
umap_data <- cbind(as.data.frame(umap_data1), as_tibble(umap_data2)) %>%
    as_tibble()
    

gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(umap_data$Cluster)))

centers <- umap_data %>%
    dplyr::select(Cluster, starts_with("UMAP_")) %>%		    
    dplyr::select(Cluster, starts_with("UMAP_")) %>%
    dplyr::group_by(Cluster) %>%		    
    dplyr::group_by(Cluster) %>%    
    dplyr::summarize_all(median) %>%
    dplyr::mutate(Lineage = NA_character_)

# Add a lineage label to final cluster
for (i in seq_along(slingLineages(sce2))) {
    last_node_in_lineage <- tail(slingLineages(sce2)[[i]], n = 1)
    centers$Lineage[which(centers$Cluster == last_node_in_lineage)] <- names(slingLineages(sce2))[i]
}
lineage_colors <- gg_color_hue(length(seq_along(slingLineages(sce2))))
lineage_names <- names(slingLineages(sce2))
    

    

p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(color = "Cluster", title = glue("Slingshot\\nRoot cluster: {root_cluster}\\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}")) +
    theme_cowplot() +
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))
#     # Adding the curves
#     for (i in seq_along(slingCurves(sce2))) {
#         curve_i <- slingCurves(sce2)[[i]]
#         curve_i <- curve_i$s[curve_i$ord, ]
#         p <- p + geom_path(data = as.data.frame(curve_i), col = lineage_colors[i], linewidth = 1)
#     }
    # Adding the lines
    for (i in seq_along(slingLineages(sce2))) {
        lineage_i <- slingLineages(sce2)[[i]]
        lineage_i2 <- centers %>%
            dplyr::filter(Cluster %in% lineage_i) %>%
            dplyr::arrange(match(Cluster, lineage_i))
        # See ?grid::arrow for options
        p <- p + geom_path(data = as.data.frame(lineage_i2), col = lineage_colors[i], linewidth = 1, position = position_jitter(w = 0.2, h = 0.2), arrow = arrow(angle = 15, ends = "last", type = "closed"))
    }
    # Add labels
    p <- p + geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black")
    p <- p + ggrepel::geom_label_repel(data = centers, mapping = aes(label = Lineage), size = 4, color = "black") +
        # Make the legend dots larger
        guides(color = guide_legend(override.aes = list(size = 8)))


pdf(glue("{prefix}slingshot_pseudotime_trajectory.pdf"), width = 8, height = 8)
print(p)
dev.off()



# ---------------------------------------------------------------------
# UMAP sce2_na_clusters
# ---------------------------------------------------------------------
# reducedDim(sce2_na_clusters, "UMAP")
umap_data1 <- sce2_na_clusters@int_colData@listData$reducedDims$UMAP
# The DFrame contains a non-standard column that cannot be coerced
umap_data2 <- subset(colData(sce2_na_clusters), select = -slingshot)
umap_data <- cbind(as.data.frame(umap_data1), as_tibble(umap_data2)) %>%
    as_tibble()

gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(umap_data$Cluster)))

centers <- umap_data %>%
    dplyr::select(Cluster, starts_with("UMAP_")) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize_all(median) %>%
    dplyr::mutate(Lineage = NA_character_)
# Add a lineage label to final cluster
for (i in seq_along(slingLineages(sce2_na_clusters))) {
    last_node_in_lineage <- tail(slingLineages(sce2_na_clusters)[[i]], n = 1)
    centers$Lineage[which(centers$Cluster == last_node_in_lineage)] <- names(slingLineages(sce2_na_clusters))[i]
}

lineage_colors <- gg_color_hue(length(seq_along(slingLineages(sce2_na_clusters))))
lineage_names <- names(slingLineages(sce2_na_clusters))

p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(color = "Cluster", title = glue("Slingshot\\nRoot cluster: {root_cluster}\\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}")) +
    theme_cowplot() +
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))
#     # Adding the curves
#     for (i in seq_along(slingCurves(sce2_na_clusters))) {
#         curve_i <- slingCurves(sce2_na_clusters)[[i]]
#         curve_i <- curve_i$s[curve_i$ord, ]
#         p <- p + geom_path(data = as.data.frame(curve_i), col = lineage_colors[i], linewidth = 1)
#     }
    # Adding the lines
    for (i in seq_along(slingLineages(sce2_na_clusters))) {
        lineage_i <- slingLineages(sce2_na_clusters)[[i]]
        lineage_i2 <- centers %>%
            dplyr::filter(Cluster %in% lineage_i) %>%
            dplyr::arrange(match(Cluster, lineage_i))
        # See ?grid::arrow for options
        p <- p + geom_path(data = as.data.frame(lineage_i2), col = lineage_colors[i], linewidth = 1, position = position_jitter(w = 0.2, h = 0.2), arrow = arrow(angle = 15, ends = "last", type = "closed"))
    }
    # Add labels
    p <- p + geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black")
    p <- p + ggrepel::geom_label_repel(data = centers, mapping = aes(label = Lineage), size = 4, color = "black") +
        # Make the legend dots larger
        guides(color = guide_legend(override.aes = list(size = 8)))

pdf(glue("{prefix}slingshot_pseudotime_trajectory_na_clusters.pdf"), width = 8, height = 8)
print(p)
dev.off()





# ---------------------------------------------------------------------
# UMAP sce2_na_clusters, facet by arm
# ---------------------------------------------------------------------
# reducedDim(sce2_na_clusters, "UMAP")
umap_data1 <- sce2_na_clusters@int_colData@listData$reducedDims$UMAP
# The DFrame contains a non-standard column that cannot be coerced
umap_data2 <- subset(colData(sce2_na_clusters), select = -slingshot)
umap_data <- cbind(as.data.frame(umap_data1), as_tibble(umap_data2)) %>%
    as_tibble() %>%
    dplyr::mutate(arm_groups = factor(case_when(
        treatment == "nilotinib_pdl1" ~ "nil_pdl1",
        TRUE  ~ "other_treatments")))


gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(umap_data$Cluster)))

centers <- umap_data %>%
    dplyr::select(Cluster, starts_with("UMAP_")) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize_all(median) %>%
    dplyr::mutate(Lineage = NA_character_)
# Add a lineage label to final cluster
for (i in seq_along(slingLineages(sce2_na_clusters))) {
    last_node_in_lineage <- tail(slingLineages(sce2_na_clusters)[[i]], n = 1)
    centers$Lineage[which(centers$Cluster == last_node_in_lineage)] <- names(slingLineages(sce2_na_clusters))[i]
}

lineage_colors <- gg_color_hue(length(seq_along(slingLineages(sce2_na_clusters))))
lineage_names <- names(slingLineages(sce2_na_clusters))

p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(color = "Cluster", title = glue("Slingshot\\nRoot cluster: {root_cluster}\\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}")) +
    facet_wrap(vars(arm_groups), nrow = 1) +
    theme_cowplot() +
    theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))
#     # Adding the curves
#     for (i in seq_along(slingCurves(sce2_na_clusters))) {
#         curve_i <- slingCurves(sce2_na_clusters)[[i]]
#         curve_i <- curve_i$s[curve_i$ord, ]
#         p <- p + geom_path(data = as.data.frame(curve_i), col = lineage_colors[i], size = 1)
#     }
    # Adding the lines
    for (i in seq_along(slingLineages(sce2_na_clusters))) {
        lineage_i <- slingLineages(sce2_na_clusters)[[i]]
        lineage_i2 <- centers %>%
            dplyr::filter(Cluster %in% lineage_i) %>%
            dplyr::arrange(match(Cluster, lineage_i))
        # See ?grid::arrow for options
        p <- p + geom_path(data = as.data.frame(lineage_i2), col = lineage_colors[i], size = 1, position = position_jitter(w = 0.2, h = 0.2), arrow = arrow(angle = 15, ends = "last", type = "closed"))
    }
    # Add labels
    p <- p + geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black")
    p <- p + ggrepel::geom_label_repel(data = centers, mapping = aes(label = Lineage), size = 4, color = "black") +
        # Make the legend dots larger
        guides(color = guide_legend(override.aes = list(size = 8)))

pdf(glue("{prefix}slingshot_pseudotime_trajectory_na_clusters_by_arm_groups.pdf"), width = 14, height = 8)
print(p)
dev.off()






# ---------------------------------------------------------------------
# Create slingshot dataset for topology
# ---------------------------------------------------------------------
sds1 <- SlingshotDataSet(sce2)
class(sds1)
sds1_na_clusters <- SlingshotDataSet(sce2_na_clusters)
class(sds1_na_clusters)
# [1] "SlingshotDataSet"
# attr(,"package")
# [1] "slingshot"

#######################################################################
# Differential Topology
#######################################################################

# Compare arms across the trajectory

# https://www.jwilber.me/permutationtest/

# Permutation test
condition <- colData(sce2)$treatment
condition <- forcats::fct_recode(condition, NP = "nilotinib_pdl1", OTHERS = "untreated_leuk", OTHERS = "nilotinib_il10r", OTHERS = "naive_hv1", OTHERS = "nilotinib")

# For each lineage
d1 <- vector()
p_values <- vector()
null_dist <- list()
for (j in seq_along(colnames(slingPseudotime(sds1, na = FALSE)))) {
    t1 <- slingPseudotime(sds1, na = FALSE)[, j]
    w1 <- slingCurveWeights(sds1)[, j]
    d1[j] <- weighted.mean(t1[condition == "OTHERS"], w1[condition == "OTHERS"]) - 
        weighted.mean(t1[condition == "NP"], w1[condition == "NP"])
     
    dist1 <- replicate(1e4, {
        condition.i <- sample(condition)
        weighted.mean(t1[condition.i == "OTHERS"], w1[condition.i == "OTHERS"]) - 
            weighted.mean(t1[condition.i == "NP"], w1[condition.i == "NP"])
    })
    
    null_dist[[j]] <- tibble(null_dist = dist1, lineage = colnames(slingPseudotime(sds1, na = FALSE))[j])
    
    p_values[j] <- paste0("Lineage ", j, ", p-value: ", mean(abs(dist1) > abs(d1[j])))
}

null_dist_tib <- bind_rows(null_dist) %>%
    dplyr::mutate(lineage = factor(lineage, labels = p_values))

d1_tib <- tibble(difference = d1, lineage = colnames(slingPseudotime(sds1, na = FALSE))) %>%
    dplyr::mutate(lineage = factor(lineage, labels = p_values))

p <- null_dist_tib %>%
    ggplot(aes(x = null_dist)) +
    geom_histogram() +
    facet_wrap(vars(lineage), labeller = label_value) +
    geom_vline(aes(xintercept = difference), data = d1_tib, colour = "red")


pdf(glue("{prefix}slingshot_diff_topology_perm_nilotinib_pdl1_vs_all_others.pdf"), width = 10, height = 10)
print(p)
dev.off()



saveRDS(sce2, glue("{prefix}sce2.rds"))





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
    # system(glue("sbatch {prefix}{curr_gemset_name}_{curr_norm_method}_{curr_res_name}.slurm"))
}





#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))





