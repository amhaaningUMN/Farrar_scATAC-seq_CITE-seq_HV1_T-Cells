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


proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "081_"
out <- glue("{prefix}adt_and_cluster")
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
#SBATCH --time=0:30:00
#SBATCH --mem=16gb
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
library(gridExtra)
library(cowplot)
library(Seurat)
library(NMF)
library(RColorBrewer)
library(openxlsx)
library(parallel)
library(gtools)
library(viridis)
library(biomaRt)
library(patchwork)





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


# ---------------------------------------------------------------------
# Capture continuous ADT expression levels from normalized counts
# ---------------------------------------------------------------------
# Get the normalized data
# Convert zero ADT counts to NAs
norm <- curr_seurat_object@assays$ADT@data %>%
    as.data.frame() %>%
    rownames_to_column(var = "adt") %>%
    as_tibble() %>%
    tidyr::pivot_longer(!adt, names_to = "barcode", values_to = "counts") %>%
    dplyr::mutate(counts_NA = case_when(
        counts == 0 ~ NA_real_,
        TRUE ~ counts))





# ---------------------------------------------------------------------
# Plot normalized expression to determine pos/neg GEMs
# ---------------------------------------------------------------------

# This is the normalized expression threshold, where cells with norm ADT
# levels higher than this, are considered "positive". This section was run interactively
# on all datasets and then decided.
adt_dsb_norm_pos_thresh <- 7


adt_features <- rownames(curr_seurat_object@assays$ADT@data)


for (i in seq_along(adt_features)) {
    curr_adt <- adt_features[i]

    p <- norm %>%
        dplyr::filter(adt == curr_adt) %>%
        dplyr::mutate(adt_pos = if_else(counts_NA > adt_dsb_norm_pos_thresh, "pos", "neg")) %>%
        ggplot(aes(x = counts_NA, fill = adt_pos)) +
        geom_histogram() +
        geom_vline(xintercept = adt_dsb_norm_pos_thresh) +
        annotate("text", x = adt_dsb_norm_pos_thresh, y = Inf, hjust = 0, vjust = 1, label = glue("ADT pos/neg threshold: {adt_dsb_norm_pos_thresh}")) +
        labs(title = glue("{curr_adt} normalized ADT expression"),
            x = "Normalized (dsb method) ADT expression",
            y = "Number of GEMs/cells per bin") +
        theme_minimal()

    pdf(glue("{prefix}{curr_adt}_histogram_pos.pdf"), height = 6, width = 6)
    print(p)
    dev.off()
}






# ---------------------------------------------------------------------
# Create ordinal levels of normalized ADT expression
# ---------------------------------------------------------------------

# NOTE: Only "ADT-positive" GEMs will be classified as low, mod, or high. Thus, the distribution
# will be based on the level of expression within positive cells only (ignoring background staining).

temp_func <- function(x) {
    factor(replace_na(as.character(x), "neg"), levels = c("neg", "low", "moderate", "high"), ordered = TRUE)
}

norm1_long <- norm %>%
    # Discretize ADT levels into low, med, high (based on Q1 and Q3)
    # The group_by is needed here to get the correct fivenum for the group
    group_by(adt) %>%
    dplyr::mutate(counts_NA_pos = if_else(counts_NA <= adt_dsb_norm_pos_thresh, NA_real_, counts_NA)) %>%
    # If the fivenum() errors out, set value to NA
    dplyr::mutate(norm_ord = tryCatch(
        cut(counts_NA_pos,
            breaks = c(-Inf, 0, fivenum(counts_NA_pos)[2], fivenum(counts_NA_pos)[4], Inf),
            labels = c(NA_character_, "low", "moderate", "high")),
        error = function(e) NA_character_)) %>%
    ungroup()

norm1 <- norm1_long %>%
    # Make wide dataset
    dplyr::select(adt, barcode, counts_NA, norm_ord) %>%
    tidyr::pivot_wider(names_from = "adt", values_from = c("counts_NA", "norm_ord")) %>%
    dplyr::rename_at(vars(starts_with("counts_NA")), function(x) {str_remove(x, "^counts_NA_")}) %>%
    dplyr::rename_at(vars(starts_with("norm_")), function(x) {str_remove(x, "^norm_")}) %>%
    purrr::modify_at(vars(starts_with("ord_")), temp_func)


for (i in seq_along(adt_features)) {
    curr_adt <- adt_features[i]

    d <- norm1_long %>%
        dplyr::filter(adt == curr_adt)

    if (!all(is.na(d$norm_ord))) {
        p <- d %>%
            ggplot(aes(x = counts_NA, fill = norm_ord)) +
            geom_histogram() +
            geom_vline(xintercept = adt_dsb_norm_pos_thresh) +
            annotate("text", x = adt_dsb_norm_pos_thresh, y = Inf, hjust = 0, vjust = 1, label = glue("ADT pos/neg threshold: {adt_dsb_norm_pos_thresh}")) +
            labs(title = glue("{curr_adt} normalized ADT expression"),
                x = "Normalized (dsb method) ADT expression",
                y = "Number of GEMs/cells per bin") +
            theme_minimal()
    } else {
        # the norm_ord column is all NAs
        warning(paste0(curr_adt, " has no ordinal levels, creating blank plot"))
        # Create an empty plot!
        p <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_void() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{curr_adt}"))
    }
    pdf(glue("{prefix}{curr_adt}_histogram_ord.pdf"), height = 6, width = 6)
    print(p)
    dev.off()
}




# ---------------------------------------------------------------------
# Add ordinal ADT levels and binary (pos/neg) ADT levels to metadata
# ---------------------------------------------------------------------




tk_adt_pos <- function(vect, thresh = NULL) {
   ifelse(vect > thresh, "pos", "neg") 
}

curr_seurat_object@meta.data <- curr_seurat_object@meta.data %>%
    dplyr::left_join(norm1, by = "barcode") %>%
    dplyr::mutate(across(starts_with("adt"), list(pos = ~tk_adt_pos(., thresh = adt_dsb_norm_pos_thresh)), .names = "pos_{.col}"))


# ---------------------------------------------------------------------
# Make ADT calls
# ---------------------------------------------------------------------

adt_pos_cells <- curr_seurat_object@meta.data %>%
    dplyr::select(barcode, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_pos_status") %>%
    dplyr::filter(adt_pos_status == "pos") %>%
    dplyr::select(-adt_pos_status) %>%
    dplyr::group_by(barcode) %>%
    tidyr::nest(data = adt_features) %>%
    dplyr::mutate(adt_features_combined = purrr::map_chr(data, ~ paste(sort(.$adt_features), collapse = ","))) %>%
    dplyr::select(barcode, adt_features_combined)


curr_seurat_object@meta.data <- curr_seurat_object@meta.data %>%
    dplyr::left_join(adt_pos_cells, by = "barcode") %>%
    tidyr::replace_na(list(adt_features_combined = "neg")) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$barcode)



# --------------------------------------------------------------------
# Save metadata
# ---------------------------------------------------------------------

saveRDS(curr_seurat_object@meta.data, glue("{prefix}cell_metadata.rds"))
write_tsv(curr_seurat_object@meta.data, glue("{prefix}cell_metadata.txt"))
openxlsx::write.xlsx(curr_seurat_object@meta.data, glue("{prefix}cell_metadata.xlsx"))

DefaultAssay(curr_seurat_object) <- "SCT"
Idents(curr_seurat_object) <- col_cluster






# ---------------------------------------------------------------------
# Capture UMAP data
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






# ---------------------------------------------------------------------
# Visualize protein levels
# ---------------------------------------------------------------------


if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_adt"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_adt"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_adt"))

for (i in seq_along(adt_features)) {
    curr_adt <- adt_features[i]

    p <- umap_data %>%
        dplyr::arrange(!!rlang::sym(curr_adt)) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(curr_adt))) +
            geom_point(size = 0.5) +
            scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::viridis(100)))) +
            labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels: {curr_adt}"), subtitle = "Cluster numbers inside plot") +
            geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
            theme_cowplot() +
            theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                strip.background = element_blank())
    pdf(glue("{prefix}umap_by_{curr_adt}.pdf"), width = 6, height = 6)
    print(p)
    dev.off()
}



setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))

p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, matches("^adt[^_]")) %>%
    tidyr::pivot_longer(matches("^adt[^_]"), names_to = "adt_features", values_to = "adt_continuous_expr") %>%
    dplyr::mutate(adt_features = str_remove(adt_features, "continuous_")) %>%
    dplyr::arrange(adt_continuous_expr) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_continuous_expr)) +
        geom_point(size = 0.5) +
        scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::viridis(100)))) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_adt.pdf"), width = 7, height = 7)
print(p)
dev.off()




# Paired UMAP with GEX and ADT
gex_features <- c("Cd4", "Entpd1", "Itga2", "Pdcd1", "Havcr2", "Lag3")
adt_features_gex_features <- c(adt_features, gex_features)

DefaultAssay(curr_seurat_object) <- "SCT"
expr_data <- Seurat::FetchData(object = curr_seurat_object, vars = gex_features, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    as_tibble()

umap_data_with_gex <- umap_data %>%
    dplyr::left_join(expr_data, by = "barcode")

p1 <- umap_data_with_gex %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, any_of(adt_features)) %>%
    tidyr::pivot_longer(any_of(adt_features), names_to = "adt_features", values_to = "adt_continuous_expr") %>%
    dplyr::mutate(adt_features = forcats::fct_relevel(adt_features, !!adt_features)) %>%
    dplyr::arrange(adt_continuous_expr) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_continuous_expr)) +
        geom_point(size = 0.5) +
        scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::viridis(100)))) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features), nrow = 1) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

p2 <- umap_data_with_gex %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, any_of(gex_features)) %>%
    tidyr::pivot_longer(any_of(gex_features), names_to = "gex_features", values_to = "gex_continuous_expr") %>%
    dplyr::mutate(gex_features = forcats::fct_relevel(gex_features, !!gex_features)) %>%
    dplyr::arrange(gex_continuous_expr) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = gex_continuous_expr)) +
        geom_point(size = 0.5) +
        scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous GEX levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(gex_features), nrow = 1) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

patchwork <- patchwork::wrap_plots(list(p1, p2), nrow = 2) +
    patchwork::plot_layout(guides = "collect") &
    patchwork::plot_annotation(title = glue("Normalized ADT and GEX levels, compared")) &
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

pdf(glue("{prefix}umap_adt_gex.pdf"), width = 15, height = 11)
print(patchwork)
dev.off()












if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/biaxial_adt"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/biaxial_adt"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/biaxial_adt"))


combos <- umap_data %>%
    dplyr::select(barcode, matches("^adt[^_]")) %>%
    tidyr::pivot_longer(matches("^adt[^_]"), names_to = "adt_features", values_to = "adt_continuous_expr") %>%
    tidyr::expand(adt_features, adt_features) %>%
    set_names("x", "y")

for (i in seq_len(nrow(combos))) {
    curr_x <- combos %>%
        dplyr::slice(i) %>%
        dplyr::pull(x)
    curr_y <- combos %>%
        dplyr::slice(i) %>%
        dplyr::pull(y)

    p <- umap_data %>%
        ggplot(aes(x = !!rlang::sym(curr_x), y = !!rlang::sym(curr_y), color = adt_features_combined)) +
        geom_point(size = 0.5) +
        # scale_color_manual(name = "ADT", values = c()) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels")) +
        theme_cowplot() +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 4))
    pdf(glue("{prefix}biaxial_adt_{curr_x}_vs_{curr_y}.pdf"), width = 20, height = 20)
    print(p)
    dev.off()
}



setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


# ---------------------------------------------------------------------
# ADT positive cells only
# ---------------------------------------------------------------------

p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_binary_expr") %>%
    dplyr::mutate(adt_binary_expr = factor(adt_binary_expr, levels = c("pos", "neg"))) %>%
    dplyr::arrange(desc(adt_binary_expr)) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_binary_expr)) +
        geom_point(size = 0.5, alpha = 1, stroke = 0) +
        scale_color_manual(name = "Exprs", values = c("blue", "gray95")) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nBinary ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        theme_cowplot() +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_adt_binary_blue.pdf"), width = 7, height = 7)
print(p)
dev.off()

p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_binary_expr") %>%
    dplyr::mutate(adt_binary_expr = factor(adt_binary_expr, levels = c("pos", "neg"))) %>%
    dplyr::arrange(desc(adt_binary_expr)) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_binary_expr)) +
        geom_point(size = 0.5, alpha = 1, stroke = 0) +
        scale_color_manual(name = "Exprs", values = c("red", "gray95")) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nBinary ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        theme_cowplot() +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_adt_binary_red.pdf"), width = 7, height = 7)
print(p)
dev.off()

p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_binary_expr") %>%
    dplyr::mutate(adt_binary_expr = factor(adt_binary_expr, levels = c("pos", "neg"))) %>%
    dplyr::arrange(desc(adt_binary_expr)) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_binary_expr)) +
        geom_point(size = 0.5, alpha = 1, stroke = 0) +
        scale_color_manual(name = "Exprs", values = c("black", "gray95")) +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nBinary ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        theme_cowplot() +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_adt_binary_black.pdf"), width = 7, height = 7)
print(p)
dev.off()


p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_binary_expr") %>%
    dplyr::mutate(adt_binary_expr_pos = if_else(adt_binary_expr == "pos", adt_features, NA_character_)) %>%
    dplyr::mutate(adt_binary_expr_pos = factor(adt_binary_expr_pos, levels = c("pos_adt2.5HP","pos_adt6.9HP", "pos_adtInsBP8E", "pos_adtInsBP8G"))) %>%
    # arrange always sorts NAs to the end, regardless of desc()
    dplyr::arrange(!is.na(adt_binary_expr_pos), adt_binary_expr_pos) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_binary_expr_pos)) +
        geom_point(size = 0.5, alpha = 1, stroke = 0) +
        scale_color_manual(name = "Exprs", values = c("firebrick1", "magenta", "blue", "forestgreen"), na.value = "grey95") +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nBinary ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

pdf(glue("{prefix}umap_adt_binary_color.pdf"), width = 7, height = 7)
print(p)
dev.off()


p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, starts_with("pos_adt")) %>%
    tidyr::pivot_longer(starts_with("pos_adt"), names_to = "adt_features", values_to = "adt_binary_expr") %>%
    dplyr::mutate(adt_binary_expr_pos = if_else(adt_binary_expr == "pos", adt_features, NA_character_)) %>%
    dplyr::mutate(adt_binary_expr_pos = factor(adt_binary_expr_pos, levels = c("pos_adt2.5HP","pos_adt6.9HP", "pos_adtInsBP8E", "pos_adtInsBP8G"))) %>%
    # arrange always sorts NAs to the end, regardless of desc()
    dplyr::arrange(!is.na(adt_binary_expr_pos), adt_binary_expr_pos) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_binary_expr_pos)) +
        geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
        scale_color_manual(name = "Exprs", values = c("firebrick1", "magenta", "blue", "forestgreen"), na.value = "grey95") +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nBinary ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        # facet_wrap(facets = vars(adt_features)) +
        guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())

pdf(glue("{prefix}umap_adt_binary_color_together.pdf"), width = 7, height = 7)
print(p)
dev.off()

p <- umap_data %>%
    dplyr::select(barcode, UMAP_1, UMAP_2, matches("^adt[^_]")) %>%
    tidyr::pivot_longer(matches("^adt[^_]"), names_to = "adt_features", values_to = "adt_continuous_expr") %>%
    dplyr::mutate(adt_continuous_expr = if_else(adt_continuous_expr > adt_dsb_norm_pos_thresh, adt_continuous_expr, NA_real_)) %>%
    # arrange always sorts NAs to the end, regardless of desc()
    dplyr::arrange(!is.na(adt_continuous_expr), adt_continuous_expr) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = adt_continuous_expr)) +
        geom_point(size = 1, stroke = 0) +
        scale_color_gradientn(name = "Exprs", colors = rev(viridis::viridis(100)), na.value = "grey95") +
        labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels"), subtitle = "Cluster numbers inside plot") +
        geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
        facet_wrap(facets = vars(adt_features)) +
        theme_cowplot() +
        theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_blank())
pdf(glue("{prefix}umap_adt_highlight_pos.pdf"), width = 7, height = 7)
print(p)
dev.off()




# ---------------------------------------------------------------------
# UMAP, ADT expression by sample_name
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name"))



for (i in seq_along(adt_features)) {
    curr_adt <- adt_features[i]

    p <- umap_data %>%
        dplyr::arrange(!!rlang::sym(curr_adt)) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(curr_adt))) +
            geom_point(size = 0.5) +
            scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::viridis(100)))) +
            labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels: {curr_adt}"), subtitle = "Cluster numbers inside plot") +
            geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
            facet_wrap(facets = vars(sample_name)) +
            theme_cowplot() +
            theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                strip.background = element_blank())
    pdf(glue("{prefix}umap_by_{curr_adt}_sample_name.pdf"), width = 12, height = 12)
    print(p)
    dev.off()
}








# ---------------------------------------------------------------------
# UMAP cont ADT, highlight by sample_name
# ---------------------------------------------------------------------


for (h in seq_along(adt_features)) {
    curr_adt <- adt_features[h]

    for (i in seq_len(length(levels(umap_data$Cluster)))) {
        curr_cluster <- levels(umap_data$Cluster)[i]

        # Individual plots
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name_cluster/{curr_adt}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name_cluster/{curr_adt}"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_sample_name_cluster/{curr_adt}"))

        p <- umap_data %>%
            dplyr::mutate(curr_adt_string = !!rlang::sym(curr_adt)) %>%
            dplyr::mutate(curr_adt_string = case_when(
                Cluster != curr_cluster ~ NA_real_,
                TRUE ~ curr_adt_string)) %>%
            dplyr::mutate(curr_adt_string = if_else(curr_adt_string == 0, NA_real_, curr_adt_string)) %>%
            # arrange always sorts NAs to the end, regardless of desc()
            dplyr::arrange(!is.na(curr_adt_string), curr_adt_string) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(x = UMAP_1, y = UMAP_2, color = curr_adt_string), size = 0.5) +
                scale_color_gradientn(name = "Exprs", colors = c(rev(viridis::viridis(100))), na.value = "gray95") +
                labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels: {curr_adt}\nCluster: {curr_cluster}"), subtitle = "Cluster numbers inside plot") +
                theme_cowplot() +
                facet_wrap(facets = vars(sample_name)) +
                theme_cowplot() +
                theme(aspect.ratio = 1,
                    plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5),
                    strip.background = element_blank())

        pdf(glue("{prefix}umap_highlight_by_{curr_adt}_sample_name_cluster_{curr_cluster}.pdf"), width = 12, height = 12)
        print(p)
        dev.off()
    }
}




# ---------------------------------------------------------------------
# UMAP, ADT expression by treatment
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment"))



for (i in seq_along(adt_features)) {
    curr_adt <- adt_features[i]

    p <- umap_data %>%
        dplyr::arrange(!!rlang::sym(curr_adt)) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(curr_adt))) +
            geom_point(size = 0.5) +
            scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::viridis(100)))) +
            labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels: {curr_adt}"), subtitle = "Cluster numbers inside plot") +
            geom_text(data = centers, mapping = aes(label = Cluster), size = 4, color = "black") +
            facet_wrap(facets = vars(treatment)) +
            theme_cowplot() +
            theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                strip.background = element_blank())
    pdf(glue("{prefix}umap_by_{curr_adt}_treatment.pdf"), width = 12, height = 12)
    print(p)
    dev.off()
}



# ---------------------------------------------------------------------
# UMAP cont ADT, highlight by treatment
# ---------------------------------------------------------------------


for (h in seq_along(adt_features)) {
    curr_adt <- adt_features[h]

    for (i in seq_len(length(levels(umap_data$Cluster)))) {
        curr_cluster <- levels(umap_data$Cluster)[i]

        # Individual plots
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment_cluster/{curr_adt}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment_cluster/{curr_adt}"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/umap_by_adt_treatment_cluster/{curr_adt}"))

        p <- umap_data %>%
            dplyr::mutate(curr_adt_string = !!rlang::sym(curr_adt)) %>%
            dplyr::mutate(curr_adt_string = case_when(
                Cluster != curr_cluster ~ NA_real_,
                TRUE ~ curr_adt_string)) %>%
            dplyr::mutate(curr_adt_string = if_else(curr_adt_string == 0, NA_real_, curr_adt_string)) %>%
            # arrange always sorts NAs to the end, regardless of desc()
            dplyr::arrange(!is.na(curr_adt_string), curr_adt_string) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(x = UMAP_1, y = UMAP_2, color = curr_adt_string), size = 0.5) +
                scale_color_gradientn(name = "Exprs", colors = c(rev(viridis::viridis(100))), na.value = "gray95") +
                labs(title = glue("{curr_gemset_name}\n{curr_res_name}\nContinuous ADT levels: {curr_adt}\nCluster: {curr_cluster}"), subtitle = "Cluster numbers inside plot") +
                theme_cowplot() +
                facet_wrap(facets = vars(treatment)) +
                theme_cowplot() +
                theme(aspect.ratio = 1,
                    plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5),
                    strip.background = element_blank())

        pdf(glue("{prefix}umap_highlight_by_{curr_adt}_treatment_cluster_{curr_cluster}.pdf"), width = 12, height = 12)
        print(p)
        dev.off()
    }
}


# ---------------------------------------------------------------------
# Export all data
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))



# Get the normalized data
# Convert zero ADT counts to NAs
norm <- curr_seurat_object@assays$ADT@data %>%
    as.data.frame() %>%
    rownames_to_column(var = "adt") %>%
    as_tibble() %>%
    tidyr::pivot_longer(!adt, names_to = "barcode", values_to = "counts") %>%
    dplyr::mutate(counts_NA = case_when(
        counts == 0 ~ NA_real_,
        TRUE ~ counts))





norm_with_meta <- norm %>%
    dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
    dplyr::select(adt, barcode, counts, counts_NA, orig.ident, sample_name, Cluster)

write_tsv(norm_with_meta, glue("{prefix}norm_adt_counts.txt"))

out_filename <- glue("{prefix}norm_adt_counts.xlsx")
wb <- openxlsx::write.xlsx(list(norm_adt_counts = norm_with_meta), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)


write_tsv(umap_data, glue("{prefix}umap_data.txt"))
saveRDS(umap_data, glue("{prefix}umap_data.rds"))








#######################################################################
# ADT pos-GEMs counts per cluster
#######################################################################

if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/counts"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/counts"), recursive = TRUE)}
setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/counts"))


p <- norm %>%
    dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
    ggplot(aes(x = adt, y = counts, fill = Cluster)) +
        geom_boxplot() +
        theme_bw() +
        labs(fill = "Cluster", title = glue("Distribution of normalized ADT counts by GEX clusters:\n{col_cluster}\n(each dot is one GEM)"),
            y = "Centered log ratio (CLR) normalized counts", x = "ADT Marker")

pdf(glue("{prefix}boxplot_norm_adt_expr_by_cluster.pdf"), width = 15, height = 5)
print(p)
dev.off()



gg_color_hue <- function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
cluster_colors <- gg_color_hue(length(levels(umap_data$Cluster)))


# Count the number of GEMs per group
for (i in seq_along(adt_features)) {
    curr_adt_feature <- adt_features[i]
    curr_colname <- glue("pos_{curr_adt_feature}")

    pos_neg_colors <- setNames(gg_color_hue(2), c("pos", "neg"))

    # ---------------------------------------------------------------------
    # Prop
    # ---------------------------------------------------------------------
    norm_curr <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum)

    filename_prefix <- glue("{prefix}table_{curr_adt_feature}_adt_by_sample_name")
    write_tsv(norm_curr, glue("{filename_prefix}.txt"))
    wb <- openxlsx::write.xlsx(list(norm_curr), file = glue("{filename_prefix}.xlsx"), rowNames = FALSE)
    openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
    openxlsx::saveWorkbook(wb, glue("{filename_prefix}.xlsx"), overwrite = TRUE)

    p <- norm_curr %>%
        ggplot(aes(x = Cluster, y = sample_name_prop, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters (all bars within a facet sum to 100%)\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Proportion of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_prop_by_sample_name.pdf"), width = 15, height = 5)
    print(p)
    dev.off()

    p <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        dplyr::filter(!!sym(curr_colname) == "pos") %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum) %>%
        ggplot(aes(x = Cluster, y = sample_name_prop, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters (all bars within a facet sum to 100%)\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Proportion of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_prop_by_sample_name_pos.pdf"), width = 15, height = 5)
    print(p)
    dev.off()

    p <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        dplyr::filter(!!sym(curr_colname) == "neg") %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum) %>%
        ggplot(aes(x = Cluster, y = sample_name_prop, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters (all bars within a facet sum to 100%)\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Proportion of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_prop_by_sample_name_neg.pdf"), width = 15, height = 5)
    print(p)
    dev.off()

    # ---------------------------------------------------------------------
    # Number
    # ---------------------------------------------------------------------
    p <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum) %>%
        ggplot(aes(x = Cluster, y = n, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Number of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_number_by_sample_name.pdf"), width = 15, height = 5)
    print(p)
    dev.off()

    p <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        dplyr::filter(!!sym(curr_colname) == "pos") %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum) %>%
        ggplot(aes(x = Cluster, y = n, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Number of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_number_by_sample_name_pos.pdf"), width = 15, height = 5)
    print(p)
    dev.off()

    p <- norm %>%
        dplyr::left_join(curr_seurat_object@meta.data, by = "barcode") %>%
        group_by(sample_name, Cluster, !!sym(curr_colname), .drop = FALSE) %>%
        tally() %>%
        dplyr::filter(!!sym(curr_colname) == "neg") %>%
        ungroup() %>%
        group_by(sample_name) %>%
        dplyr::mutate(sample_name_sum = sum(n)) %>%
        dplyr::mutate(sample_name_prop = n / sample_name_sum) %>%
        ggplot(aes(x = Cluster, y = n, fill = !!sym(curr_colname))) +
            geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
            facet_wrap(vars(sample_name), scales = "free_y", nrow = 1) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
            scale_fill_manual(values = pos_neg_colors) +
            coord_flip() +
            theme_bw() +
            labs(fill = "adt-pos/neg", title = glue("Distribution of {curr_adt_feature} GEMs by GEX clusters\n{curr_gemset_name}, {curr_norm_method}, {curr_res_name}"),
                y = "Number of GEMs in each cluster", x = "Cluster")

    pdf(glue("{prefix}barplot_{curr_adt_feature}_adt_number_by_sample_name_neg.pdf"), width = 15, height = 5)
    print(p)
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
for (m  in seq_along(long_string_r_vect)) {
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
