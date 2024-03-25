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



# Set R working directory
proj_dir <- "/home/fifeb/shared/ris/knut0297/cd4_ins_pd1_scrnaseq_20230501"



# Set up output directory and filename prefix
out_dir <- glue("{proj_dir}/code_out/101_tetramer_de")
if (!dir.exists(glue("{out_dir}"))) {dir.create(glue("{out_dir}"), recursive = TRUE)}
setwd(glue("{out_dir}"))

prefix <- "101_"



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





# ---------------------------------------------------------------------
# Find DE genes between Tetramer pos cell groups
# ---------------------------------------------------------------------

# This file contains the tk_adt_de() function
source(glue("{proj_dir}/code/001_r_functions.R"))

curr_out <- list()

for (h in seq_along(dataset_list)) {
    curr_gemset_name <- dataset_list[[h]]$gemset_name
    curr_norm_method <- dataset_list[[h]]$norm_method
    curr_resolution <- dataset_list[[h]]$resolutions
    curr_res_name <- glue("SCT_snn_res_{curr_resolution}")


    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


    # ---------------------------------------------------------------------
    # Input data
    # ---------------------------------------------------------------------

    curr_seurat_object <- readRDS(glue("{proj_dir}/code_out/041_pca_umap/{curr_gemset_name}/{curr_norm_method}/041_seurat_object.rds"))
    curr_seurat_object@meta.data <- readRDS(glue("{proj_dir}/code_out/081_tetramer_and_cluster/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/081_cell_metadata.rds"))



    # ---------------------------------------------------------------------
    # Summarize the number of GEMs as TET pos/neg
    # ---------------------------------------------------------------------

    tet <- curr_seurat_object@meta.data %>%
        dplyr::select(barcode, starts_with("pos_")) %>%
        tidyr::pivot_longer(-barcode) %>%
        dplyr::filter(value == "pos") %>%
        dplyr::select(-value) %>%
        dplyr::group_by(barcode) %>%
        dplyr::summarize(data = list(name)) %>%
        dplyr::mutate(data_collapsed = sapply(data, function(x) paste0(sort(x), collapse = ",")))




    tet_tally <- tet %>%
        dplyr::group_by(data_collapsed) %>%
        tally()

    write_tsv(tet_tally, glue("{prefix}tetramer_positive_tally.txt"))

    p <- tet %>%
        ggplot(aes(x = data_collapsed)) +
        geom_bar() +
        geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
        ggupset::axis_combmatrix(sep = ",") +
        ggupset::theme_combmatrix()


    pdf(glue("{prefix}tetramer_positive_tally.pdf"))
    print(p)
    dev.off()



    # ---------------------------------------------------------------------
    # Test DE genes between treatments, within certain tetramer-pos cells
    # ---------------------------------------------------------------------

    tet_treat <- curr_seurat_object@meta.data %>%
        dplyr::select(barcode, treatment, starts_with("pos_")) %>%
        tidyr::pivot_longer(-c(barcode, treatment)) %>%
        dplyr::filter(value == "pos") %>%
        dplyr::select(-value) %>%
        dplyr::group_by(treatment, barcode) %>%
        dplyr::summarize(data = list(name)) %>%
        dplyr::mutate(data_collapsed = sapply(data, function(x) paste0(sort(x), collapse = ",")))


    tet_sets_raw <- tet_treat %>%
        dplyr::pull(data_collapsed) %>%
        unique()
    tet_sets <- tet_sets_raw %>%
        str_replace_all(",", "__")

    # Run only for datasets that have more than one treatment
    if (str_detect(curr_gemset_name, "all$")) {
        for (g in seq_along(tet_sets)) {
            curr_tet_set_raw <- tet_sets_raw[g]
            curr_tet_set <- tet_sets[g]

            if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}"), recursive = TRUE)}
            setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}"))


            comparisons <- list(
                pd1_vs_ctrl = list(
                    description = glue("pd1 vs control, in only {curr_tet_set} cells)"),
                    name_a = "pd1",
                    name_b = "control",
                    barcodes_a = (tet_treat %>%
                        dplyr::filter(treatment == "pd1" & data_collapsed == glue("{curr_tet_set_raw}")) %>%
                        pull(barcode)),
                    barcodes_b = (tet_treat %>%
                        dplyr::filter(treatment == "control" & data_collapsed == glue("{curr_tet_set_raw}")) %>%
                        pull(barcode))
                ),
                pdl1_vs_ctrl = list(
                    description = glue("pdl1 vs control, in only {curr_tet_set} cells)"),
                    name_a = "pdl1",
                    name_b = "control",
                    barcodes_a = (tet_treat %>%
                        dplyr::filter(treatment == "pdl1" & data_collapsed == glue("{curr_tet_set_raw}")) %>%
                        pull(barcode)),
                    barcodes_b = (tet_treat %>%
                        dplyr::filter(treatment == "control" & data_collapsed == glue("{curr_tet_set_raw}")) %>%
                        pull(barcode))
                )
            )




            for (m in seq_along(comparisons)) {
                curr_description <- comparisons[[m]]$description
                curr_name_a <- comparisons[[m]]$name_a
                curr_name_b <- comparisons[[m]]$name_b
                curr_barcodes_a <- comparisons[[m]]$barcodes_a
                curr_barcodes_b <- comparisons[[m]]$barcodes_b
                curr_comparison_name <- glue("{curr_tet_set}_{names(comparisons)[m]}")


                if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}"), recursive = TRUE)}
                setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}"))

                input_seurat_ob <- curr_seurat_object
                input_seurat_ob@meta.data <- input_seurat_ob@meta.data %>%
                    dplyr::mutate(comp_fct = case_when(
                        barcode %in% curr_barcodes_a ~ "group_a",
                        barcode %in% curr_barcodes_b ~ "group_b",
                        TRUE ~ NA_character_
                        )) %>%
                    dplyr::mutate(comp_fct = factor(comp_fct, levels = c("group_a", "group_b")))

                comp_sizes <- input_seurat_ob@meta.data %>%
                    dplyr::count(comp_fct) %>%
                    dplyr::filter(!is.na(comp_fct))

                # Run function
                # Make sure there are two groups and more than three cells/group (required by Seurat)
                if (nrow(comp_sizes) == 2 & all(comp_sizes$n >= 3)) {
                    curr_comparison_fct <- input_seurat_ob@meta.data$comp_fct
                    names(curr_comparison_fct) <- input_seurat_ob@meta.data$barcode
                    input_seurat_ob@active.ident <- curr_comparison_fct


                    # Run function
                    curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]] <- tk_adt_de(
                        seurat_object = input_seurat_ob,
                        comparison_name = curr_comparison_name,
                        comparison_fct = curr_comparison_fct,
                        lfc_threshold = 0.1,
                        description = curr_description
                    )
                    # Test within a cluster
                    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster"), recursive = TRUE)}
                    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster"))

                    input_seurat_ob <- curr_seurat_object

                    curr_cluster_list <-input_seurat_ob@meta.data$Cluster %>%
                        levels()

                    for (o in seq_along(curr_cluster_list)) {
                        curr_cluster <- curr_cluster_list[o]

                    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"), recursive = TRUE)}
                    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_tet_set}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"))

                        input_seurat_ob@meta.data <- input_seurat_ob@meta.data %>%
                            dplyr::mutate(comp_fct = case_when(
                                barcode %in% curr_barcodes_a & Cluster == curr_cluster ~ "group_a",
                                barcode %in% curr_barcodes_b & Cluster == curr_cluster ~ "group_b",
                                TRUE ~ NA_character_
                                )) %>%
                            dplyr::mutate(comp_fct = factor(comp_fct, levels = c("group_a", "group_b")))

                        comp_sizes <- input_seurat_ob@meta.data %>%
                            dplyr::count(comp_fct) %>%
                            dplyr::filter(!is.na(comp_fct))

                        # Run function
                        # Make sure there are two groups and more than three cells/group (required by Seurat)
                        if (nrow(comp_sizes) == 2 & all(comp_sizes$n >= 3)) {
                            curr_comparison_fct <- input_seurat_ob@meta.data$comp_fct
                            names(curr_comparison_fct) <- input_seurat_ob@meta.data$barcode
                            input_seurat_ob@active.ident <- curr_comparison_fct

                            curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_cluster"]][[curr_cluster]] <- tk_adt_de(
                                seurat_object = input_seurat_ob,
                                comparison_name = glue("Cluster_{curr_cluster}_{curr_comparison_name}"),
                                comparison_fct = curr_comparison_fct,
                                lfc_threshold = 0.1,
                                description = curr_description
                            )
                        } else {
                            curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_cluster"]][[curr_cluster]] <- "not_enough_cells_per_group"
                            write_tsv(tibble(not_enough_cells_per_group = ""), "not_enough_cells_per_group.txt")
                        }
                    }
                } else {
                    curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]] <- "not_enough_cells_per_group"
                    write_tsv(tibble(not_enough_cells_per_group = ""), "not_enough_cells_per_group.txt")
                }
            }
        } 
    }












    # ---------------------------------------------------------------------
    # get comparisons
    # ---------------------------------------------------------------------


    comparisons <- list(
        tet2.5HP_vs_tetInsBP8G = list(
            description = "2.5HP+, 6.9HPneg, P8Eneg, P8Gneg vs P8G+, 2.5HPneg, 6.9HPneg, P8Eneg",
            name_a = "tet2.5HP",
            name_b = "tetInsBP8G",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet2.5HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8G") %>%
                pull(barcode))
        ),
        tet2.5HP_vs_tetInsBP8E = list(
            description = "2.5HP+, 6.9HPneg, P8Eneg, P8Gneg vs P8E+, 2.5HPneg, 6.9HPneg, P8Gneg",
            name_a = "tet2.5HP",
            name_b = "tetInsBP8E",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet2.5HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E") %>%
                pull(barcode))
        ),
        tet2.5HP_vs_tetInsBP8G_or_tetInsBP8E = list(
            description = "2.5HP+, 6.9HPneg, P8Eneg, P8Gneg vs  P8G+, 2.5HPneg, 6.9HPneg, P8Eneg and P8E+, 2.5HPneg, 6.9HPneg, P8Gneg",
            name_a = "tet2.5HP",
            name_b = "tetInsBP8G_or_tetInsBP8E",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet2.5HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8G" | data_collapsed == "pos_tetInsBP8E") %>%
                pull(barcode))
        ),
        tetInsBP8E_vs_tetInsBP8G = list(
            description = "P8E+, 6.9HPneg, 2.5HPneg, P8Gneg vs P8G+, 6.9HPneg, 2.5HPneg, P8Eneg",
            name_a = "tetInsBP8E",
            name_b = "tetInsBP8G",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8G") %>%
                pull(barcode))
        ),
        tet2.5HP_vs_tet6.9HP = list(
            description = "2.5HP+, 6.9HPneg, P8Eneg, P8Gneg vs 6.9HP+, 2.5HPneg, P8Eneg, P8Gneg",
            name_a = "tet2.5HP",
            name_b = "tet6.9HP",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet2.5HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet6.9HP") %>%
                pull(barcode))
        ),
        tet6.9HP_vs_tetInsBP8E = list(
            description = "6.9HP+, 2.5HPneg, P8Eneg, P8Gneg vs P8E+, 2.5HPneg, 6.9HPneg, P8Gneg",
            name_a = "tet6.9HP",
            name_b = "tetInsBP8E",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet6.9HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E") %>%
                pull(barcode))
        ),
        tet6.9HP_vs_tetInsBP8G = list(
            description = "6.9HP+, 2.5HPneg, P8Eneg, P8Gneg vs P8G+, 2.5HPneg, 6.9HPneg, P8Eneg",
            name_a = "tet6.9HP",
            name_b = "tetInsBP8G",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tet6.9HP") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8G") %>%
                pull(barcode))
        ),
        tetInsBP8G_vs_tetInsBP8E_and_tetInsBP8G = list(
            description = "P8G+, 6.9HPneg, 2.5HPneg, P8Gneg vs P8G+, 6.9HPneg, 2.5HPneg, P8E+",
            name_a = "tetInsBP8G",
            name_b = "tetInsBP8E_and_tetInsBP8G",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8G") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E,pos_tetInsBP8G") %>%
                pull(barcode))
        ),
        tetInsBP8E_vs_tetInsBP8E_and_tetInsBP8G = list(
            description = "P8E+, 6.9HPneg, 2.5HPneg, P8Gneg vs P8G+, 6.9HPneg, 2.5HPneg, P8E+",
            name_a = "tetInsBP8E",
            name_b = "tetInsBP8E_and_tetInsBP8G",
            barcodes_a = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E") %>%
                pull(barcode)),
            barcodes_b = (tet %>%
                dplyr::filter(data_collapsed == "pos_tetInsBP8E,pos_tetInsBP8G") %>%
                pull(barcode))
        )
    )





    for (m in seq_along(comparisons)) {
        curr_description <- comparisons[[m]]$description
        curr_name_a <- comparisons[[m]]$name_a
        curr_name_b <- comparisons[[m]]$name_b
        curr_barcodes_a <- comparisons[[m]]$barcodes_a
        curr_barcodes_b <- comparisons[[m]]$barcodes_b
        curr_comparison_name <- names(comparisons)[m]


        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}"))

        input_seurat_ob <- curr_seurat_object
        input_seurat_ob@meta.data <- input_seurat_ob@meta.data %>%
            dplyr::mutate(comp_fct = case_when(
                barcode %in% curr_barcodes_a ~ "group_a",
                barcode %in% curr_barcodes_b ~ "group_b",
                TRUE ~ NA_character_
                )) %>%
            dplyr::mutate(comp_fct = factor(comp_fct, levels = c("group_a", "group_b")))

        curr_comparison_fct <- input_seurat_ob@meta.data$comp_fct
        names(curr_comparison_fct) <- input_seurat_ob@meta.data$barcode
        input_seurat_ob@active.ident <- curr_comparison_fct


        # Run function
        curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]] <- tk_adt_de(
            seurat_object = input_seurat_ob,
            comparison_name = curr_comparison_name,
            comparison_fct = curr_comparison_fct,
            lfc_threshold = 0.1,
            description = curr_description
        )
        # Test within a cluster
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster"))

        input_seurat_ob <- curr_seurat_object

        curr_cluster_list <-input_seurat_ob@meta.data$Cluster %>%
            levels()

        for (o in seq_along(curr_cluster_list)) {
            curr_cluster <- curr_cluster_list[o]

            if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"), recursive = TRUE)}
            setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_cluster/cluster_{curr_cluster}"))

            input_seurat_ob@meta.data <- input_seurat_ob@meta.data %>%
                dplyr::mutate(comp_fct = case_when(
                    barcode %in% curr_barcodes_a & Cluster == curr_cluster ~ "group_a",
                    barcode %in% curr_barcodes_b & Cluster == curr_cluster ~ "group_b",
                    TRUE ~ NA_character_
                    )) %>%
                dplyr::mutate(comp_fct = factor(comp_fct, levels = c("group_a", "group_b")))

            comp_sizes <- input_seurat_ob@meta.data %>%
                dplyr::count(comp_fct) %>%
                dplyr::filter(!is.na(comp_fct))

            # Run function
            # Make sure there are two groups and more than three cells/group (required by Seurat)
            if (nrow(comp_sizes) == 2 & all(comp_sizes$n >= 3)) {
                curr_comparison_fct <- input_seurat_ob@meta.data$comp_fct
                names(curr_comparison_fct) <- input_seurat_ob@meta.data$barcode
                input_seurat_ob@active.ident <- curr_comparison_fct

                curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_cluster"]][[curr_cluster]] <- tk_adt_de(
                    seurat_object = input_seurat_ob,
                    comparison_name = glue("Cluster_{curr_cluster}_{curr_comparison_name}"),
                    comparison_fct = curr_comparison_fct,
                    lfc_threshold = 0.1,
                    description = curr_description
                )
            } else {
                curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_cluster"]][[curr_cluster]] <- "not_enough_cells_per_group"
                write_tsv(tibble(not_enough_cells_per_group = ""), "not_enough_cells_per_group.txt")
            }
        }
        # Test within a treatment
        if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment"), recursive = TRUE)}
        setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment"))

        input_seurat_ob <- curr_seurat_object

        curr_cluster_list <-input_seurat_ob@meta.data$treatment %>%
            factor() %>%
            levels()

        for (o in seq_along(curr_cluster_list)) {
            curr_cluster <- curr_cluster_list[o]

            if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment/treatment_{curr_cluster}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment/treatment_{curr_cluster}"), recursive = TRUE)}
            setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/{curr_comparison_name}/by_treatment/treatment_{curr_cluster}"))

            input_seurat_ob@meta.data <- input_seurat_ob@meta.data %>%
                dplyr::mutate(comp_fct = case_when(
                    barcode %in% curr_barcodes_a & treatment == curr_cluster ~ "group_a",
                    barcode %in% curr_barcodes_b & treatment == curr_cluster ~ "group_b",
                    TRUE ~ NA_character_
                    )) %>%
                dplyr::mutate(comp_fct = factor(comp_fct, levels = c("group_a", "group_b")))

            comp_sizes <- input_seurat_ob@meta.data %>%
                dplyr::count(comp_fct) %>%
                dplyr::filter(!is.na(comp_fct))

            # Run function
            # Make sure there are two groups and more than three cells/group (required by Seurat)
            if (nrow(comp_sizes) == 2 & all(comp_sizes$n >= 3)) {
                curr_comparison_fct <- input_seurat_ob@meta.data$comp_fct
                names(curr_comparison_fct) <- input_seurat_ob@meta.data$barcode
                input_seurat_ob@active.ident <- curr_comparison_fct

                curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_treatment"]][[curr_cluster]] <- tk_adt_de(
                    seurat_object = input_seurat_ob,
                    comparison_name = glue("treatment_{curr_cluster}_{curr_comparison_name}"),
                    comparison_fct = curr_comparison_fct,
                    lfc_threshold = 0.1,
                    description = curr_description
                )
            } else {
                curr_out[[curr_gemset_name]][[curr_norm_method]][[curr_res_name]][[curr_comparison_name]][["by_treatment"]][[curr_cluster]] <- "not_enough_cells_per_group"
                write_tsv(tibble(not_enough_cells_per_group = ""), "not_enough_cells_per_group.txt")
            }
        }
    }
}

if (!dir.exists(glue("{out_dir}"))) {dir.create(glue("{out_dir}"), recursive = TRUE)}
setwd(glue("{out_dir}"))

saveRDS(comparisons, glue("{prefix}comparisons_list.rds"))
saveRDS(curr_out, glue("{prefix}tet_de_list.rds"))



#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))










