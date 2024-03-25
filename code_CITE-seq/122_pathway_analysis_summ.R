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
prefix <- "122_"
out <- glue("{prefix}pathway_analysis_summ")
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
    list(gemset_name = "adtCD4neg", norm_method = "cc_none_gene_none", resolutions = "0.5"),
    list(gemset_name = "adtCD4pos", norm_method = "cc_none_gene_none", resolutions = "0.5"),
    list(gemset_name = "all", norm_method = "cc_none_gene_none", resolutions = "0.5")
    )


# ---------------------------------------------------------------------
# Loop through each dataset
# ---------------------------------------------------------------------

for (h in seq_along(dataset_list)) {
    curr_gemset_name <- dataset_list[[h]]$gemset_name
    curr_norm_method <- dataset_list[[h]]$norm_method
    curr_resolution <- dataset_list[[h]]$resolutions
    curr_res_name <- glue("SCT_snn_res_{curr_resolution}")


    if (!dir.exists(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))) {dir.create(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}"))


    input <- Sys.glob(glue("{proj_dir}/code_out/121_pathway_analysis/{curr_gemset_name}/{curr_norm_method}/{curr_res_name}/cluster_*/suppl_files/gsea_msigdb.tsv"))

    info_list <- list()
    for (i in seq_along(input)) {
        curr_gemset <- str_split(input[i], "/", simplify = TRUE)[[10]]
        curr_norm <- str_split(input[i], "/", simplify = TRUE)[[11]]
        curr_res <- str_split(input[i], "/", simplify = TRUE)[[12]]
        curr_comparison <- str_split(input[i], "/", simplify = TRUE)[[13]]
        curr_comparison_type <- if_else(str_detect(curr_comparison, "all$"), "one_vs_all", "pairwise")
        curr_path <- input[i]
        curr_comparison_split <- str_remove(curr_comparison, "cluster_")
        curr_comp1 <- str_split(curr_comparison_split, "_vs_", simplify = TRUE)[[1]]
        curr_comp2 <- str_split(curr_comparison_split, "_vs_", simplify = TRUE)[[2]]

        info_list[[i]] <- tibble(
            gemset = curr_gemset,
            norm = curr_norm,
            res = curr_res,
            comparison = curr_comparison,
            comp1 = curr_comp1,
            comp2 = curr_comp2,
            comparison_type = curr_comparison_type,
            path = curr_path
        )
    }
    info <- bind_rows(info_list)


    info2 <- info %>%
        dplyr::filter(comparison_type == "pairwise")

    comp2 <- info2 %>%
        dplyr::pull(comparison)

    gsea_list <- list()
    for (i in seq_len(nrow(info2))) {
        curr_info2 <- info2 %>%
            slice(i)

        # Add metadata to GSEA output
        curr_data <- read_tsv(info2$path[i]) %>%
            dplyr::full_join(curr_info2, by = character())

        gsea_list[[i]] <- curr_data
        names(gsea_list)[i] <- info2$comparison[i]
    }
    gsea <- bind_rows(gsea_list)

    # Find a unique set of gene sets with results (gene set names in column 1 of each list)
    sets_all <- gsea %>%
        dplyr::select(ID, gene_set_database) %>%
        distinct() %>%
        dplyr::mutate(database_geneset = glue("{gene_set_database}_{ID}"))

    # Remove some gene sets
    sets <- sets_all %>%
        dplyr::filter(!str_detect(gene_set_database, "^c7"))

    my_colors <- colorRampPalette(c("royalblue1", "white", "firebrick1"))(300)

    # for (i in str_which(sets$database_geneset, "mueller_AnergicSignature")) {
    for (i in seq_len(nrow(sets))) {
        curr_geneset <- sets$ID[i]
        curr_database <- sets$gene_set_database[i]

        filter_gsea <- gsea %>%
            dplyr::filter(ID == curr_geneset) %>%
            dplyr::select(comp1, comp2, NES)

        reverse_comparison <- filter_gsea %>%
            dplyr::rename(
                comp1_orig = "comp1",
                comp2_orig = "comp2"
            ) %>%
            dplyr::rename(
                comp1 = "comp2_orig",
                comp2 = "comp1_orig"
            ) %>%
            dplyr::mutate(NES = -1 * NES)

        combined <- dplyr::full_join(filter_gsea, reverse_comparison)

        all <- combined %>%
            tidyr::expand(comp1, comp2)

        all_combined <- full_join(combined, all)

        curr_gsea <- all_combined %>%
            dplyr::mutate(comp1 = forcats::fct_expand(comp1, unique(c(comp1, comp2)))) %>%
            dplyr::mutate(comp2 = forcats::fct_expand(comp2, unique(c(comp1, comp2)))) %>%
            dplyr::mutate(comp1 = forcats::fct_relevel(comp1, gtools::mixedsort)) %>%
            dplyr::mutate(comp2 = fct_rev(forcats::fct_relevel(comp2, gtools::mixedsort)))

        my_limits <- c(-abs(max(curr_gsea$NES)), abs(max(curr_gsea$NES)))

        p <- curr_gsea %>%
            ggplot(aes(comp1, comp2, fill = NES)) +
            geom_tile() +
            scale_x_discrete(position = "top") +
            scale_fill_gradientn(name = "Normalized\nEnrichment Score", colors = my_colors, limits = my_limits, breaks = scales::breaks_extended(n = 8), na.value = "grey75") +
            labs(
                title = glue("GSEA\n{curr_geneset}\ncluster_X_vs_cluster_Y"),
                x = glue("Cluster X"),
                y = glue("Cluster Y")
            ) +
            coord_fixed() +
            theme_minimal()

        if (!dir.exists(glue("{curr_database}"))) {
            dir.create(glue("{curr_database}"), recursive = TRUE)
        }
        pdf(glue("{curr_database}/{prefix}{curr_geneset}.pdf"))
        print(p)
        dev.off()
    }
}

#
# x <- gsea %>%
#     dplyr::filter(comparison == "cluster_7_vs_11") %>%
#     dplyr::arrange(NES) %>%
#     dplyr::mutate(ID = forcats::fct_reorder(factor(ID), NES))
#
#
# label <-  ifelse(x$ID %in% levels(x$ID)[1:5], x$ID, "")
#
#
# p <- x %>%
#     ggplot(aes(ID, NES, fill = qvalue)) +
#     geom_point()
#
#



#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
