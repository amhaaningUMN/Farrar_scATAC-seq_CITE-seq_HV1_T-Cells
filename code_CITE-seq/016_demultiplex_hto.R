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
prefix <- "016_"
out <- glue("{prefix}demultiplex_hto")
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
# HTO demultiplexing with GMM-demux
#######################################################################

# Seurat CLR function
# https://github.com/satijalab/seurat/blob/a1294c4d363780548dbf9cc4a4abb3a6078a6d64/R/preprocessing.R#L2484-L2486
clr_function <- function(x) {
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}






# Get data
hto_raw <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_hto_raw.rds"))
gex_filtered <- readRDS(glue("{proj_dir}/code_out/015_tidy_count_tables/015_gex_filtered.rds"))

barcode_subset <- list(TRUE,
    which(colSums(hto_raw$matrix) > 100),
    which(colnames(hto_raw$matrix) %in% colnames(gex_filtered$matrix)))
names(barcode_subset) <- c("all", "greater_than_100umi", "match_gex_filtered")

for (j in seq_along(barcode_subset)) {

    curr_barcode_subset <- barcode_subset[[j]]
    curr_barcode_subset_name <- names(barcode_subset)[j]

    if (!dir.exists(glue("{out_dir}/hto/{curr_barcode_subset_name}"))) {dir.create(glue("{out_dir}/hto/{curr_barcode_subset_name}"), recursive = TRUE)}
    setwd(glue("{out_dir}/hto/{curr_barcode_subset_name}"))

    hto <- samples %>%
        dplyr::filter(application == "hashtag") %>%
        dplyr::pull(sample_name) %>%
        as.character(.)

    hto_colnames <- paste0(hto, "_hto")


    # Export HTOs
    curr_features <- hto_raw$features
    curr_barcodes <- colnames(hto_raw$matrix)[curr_barcode_subset]
    curr_matrix <- hto_raw$matrix[, curr_barcode_subset]

    curr_matrix_clr <- apply(
        X = curr_matrix,
        MARGIN = 1,
        FUN = clr_function)
    curr_matrix_clr <- Matrix::t(curr_matrix_clr)


    # Histogram
    curr_matrix_long <- curr_matrix %>%
        as.matrix(.) %>%
        tibble::as_tibble(rownames = "unique_id") %>%
        tidyr::pivot_longer(!unique_id, values_to = "hto_counts")

    curr_matrix_clr_long <- curr_matrix_clr %>%
        tibble::as_tibble(rownames = "unique_id") %>%
        tidyr::pivot_longer(!unique_id, values_to = "hto_counts_clr")

    curr_matrix_tbl <- curr_matrix_long %>%
        dplyr::full_join(curr_matrix_clr_long, by = c("unique_id", "name"))

    saveRDS(curr_matrix_tbl, glue("{prefix}hto_counts_tbl.rds"))
    write_tsv(curr_matrix_tbl, glue("{prefix}hto_counts_tbl.txt"))

    curr_matrix_tbl_wide <- curr_matrix_tbl %>%
        tidyr::pivot_wider(names_from = "unique_id", values_from = c("hto_counts", "hto_counts_clr")) %>%
        dplyr::rename(barcode = "name")


    p <- curr_matrix_tbl %>%
        ggplot(aes(x = hto_counts_clr)) +
        geom_histogram() +
        facet_wrap(~ unique_id) +
        labs(title = glue("{curr_barcode_subset_name}\nDistribution of centered log ratio (CLR)\ntransformed UMI counts split by HTO feature"),
            x = "CLR transformed expression (30 bins)",
            y = "Number of GEMs with expression value")

    pdf(glue("{prefix}clr_transformed_counts_histogram.pdf"), width = 8, height = 8)
    print(p)
    dev.off()


    # GMM-demux (export raw counts, not CLR transformed)
    write_tsv(tibble(curr_features), glue("features.tsv.gz"), col_names = FALSE)
    write_tsv(tibble(curr_barcodes), glue("barcodes.tsv.gz"), col_names = FALSE)
    sparse_mtx <- Matrix::Matrix(curr_matrix, sparse = TRUE)
    Matrix::writeMM(sparse_mtx, file = glue("matrix.mtx"))
    system(glue("gzip matrix.mtx"))

    matrix_dir <- "."
    hto_names <- read_tsv(glue("{matrix_dir}/features.tsv.gz"), col_names = FALSE) %>%
        dplyr::pull() %>%
        paste(., collapse = ",")
    system_gmm_demux <- toddr::robust_system(glue("GMM-demux -o 'SSD_mtx' -f 'SSD_mtx_report' {matrix_dir} {hto_names}"))
    toddr::robust_system_check(system_gmm_demux)

    gmm_demux_config <- read_csv(glue("SSD_mtx_report/GMM_full.config"), col_names = FALSE) %>%
        dplyr::rename(gmm_demux_hto_cluster_id = "X1", hto_class1 = "X2") %>%
        # Add columns for each of the possible HTO tags
        tibble::add_column(!!!magrittr::set_names(as.list(rep("", length(hto_colnames))), nm = hto_colnames))

    # Label whether a barcode includes a particular HTO tag, in individual columns
    gmm_demux_config2 <- gmm_demux_config
    for (i in seq_len(dim(gmm_demux_config)[1])) {
        for (j in seq_along(hto_colnames)) {
            if (str_detect(gmm_demux_config$hto_class1[i], fixed(hto[j]))) {
                gmm_demux_config2[i, hto_colnames[j]] <- "yes"
            } else {
                gmm_demux_config2[i, hto_colnames[j]] <- NA
            }
        }
    }


    # Read the gmm_demux results
    gmm_demux <- read_csv(glue("SSD_mtx_report/GMM_full.csv")) %>%
        dplyr::rename(barcode = "...1", gmm_demux_hto_cluster_id = "Cluster_id", gmm_demux_confidence = "Confidence") %>%
        dplyr::left_join(gmm_demux_config2, by = "gmm_demux_hto_cluster_id") %>%
        dplyr::mutate(gem_type = case_when(
            hto_class1 == "negative" ~ "negative",
            !str_detect(hto_class1, "-") & hto_class1 != "negative" ~ "same_sample_gem",
            str_detect(hto_class1, "-") ~ "multi_sample_gem",
            TRUE ~ NA_character_
        )) %>%
        dplyr::mutate(hto_class2 = case_when(
            gem_type == "same_sample_gem" ~ hto_class1,
            gem_type == "negative" ~ "negative",
            gem_type == "multi_sample_gem" ~ "multi_sample_gem",
            TRUE ~ NA_character_
        )) %>%
        dplyr::left_join(samples, by = c("hto_class1" = "sample_name"), keep = TRUE) %>%
        purrr::modify_at(vars(!any_of("gmm_demux_confidence")), factor) %>%
        dplyr::mutate(hto_class1 = fct_relevel(hto_class1, gtools::mixedsort)) %>%
        dplyr::mutate(hto_class2 = fct_relevel(hto_class2, gtools::mixedsort)) %>%
        # Start with all possible HTO names, then lump all levels except for the 'n' most frequent (i.e. most number of GEMs labeled with the HTO)
        dplyr::mutate(hto_class3 = forcats::fct_rev(forcats::fct_infreq(forcats::fct_lump_n(hto_class1, n = 50, other_level = "other", ties.method = "min")))) %>%
        dplyr::left_join(curr_matrix_tbl_wide, by = "barcode")


    saveRDS(gmm_demux, glue("{prefix}gmm_demux.rds"))
    write_tsv(gmm_demux, glue("{prefix}gmm_demux.txt"))


    # ---------------------------------------------------------------------
    # UMI counts box plots
    # ---------------------------------------------------------------------

    p <- gmm_demux %>%
        dplyr::select(barcode, hto_class2, starts_with("hto_counts_clr")) %>%
        tidyr::pivot_longer(!c(barcode, hto_class2)) %>%
        dplyr::mutate(name = fct_relevel(name, gtools::mixedsort)) %>%
        ggplot(aes(x = name, y = value, fill = name)) +
        geom_boxplot() +
        facet_wrap(~ hto_class2) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = glue("{curr_barcode_subset_name}\nDistribution of centered log ratio (CLR)\ntransformed UMI counts split by HTO classification call (via GMM-demux)"),
            x = "GEMs grouped by HTO classification call (GMM-demux)",
            y = "CLR transformed expression")

    pdf(glue("{prefix}clr_transformed_counts_boxplot.pdf"), width = 12, height = 12)
    print(p)
    dev.off()


    # ---------------------------------------------------------------------
    # Summary stats
    # ---------------------------------------------------------------------


    gmm_demux_tally1 <- gmm_demux %>%
        group_by(hto_class2) %>%
        tally(name = "number_of_gems") %>%
        arrange(desc(number_of_gems))

    write_tsv(gmm_demux_tally1, glue("{prefix}gem_tally_by_hto_class2.txt"))



    p <- gmm_demux %>%
        ggplot(aes(x = forcats::fct_infreq(hto_class2))) +
            geom_bar(fill = "dodgerblue") +
            ggrepel::geom_text_repel(stat = "count", aes(label = after_stat(count)), size = 3, position = position_stack(vjust = 1), point.size = NA) +
            coord_flip() +
            labs(title = "Number of barcodes (GEMs)\nclassified by HTO expression", x = "HTO Classification", y = "Number of barcodes")
    pdf(glue("{prefix}gem_tally_by_hto_class2_single_sample_gems_barplot.pdf"), width = 7, height = 7)
    print(p)
    dev.off()


    gmm_demux_tally4 <- gmm_demux %>%
        group_by(gem_type) %>%
        tally(name = "number_of_gems") %>%
        arrange(desc(number_of_gems))

    p <- gmm_demux %>%
        ggplot(aes(x = forcats::fct_infreq(gem_type), fill = gem_type)) +
            geom_bar() +
            geom_text(stat = "count", aes(label = after_stat(count)), size = 3, position = position_stack(vjust = 0.5)) +
            theme(legend.position = "none") +
            labs(title = "Number of barcodes (GEMs)\nclassified by HTO expression", x = "HTO Classification", y = "Number of barcodes")
    pdf(glue("{prefix}gem_tally_by_gem_type_barplot.pdf"), width = 7, height = 7)
    print(p)
    dev.off()




    # ---------------------------------------------------------------------
    # Plot GMM-demux confidence values
    # ---------------------------------------------------------------------



    # GMM-demux confidence values
    p <- gmm_demux %>%
        ggplot(aes(x = gmm_demux_confidence, y = hto_class1, fill = after_stat(x))) +
            ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
            scale_fill_viridis_c(option = "viridis") +
            # ggridges theme modifications
            scale_x_continuous(expand = c(0, 1)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_cartesian(clip = "off") +
            theme_minimal(base_size = 14) +
            theme(axis.text.y = element_text(vjust = 0)) +
            theme(legend.position = "none") +
            labs(title = "Density of GMM-demux HTO-Calling\nConfidence Values by HTO class", x = "GMM-Demux Confidence [0, 1.0]", y = "Proportion of barcodes, split by HTO class\n(sorted by most frequent HTO class)")
    pdf(glue("{prefix}gmm_demux_confidence_by_hto_class1_ridgeplot.pdf"), width = 10, height = 8)
    print(p)
    dev.off()



    # GMM-demux confidence values
    p <- gmm_demux %>%
        ggplot(aes(x = gmm_demux_confidence, y = hto_class2, fill = after_stat(x))) +
            ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
            scale_fill_viridis_c(option = "viridis") +
            # ggridges theme modifications
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_cartesian(clip = "off") +
            theme_minimal(base_size = 14) +
            theme(axis.text.y = element_text(vjust = 0)) +
            theme(legend.position = "none") +
            labs(title = "Density of GMM-demux HTO-Calling\nConfidence Values by HTO class", x = "GMM-Demux Confidence [0, 1.0]", y = "Proportion of barcodes, split by HTO class\n(sorted by most frequent HTO class)")
    pdf(glue("{prefix}gmm_demux_confidence_by_hto_class2_ridgeplot.pdf"), width = 10, height = 8)
    print(p)
    dev.off()
}







#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
