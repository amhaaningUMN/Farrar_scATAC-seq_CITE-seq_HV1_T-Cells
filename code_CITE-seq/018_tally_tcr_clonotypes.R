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
library(gtools)
library(openxlsx)



#######################################################################
# Script parameters
#######################################################################


proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "018_"
out <- glue("{prefix}tally_tcr_clonotypes")
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
# Get data
# ---------------------------------------------------------------------


hto <- readRDS(glue("{proj_dir}/code_out/016_demultiplex_hto/hto/match_gex_filtered/016_gmm_demux.rds")) %>%
    dplyr::select(barcode, hto_class1, hto_class2, hto_class3, gem_type, treatment, sample_name)

adt <- readRDS(glue("{proj_dir}/code_out/017_normalize_adt/017_adt_dsb_normalized.rds")) %>%
    # Transpose the matrix
    t(.) %>%
    tibble::as_tibble(rownames = "barcode") %>%
    # Add prefix to most colnames
    dplyr::rename_with(.fn = ~ paste0("adt", .x), .cols = -barcode)

contigs <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_contigs.rds")) %>%
    dplyr::left_join(hto, by = "barcode") %>%
    dplyr::mutate(across(c(hto_class1, hto_class2, hto_class3, gem_type, sample_name), as.character)) %>%
    dplyr::mutate(across(c(hto_class1, hto_class2, hto_class3, gem_type, sample_name), ~ replace_na(.x, "negative")))

clonotypes <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_clonotypes.rds"))
clonotype_id_cdr3aa_unique <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_clonotype_id_cdr3aa_unique.rds"))
clonotype_id_cdr3nt_unique <- readRDS(glue("{proj_dir}/code_out/014_tcr_clonotypes/014_clonotype_id_cdr3nt_unique.rds"))




# ---------------------------------------------------------------------
# Custom function
# ---------------------------------------------------------------------


tk_combine_counts_percentages <- function(counts, percentages) {
    c_and_p <- glue("{counts}\n({percentages})")
}


# ---------------------------------------------------------------------
# Tally clonotypes by raw_clonotype_id, split by sample
# ---------------------------------------------------------------------

adt_by_clonotype_id_raw <- contigs %>%
    dplyr::left_join(adt, by = "barcode") %>%
    group_by(raw_clonotype_id) %>%
    dplyr::summarize(across(starts_with("adt"), list(median_clonotype_id_raw = ~median(., na.rm = TRUE))))


# This mimics the web_summary.html
gem_tally_by_clonotype_id_raw <- contigs %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(raw_clonotype_id, hto_class2) %>%
    tally() %>%
    tidyr::pivot_wider(names_from = "hto_class2", names_glue = "{hto_class2}", values_from = "n") %>%
    ungroup() %>%
    dplyr::mutate(total_n_of_gems = rowSums(across(where(is.integer)), na.rm = TRUE)) %>%
    dplyr::mutate(prop_of_total_gems = total_n_of_gems / sum(total_n_of_gems)) %>%
    dplyr::left_join(dplyr::select(.data = clonotypes, clonotype_id, clonotype_id_cdr3aa, clonotype_id_cdr3nt, cdr3s_aa, cdr3s_nt), by = c("raw_clonotype_id" = "clonotype_id")) %>%
    dplyr::mutate(count_tra = str_count(cdr3s_aa, "TRA:")) %>%
    dplyr::mutate(count_trb = str_count(cdr3s_aa, "TRB:")) %>%
    dplyr::mutate(count_tra_trb = factor(glue("tra_{count_tra}_trb_{count_trb}"))) %>%
    dplyr::mutate(raw_clonotype_id_number = as.integer(str_extract(raw_clonotype_id, "\\d+"))) %>%
    dplyr::arrange(order(gtools::mixedorder(raw_clonotype_id))) %>%
    dplyr::mutate(clonotype_id_cdr3aa_number = as.integer(str_remove(clonotype_id_cdr3aa, fixed("clonotype_cdr3aa")))) %>%
    dplyr::mutate(clonotype_id_cdr3nt_number = as.integer(str_remove(clonotype_id_cdr3nt, fixed("clonotype_cdr3nt")))) %>%
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::relocate(raw_clonotype_id, raw_clonotype_id_number, clonotype_id_cdr3aa, clonotype_id_cdr3aa_number,
        clonotype_id_cdr3nt, clonotype_id_cdr3nt_number, negative, multi_sample_gem) %>%
    dplyr::relocate(c(total_n_of_gems, prop_of_total_gems, count_tra, count_trb, count_tra_trb, cdr3s_aa, cdr3s_nt), .after = last_col()) %>%
    dplyr::left_join(adt_by_clonotype_id_raw, by = "raw_clonotype_id")




write_tsv(gem_tally_by_clonotype_id_raw, glue("{prefix}gem_tally_by_clonotype_id_raw.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_raw.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_clonotype_id_raw), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 2, widths = 6)
openxlsx::setColWidths(wb, sheet = 1, cols = 3, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 4, widths = 6)
openxlsx::setColWidths(wb, sheet = 1, cols = 5, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 6, widths = 6)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)



# Tally clonotypes by number of TRA or TRB seqs
clonotype_id_raw_tally_by_tra_trb <- gem_tally_by_clonotype_id_raw %>%
    group_by(count_tra_trb) %>%
    tally()

write_tsv(clonotype_id_raw_tally_by_tra_trb, glue("{prefix}clonotype_id_raw_tally_by_tra_trb.txt"))
out_filename <- glue("{prefix}clonotype_id_raw_tally_by_tra_trb.xlsx")
wb <- openxlsx::write.xlsx(list(tally = clonotype_id_raw_tally_by_tra_trb), file = out_filename, rowNames = FALSE)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)




p <- gem_tally_by_clonotype_id_raw %>%
    ggplot(aes(fct_infreq(count_tra_trb))) +
    geom_bar(fill = "dodgerblue") +
    geom_text(stat = "count", aes(label = after_stat(tk_combine_counts_percentages(count, scales::percent(prop))), group = 1), size = 3, position = position_stack(vjust = 0.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
    labs(title = "Distribution of TRA and TRB seqs found in each 'clonotype_id_raw'", x = "", y = "Number of clonotypes that contain TRA and TRB seq profiles")
pdf(glue("{prefix}clonotype_id_raw_tally_by_tra_trb.pdf"), width = 7, height = 7)
print(p)
dev.off()





# ---------------------------------------------------------------------
# Tally clonotypes by clonotype_id_cdr3aa, split by sample
# ---------------------------------------------------------------------

adt_by_clonotype_id_cdr3aa <- contigs %>%
    dplyr::left_join(adt, by = "barcode") %>%
    group_by(clonotype_id_cdr3aa) %>%
    dplyr::summarize(across(starts_with("adt"), list(median_clonotype_id_raw = ~median(., na.rm = TRUE))))



# This mimics the web_summary.html
gem_tally_by_clonotype_id_cdr3aa <- contigs %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(clonotype_id_cdr3aa, hto_class2) %>%
    tally() %>%
    tidyr::pivot_wider(names_from = "hto_class2", names_glue = "{hto_class2}", values_from = "n") %>%
    ungroup() %>%
    dplyr::mutate(total_n_of_gems = rowSums(across(where(is.integer)), na.rm = TRUE)) %>%
    dplyr::mutate(prop_of_total_gems = total_n_of_gems / sum(total_n_of_gems)) %>%
    dplyr::left_join(clonotype_id_cdr3aa_unique, by = "clonotype_id_cdr3aa") %>%
    dplyr::mutate(count_tra = str_count(cdr3s_aa, "TRA:")) %>%
    dplyr::mutate(count_trb = str_count(cdr3s_aa, "TRB:")) %>%
    dplyr::mutate(count_tra_trb = factor(glue("tra_{count_tra}_trb_{count_trb}"))) %>%
    dplyr::arrange(order(gtools::mixedorder(clonotype_id_cdr3aa))) %>%
    dplyr::mutate(clonotype_id_cdr3aa_number = as.integer(str_remove(clonotype_id_cdr3aa, fixed("clonotype_cdr3aa")))) %>%
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::relocate(clonotype_id_cdr3aa, clonotype_id_cdr3aa_number,
        negative, multi_sample_gem) %>%
    dplyr::relocate(c(total_n_of_gems, prop_of_total_gems, count_tra, count_trb, count_tra_trb, cdr3s_aa), .after = last_col()) %>%
    dplyr::left_join(adt_by_clonotype_id_cdr3aa, by = "clonotype_id_cdr3aa")






write_tsv(gem_tally_by_clonotype_id_cdr3aa, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_cdr3aa.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_clonotype_id_cdr3aa), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 2, widths = 6)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)



# Tally clonotypes by number of TRA or TRB seqs
clonotype_id_cdr3aa_tally_by_tra_trb <- gem_tally_by_clonotype_id_cdr3aa %>%
    group_by(count_tra_trb) %>%
    tally()

write_tsv(clonotype_id_cdr3aa_tally_by_tra_trb, glue("{prefix}clonotype_id_cdr3aa_tally_by_tra_trb.txt"))
out_filename <- glue("{prefix}clonotype_id_cdr3aa_tally_by_tra_trb.xlsx")
wb <- openxlsx::write.xlsx(list(tally = clonotype_id_cdr3aa_tally_by_tra_trb), file = out_filename, rowNames = FALSE)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)




p <- gem_tally_by_clonotype_id_cdr3aa %>%
    ggplot(aes(fct_infreq(count_tra_trb))) +
    geom_bar(fill = "dodgerblue") +
    geom_text(stat = "count", aes(label = after_stat(tk_combine_counts_percentages(count, scales::percent(prop))), group = 1), size = 3, position = position_stack(vjust = 0.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
    labs(title = "Distribution of TRA and TRB seqs found in each 'clonotype_id_cdr3aa'", x = "", y = "Number of clonotypes that contain TRA and TRB seq profiles")
pdf(glue("{prefix}clonotype_id_cdr3aa_tally_by_tra_trb.pdf"), width = 7, height = 7)
print(p)
dev.off()








# ---------------------------------------------------------------------
# Tally clonotypes by clonotype_id_cdr3nt, split by sample
# ---------------------------------------------------------------------

adt_by_clonotype_id_cdr3nt <- contigs %>%
    dplyr::left_join(adt, by = "barcode") %>%
    group_by(clonotype_id_cdr3nt) %>%
    dplyr::summarize(across(starts_with("adt"), list(median_clonotype_id_raw = ~median(., na.rm = TRUE))))


# This mimics the web_summary.html
gem_tally_by_clonotype_id_cdr3nt <- contigs %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(clonotype_id_cdr3nt, hto_class2) %>%
    tally() %>%
    tidyr::pivot_wider(names_from = "hto_class2", names_glue = "{hto_class2}", values_from = "n") %>%
    ungroup() %>%
    dplyr::mutate(total_n_of_gems = rowSums(across(where(is.integer)), na.rm = TRUE)) %>%
    dplyr::mutate(prop_of_total_gems = total_n_of_gems / sum(total_n_of_gems)) %>%
    dplyr::left_join(clonotype_id_cdr3nt_unique, by = "clonotype_id_cdr3nt") %>%
    dplyr::mutate(count_tra = str_count(cdr3s_nt, "TRA:")) %>%
    dplyr::mutate(count_trb = str_count(cdr3s_nt, "TRB:")) %>%
    dplyr::mutate(count_tra_trb = factor(glue("tra_{count_tra}_trb_{count_trb}"))) %>%
    dplyr::arrange(order(gtools::mixedorder(clonotype_id_cdr3nt))) %>%
    dplyr::mutate(clonotype_id_cdr3nt_number = as.integer(str_remove(clonotype_id_cdr3nt, fixed("clonotype_cdr3nt")))) %>%
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::relocate(clonotype_id_cdr3nt, clonotype_id_cdr3nt_number,
        negative, multi_sample_gem) %>%
    dplyr::relocate(c(total_n_of_gems, prop_of_total_gems, count_tra, count_trb, count_tra_trb, cdr3s_nt), .after = last_col()) %>%
    dplyr::left_join(adt_by_clonotype_id_cdr3nt, by = "clonotype_id_cdr3nt")





write_tsv(gem_tally_by_clonotype_id_cdr3nt, glue("{prefix}gem_tally_by_clonotype_id_cdr3nt.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_cdr3nt.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_clonotype_id_cdr3nt), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 2, widths = 6)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)



# Tally clonotypes by number of TRA or TRB seqs
clonotype_id_cdr3nt_tally_by_tra_trb <- gem_tally_by_clonotype_id_cdr3nt %>%
    group_by(count_tra_trb) %>%
    tally()

write_tsv(clonotype_id_cdr3nt_tally_by_tra_trb, glue("{prefix}clonotype_id_cdr3nt_tally_by_tra_trb.txt"))
out_filename <- glue("{prefix}clonotype_id_cdr3nt_tally_by_tra_trb.xlsx")
wb <- openxlsx::write.xlsx(list(tally = clonotype_id_cdr3nt_tally_by_tra_trb), file = out_filename, rowNames = FALSE)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)




p <- gem_tally_by_clonotype_id_cdr3nt %>%
    ggplot(aes(fct_infreq(count_tra_trb))) +
    geom_bar(fill = "dodgerblue") +
    geom_text(stat = "count", aes(label = after_stat(tk_combine_counts_percentages(count, scales::percent(prop))), group = 1), size = 3, position = position_stack(vjust = 0.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
    labs(title = "Distribution of TRA and TRB seqs found in each 'clonotype_id_cdr3nt'", x = "", y = "Number of clonotypes that contain TRA and TRB seq profiles")
pdf(glue("{prefix}clonotype_id_cdr3nt_tally_by_tra_trb.pdf"), width = 7, height = 7)
print(p)
dev.off()






# ---------------------------------------------------------------------
# Find top clonotypes by clonotype_id_cdr3aa, split by sample
# ---------------------------------------------------------------------



# Find the top number of clonotypes by sample
gem_tally_by_clonotype_id_cdr3aa_sample_name_top <- contigs %>%
    distinct(barcode, .keep_all = TRUE) %>%
    group_by(hto_class2, clonotype_id_cdr3aa) %>%
    tally() %>%
    arrange(desc(n)) %>%
    slice_head(n = 5) %>%
    dplyr::rename(n_gems_with_clonotype = "n")


saveRDS(gem_tally_by_clonotype_id_cdr3aa_sample_name_top, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.rds"))

write_tsv(gem_tally_by_clonotype_id_cdr3aa_sample_name_top, glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.txt"))
out_filename <- glue("{prefix}gem_tally_by_clonotype_id_cdr3aa_sample_name_top.xlsx")
wb <- openxlsx::write.xlsx(list(tally = gem_tally_by_clonotype_id_cdr3aa_sample_name_top), file = out_filename, rowNames = FALSE)
openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
openxlsx::setColWidths(wb, sheet = 1, cols = 2, widths = 30)
openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)




#######################################################################
# Heatmaps
#######################################################################

# ---------------------------------------------------------------------
# Number of GEMs per clonotype with ADT: clonotype_id_raw
# ---------------------------------------------------------------------

no_pivot_cols1 <- gem_tally_by_clonotype_id_raw %>%
    dplyr::select(starts_with("adt")) %>%
    colnames()

no_pivot_cols2 <- c("raw_clonotype_id", "raw_clonotype_id_number", "clonotype_id_cdr3aa", "clonotype_id_cdr3aa_number", "clonotype_id_cdr3nt", "clonotype_id_cdr3nt_number",
    "total_n_of_gems", "prop_of_total_gems", "count_tra", "count_trb", "count_tra_trb", "cdr3s_aa", "cdr3s_nt")

no_pivot_cols <- c(no_pivot_cols1, no_pivot_cols2)

p1_data <- gem_tally_by_clonotype_id_raw %>%
    tidyr::pivot_longer(cols = !all_of(no_pivot_cols)) %>%
    # Show only clonotypes with at least 1 GEM
    dplyr::filter(value > 1 & !is.na(value)) %>%
    dplyr::mutate(raw_clonotype_id = fct_rev(fct_relevel(raw_clonotype_id, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort))

p1 <- p1_data %>%
    ggplot(aes(x = name, y = raw_clonotype_id, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Number of\nGEMs", colors = c(viridis::viridis(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "Original cellranger based clonotype\n(clonotype_id_raw)")




p2 <- gem_tally_by_clonotype_id_raw %>%
    # Rename cols, by removing suffix
    dplyr::rename_with(.fn = ~ str_remove(.x, "_median_clonotype_id_raw"), .cols = starts_with("adt")) %>%
    tidyr::pivot_longer(cols = starts_with("adt")) %>%
    # Keep only same clonotypes as plotted above
    dplyr::filter(raw_clonotype_id %in% p1_data$raw_clonotype_id) %>%
    dplyr::mutate(raw_clonotype_id = fct_rev(fct_relevel(raw_clonotype_id, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort)) %>%
    ggplot(aes(x = name, y = raw_clonotype_id, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Within clonotype\nmedian of dsb\nnormalized ADT values", colors = c(viridis::plasma(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "")

patchwork <- patchwork::wrap_plots(list(p1, p2), ncol = 2) +
    patchwork::plot_layout(widths = c(4, 1),
        guides = "collect") &
    patchwork::plot_annotation(title = glue("clonotype_id_raw"),
        subtitle = "Only clonotypes with >1 GEM are shown") &
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

pdf(glue("{prefix}gem_tally_with_adt_heatmap_clonotype_id_raw.pdf"), width = 7, height = 11)
print(patchwork)
dev.off()





# ---------------------------------------------------------------------
# Number of GEMs per clonotype with ADT: clonotype_id_cdr3aa
# ---------------------------------------------------------------------


no_pivot_cols1 <- gem_tally_by_clonotype_id_raw %>%
    dplyr::select(starts_with("adt")) %>%
    colnames()

no_pivot_cols2 <- c("clonotype_id_cdr3aa", "clonotype_id_cdr3aa_number",
    "total_n_of_gems", "prop_of_total_gems", "count_tra", "count_trb", "count_tra_trb", "cdr3s_aa")

no_pivot_cols <- c(no_pivot_cols1, no_pivot_cols2)

p1_data <- gem_tally_by_clonotype_id_cdr3aa %>%
    tidyr::pivot_longer(cols = !all_of(no_pivot_cols)) %>%
    # Show only clonotypes with at least 1 GEM
    dplyr::filter(value > 1 & !is.na(value)) %>%
    dplyr::mutate(clonotype_id_cdr3aa = fct_rev(fct_relevel(clonotype_id_cdr3aa, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort))
p1 <- p1_data %>%
    ggplot(aes(x = name, y = clonotype_id_cdr3aa, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Number of\nGEMs", colors = c(viridis::viridis(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "Clonotype based on CDR3 amino acid seq\n(clonotype_id_cdr3aa)")



p2 <- gem_tally_by_clonotype_id_cdr3aa %>%
    # Rename cols, by removing suffix
    dplyr::rename_with(.fn = ~ str_remove(.x, "_median_clonotype_id_raw"), .cols = starts_with("adt")) %>%
    tidyr::pivot_longer(cols = starts_with("adt")) %>%
    # Keep only same clonotypes as plotted above
    dplyr::filter(clonotype_id_cdr3aa %in% p1_data$clonotype_id_cdr3aa) %>%
    dplyr::mutate(clonotype_id_cdr3aa = fct_rev(fct_relevel(clonotype_id_cdr3aa, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort)) %>%
    ggplot(aes(x = name, y = clonotype_id_cdr3aa, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Within clonotype\nmedian of dsb\nnormalized ADT values", colors = c(viridis::plasma(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "")

patchwork <- patchwork::wrap_plots(list(p1, p2), ncol = 2) +
    patchwork::plot_layout(widths = c(4, 1),
        guides = "collect") &
    patchwork::plot_annotation(title = glue("clonotype_id_cdr3aa"),
        subtitle = "Only clonotypes with >1 GEM are shown") &
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

pdf(glue("{prefix}gem_tally_with_adt_heatmap_clonotype_id_cdr3aa.pdf"), width = 7, height = 11)
print(patchwork)
dev.off()






# ---------------------------------------------------------------------
# Number of GEMs per clonotype with ADT: clonotype_id_cdr3nt
# ---------------------------------------------------------------------


no_pivot_cols1 <- gem_tally_by_clonotype_id_raw %>%
    dplyr::select(starts_with("adt")) %>%
    colnames()

no_pivot_cols2 <- c("clonotype_id_cdr3nt", "clonotype_id_cdr3nt_number",
    "total_n_of_gems", "prop_of_total_gems", "count_tra", "count_trb", "count_tra_trb", "cdr3s_nt")

no_pivot_cols <- c(no_pivot_cols1, no_pivot_cols2)

p1_data <- gem_tally_by_clonotype_id_cdr3nt %>%
    tidyr::pivot_longer(cols = !all_of(no_pivot_cols)) %>%
    # Show only clonotypes with at least 1 GEM
    dplyr::filter(value > 1 & !is.na(value)) %>%
    dplyr::mutate(clonotype_id_cdr3nt = fct_rev(fct_relevel(clonotype_id_cdr3nt, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort))
p1 <- p1_data %>%
    ggplot(aes(x = name, y = clonotype_id_cdr3nt, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Number of\nGEMs", colors = c(viridis::viridis(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "Clonotype based on CDR3 nucleotide seq\n(clonotype_id_cdr3nt)")



p2 <- gem_tally_by_clonotype_id_cdr3nt %>%
    # Rename cols, by removing suffix
    dplyr::rename_with(.fn = ~ str_remove(.x, "_median_clonotype_id_raw"), .cols = starts_with("adt")) %>%
    tidyr::pivot_longer(cols = starts_with("adt")) %>%
    # Keep only same clonotypes as plotted above
    dplyr::filter(clonotype_id_cdr3nt %in% p1_data$clonotype_id_cdr3nt) %>%
    dplyr::mutate(clonotype_id_cdr3nt = fct_rev(fct_relevel(clonotype_id_cdr3nt, gtools::mixedsort))) %>%
    dplyr::mutate(name = fct_relevel(name, gtools::mixedsort)) %>%
    ggplot(aes(x = name, y = clonotype_id_cdr3nt, fill = value)) +
    scale_x_discrete(position = "top") +
    geom_tile() +
    scale_fill_gradientn(name = "Within clonotype\nmedian of dsb\nnormalized ADT values", colors = c(viridis::plasma(100))) +
    theme_minimal() +
    # hjust = 0 makes the labels align at the top, next to the plot
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        axis.text.y = element_text(size = 0.5, vjust = 0.5, hjust = 1,
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.y = element_line(colour = "grey92")) +
    labs(title = "",
        x = "",
        y = "")

patchwork <- patchwork::wrap_plots(list(p1, p2), ncol = 2) +
    patchwork::plot_layout(widths = c(4, 1),
        guides = "collect") &
    patchwork::plot_annotation(title = glue("clonotype_id_cdr3nt"),
        subtitle = "Only clonotypes with >1 GEM are shown") &
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

pdf(glue("{prefix}gem_tally_with_adt_heatmap_clonotype_id_cdr3nt.pdf"), width = 7, height = 11)
print(patchwork)
dev.off()






#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))
