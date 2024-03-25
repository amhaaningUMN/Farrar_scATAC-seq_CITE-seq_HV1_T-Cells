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
prefix <- "012_"
out <- glue("{prefix}count_fb")
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
# Use kallisto and bustools to count HTO and ADT barcodes
#######################################################################


# ---------------------------------------------------------------------
# Count the number of HTO and ADT tags
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/kallisto_bustools"))) {dir.create(glue("{out_dir}/kallisto_bustools"), recursive = TRUE)}
setwd(glue("{out_dir}/kallisto_bustools"))


all_feature_ref <- samples %>%
    dplyr::filter(seq_library == "FB") %>%
    dplyr::select(sample_name, biolegend_barcode)

write_csv(all_feature_ref, glue("{prefix}all_feature_ref.csv"))

# Run "kite" to create an include list of HTO/ADT barcodes allowing for 1 bp errors
file_t2g <- glue("{prefix}features_mismatch.t2g")
file_fa <- glue("{prefix}features_mismatch.fa")
file_feature <- glue("{prefix}all_feature_ref.csv")
system_kite_featuremap <- toddr::robust_system(glue("featuremap.py --t2g {file_t2g} --fa {file_fa} {file_feature} --header"))
toddr::robust_system_check(system_kite_featuremap)



# Create kallisto index
barcode_nt_length <- 15
file_idx <- glue("{prefix}features_mismatch.idx")
system_kallisto_index <- toddr::robust_system(glue("kallisto index -i {file_idx} -k {barcode_nt_length} {file_fa}"))
toddr::robust_system_check(system_kallisto_index)

system_kallisto_inspect <- toddr::robust_system(glue("kallisto inspect {file_idx}"))
toddr::robust_system_check(system_kallisto_inspect)


# https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-
# The 10X Genomics kit used, was:
# "Single Cell 5' R2-only", "SC5P-R2",
# "whitelist": "737K-august-2016"
# Version 2 chemistry
# # 737k-august-2016.txt:	Single Cell 3' v2, Single Cell 5' v1 and v2, Single Cell 5' HT v2

fastqs <- samples %>%
    dplyr::filter(seq_library == "FB") %>%
    dplyr::select(fastq_r1, fastq_r2) %>%
    dplyr::distinct()

# Create BUS file
bus_out_dir <- "."
chem_ver <- "10xv2"
system_kallisto_bus <- toddr::robust_system(glue("kallisto bus -i {file_idx} -o {bus_out_dir} -x {chem_ver} -t {curr_threads} --verbose {fastqs$fastq_r1} {fastqs$fastq_r2}"))
toddr::robust_system_check(system_kallisto_bus)



# Correct and limit barcodes to those provided in whitelist
# Use the whitelist provided by 10X genomics listing all possible barcodes for this chemistry
barcode_include_list <- glue("{out_dir}/737K-august-2016.txt")
system_bustools_correct <- toddr::robust_system(glue("bustools correct -w {barcode_include_list} -o output_corrected.bus output.bus"))
toddr::robust_system_check(system_bustools_correct)

# Sort
system_bustools_sort <- toddr::robust_system(glue("bustools sort -t {curr_threads} -o output_sorted.bus output_corrected.bus"))
toddr::robust_system_check(system_bustools_sort)

# Count
system_bustools_count <- toddr::robust_system(glue("bustools count -o {prefix}raw_fb --genecounts -g {file_t2g} -e matrix.ec -t transcripts.txt output_sorted.bus"))
toddr::robust_system_check(system_bustools_count)

# Delete extraneous empty dir that bustools count creates
unlink(glue("{prefix}raw_fb"), recursive = TRUE)




#######################################################################
# Explore the HTO ADT counts data generated by kallisto bustools
#######################################################################

# ---------------------------------------------------------------------
# Import the ADT HTO counts tables into R list
# ---------------------------------------------------------------------

fb_raw_allfeatures <- list()
# Import HTO and ADT counts
fb_raw_allfeatures$features <- read_tsv(glue("{out_dir}/kallisto_bustools/{prefix}raw_fb.genes.txt"), col_names = FALSE) %>%
    dplyr::rename(unique_id = "X1")

fb_raw_allfeatures$barcodes <- read_tsv(glue("{out_dir}/kallisto_bustools/{prefix}raw_fb.barcodes.txt"), col_names = FALSE) %>%
    dplyr::rename(barcode = "X1")

fb_raw_allfeatures$matrix <- Matrix::readMM(glue("{out_dir}/kallisto_bustools/{prefix}raw_fb.mtx")) %>%
    as.matrix(.) %>%
    t(.) %>%
    as(object = ., Class = "dgCMatrix") %>%
    magrittr::set_rownames(fb_raw_allfeatures$features$unique_id) %>%
    magrittr::set_colnames(fb_raw_allfeatures$barcodes$barcode)






# ---------------------------------------------------------------------
# Export raw counts list
# ---------------------------------------------------------------------

if (!dir.exists(glue("{out_dir}/explore_counts"))) {dir.create(glue("{out_dir}/explore_counts"), recursive = TRUE)}
setwd(glue("{out_dir}/explore_counts"))

saveRDS(fb_raw_allfeatures, glue("{prefix}fb_raw_allfeatures.rds"))





# ---------------------------------------------------------------------
# Plot raw count totals
# ---------------------------------------------------------------------




# For each HTO/ADT, sum the UMI counts across all possible GEMs
fb_rowsums <- rowSums(fb_raw_allfeatures$matrix)

# For each HTO/ADT, sum the GEMs that contain more than 10 UMIs. (This value represents
# the number of unique barcode GEMs that contain some useful data).
fb_rowsums_gt10 <- rowSums(fb_raw_allfeatures$matrix > 10)



raw_fb_counts <- tibble(feature = factor(names(fb_rowsums), levels = names(fb_rowsums)),
    fb = fb_rowsums,
    fb_rowsums_gt10 = fb_rowsums_gt10)


write_tsv(raw_fb_counts, glue("{prefix}sum_total_fb.txt"))


p <- raw_fb_counts %>%
    tidyr::pivot_longer(!feature) %>%
    dplyr::filter(!str_detect(name, "gt10")) %>%
    ggplot(aes(x = forcats::fct_rev(feature), y = value, fill = name)) +
    geom_col() +
    facet_grid(cols = vars(name)) +
    coord_flip() +
    guides(fill = guide_legend("Seq Lib")) +
    labs(title = "Total UMI counts for each ADT/HTO tag",
        subtitle = "Data colored by sequencing library name",
        x = "ADT or HTO feature name",
        y = "Sum total of all UMI counts for each ADT/HTO tag, across all GEMs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(glue("{prefix}sum_total_fb_umi.pdf"))
print(p)
dev.off()



p <- raw_fb_counts %>%
    tidyr::pivot_longer(!feature) %>%
    dplyr::filter(str_detect(name, "gt10")) %>%
    ggplot(aes(x = forcats::fct_rev(feature), y = value, fill = name)) +
    geom_col() +
    facet_grid(cols = vars(name)) +
    coord_flip() +
    guides(fill = guide_legend("Seq Lib")) +
    labs(title = "Number of barcodes with >10 ADT/HTO counts",
        subtitle = "Data colored by sequencing library name",
        x = "ADT or HTO feature name",
        y = "Sum total of unique barcodes with >10 ADT/HTO counts") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(glue("{prefix}sum_total_fb_gems.pdf"))
print(p)
dev.off()














#######################################################################
# Save session info
#######################################################################


toddr::write_session_info(glue("{out_dir}/{prefix}"))










