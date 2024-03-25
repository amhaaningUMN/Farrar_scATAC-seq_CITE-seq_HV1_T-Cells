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
library(openxlsx)
library(glue)




#######################################################################
# Script parameters
#######################################################################



proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "010_"
out <- glue("{prefix}samples")
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



samples <- read_tsv(glue("{proj_dir}/input/metadata.tsv")) %>%
    purrr::modify_at(vars(everything()), factor)


write_tsv(samples, file = glue("{out_dir}/{prefix}samples.txt"))
openxlsx::write.xlsx(list(samples = samples), glue("{out_dir}/{prefix}samples.xlsx"))
saveRDS(samples, file = glue("{out_dir}/{prefix}samples.rds"))



#######################################################################
# Save session info
#######################################################################



# ---------------------------------------------------------------------
# Write out session info
# ---------------------------------------------------------------------


toddr::write_session_info(glue("{out_dir}/{prefix}"))

