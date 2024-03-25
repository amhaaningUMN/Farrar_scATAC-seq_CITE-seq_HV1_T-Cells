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
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)
library(msigdbr)
library(VISION)
library(pathwayPCA)





#######################################################################
# Script parameters
#######################################################################



proj <- "cd4_hv1_nilotinib_pdl1_il10_citeseq_20231201"
prefix <- "026_"
out <- glue("{prefix}clean_gene_sets_and_msigdb")
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
# msigdbr
# ---------------------------------------------------------------------
# Display available species
msigdbr::msigdbr_species()
msigdbr::msigdbr_collections()

# Subset by collection
# m_df = msigdbr(species = "Mus musculus", category = "H")
# head(m_df)
#
# m_df = msigdbr(species = "Mus musculus", category = "C7")
# head(m_df)
#
# m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# head(m_df)
#



mm_msigdb1 <- msigdbr(species = "Mus musculus") %>%
    # Combine columns
    dplyr::mutate(colon_to_underscore = str_replace_all(gs_subcat, ":", "_")) %>%
    tidyr::unite(col = "database", c(gs_cat, colon_to_underscore)) %>%
    dplyr::mutate(database = tolower(database)) %>%
    # Find strings that end with undescore, remove trailing underscore
    dplyr::mutate(database = case_when(
        str_detect(database, "_$") ~ str_sub(database, 1, -2),
        TRUE ~ database
        )) %>%
    # Add the MSigDB web links
    dplyr::mutate(description = glue("http://www.broadinstitute.org/gsea/msigdb/cards/{gs_name}"))


mm_msigdb1_small <- mm_msigdb1 %>% 
    dplyr::select(c(gene_symbol, human_ensembl_gene)) %>%
    dplyr::distinct()

# ---------------------------------------------------------------------
# Get responder nonresponder genes
# ---------------------------------------------------------------------




resp <- read_tsv(glue("{proj_dir}/input/upregulated_in_responders.txt"), col_names = FALSE) %>%
    dplyr::rename(human_ensembl_gene = "X1") %>%
    left_join(mm_msigdb1_small, by = "human_ensembl_gene") %>%
    tidyr::drop_na(gene_symbol) %>%
    dplyr::pull(gene_symbol)


nonresp <- read_tsv(glue("{proj_dir}/input/upregulated_in_nonresponders.txt"), col_names = FALSE) %>%
    dplyr::rename(human_ensembl_gene = "X1") %>%
    left_join(mm_msigdb1_small, by = "human_ensembl_gene") %>%
    tidyr::drop_na(gene_symbol) %>%
    dplyr::pull(gene_symbol)






# ---------------------------------------------------------------------
# Get custom gene lists
# ---------------------------------------------------------------------


gmt <- pathwayPCA::read_gmt(glue("{proj_dir}/input/raw_custom_gene_sets.gmt"), description = TRUE)
names(gmt$pathways) <- gmt$TERMS

# Add new gene sets
gmt$pathways[[length(gmt$pathways) + 1]] <- resp
names(gmt$pathways)[length(gmt$pathways)] <- "Zhao_2021_responders_UP"
gmt$TERMS[length(gmt$TERMS) + 1] <- "Zhao_2021_responders_UP"
gmt$description[length(gmt$description) + 1] <- "Gene upregulated in blinatumomab treated patients that responded to drug"


gmt$pathways[[length(gmt$pathways) + 1]] <- nonresp
names(gmt$pathways)[length(gmt$pathways)] <- "Zhao_2021_nonresponders_UP"
gmt$TERMS[length(gmt$TERMS) + 1] <- "Zhao_2021_nonresponders_UP"
gmt$description[length(gmt$description) + 1] <- "Gene upregulated in blinatumomab treated patients that did not respond to drug"


# Crete long table of data
gmt_list <- list()
for (i in seq_along(gmt$pathways)) {
    curr_vect <- gmt$pathways[[i]]
    curr_ont <- gmt$TERMS[[i]]
    curr_desc <- gmt$description[[i]]
    curr_database <- "custom"
    # Clean up list
    curr_vect <- curr_vect[!is.na(curr_vect)]
    curr_vect <- curr_vect[curr_vect != ""]

    gmt_list[[i]] <- tibble(ont = curr_ont, input_symbol = curr_vect, database = curr_database, description = curr_desc) 
}
gmt_tbl <- bind_rows(gmt_list) %>%
    dplyr::distinct()




# ---------------------------------------------------------------------
# Fix any gene names
# ---------------------------------------------------------------------

gene_db_url <- "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
download.file(gene_db_url, basename(gene_db_url), method = "wget")
gene_db <- read_tsv(gzfile(basename(gene_db_url)))



genes_clean <- toddr::check_gene_db(gmt_tbl$input_symbol, gene_db)



# Print data that has a fuzzy match, but still no final symbol
genes_clean_manual_review <- genes_clean %>%
    dplyr::distinct() %>%
    dplyr::filter(is.na(final_symbol) & is.na(possible_symbol) & note != "unknown")

sapply(genes_clean_manual_review$possible_df, print)



# GNLY is human only
# PZLF is human only (official name is ZBTB16)
# gzmh is the human gene, ortholog in mouse is Gzmc
# PDCD1 (human) is "Pdcd1" in mouse, which is CD279
# CDNK1a -- I think is a misspelling and should be Cdkn1a
# Rorgt is an isoform of "Rorc"
# FASLG is the human gene, Fasl is the mouse gene
# IFNAR should be Ifnar1
# IL17 should be Il17rb
# TGFBeta should be Tgfb1

genes_clean_manual_update <- genes_clean %>%
    dplyr::mutate(final_symbol = possible_symbol) %>%
    dplyr::mutate(final_symbol = case_when(
        input_symbol == "Atp5l" ~ "Atp5mg",
        input_symbol == "TNFa" ~ "Tnf",
        input_symbol == "ctla2" ~ "Ctla2a",
        input_symbol == "IL2R" ~ "Il2ra",
        input_symbol == "Itgbp7" ~ "Itgb7",
        input_symbol == "Ifgbp4" ~ "Igfbp4",
        input_symbol == "HACVR2" ~ "Havcr2",
        input_symbol == "gzmh" ~ "Gzmc",
        input_symbol == "CD279" ~ "Pdcd1",
        input_symbol == "CDNK1a" ~ "Cdkn1a",
        input_symbol == "FasLG" ~ "Fasl",
        input_symbol == "IFNAR" ~ "Ifnar1",
        input_symbol == "IL17" ~ "Il17rb",
        input_symbol == "Rorgt" ~ "Rorc",
        input_symbol == "H2afy" ~ "Macroh2a1",
        input_symbol == "Tceb2" ~ "Elob",
        input_symbol == "Ero1l" ~ "Ero1a",
        input_symbol == "Tmem2" ~ "Cemip2",
        input_symbol == "Gm440" ~ "Rundc3b",
        input_symbol == "Rnf12" ~ "Rlim",
        input_symbol == "Trim34" ~ "Trim34a",
        input_symbol == "Brp44l" ~ "Mpc1",
        input_symbol == "Gm1752" ~ "Nxpe3",
        input_symbol == "Nelf" ~ "Nsmf",
        input_symbol == "Ptk9" ~ "Twf1",
        input_symbol == "Rnf19" ~ "Rnf19a",
        input_symbol == "Sdccag1" ~ "Nemf",
        input_symbol == "Adfp" ~ "Plin2",
        input_symbol == "Bhlhb2" ~ "Bhlhe40",
        input_symbol == "Ccdc5" ~ "Haus1",
        input_symbol == "Tmem20" ~ "Slc35g1",
        input_symbol == "Tubb2c" ~ "Tubb4b",
        input_symbol == "Ctgf" ~ "Ccn2",
        input_symbol == "Mtap2" ~ "Map2",
        input_symbol == "Ptpla" ~ "Hacd1",
        input_symbol == "Lyz" ~ "Lyz1",
        input_symbol == "Ccdc11" ~ "Cfap53",
        input_symbol == "Gm129" ~ "Ciart",
        input_symbol == "Tmem8" ~ "Pgap6",
        input_symbol == "Sgol2" ~ "Sgo2a",
        input_symbol == "Whsc1" ~ "Nsd2",
        input_symbol == "Olfr60" ~ "Or13a27",
        input_symbol == "Olfr46" ~ "Or13a18",
        input_symbol == "Arntl" ~ "Bmal1",
        input_symbol == "Lrrc6" ~ "Dnaaf11",
        TRUE ~ final_symbol)) %>%
    dplyr::rename(gene_description = "description") %>%
    dplyr::select(-possible_df) %>%
    dplyr::distinct()







# Merge back with original table
genes_updated <- gmt_tbl %>%
    dplyr::left_join(genes_clean_manual_update, by = "input_symbol") %>%
    dplyr::select(input_symbol, final_symbol, everything()) %>%
    dplyr::distinct()

genes_final <- genes_updated %>%
    tidyr::drop_na(final_symbol)
genes_dropped <- genes_updated %>%
    dplyr::filter(is.na(final_symbol))
    
write_tsv(genes_final, glue("{prefix}custom_gene_sets_included_genes.txt"))
openxlsx::write.xlsx(as.data.frame(genes_final), file = glue("{prefix}custom_gene_sets_included_genes.xlsx"))
write_tsv(genes_dropped, glue("{prefix}custom_gene_sets_removed_genes.txt"))
openxlsx::write.xlsx(as.data.frame(genes_dropped), file = glue("{prefix}custom_gene_sets_removed_genes.xlsx"))




# Write a properly formatted GMT file (with differing column numbers depending
# on number of genes in the gene set). Do not allow NA as a placeholder.
gene_set_names <- unique(genes_final$ont)
con <- file(glue("{prefix}custom_gene_sets_included_genes.gmt"), "w")
for (i in seq_along(gene_set_names)) {
    curr_name <- gene_set_names[i]
    curr <- genes_final %>%
        dplyr::filter(ont == curr_name)
    curr_genes <- curr %>%
        dplyr::arrange(final_symbol) %>%
        dplyr::pull(final_symbol)
    curr_desc <- unique(curr$description)
    vect <- c(curr_name, curr_desc, curr_genes)
    line <- glue::glue_collapse(vect, sep = "\t")
    writeLines(line, con = con)
}
close(con)




gmt_custom <- genes_final %>%
    dplyr::select(
        ont,
        symbol = "final_symbol",
        database,
        description
    )


# ---------------------------------------------------------------------
# Create long custom gene set table
# ---------------------------------------------------------------------



# add in gene name for each entrez gene id
entrez <- AnnotationDbi::mapIds(org.Mm.eg.db,
                    keys = unique(gmt_custom$symbol),
                    keytype = "SYMBOL",
                    column = "ENTREZID",
                    multiVals = "first") %>%
    tibble::enframe(name = "symbol", value = "entrez_gene") %>%
    purrr::modify_at(vars(entrez_gene), as.integer)

gmt_custom2 <- gmt_custom %>%
    dplyr::left_join(entrez, by = "symbol") %>%
    dplyr::rename(gene = "entrez_gene") %>%
    dplyr::select(ont, gene, symbol, database, description) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::distinct()

dups <- gmt_custom2 %>%
    tibble::as_tibble() %>%
    dplyr::group_by_all() %>%
    dplyr::add_tally() %>%
    dplyr::filter(n > 1) %>%
    dplyr::arrange(desc(n))





# ---------------------------------------------------------------------
# Get mouse specific nomenclature for gene names
# ---------------------------------------------------------------------


# The msigdbr package provides mouse gene symbols. But some of those symbols are different than the symbols
# associated with a particular entrez id.



# Find unique entrez ids (which are mouse specific)
unique_entrez_gene <- unique(mm_msigdb1$entrez_gene)

# add in gene name for each entrez gene id
offical_mouse_gene_symbol <- AnnotationDbi::mapIds(org.Mm.eg.db,
                    keys = as.character(unique_entrez_gene),
                    keytype = "ENTREZID",
                    column = "SYMBOL",
                    multiVals = "first") %>%
    tibble::enframe(name = "entrez_gene", value = "offical_mouse_gene_symbol") %>%
    purrr::modify_at(vars(entrez_gene), as.integer)

mm_msigdb <- mm_msigdb1 %>%
    dplyr::left_join(offical_mouse_gene_symbol, by = "entrez_gene")



# Find all the symbols that differ between human and mouse, but are the same ortholog.
mm_msigdb_symbol_diff <- mm_msigdb %>%
	dplyr::select(one_of("entrez_gene", "human_gene_symbol", "offical_mouse_gene_symbol", "gene_symbol")) %>%
	filter(!tolower(gene_symbol) == tolower(offical_mouse_gene_symbol))



# However, the symbols provided by msigdbr ARE the symbols used by the GTF file when mapping -- so use
# those msigdbr provided symbols downstream.

# For example, the entrez id for Dpcd is 226162
# The offical symbol is "Dpcd", but the symbol provided by msigdbr in the "gene_symbol" column is "Gm17018"
# and "Gm17018" is found in the cellranger GTF
# grep "Dpcd" /home/lmnp/knut0297/software/modules/cellranger/5.0.1/ref_downloads/refdata-gex-mm10-2020-A/genes/genes.gtf
# NO HITS
# grep "Gm17018" /home/lmnp/knut0297/software/modules/cellranger/5.0.1/ref_downloads/refdata-gex-mm10-2020-A/genes/genes.gtf
# LOTS OF HITS against transcript ids



write_tsv(mm_msigdb, glue("{out_dir}/{prefix}mm_msigdb_long_all_columns.tsv"))




# ---------------------------------------------------------------------
# Create a data frame for clusterProfiler (for genes as Entrez Gene IDs).
# ---------------------------------------------------------------------

# Use the "gene_symbol" column, which is associated with the GTF used for mapping.



mm_msigdb_long <- mm_msigdb %>%
	dplyr::select(gs_name, entrez_gene, gene_symbol, database, description) %>%
	as.data.frame()

colnames(mm_msigdb_long) <- c("ont", "gene", "symbol", "database", "description")

mm_msigdb_long <- mm_msigdb_long %>%
    bind_rows(gmt_custom2)
# Use with clusterProfiler
# enricher(gene = genes_entrez, TERM2GENE = m_t2g, ...)



write_tsv(mm_msigdb_long, glue("{out_dir}/{prefix}mm_msigdb_long.tsv"))









# ---------------------------------------------------------------------
# Data description only
# ---------------------------------------------------------------------



# Make a "wide" database file that contains only the descriptions of the gene sets
un <- unique(mm_msigdb_long$ont)
unique_index <- !duplicated(mm_msigdb_long$ont)
uniq_ont <- mm_msigdb_long$ont[unique_index]
all(uniq_ont == un)
# TRUE

mm_msigdb_description <- mm_msigdb_long[unique_index, ]
mm_msigdb_description$gene <- NULL
mm_msigdb_description$symbol <- NULL


mm_msigdb_description <- mm_msigdb_description[order(mm_msigdb_description$database), ]
class(mm_msigdb_description$description) <- "hyperlink"

write_tsv(mm_msigdb_description, glue("{out_dir}/{prefix}mm_msigdb_description.tsv"))

wb <- openxlsx::write.xlsx(list(mm_msigdb_description = mm_msigdb_description), file = paste0(out_dir, "/", prefix, "mm_msigdb_description.xlsx"))
setColWidths(wb, sheet = 1, cols = 1, widths = 76)
setColWidths(wb, sheet = 1, cols = 2, widths = 12)
setColWidths(wb, sheet = 1, cols = 3, widths = 50)
saveWorkbook(wb, paste0(out_dir, "/", prefix, "mm_msigdb_description.xlsx"), overwrite = TRUE)






# ---------------------------------------------------------------------
# Create gene signatures objects for use with VISION
# ---------------------------------------------------------------------

# For the sigData vector, the names represent gene names and the values (1 or -1) represent the
# 'sign' of the gene.  For an unsigned signature, just use 1 for the value of every gene.
# sigData <- c(
#     Gene.A = 1, Gene.B = 1, Gene.C = 1,
#     Gene.D = -1, Gene.B = -1, Gene.C = -1
# )
# sig <- createGeneSignature(name = "Interesting Process", sigData = sigData)
# sig1 <- createGeneSignature(name = "Interesting Process", ... )
# sig2 <- createGeneSignature(name = "Another Interesting Process", ... )
#
# mySignatures <- c(sig1, sig2)
#
# vis <- Vision(data = expressionMat, signatures = mySignatures)
#
#







ont_symbol_list <- mm_msigdb_long %>%
	dplyr::select(one_of("ont", "symbol", "database")) %>%
	nest(data = one_of("symbol"))

ont_entrez_list <- mm_msigdb_long %>%
	dplyr::select(one_of("ont", "gene", "database")) %>%
	nest(data = one_of("gene"))




f <- function(i, ont_list) {
	curr_gs_name <- ont_list[i, ]$ont[[1]]
	gs <- ont_list[i, ]$data[[1]]  %>% pull(1)
	len <- length(gs)
	gs_vect <- rep(1, len)
	names(gs_vect) <- gs
	curr_database <- ont_list[i, ]$database[[1]]
	sig <- VISION::createGeneSignature(name = curr_gs_name, sigData = gs_vect, metadata = curr_database)
	return(sig)
}

ont_symbol_list_sig <- purrr::map(seq_along(ont_symbol_list$ont), f, ont_list = ont_symbol_list)
ont_entrez_list_sig <- purrr::map(seq_along(ont_entrez_list$ont), f, ont_list = ont_entrez_list)





# These lists can be used with VISION like this:
# vis <- Vision(data = expressionMat, signatures = ont_symbol_list_sig)
# Subset the signature based on database
ont_symbol_list_sig[[1]]@metaData == "c3_mir_mir_legacy"
# [1] TRUE
ont_entrez_list_sig[[1]]@metaData == "c3_mir_mir_legacy"
# [1] TRUE





#######################################################################
# Save session info
#######################################################################

# ---------------------------------------------------------------------
# Save just the VISION gene signatures
# ---------------------------------------------------------------------

saveRDS(ont_symbol_list_sig, glue("{out_dir}/{prefix}ont_symbol_list_sig.rds"))
saveRDS(ont_entrez_list_sig, glue("{out_dir}/{prefix}ont_entrez_list_sig.rds"))




# ---------------------------------------------------------------------
# Save all R objects
# ---------------------------------------------------------------------

# save.image(file = glue("{out_dir}/{prefix}r_objects.RData"))
write_tsv(data.frame(RData_objects = ls()), glue("{out_dir}/{prefix}r_objects.txt"))


# ---------------------------------------------------------------------
# Write out session info
# ---------------------------------------------------------------------


toddr::write_session_info(glue("{out_dir}/{prefix}"))
