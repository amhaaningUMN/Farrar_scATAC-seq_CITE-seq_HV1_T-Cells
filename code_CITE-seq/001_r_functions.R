#!/usr/bin/env Rscript


# ---------------------------------------------------------------------
# adt de
# ---------------------------------------------------------------------


tk_adt_de <- function(seurat_object, comparison_name, comparison_fct, lfc_threshold, description) {
    ident1 <- levels(comparison_fct)[1]
    ident2 <- levels(comparison_fct)[2]

    cells_per_group <- comparison_fct %>%
        table()

# Run DE analysis
    de <- Seurat::FindMarkers(object = seurat_object, assay = "SCT", slot = "data", ident.1 = ident1, ident.2 = ident2, min.pct = 0, test.use = "wilcox", logfc.threshold = lfc_threshold)
# Extract significant DE genes
    lfc_col <- which(grepl("log2FC", colnames(de)))
    adj_p_col <- which(grepl("p_val_adj", colnames(de)))
    de <- de[!is.na(de[, adj_p_col]), ]
    de_signif <- de[abs(de[, lfc_col]) >= lfc_threshold & de[, adj_p_col] <= 0.01, ]
    de_signif <- de_signif[order(de_signif[, lfc_col], decreasing = TRUE), ]
    de <- rownames_to_column(as.data.frame(de), var = "symbol") %>%
        as_tibble() %>%
        arrange(desc(avg_log2FC))
    de_signif <- rownames_to_column(as.data.frame(de_signif), var = "symbol") %>%
        as_tibble() 
# How many genes are significant?
    number_of_signif_genes <- dim(de_signif)[1]
    results <- list(comparison_name = comparison_name, group_a = ident1, group_b = ident2, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)



# Write this file out
    write_tsv(de, glue("{comparison_name}_de.txt"))
    write_tsv(de_signif, glue("{comparison_name}_de_signif.txt"))

    out_filename <- glue("{comparison_name}_de.xlsx")
    wb <- openxlsx::write.xlsx(list(de = de), file = out_filename, rowNames = FALSE)
    openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
    openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)	

    out_filename <- glue("{comparison_name}_de_signif.xlsx")
    wb <- openxlsx::write.xlsx(list(de_signif = de_signif), file = out_filename, rowNames = FALSE)
    openxlsx::setColWidths(wb, sheet = 1, cols = 1, widths = 20)
    openxlsx::saveWorkbook(wb, out_filename, overwrite = TRUE)	



# Volcano

    top_and_bottom_n <- 10
    my_lfc_threshold <- lfc_threshold
    my_p_val_adj_threshold <- 0.01
    data1 <- de %>%
        dplyr::mutate(lfc_threshold = if_else(abs(avg_log2FC) >= my_lfc_threshold, TRUE, FALSE)) %>%
        dplyr::mutate(p_val_adj_threshold = if_else(p_val_adj <= my_p_val_adj_threshold, TRUE, FALSE)) %>%
        dplyr::mutate(signif_threshold = if_else(lfc_threshold == TRUE & p_val_adj_threshold == TRUE, TRUE, FALSE))

    # Figure out the top and bottom significant genes to include as labels
    up_genes <- data1 %>%
        dplyr::filter(avg_log2FC >= 0 & signif_threshold == TRUE) %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        slice_head(n = top_and_bottom_n) %>%
        pull(symbol)
    dn_genes <- data1 %>%
        dplyr::filter(avg_log2FC < 0 & signif_threshold == TRUE) %>%
        dplyr::arrange(avg_log2FC) %>%
        slice_head(n = top_and_bottom_n) %>%
        pull(symbol)
    my_gene_labels <- c(up_genes, dn_genes)


    p1 <- data1 %>%
        dplyr::mutate(signif_threshold = factor(as.character(signif_threshold), levels = c("TRUE", "FALSE"))) %>%
        dplyr::mutate(gene_labels = if_else(symbol %in% my_gene_labels, symbol, "")) %>%
        ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = signif_threshold)) +
        geom_vline(xintercept = c(-my_lfc_threshold, my_lfc_threshold), linetype = "dashed", color = "grey75") +
        geom_hline(yintercept = -log10(my_p_val_adj_threshold), linetype = "dashed", color = "grey75") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        scale_color_manual(values = c("firebrick1", "gray95"), drop = FALSE) +
        ggrepel::geom_text_repel(aes(label = gene_labels), max.overlaps = 30) +
        coord_cartesian(clip = "off") +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        labs(title = glue("Signifcant DE Genes"),
            subtitle = glue("Top 10 up and down-regulated genes labeled"),
            x = "log2(fold change)",
            y = "-log10(adjusted P value)") +
        theme_minimal() +
        theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8))
    
    pdf(glue("{comparison_name}_volcano.pdf"), width = 6, height = 6)
    print(p1)
    dev.off()

# UMAP
    umap_data <- cbind(seurat_object@reductions$umap@cell.embeddings, seurat_object@meta.data) %>%
        as_tibble()

    centers <- umap_data %>%
        dplyr::select(Cluster, starts_with("UMAP_")) %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarize_all(median)

    p2 <- umap_data %>%
        # Order by facet_value
        # arrange always sorts NAs to the end, regardless of desc()
        dplyr::arrange(!is.na(comp_fct), comp_fct) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = comp_fct)) +
            geom_point(size = 0.5) +
            scale_color_manual(name = "Comparison\nGroups", na.value = "grey95", values = c("firebrick1", "dodgerblue")) +
            labs(title = "Cells in comparison", subtitle = "Cluster numbers inside plot") +
            geom_text(data = centers, mapping = aes(label = Cluster), size = 2, color = "black") +
            theme_minimal() +
            guides(color = guide_legend(override.aes = list(size = 6))) +
            theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5, size = 8),
                legend.position = "right",
                strip.background = element_blank(),
                strip.text = element_text(size = 6))

    pdf(glue("{comparison_name}_umap.pdf"), width = 6, height = 6)
    print(p2)
    dev.off()
    
    patchwork <- patchwork::wrap_plots(list(p1, p2), nrow = 1) +
        patchwork::plot_layout(guides = "collect") &
        patchwork::plot_annotation(
            title = glue("{comparison_name}"),
            subtitle = glue("Number of cells:\n{names(cells_per_group)[1]}: {cells_per_group[1]}\n{names(cells_per_group)[2]}: {cells_per_group[2]}\n({description})")) &
            theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                legend.position = "right")

    pdf(glue("{comparison_name}_umap_volcano.pdf"), width = 12, height = 6)
    print(patchwork)
    dev.off()







    return(results)

}







#######################################################################
# Enrichment Testing
#######################################################################


# Wrapper function for the core enricher function
# Takes a de list, runs the hypergeometric test for all, up, and down genes
tk_cluster_profiler <- function(de, log2fc_cutoff = 0.2, padj_cutoff = 0.05, logFC_col = NULL, padj_col = NULL, ensembl_id_col = NULL, entrez_id_col = NULL, geneset_to_gene_db_long, database_types = c("h", "c6"), plot_title = "Plot Title", run_ora = TRUE, run_go = TRUE, run_kegg = TRUE, run_gsea = TRUE) {
	
	# Create dir for text file output
	if (!dir.exists("suppl_files")) {dir.create("suppl_files")}


	# Make the de table nicer -- add links to outside databases
	tk_print_ncbi <- function(entrez_id) {
		if (str_detect(entrez_id, "[[:digit:]]") & !is.na(str_detect(entrez_id, "[[:digit:]]"))) {
			x <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", entrez_id)
		} else {
			x <- NA
		}
		return(x)
	}
	de$ncbi_gene_db <- sapply(de[[entrez_id_col]], function(x) {tk_print_ncbi(x)})
	class(de$ncbi_gene_db) <- "hyperlink"
	de$ensembl_gene_db <- sapply(de[[ensembl_id_col]], function(x) {paste0("http://www.ensembl.org/id/", x)})
	class(de$ensembl_gene_db) <- "hyperlink"


	# Get only significant/substantial genes 
	de_signif <- de[abs(de[[logFC_col]]) > log2fc_cutoff & (!is.na(de[[padj_col]]) & de[[padj_col]] < padj_cutoff), ]
	
	
	de_signif_list <- list()
	# Get all entrez ids that are known
	de_signif_list[["all"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]")]
	de_signif_list[["up"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]") & de_signif[[logFC_col]] > 0]
	de_signif_list[["dn"]] <- de_signif[[entrez_id_col]][str_detect(de_signif[[entrez_id_col]], "[[:digit:]]") & de_signif[[logFC_col]] < 0]
	

	
	# OVER-REPRESENTATION ANALYSIS
	# For all, up, down lists
	if (run_ora == TRUE) {
		enricher_list <- list()
		enricher_list_df <- list()
		enricher_list_list <- list()
		for (i in seq_along(names(de_signif_list))) {
			if (length(de_signif_list[[i]]) > 2) {
				print(names(de_signif_list)[i])
				enricher_list[[names(de_signif_list)[i]]] <- tk_enricher_core(de_signif_list[[i]], geneset_to_gene_db_long, database_types)
				enricher_list_df[[names(de_signif_list)[i]]] <- enricher_list[[names(de_signif_list)[i]]][[1]]
				enricher_list_list[[names(de_signif_list)[i]]] <- enricher_list[[names(de_signif_list)[i]]][[2]]
			} else {
				print(paste("The length of de_signif_list is less than 2 genes:", names(de_signif_list[1])))
				enricher_list_df[[names(de_signif_list)[i]]] <- data.frame()
				enricher_list_list[[names(de_signif_list)[i]]] <- list()
			}
		}
		# Clean up and sort the enrichment results
		for (j in seq_along(names(enricher_list_df))) {
			if (dim(enricher_list_df[[j]])[1] == 0) {
				enricher_list_df[[j]] <- data.frame(ID = "Zero pathways enriched.")
			} else {
				enricher_list_df[[j]] <- enricher_list_df[[j]][order(enricher_list_df[[j]]$qvalue, decreasing = FALSE), ]
				enricher_list_df[[j]]$gene_set_description <- geneset_to_gene_db_long$description[match(enricher_list_df[[j]]$ID, geneset_to_gene_db_long$ont)]
				class(enricher_list_df[[j]]$gene_set_description) <- "hyperlink"
			}
		}
		# Write out text files
		for (i in seq_along(enricher_list_df)) {
			write_tsv(enricher_list_df[[i]], path = paste0("suppl_files/msigdb_enrichment_", names(enricher_list_df)[i], ".tsv"))
		}
		# Write out enrichment excel file
		excel_list <- enricher_list_df
		names(excel_list) <- paste0("msigdb_enrichment_", names(enricher_list_df))
		wb <- openxlsx::write.xlsx(excel_list, file = paste0("msigdb_enrichment.xlsx"))
		setColWidths(wb, sheet = 1, cols = 1, widths = 60)
		setColWidths(wb, sheet = 2, cols = 1, widths = 60)
		setColWidths(wb, sheet = 3, cols = 1, widths = 60)
		saveWorkbook(wb, paste0("msigdb_enrichment.xlsx"), overwrite = TRUE)		

		# ORA Plots????
	}


	
	
	
	# GO ANALYSIS
	if (run_go == TRUE) {
		print("running GO Group Analysis")
		go_group_list <- list()
		for (i in seq_along(de_signif_list)) {
			go_group_list[[names(de_signif_list)[i]]] <- tk_go_group(geneset = de_signif_list[[i]], go_level = c(4))
		}
		# Write out multi-tab excel
		for (i in seq_along(go_group_list)) {		
			wb <- openxlsx::write.xlsx(go_group_list[[i]], file = paste0("go_groups_", names(go_group_list)[i], ".xlsx"))
			for (j in seq_along(go_group_list[[i]])) {
				setColWidths(wb, sheet = j, cols = 1, widths = 13)
				setColWidths(wb, sheet = j, cols = 2, widths = 40)
				saveWorkbook(wb, paste0("go_groups_", names(go_group_list)[i], ".xlsx"), overwrite = TRUE)
			}
		}
		# Write out text files
		for (i in seq_along(go_group_list)) {
			for (j in seq_along(go_group_list[[i]])) {
				write_tsv(go_group_list[[i]][[j]], path = paste0("suppl_files/go_groups_", names(go_group_list[[i]][j]), "_", names(go_group_list)[i], ".tsv"))
			}
		}
		print("running GO Enrichment Analysis")
		go_enrich_list <- list()
		for (i in seq_along(de_signif_list)) {
			go_enrich_list[[names(de_signif_list)[i]]] <- tk_go_enrich(geneset = de_signif_list[[i]], universe = de[[entrez_id_col]])
		}
		# Write out multi-tab excel
		for (i in seq_along(go_enrich_list)) {		
			wb <- openxlsx::write.xlsx(go_enrich_list[[i]], file = paste0("go_enrich_", names(go_enrich_list)[i], ".xlsx"))
			for (j in seq_along(go_enrich_list[[i]])) {
				setColWidths(wb, sheet = j, cols = 1, widths = 13)
				setColWidths(wb, sheet = j, cols = 2, widths = 40)
				saveWorkbook(wb, paste0("go_enrich_", names(go_enrich_list)[i], ".xlsx"), overwrite = TRUE)
			}
		}
		# Write out text files
		for (i in seq_along(go_enrich_list)) {
			for (j in seq_along(go_enrich_list[[i]])) {
				write_tsv(go_enrich_list[[i]][[j]], path = paste0("suppl_files/go_enrich_", names(go_enrich_list[[i]][j]), "_", names(go_enrich_list)[i], ".tsv"))
			}
		}
	}

	prerank <- tk_create_prerank(de, method = "lfc", logFC_col = logFC_col, padj_col = padj_col, sort = TRUE)


	# KEGG PATHWAYS ANALYSIS
	if (run_kegg == TRUE) {
		print("running KEGG")
		# Run function
		# CHOOSE BOTH UP AND DOWNREGULATED GENES
		kegg_results <- tk_kegg_enrich(geneset = de_signif_list[["all"]], prerank)
		#Excel list

		if (!is.null(kegg_results[["kegg_enrich_all"]])) {
			kegg_enrich_all <- kegg_results[["kegg_enrich_all"]]
			class(kegg_enrich_all$database_link) <- "hyperlink"
		} else {
			kegg_enrich_all <- data.frame(ID = "Zero pathways enriched.")
		}
		if (!is.null(kegg_results[["kegg_enrich_signif"]])) {
			kegg_enrich_signif <- kegg_results[["kegg_enrich_signif"]]
			class(kegg_enrich_signif$database_link) <- "hyperlink"
		} else {
			kegg_enrich_signif <- data.frame(ID = "Zero pathways significantly enriched.")
		}
        kegg_enrich_list <- list(kegg_enrich_signif = kegg_enrich_signif, kegg_enrich_all = kegg_enrich_all)
		# Write out  excel
		wb <- openxlsx::write.xlsx(kegg_enrich_list, file = paste0("kegg_enrich.xlsx"))
		setColWidths(wb, sheet = 1, cols = 1, widths = 11)
		setColWidths(wb, sheet = 1, cols = 2, widths = 40)
		setColWidths(wb, sheet = 2, cols = 1, widths = 11)
		setColWidths(wb, sheet = 2, cols = 2, widths = 40)
		saveWorkbook(wb, paste0("kegg_enrich.xlsx"), overwrite = TRUE)
		# Write out text files
		for (i in seq_along(kegg_enrich_list)) {
			write_tsv(kegg_enrich_list[[i]], path = paste0("suppl_files/", names(kegg_enrich_list)[i], ".tsv"))
		}
	}
	
	
	
	# GSEA
	if (run_gsea == TRUE) {
		print("running GSEA")
		enricher_gsea_output <- tk_enricher_gsea(sorted_prerank_df = prerank, geneset_to_gene_db_long = geneset_to_gene_db_long, plot_title = plot_title, qvalue_threshold = padj_cutoff, database_types = database_types)
		#universe = sorted_prerank_df[[entrez_id_col]]
	}
}









# ---------------------------------------------------------------------
# Run enricher hypergeometric test
# ---------------------------------------------------------------------

tk_enricher_core <- function(geneset, geneset_to_gene_db_long, database_types = levels(geneset_to_gene_db_long$database)) {
	# all databases
	# database_types <- levels(geneset_to_gene_db_long$database)
	# OR MANUAL EDIT TO CHOOSE FOR YOUR FAVORITE DATABASES
	# database_types <- c("c2_cgp", "c2_cp_reactome", "c2_cp_biocarta", "c2_cp_kegg", "h", "c6")
	# database_types <- levels(geneset_to_gene_db_long$database)[grepl("_all", levels(geneset_to_gene_db_long$database))]
	hypergeo_msigdb_df <- data.frame()
	hypergeo_msigdb_list <- list()
	for (i in seq_along(database_types)) {
		# get only one database type at a time
		# Print database name
		print(database_types[i])
		geneset_to_gene_db_long_curr <- filter(geneset_to_gene_db_long, database == database_types[i])
		# Print the number of gene sets in database
		print(length(levels(factor(geneset_to_gene_db_long_curr$ont))))
		hypergeo_msigdb_curr <- enricher(geneset, TERM2GENE = geneset_to_gene_db_long_curr, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)
		if (!is.null(hypergeo_msigdb_curr)) {
			hypergeo_msigdb_curr <- setReadable(hypergeo_msigdb_curr, org.Mm.eg.db, keyType = "ENTREZID")
		}
		hypergeo_msigdb_list[[database_types[i]]] <- hypergeo_msigdb_curr
		# Clean up results
		if (!is.null(hypergeo_msigdb_curr)) {
			print("appending all list")
			hypergeo_msigdb_curr_df <- hypergeo_msigdb_curr@result
			hypergeo_msigdb_curr_df <- hypergeo_msigdb_curr_df[hypergeo_msigdb_curr_df$Count > 1, , drop = FALSE]
			if (dim(hypergeo_msigdb_curr_df)[1] == 0) {
				# Zero enriched genes after filtering
				print("NOT appending to list") 
			} else {
				hypergeo_msigdb_curr_df$Description <- NULL
				hypergeo_msigdb_curr_df$gene_set_database <- database_types[i]
				hypergeo_msigdb_df <- rbind(hypergeo_msigdb_df, hypergeo_msigdb_curr_df)			
			}
		} else { 
			print("NOT appending to list") 
		}
		# NOTE:
		# if you get some error like this:
		# --> No gene can be mapped....
		# --> Expected input gene ID: 235293,17125,18710,225326,225010,104709
		# --> return NULL...

		# All your genes in your geneset of interest (upregulated) must be present in the gene sets of interest. If you select
		# a list of gene sets, find the total uniqe set of genes available, the genes you're interested in testing ALL must be present
		# in that gene set list.
		# up %in% unique(geneset_to_gene_db_long_curr$gene)
	}
	core_list <- list(hypergeo_msigdb_df = hypergeo_msigdb_df, hypergeo_msigdb_list = hypergeo_msigdb_list)
	return(core_list)
}





# ---------------------------------------------------------------------
# Gene Ontology - Group info
# ---------------------------------------------------------------------

tk_go_group <- function(geneset, go_level = c(2, 3, 4)) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		go_list <- list()
		return(go_list)
	} else {
		ontology <- c("CC", "MF", "BP")
		go_list <- list()
		for (i in seq_along(ontology)) {
			for (j in seq_along(go_level)) {
			name <- paste0("go_", ontology[i], "_", go_level[j])
			go <- groupGO(gene = geneset, OrgDb = org.Mm.eg.db, ont = ontology[i], level = go_level[j], readable = TRUE)
			# Sort results
			go_sorted <- go@result[order(go@result$Count, decreasing = TRUE), ]
			go_list[[name]] <- go_sorted
			}
		}
		return(go_list)
	}
}


# ---------------------------------------------------------------------
# Gene Ontology - Group enrichment
# ---------------------------------------------------------------------

tk_go_enrich <- function(geneset, universe) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		go_list <- list()
		return(go_list)
	} else {
		ontology <- c("CC", "MF", "BP")
		go_list <- list()
		for (i in seq_along(ontology)) {
			name <- paste0("go_", ontology[i])
			go <- enrichGO(gene = geneset, OrgDb = org.Mm.eg.db, ont = ontology[i], readable = TRUE, universe = universe, pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
			# Sort the results
			go_sorted <- go@result[go@result$qvalue <= 0.05, ]
			go_sorted <- go_sorted[order(go_sorted$qvalue, decreasing = FALSE), ]
			go_list[[name]] <- go_sorted
		}
		return(go_list)
	}
}

###### bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
###### enrichMap(bp)






# ---------------------------------------------------------------------
# KEGG
# ---------------------------------------------------------------------


tk_kegg_enrich <- function(geneset, prerank) {
	if (length(geneset) == 0) {
		print("geneset is empty")
		kegg_enrich_list <- list()
		return(kegg_enrich_list)
	} else {
		# Make sure gene set is only entrez_ids, (digits)
		geneset <- geneset[str_detect(geneset, "\\d")]
		# Run KEGG enrichment test
		kegg_enrich <- enrichKEGG(gene = geneset,
							organism = "hsa",
							keyType = "kegg",
       						pvalueCutoff = 0, 
       						pAdjustMethod = "BH",
       						minGSSize = 10, 
       						maxGSSize = 500, 
       						qvalueCutoff = 0)
       	if (is.null(kegg_enrich)) {
       		# There are no enrichment results!
       		kegg_enrich_all <- NULL
       		kegg_enrich_signif <- NULL
       	} else {
			kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
			kegg_enrich_all <- kegg_enrich@result
			kegg_enrich_all$database_link <- paste0("https://www.genome.jp/dbget-bin/www_bget?", kegg_enrich_all$ID)
			kegg_enrich_signif <- kegg_enrich_all[kegg_enrich_all$p.adjust < 0.05 & kegg_enrich_all$qvalue < 0.2, ]
			# Make nice names for filenames
			description_clean <- str_replace_all(kegg_enrich_signif$Description, " ", "_")
			# Create entrez_id named vector of logFC values for all genes
			gene_list_tb <- prerank %>%
							dplyr::select(entrez_id, logFC) %>%
							filter(str_detect(entrez_id, "\\d")) %>%
							arrange(desc(logFC))
			gene_list <- gene_list_tb[["logFC"]]
			names(gene_list) <- gene_list_tb[["entrez_id"]]
			# Print out all the figs
			orig_wd <- getwd()
			results_out_dir <- glue("{orig_wd}/kegg_enrich")
			if (!dir.exists(results_out_dir)) {dir.create(results_out_dir, recursive = TRUE)}
			setwd(results_out_dir)
			# TEMP TEST
			options(bitmapType = 'cairo')
			for (i in seq_along(kegg_enrich_signif$ID)) {
				# Print KEGG plots
				p1 <- pathview(gene.data  = gene_list,
							pathway.id = kegg_enrich_signif$ID[i],
							species    = "hsa",
							low = list(gene = "blue", cpd = "green"),
							mid = list(gene = "white", cpd = "white"),
							high = list(gene = "firebrick1", cpd = "yellow"),
							out.suffix = description_clean[i])
				# Print pathview PDFs
				p2 <- pathview(gene.data  = gene_list,
							pathway.id = kegg_enrich_signif$ID[i],
							species    = "hsa",
							low = list(gene = "blue", cpd = "green"),
							mid = list(gene = "white", cpd = "white"),
							high = list(gene = "firebrick1", cpd = "yellow"),
							out.suffix = description_clean[i],
							kegg.native = FALSE)
			}
			# Delete intermediate files
			unlink("*.xml")
			unlink("*[[:digit:]].png")
			setwd(orig_wd)
		}	
	}
	return(list(kegg_enrich_all = kegg_enrich_all, kegg_enrich_signif = kegg_enrich_signif))
}






# ---------------------------------------------------------------------
# Create pre-ranked gene list for GSEA
# ---------------------------------------------------------------------

tk_create_prerank <- function(tibble, method = NULL, logFC_col = NULL, padj_col = NULL, sort = TRUE) {
	# method must be one of: "pval_sign" or "lfc"
	# Check input parameters
	if (is.null(logFC_col)) {
		print("You need to provide the column index for log fold change values")
		stop()
	}
	if (is.null(padj_col)) {
		print("You need to provide the column index for adjusted p-values")
		stop()
	}
	# Create prerank table
	if (method == "pval_sign") {
		tibble <- add_column(tibble, "log2FoldChange_sign" = sign(tibble[[logFC_col]]))
		tibble <- add_column(tibble, "minuslog10pval" = -log10(tibble[[padj_col]]))
		# If any are Inf, set to very high value (e.g. 300)
		tibble$minuslog10pval[which(tibble$minuslog10pval == "Inf")] <- 300
		# Create the prerank_metric: -log10pval * fold change sign
		tibble <- add_column(tibble, "prerank_metric" = (tibble$minuslog10pval * tibble$log2FoldChange_sign))
		# Remove genes that have p values with "NA"
		tibble <- tibble[!(is.na(tibble$minuslog10pval)), ]
		if (sort == TRUE) {
			tibble <- tibble[order(tibble$prerank_metric, decreasing = TRUE), ]
		}
		return(tibble)
	} else if (method == "lfc") {
		tibble <- tibble %>% 
					add_column(prerank_metric = tibble[[logFC_col]]) %>%
					drop_na(all_of(padj_col)) %>%
					filter(entrez_id != "unknown")
		if (sort == TRUE) {
			tibble <- tibble %>% 
					arrange(desc(tibble[[logFC_col]]))
		}	
		return(tibble)
	} else {
		print("You need to provide a method argument that is either pval_sign or lfc")
		stop()
	}
}





# ---------------------------------------------------------------------
# GSEA function
# ---------------------------------------------------------------------

tk_enricher_gsea <- function(sorted_prerank_df, geneset_to_gene_db_long, plot_title = "Plot Title", qvalue_threshold = 0.05, database_types = levels(geneset_to_gene_db_long$database)) {

	prerank_gene_list <- sorted_prerank_df$prerank_metric
	names(prerank_gene_list) <- as.character(sorted_prerank_df$entrez_id)

	# DATABASES
	# all:
	# database_types <- levels(geneset_to_gene_db_long$database)
	# OR MANUAL EDIT TO CHOOSE FOR YOUR FAVORITE DATABASES
	# database_types <- levels(geneset_to_gene_db_long$database)[grepl("_all", levels(geneset_to_gene_db_long$database))]
	# database_types <- c("c2_cgp", "c2_cp_reactome", "c2_cp_biocarta", "c2_cp_kegg", "h", "c6")
	
	
	gsea_msigdb_df <- data.frame()
	gsea_msigdb_signif_df <- data.frame()
	gsea_msigdb_list <- list()
	for (i in seq_along(database_types)) {
		# get only one database type at a time
		# Print database name
		print(database_types[i])
		geneset_to_gene_db_long_curr <- dplyr::filter(geneset_to_gene_db_long, database == database_types[i])
	    print(as_tibble(geneset_to_gene_db_long_curr))	
        geneset_to_gene_db_long_curr %>%
            group_by(ont) %>%
            count()

        # Print the number of gene sets in database
		print(length(levels(factor(geneset_to_gene_db_long_curr$ont))))
		gsea_msigdb_curr <- clusterProfiler::GSEA(prerank_gene_list, TERM2GENE = geneset_to_gene_db_long_curr, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH")
		if (length(gsea_msigdb_curr$ID) > 0) {
			gsea_msigdb_curr <- setReadable(gsea_msigdb_curr, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
			gsea_msigdb_curr@result$gene_set_database <- database_types[i]
			gsea_msigdb_curr@result$gene_set_description <- geneset_to_gene_db_long$description[match(gsea_msigdb_curr$ID, geneset_to_gene_db_long$ont)]
			# Sometimes with only 1 signif geneset, the qvalues are logical with NA. This causes problems when plotting
			if (!is.logical(gsea_msigdb_curr@result$qvalue)) {
				gsea_msigdb_curr_df <- gsea_msigdb_curr@result %>%
                    dplyr::arrange(qvalue)
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr@result %>%
                    dplyr::filter(qvalue <= qvalue_threshold) %>%
                    dplyr::arrange(qvalue)
			} else {
				gsea_msigdb_curr_df <- gsea_msigdb_curr@result %>%
                    dplyr::arrange(p.adjust)
				# Else, filter on adjusted p values (which are similar)
				gsea_msigdb_curr_signif_df <- gsea_msigdb_curr@result %>%
                    dplyr::filter(p.adjust <= qvalue_threshold) %>%
                    dplyr::arrange(p.adjust)
			}
			# Combine results from curr df into full df
			gsea_msigdb_df <- rbind(gsea_msigdb_df, gsea_msigdb_curr_df)
			gsea_msigdb_signif_df <- rbind(gsea_msigdb_signif_df, gsea_msigdb_curr_signif_df)
			# Combine results from curr database to all others
			gsea_msigdb_list[[database_types[i]]] <- gsea_msigdb_curr
			# Make individual plots for all signif gene sets
			for (j in seq_along(gsea_msigdb_curr_signif_df$ID)) {
				if (gsea_msigdb_curr_signif_df$NES[j] > 0) {
					nes_direction <- "pos"
				} else {
					nes_direction <- "neg"
				}
				if (!dir.exists(paste0("msigdb_gsea/", database_types[i]))) {dir.create(paste0("msigdb_gsea/", database_types[i]), recursive = TRUE)}
				pdf(paste0("msigdb_gsea/", database_types[i], "/", str_pad(j, 4, pad = "0"), "_", nes_direction, "_", gsea_msigdb_curr_signif_df$ID[j], ".pdf"))
				print(
					enrichplot::gseaplot2(gsea_msigdb_curr, geneSetID = gsea_msigdb_curr_signif_df$ID[j], title = paste0(plot_title, "\n", gsea_msigdb_curr_signif_df$ID[j], "\n", "NES = ", round(gsea_msigdb_curr_signif_df$NES[j], 4), "\n", "qvalue = ", round(gsea_msigdb_curr_signif_df$qvalue[j], 4)))
				)
				dev.off()
			}
		}
	}
	# Clean up and sort data frame, export to Excel file
	gsea_msigdb_df <- gsea_msigdb_df %>%
        dplyr::arrange(qvalue)
	gsea_msigdb_df$Description <- NULL
	class(gsea_msigdb_df$gene_set_description) <- "hyperlink"

	gsea_msigdb_signif_df <- gsea_msigdb_signif_df %>%
        dplyr::arrange(qvalue)
	gsea_msigdb_signif_df$Description <- NULL
	class(gsea_msigdb_signif_df$gene_set_description) <- "hyperlink"
	
	prerank_excel <- as.data.frame(sorted_prerank_df)
	class(prerank_excel$ncbi_gene_db) <- "hyperlink"
	class(prerank_excel$ensembl_gene_db) <- "hyperlink"
	prerank_excel_short <- prerank_excel[, c("symbol", "prerank_metric")]
		
		
	excel_list <- list(gsea_msigdb_signif_df, gsea_msigdb_df, prerank_excel_short, prerank_excel)
	names(excel_list) <- c("gsea_signif", "gsea_all", "preranked_gene_order", "preranked_all_columns")
	wb <- openxlsx::write.xlsx(excel_list, file = paste0("msigdb_gsea.xlsx"))
	setColWidths(wb, sheet = 1, cols = 1, widths = 60)
	saveWorkbook(wb, paste0("msigdb_gsea.xlsx"), overwrite = TRUE)
	
	# Create dir for text file output
	if (!dir.exists("suppl_files")) {dir.create("suppl_files")}
	write_tsv(gsea_msigdb_df, "suppl_files/gsea_msigdb.tsv")
	write_tsv(gsea_msigdb_signif_df, "suppl_files/gsea_msigdb_signif.tsv")
	write_tsv(prerank_excel_short, "suppl_files/prerank_less_cols.tsv")
	write_tsv(prerank_excel, "suppl_files/prerank.tsv")
	return(gsea_msigdb_list)
}




#######################################################################
# Special enrichment testing
#######################################################################



tk_parallel_cluster_profiler <- function(k, all_de_data, gene_info, mm_msigdb_long, out_dir, curr_sample_name, curr_res_name, curr_de_comparison) {
    curr_de_comparison <- names(all_de_data)[k]
    print(curr_de_comparison)

    if (!dir.exists(glue("{out_dir}/{curr_sample_name}/{curr_res_name}/{curr_de_comparison}"))) {dir.create(glue("{out_dir}/{curr_sample_name}/{curr_res_name}/{curr_de_comparison}"), recursive = TRUE)}
    setwd(glue("{out_dir}/{curr_sample_name}/{curr_res_name}/{curr_de_comparison}"))

    curr_de_data <- all_de_data[k]

    # Cluster profiler parameters
    de <- curr_de_data[[1]]$de %>%
        # Add ensembl and entrez ids
        dplyr::left_join(dplyr::select(.data = gene_info, ensembl_id, entrez_id, symbol), by = "symbol")

    log2fc_cutoff <- 0.25
    padj_cutoff <- 0.05
    logFC_col <- 3
    padj_col <- 6
    ensembl_id_col <- 7
    entrez_id_col <- 8
    geneset_to_gene_db_long <- mm_msigdb_long
    database_types <- c("h", "c2_cgp", "c2_cp", "c7")
    #database_types <- levels(mm_msigdb_long$database)
    plot_title <- paste(curr_res_name, curr_de_comparison, sep = "")
    run_ora <- FALSE
    run_go <- FALSE
    run_kegg <- FALSE
    run_gsea <- TRUE


    results <- tk_cluster_profiler(de = de, log2fc_cutoff = log2fc_cutoff, padj_cutoff = padj_cutoff, logFC_col = logFC_col, padj_col = padj_col, ensembl_id_col = ensembl_id_col, entrez_id_col = entrez_id_col, geneset_to_gene_db_long = geneset_to_gene_db_long, database_types = database_types, plot_title = plot_title, run_ora = run_ora, run_go = run_go, run_kegg = run_kegg, run_gsea = run_gsea)
    names(results) <- curr_de_comparison
    return(results)
}


