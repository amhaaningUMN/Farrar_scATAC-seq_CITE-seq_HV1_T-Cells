
# seurat_object = curr_seurat_object
# seurat_object_name = curr_seurat_object@project.name
# lfc_threshold = 0.25
# num_genes_diff_between_clusters_threshold = 0



tk_citeseq_cluster_merge2 <- function(seurat_object, seurat_object_name, lfc_threshold, num_genes_diff_between_clusters_threshold) {
	print(seurat_object_name)
	# Get the matrix and labels from the clustering results
	cluster_cell_class <- seurat_object@meta.data[, grepl("SCT_snn_res", colnames(seurat_object@meta.data)), drop = FALSE]
	cluster_res_name <- str_replace_all(colnames(cluster_cell_class), fixed("res."), "res_")
	cell_ids <- rownames(cluster_cell_class)
	orig_assay <- DefaultAssay(object = seurat_object)
	orig_wd <- getwd()
	return_list <- list()
	# For each clustering resolution
	for (j in seq_along(cluster_res_name))  {
	    numb_of_samples <- length(levels(factor(seurat_object@meta.data$orig.ident)))
		# Create dir for all results
		out_dir <- cluster_res_name[j]
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}"))
		print(cluster_res_name[j])
		# Factor the current clustering resolution cell classes
		class_orig <- factor(cluster_cell_class[, j])
		names(class_orig) <- cell_ids
		# Re-assign the class names to @ident slot with current cell classes
		seurat_object@active.ident <- class_orig
		
		
		# Figure out the right assay to use
		DefaultAssay(seurat_object) <- "SCT"
		print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
		
		
		
		
		print("Starting DE set1, PAIRWISE DE and MERGE CLUSTERS BASED ON SIGNIF DE GENES")
		# set1: PAIRWISE DE and MERGE CLUSTERS BASED ON SIGNIF DE GENES
		# Use all cells in a cluster, regardless of their orig.ident (i.e. sample)
		tk_merge_clusters_by_de_results <- tk_merge_clusters_by_de(seurat_object, cluster_res_name = cluster_res_name[j], lfc_threshold, num_genes_diff_between_clusters_threshold)
		# Extract results
		merge_data_list <- tk_merge_clusters_by_de_results[["merge_data_list"]]
		iteration_data_list <- tk_merge_clusters_by_de_results[["iteration_data_list"]]
		final_classes <- tk_merge_clusters_by_de_results[["final_classes"]]
	
	
		# MERGING HAS STOPPED, The clusters are now finalized
		# Re-assign the identity factor
		seurat_object@active.ident <- final_classes		
		curr_number_of_clusters_final <- length(levels(factor(seurat_object@active.ident)))		






		# RE-WRITE CLUSTER NAMES (e.g. convert the 1.2.3.4 names to 1) 
		# Rename the final clusters to be regular integers, for plotting and file output
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_orig")]] <- as.character(final_classes)
		# Rename the final clusters to be regular integers, for plotting and file output
		# Make them zero-based index and sorted by size
		nclasses <- length(levels(final_classes))
		cluster_cell_class_factor_final <- factor(forcats::fct_infreq(final_classes), labels = c(1:nclasses - 1))
		names(cluster_cell_class_factor_final) <- names(final_classes)
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final")]] <- cluster_cell_class_factor_final	
		# Map the ugly names to final names
		mapping_full <- data.frame(de_merge_orig = final_classes, de_merge_final = cluster_cell_class_factor_final)
		mapping <- unique(mapping_full)
		mapping <- mapping %>% dplyr::arrange(de_merge_final)
		write_tsv(mapping, glue("{orig_wd}/{out_dir}/mapping_cluster_labels.txt"), col_names = TRUE)
		rownames(mapping) <- NULL
		# Export data here
		iteration_data_list <- iteration_data_list[!sapply(iteration_data_list, is.null)]
		final_cluster_labels <- data.frame(cell_id = names(cluster_cell_class_factor_final), clust_orig = class_orig, clust_de_merge_orig = final_classes, clust_de_merge_final = cluster_cell_class_factor_final)
		write_tsv(final_cluster_labels, glue("{orig_wd}/{out_dir}/final_cluster_labels.txt"), col_names = TRUE)
        seurat_object@meta.data$orig.ident <- factor(as.character(seurat_object@meta.data$orig.ident), levels = gtools::mixedsort(as.character(seurat_object@meta.data$orig.ident)), labels = gtools::mixedsort(as.character(seurat_object@meta.data$orig.ident)))


		
		# UPDATE MERGE DATA LIST
		# Get only the final clusters info
		de_merge_orig_names <- names(merge_data_list)
		de_merge_orig_names2 <- str_replace(de_merge_orig_names, "cluster_", "")
		de_merge_orig_names3 <- str_split(de_merge_orig_names2, "_vs_", simplify = FALSE)
		de_data_list <- list()
		for (i in seq_along(de_merge_orig_names3)) {
			if (all(de_merge_orig_names3[[i]] %in% mapping$de_merge_orig)) {
				curr_idx <- length(de_data_list) + 1
				de_data_list[[curr_idx]] <- merge_data_list[[i]]
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "name"] <- "name_orig"
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "cluster_a"] <- "cluster_a_orig"
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "cluster_b"] <- "cluster_b_orig"
				de_data_list[[curr_idx]]$cluster_a <- as.character(mapping$de_merge_final[mapping$de_merge_orig %in% merge_data_list[[i]]$cluster_a])
				de_data_list[[curr_idx]]$cluster_b <- as.character(mapping$de_merge_final[mapping$de_merge_orig %in% merge_data_list[[i]]$cluster_b])
				de_data_list[[curr_idx]]$name <- paste0("cluster_", de_data_list[[curr_idx]]$cluster_a, "_vs_", de_data_list[[curr_idx]]$cluster_b)
				de_data_list[[curr_idx]]$de_type <- "de_pairwise"
				names(de_data_list)[curr_idx] <- de_data_list[[curr_idx]]$name
			}
		}
		#lapply(de_data_list, function(x) {print(paste0("final: ", x$name)); print(paste0("orig: ", x$name_orig))})
		
		
		
		
		print("Starting DE set2, DE ONE VS ALL")
		# set2: DE ONE VS ALL
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# testing for cluster "a" vs all other cells (all other clusters combined)
		# Use all cells in a cluster, regardless of their orig.ident (i.e. sample)
		# Add these results to the de_data_list
		print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
		m_vect <- 1:curr_number_of_clusters_final
		out <- mclapply(X = m_vect, FUN = tk_parallel_de_one_vs_all, seurat_object = seurat_object, lfc_threshold = lfc_threshold,
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		# name the list elements
		names(out) <- base::sapply(out, function(x) x[["name"]])
		# update de_data_list	
		# append lists
		length_de_data_list <- length(de_data_list)
		for (i in seq_along(out)) {
			out[[i]]$de_type <- "de_one_vs_all"
			idx <- length_de_data_list + i
			de_data_list[[idx]] <- out[[i]]
			names(de_data_list)[idx] <- out[[i]]$name
		}
		
		
		
		
		
        if (numb_of_samples > 1) {
            print("Multiple samples in analysis. Starting across sample analyses.")


            print("Starting DE set3, PAIRWISE CONSERVED DE LISTS")
            # set3: PAIRWISE CONSERVED DE LISTS
            seurat_object@active.ident <- cluster_cell_class_factor_final
            # Generate pairwise cluster numbering for all cluster comparisons
            combos <- combn(levels(cluster_cell_class_factor_final), 2)
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Run all pairwise cluster comparisons, find DE genes
            k_vect <- 1:dim(combos)[2]
            # Help with parallel lapply:
            # https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
            curr_cluster_res_name <- glue("{cluster_res_name[j]}_de_merge_final")
            out <- mclapply(X = k_vect, FUN = tk_int_parallel_de_pairwise, seurat_object = seurat_object, curr_clusters = cluster_cell_class_factor_final, curr_number_of_clusters = curr_number_of_clusters_final, combos = combos, lfc_threshold = lfc_threshold, de_data_list = de_data_list, meta_data_col = curr_cluster_res_name,
                mc.cores = min(16, as.integer(system("echo $THREADS", intern = TRUE))))
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_pairwise_conserved"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		
        
            print("Starting DE set4, ONE VS ALL CONSERVED DE LISTS")
            # set4: ONE VS ALL CONSERVED DE LISTS
            seurat_object@active.ident <- cluster_cell_class_factor_final
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Run all pairwise cluster comparisons, find DE genes
            m_vect <- 1:curr_number_of_clusters_final
            # Help with parallel lapply:
            # https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
            curr_cluster_res_name <- glue("{cluster_res_name[j]}_de_merge_final")
            out <- mclapply(X = m_vect, FUN = tk_int_parallel_de_one_vs_all, seurat_object = seurat_object, lfc_threshold = lfc_threshold, meta_data_col = curr_cluster_res_name,
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_one_vs_all_conserved"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		


        




            print("Starting DE set5, DE LISTS ACROSS SAMPLES -- WITHIN SAME CLUSTER")
            # set5: DE LISTS ACROSS SAMPLES -- WITHIN SAME CLUSTER
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Create a new identity labeling scheme that combines the original sample name and final cluster name
            col_cluster <- colnames(seurat_object@meta.data)[str_detect(colnames(seurat_object@meta.data), "_de_merge_final$")]
            col_sample <- "orig.ident"
            new_levels <- tidyr::crossing(seurat_object@meta.data[[col_cluster]], seurat_object@meta.data[[col_sample]]) %>%
                dplyr::rename(col_cluster = `seurat_object@meta.data[[col_cluster]]`, col_sample = `seurat_object@meta.data[[col_sample]]`) %>%
                dplyr::mutate(new_levels = as.character(glue("{col_cluster}.{col_sample}")))
            cluster_sample <- paste0(seurat_object@meta.data[[col_cluster]], ".", seurat_object@meta.data[[col_sample]])
            cluster_sample <- factor(cluster_sample, levels = new_levels$new_levels)
            cluster_sample <- fct_relevel(cluster_sample, gtools::mixedsort)
            seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final_cluster_sample")]] <- cluster_sample
            names(cluster_sample) <- rownames(seurat_object@meta.data)
            seurat_object@active.ident <- cluster_sample
            cluster_sample_tbl <- tibble(cluster = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[1]),
                sample = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[2]))
            uniq_clust <- unique(cluster_sample_tbl$cluster)
            cluster_sample_list <- list()
            for (i in seq_along(uniq_clust)) {
                curr <- cluster_sample_tbl %>% filter(cluster == uniq_clust[i])
                cluster_sample_list[[i]] <- list(cluster = uniq_clust[i], sample = curr$sample)
            }
            t_vect <- 1:curr_number_of_clusters_final
            out1 <- mclapply(X = t_vect, FUN = tk_int_parallel_de_across_samples, seurat_object = seurat_object, cluster_sample_list = cluster_sample_list, lfc_threshold = lfc_threshold,
                mc.cores = min(16, as.integer(system("echo $THREADS", intern = TRUE))))
            # Flatten one level of the lists
            out <- unlist(out1, recursive = FALSE)
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_across_samples"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		












            print("Starting DE set6, DE LISTS ACROSS SAMPLES and ACROSS CLUSTERS")
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Create a new identity labeling scheme that combines the original sample name and final cluster name



            # col_cluster <- colnames(seurat_object@meta.data)[str_detect(colnames(seurat_object@meta.data), "_de_merge_final$")]
            # col_sample <- "orig.ident"
            # new_levels <- tidyr::crossing(seurat_object@meta.data[[col_cluster]], seurat_object@meta.data[[col_sample]]) %>%
            #     dplyr::rename(col_cluster = `seurat_object@meta.data[[col_cluster]]`, col_sample = `seurat_object@meta.data[[col_sample]]`) %>%
            #     dplyr::mutate(new_levels = as.character(glue("{col_cluster}.{col_sample}")))
            #
            #
            # cluster_sample <- paste0(seurat_object@meta.data[[col_cluster]], ".", seurat_object@meta.data[[col_sample]])
            # cluster_sample <- factor(cluster_sample, levels = new_levels$new_levels)
            # cluster_sample <- fct_relevel(cluster_sample, gtools::mixedsort)
            # seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final_cluster_sample")]] <- cluster_sample
            # names(cluster_sample) <- rownames(seurat_object@meta.data)
            # seurat_object@active.ident <- cluster_sample
            # cluster_sample_tbl <- tibble(cluster = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[1]),
            #     sample = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[2]))
            # uniq_clust <- unique(cluster_sample_tbl$cluster)
            # cluster_sample_list <- list()
            # for (i in seq_along(uniq_clust)) {
            #     curr <- cluster_sample_tbl %>% filter(cluster == uniq_clust[i])
            #     cluster_sample_list[[i]] <- list(cluster = uniq_clust[i], sample = curr$sample)
            # }

            cluster_sample_matrix <- combn(new_levels$new_levels, 2)

            t_vect <- seq_len(ncol(cluster_sample_matrix))
            out1 <- mclapply(X = t_vect, FUN = tk_int_parallel_de_across_samples_and_clusters, seurat_object = seurat_object, cluster_sample_matrix = cluster_sample_matrix, lfc_threshold = lfc_threshold,
                mc.cores = min(44, as.integer(system("echo $THREADS", intern = TRUE))))
            # Flatten one level of the lists
            out <- unlist(out1, recursive = FALSE)
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_across_samples_and_clusters"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		
        
        } else if (numb_of_samples == 1) {
            print("One sample in analysis. Skipping across sample analyses.")
        }










		print("Making UMAP plots")
		# UMAP PLOTS
		# Plot original UMAP
		umap_plot_list <- list()
		# Re-assign the identity factor
		# Make this a nicely sorted factor of integers, zero-based
		class_orig_factor <- factor(as.numeric(class_orig) - 1)
		names(class_orig_factor) <- names(class_orig)
		seurat_object@active.ident <- class_orig_factor
		plot_data <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
		print("Using seurat_object@reductions$umap@cell.embeddings for UMAP")
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[1]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nInitial pre-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[1] <- paste0(cluster_res_name[j], "_before_de_merge")
		# Plot final UMAP
		# Re-assign the identity factor
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# get centers for labels
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[2]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nFinal post-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[2] <- paste0(cluster_res_name[j], "_final")
		pdf("umap_final_de_merge.pdf", width = 6, height = 6)
		print(umap_plot_list[[2]])
		dev.off()
		pdf("umap_original.pdf", width = 6, height = 6)
		print(umap_plot_list[[1]])
		dev.off()
		# Get the final iteration counter for plotting
		iteration_counter <- length(names(iteration_data_list))
		iteration_counter_final <- iteration_counter - 1
		pdf("umap_original_and_final_de_merge.pdf", width = 10, height = 6)
		grid.arrange(
			grobs = umap_plot_list,
			nrow = 1,
			top = paste0("UMAP plots, Pre- and Post-DE merging.\n", iteration_counter_final , " cluster merging steps required."))
		dev.off()
		
		
		
# 		# Print the names of de_data_list
# 		lapply(de_data_list, function(x) {print(x$name)})
# 		lapply(de_data_list, function(x) {names(x)})
# 		lapply(de_data_list, function(x) {print(head(x$de_signif))})
# 		# Not tibbles
# 		cluster_5_vs_all
# 		cluster_4_vs_5
		
		


        print("Starting most SCT variable gene plots")
        if (!is.null(seurat_object@assays$SCT)) {
            # Requires SCT assay
            most_variable_features_sorted_n50 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 50) %>%
                base::rownames(.) 
            most_variable_features_sorted_n100 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 100) %>%
                base::rownames(.) 
            most_variable_features_sorted_n200 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 200) %>%
                base::rownames(.)
            # Heatmap (top most variable genes across all cells)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"))) {dir.create(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"))
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n50, filename_prefix = "heatmap_top50_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 50 most variable genes across all cells\nEach row independently scaled")
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n100, filename_prefix = "heatmap_top100_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 100 most variable genes across all cells\nEach row independently scaled")
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n200, filename_prefix = "heatmap_top200_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 200 most variable genes across all cells\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
        }










		print("Starting set1 plots: de_pairwise")
		# Set1: GET PAIRWISE PLOTS
		r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_pairwise")
		out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_pairwise_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		names(out) <- base::sapply(out, function(x) x[["name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])       
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		tk_int_parallel_de_pairwise_plots_results <- out
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# The openxlsx package is using an old zip method, so you might get warnings (but should still work)
		# Note: zip::zip() is deprecated, please use zip::zipr() instead
		# https://github.com/awalker89/openxlsx/issues/454
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_pairwise.xlsx")
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_pairwise_signif.xlsx")
        # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
        comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
        comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_pairwise/umap_gene_expression_flattened"))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_pairwise/umap_gene_expression_flattened"))
		# Heatmap
		# Get only top DE genes from each comparison
		top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
        top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
        # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
        top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
        top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
		# Heatmap (top up and and down 5 genes)
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))    
		# Heatmap (top up or down genes, farthest from the mean)       
 		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))             
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_pairwise_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))
		# Dotplots
		seurat_object@active.ident <- cluster_cell_class_factor_final
		if (!file.exists(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))
		tk_int_dotplots(seurat_object, tk_int_parallel_de_pairwise_plots_results)
		setwd(glue("{orig_wd}/{out_dir}"))
		










		print("Starting set2 plots: de_one_vs_all")
		# Set2: GET ONE_VS_ALL PLOTS AND LISTS
		r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_one_vs_all")
		out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_one_vs_all_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		names(out) <- base::sapply(out, function(x) x[["name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		tk_int_parallel_de_one_vs_all_plots_results <- out
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# The openxlsx package is using an old zip method, so you might get warnings (but should still work)
		# Note: zip::zip() is deprecated, please use zip::zipr() instead
		# https://github.com/awalker89/openxlsx/issues/454
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_one_vs_all.xlsx")
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_one_vs_all_signif.xlsx")
        # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
        comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
        comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_one_vs_all/umap_gene_expression_flattened"))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_one_vs_all/umap_gene_expression_flattened"))		
		# Heatmap
		# Get only top DE genes from each comparison
		top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
        top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
        # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
        top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
        top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
		# Heatmap (top up and and down 5 genes)
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))    
		# Heatmap (top up or down genes, farthest from the mean)       
 		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))             
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_one_vs_all_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))
		# Dotplots
		seurat_object@active.ident <- cluster_cell_class_factor_final
		if (!file.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))
		tk_int_dotplots(seurat_object, tk_int_parallel_de_one_vs_all_plots_results)
		setwd(glue("{orig_wd}/{out_dir}"))
			








		
        if (numb_of_samples > 1) {
            print("Multiple samples in analysis. Starting across sample plots.")


        
        
            print("Starting set3 plots: de_pairwise_conserved")
            # Set3: GET PAIRWISE CONSERVED PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_pairwise_conserved")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_pairwise_conserved_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_pairwise_conserved_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_pairwise_cons.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_pairwise_cons_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_pairwise_cons/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_pairwise_cons/umap_gene_expression_flattened"))		
            # Heatmap 
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg_cons, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg_cons, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_pairwise_cons_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_pairwise_conserved_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))


                
        
        
        








            print("Starting set4 plots: de_one_vs_all_conserved")
            # Set4: GET ONE_VS_ALL CONSERVED PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_one_vs_all_conserved")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_one_vs_all_conserved_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_one_vs_all_conserved_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_one_vs_all_cons.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_one_vs_all_cons_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_one_vs_all_cons/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_one_vs_all_cons/umap_gene_expression_flattened"))
            # Heatmap
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg_cons, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg_cons, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))             
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_one_vs_all_cons_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_one_vs_all_conserved_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))


                
        
        
        
        
            print("Starting set5 plots: de_across_samples")
            # Set5: DE ACROSS SAMPLES -- within same cluster -- PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_across_samples")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_across_samples_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_across_samples_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_across_samples.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_across_samples_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_across_samples/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_across_samples/umap_gene_expression_flattened"))            
            # Heatmap
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))             
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_across_samples_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_across_samples_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))


















            print("Starting set6 plots: de_across_samples_and_clusters")
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_across_samples_and_clusters")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_across_samples_and_clusters_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_across_samples_and_clusters_plots_results <- out
                # ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
                # excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
                # ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
                # excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
                # # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
                # # Note: zip::zip() is deprecated, please use zip::zipr() instead
                # # https://github.com/awalker89/openxlsx/issues/454
                # # xlsx tab names cannot be longer than 31 characters, clip them
                # excel_list_of_df_de_XLSX <- excel_list_of_df_de
                # names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
                # excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
                # names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
                # openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_across_samples_and_clusters.xlsx")
                # openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_across_samples_and_clusters_signif.xlsx")
                # # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
                # comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples_and_clusters/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
                # comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples_and_clusters/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
                # base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_across_samples_and_clusters/umap_gene_expression_flattened"))
                # base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_across_samples_and_clusters/umap_gene_expression_flattened"))            
            # Heatmap
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))             
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_across_samples_and_clusters_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples_and_clusters/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_across_samples_and_clusters_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))




        } else if (numb_of_samples == 1) {
            print("One sample in analysis. Skipping across sample plots.")
        }






        ## EXPORT DATA
		DefaultAssay(seurat_object) <- orig_assay


        saveRDS(seurat_object@meta.data, "cell_metadata.rds")
        write_tsv(as_tibble(seurat_object@meta.data, rownames = "cell_id"), "cell_metadata.txt")
		# Export results to list of list
		curr_return_list <- list(cluster_cell_class_factor_final = cluster_cell_class_factor_final,
			final_cluster_labels = final_cluster_labels, 
			iteration_data_list = iteration_data_list, 
			merge_data_list = merge_data_list, 
			de_data_list = de_data_list,
			mapping = mapping, 
			mapping_full = mapping_full, 
			lfc_threshold = lfc_threshold, 
			num_genes_diff_between_clusters_threshold = num_genes_diff_between_clusters_threshold,
			umap_plot_list = umap_plot_list,
			seurat_object_meta_data = seurat_object@meta.data)
		return_list[[j]] <- curr_return_list
		names(return_list)[j] <- cluster_res_name[j]
	}
	return(return_list)
}





tk_int_parallel_de_across_samples_and_clusters <- function(t, seurat_object = seurat_object, cluster_sample_matrix = cluster_sample_matrix, lfc_threshold = lfc_threshold) {
    curr_sample_a <- cluster_sample_matrix[1, t]
    curr_sample_b <- cluster_sample_matrix[2, t]
    ident_1 <- curr_sample_a
    ident_2 <- curr_sample_b
    curr_name <- paste0(curr_sample_a,"_vs_", curr_sample_b)
    # Initialize to "not" skip 
    sample_specific_cells <- seurat_object@active.ident == ident_1
    ident_1_cells_num <- dim(GetAssayData(seurat_object, slot = "data")[, sample_specific_cells, drop = FALSE])[2]
    sample_specific_cells <- seurat_object@active.ident == ident_2
    ident_2_cells_num <- dim(GetAssayData(seurat_object, slot = "data")[, sample_specific_cells, drop = FALSE])[2]
    if (ident_1_cells_num < 3 | ident_2_cells_num < 3) {
        # Too few files for DE analysis, skip DE analysis
        print(paste0("Too few files for DE analysis, skipping DE analysis for ident: ", curr_name))
        # Create empty tibble
        de_blank <- matrix(0, 0, ncol = 6)
        colnames(de_blank) <- c("symbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
        de_blank <- as_tibble(de_blank)
        de_blank <- de_blank %>%
            mutate_all(as.character)
        de <- de_blank
        de_signif <- de_blank
        number_of_signif_genes <- 0
        results <- list(name = curr_name, sample_a = curr_sample_a, sample_b = curr_sample_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
    } else {
        print(paste0(curr_name, " DE analysis."))
        # Find DE genes
        de <- FindMarkers(object = seurat_object, slot = "data", ident.1 = ident_1, ident.2 = ident_2, min.pct = 0, test.use = "wilcox", logfc.threshold = 0)
        de <- de %>%
            tibble::rownames_to_column(var = "symbol") %>%
            as_tibble()
        # Pass thresholds
        # Sort the table by the first log2FC column group
        de_signif <- de %>%
            dplyr::filter(abs(avg_log2FC) >= lfc_threshold & p_val_adj <= 0.01) %>%
            dplyr::arrange(desc(avg_log2FC))	
        # How many genes are significant?
        number_of_signif_genes <- dim(de_signif)[1]
        results <- list(name = curr_name, sample_a = curr_sample_a, sample_b = curr_sample_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
    }
    all_results <- list(results)
    names(all_results)[1] <- curr_name
	return(all_results)
}













tk_int_parallel_de_across_samples_and_clusters_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	if (!dir.exists("de_across_samples_and_clusters/de_text_files")) {dir.create("de_across_samples_and_clusters/de_text_files", recursive = TRUE)}
	# if (!dir.exists("de_across_samples_and_clusters/umap_gene_expression")) {dir.create("de_across_samples_and_clusters/umap_gene_expression", recursive = TRUE)}
	# DE LISTS ALL_VS_ONE
	write_tsv(de_data_list[[r]]$de, paste0("de_across_samples_and_clusters/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_across_samples_and_clusters/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS ALL_VS_ONE
	starting_dir <- getwd()
	# setwd("de_across_samples_and_clusters/umap_gene_expression")
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC > 0]
	# tk_int_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot)
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	# tk_int_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot)
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(starting_dir)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}
