library(tidyverse)
library(ggpubr)

# plot_gene_count() -------------------------------------------------------

# => Plots n counts by genotype and treatment (or exposition) for a given transcript
# => Uses DESeq2::plotCounts() which normalizes counts by the estimated size factors and adds a pseudocount of 1/2 to allow for log scale plotting
plot_gene_count <- function(dds, 
                            res,
                            gene_metadata, 
                            dds_gene_id_type = "ensembl_gene_id_version", 
                            user_gene_id_type = "gene_name", 
                            gene,
                            dataset = "pff",
                            pr = F) {
  
  dds_gene <- gene_metadata[gene_metadata[[user_gene_id_type]] == gene, dds_gene_id_type][1,][[1]]
  
  if(dataset == "pff") {
    
    p <- ggplot(data = DESeq2::plotCounts(dds = dds, 
                                          gene = which(res@rownames == dds_gene),
                                          intgroup = c("genotype", "treatment", "line_name"),
                                          returnData = T) %>% 
                  dplyr::mutate(genotype_treatment = paste0(genotype, " ", treatment)),
                mapping = aes(x = genotype_treatment, y = count, fill = line_name))
    
    if(isFALSE(pr)) {
      p <- p +
        geom_line(aes(group = line_name), color = "orange", linetype = "dashed") 
    }
    p <- p +
      geom_point(position = position_dodge(width = 0.2), size = 6, shape = 21) +
      theme_pubr(legend = "right") +
      scale_fill_manual(values = c("#00b4d8", "#0096c7", "#0077b6", "#ff4d6d", "#c9184a", "#a4133c", "#29a655")) +
      scale_y_log10(n.breaks = 10) +
      rotate_x_text(45) +
      labs(x = " genotype / treatment",
           y = paste0(gene, " normalized count"),
           fill = "line ID")
    
  } else if(dataset == "coculture") {
    
    p <-  ggplot(data = DESeq2::plotCounts(dds = dds, 
                                           gene = which(res@rownames == dds_gene),
                                           intgroup = c("genotype", "exposition", "line_name"),
                                           returnData = T) %>%
                   dplyr::mutate(genotype_exposition = paste0(genotype, " ", exposition),
                                 order = case_when(genotype_exposition == "CTRL unexposed" ~ 1L, 
                                                   genotype_exposition == "CTRL TRIP-exposed" ~ 2L, 
                                                   genotype_exposition == "KO unexposed" ~ 3L, 
                                                   genotype_exposition == "KO TRIP-exposed" ~ 4L,
                                                   genotype_exposition == "TRIP-exposed" ~ 5L,
                                                   TRUE ~ NA_integer_),
                                 genotype_exposition = reorder(genotype_exposition, order)),
                 mapping = aes(x = genotype_exposition,
                               y = count, 
                               fill = line_name))
    if(isFALSE(pr)) {
      p <- p +
        geom_line(aes(group = line_name), color = "orange", linetype = "dashed") 
    }
    p <- p +
      geom_point(position = position_dodge(width = 0.2), size = 6, shape = 21) +
      theme_pubr(legend = "right") +
      scale_fill_manual(values = c("#00b4d8", "#0096c7", "#0077b6", "#29a655", "#a4133c")) +
      scale_y_log10(n.breaks = 10) +
      rotate_x_text(45) +
      labs(x = " genotype / exposition",
           y = paste0(gene, " normalized count"),
           fill = "line ID")
    
  }
  p
}

# compute_dge() -------------------------------------------------------------------------

compute_dge <- function(files,
                        tx2gene,
                        colData,
                        min_counts = 0,
                        min_samples = 0,
                        design,
                        alpha = 0.1,
                        lfcThreshold = 0,
                        altHypothesis = "greaterAbs",
                        effect_name = NA_character_,
                        filterFun = NULL,
                        feature_id_var_name = "ensembl_gene_id_version",
                        lfc_shrinkage_method = "apeglm",
                        gene_metadata = NULL,
                        parallel = F,
                        n_workers = 8) {
  
  colData <- colData %>%
    mutate_if(.predicate = is.factor,
              .funs = droplevels)
  
  txi <- tximport::tximport(files = files,
                            type = "salmon",
                            txOut = FALSE,
                            countsFromAbundance = "no",
                            tx2gene = tx2gene)
  
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi,
                                          colData = colData,
                                          design = ~ 1)
  
  dds <- dds[rowSums(BiocGenerics::counts(dds) >= min_counts) >= min_samples,]
  
  user_design <- model.matrix(object = formula(design), 
                              data = colData)
  attributes(user_design) <- attributes(user_design)[c("dim", "dimnames")]
  
  reduced_user_design <- user_design[, which(colSums(user_design) != 0)]
  
  BiocGenerics::design(dds) <- reduced_user_design
  
  if (parallel) {
    
    dds <- DESeq2::DESeq(object = dds, 
                         parallel = parallel,
                         BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))
    
  } else {
    
    dds <- DESeq2::DESeq(object = dds,
                         parallel = parallel)
    
  }
  
  normalised_counts <- as.data.frame(BiocGenerics::counts(dds, normalized = T)) %>% 
    rownames_to_column(var = feature_id_var_name)
  
  calculated_vals <- as.data.frame(S4Vectors::mcols(dds, use.names = T))
  
  calculated_vals_description <- as.data.frame(S4Vectors::mcols(S4Vectors::mcols(dds, use.names = T)))
  
  if(is.na(effect_name)) {
    
    effect_name <- DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))]
    
  }
  
  message("\nSelected design: ", design)
  
  message("\nAvailable effect names:\n", 
          paste0(DESeq2::resultsNames(dds), collapse = "\n"))
  
  message("\nSelected effect name: ", effect_name, "\n")
  
  if (!effect_name %in% DESeq2::resultsNames(dds)) {
    
    stop("\nSelected effect name is not an available individual effect name for the current design\n")
    
  }
  
  if (is.null(filterFun)) {
    
    res <- DESeq2::results(object = dds,
                           alpha = alpha,
                           name = effect_name,
                           lfcThreshold = lfcThreshold,
                           altHypothesis = altHypothesis)
    
  } else {
    
    res <- DESeq2::results(object = dds,
                           alpha = alpha,
                           name = effect_name,
                           lfcThreshold = lfcThreshold,
                           altHypothesis = altHypothesis,
                           filterFun = filterFun)
    
  }
  
  res_df <- as.data.frame(res) %>% 
    rownames_to_column(var = feature_id_var_name)
  
  if (parallel) {
    
    res_shrunken_lfc <- DESeq2::lfcShrink(dds = dds,
                                          coef = effect_name,
                                          type = lfc_shrinkage_method,
                                          parallel = parallel,
                                          BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))
    
  } else {
    
    res_shrunken_lfc <- DESeq2::lfcShrink(dds = dds,
                                          coef = effect_name,
                                          type = lfc_shrinkage_method,
                                          parallel = parallel)
    
  }
  
  res_df <- res_df %>% 
    full_join(as.data.frame(res_shrunken_lfc) %>%
                rownames_to_column(var = feature_id_var_name) %>% 
                dplyr::rename(shrunken_log2FoldChange = log2FoldChange,
                              shrunken_lfcSE = lfcSE) %>% 
                dplyr::select(all_of(feature_id_var_name), shrunken_log2FoldChange, shrunken_lfcSE),
              by = feature_id_var_name) %>% 
    dplyr::mutate(sign_log2fc_times_minus_log10pvalue = sign(log2FoldChange) * -log(x = pvalue, base = 10)) %>% 
    dplyr::arrange(desc(sign_log2fc_times_minus_log10pvalue))
  
  deg <- with(res_df, res_df[!is.na(padj) & padj < alpha,]) %>% 
    dplyr::arrange(desc(log2FoldChange))
  
  if(!is.null(gene_metadata)) {
    
    if (!feature_id_var_name %in% colnames(gene_metadata)) {
      
      message("Column ", feature_id_var_name, " not present in gene_metadata: not joining res_df and deg to gene_metadata\n")
      
    } else {
      
      res_df <- res_df %>% 
        left_join(gene_metadata, by = feature_id_var_name) 
      
      deg <- deg %>% 
        left_join(gene_metadata, by = feature_id_var_name) 
      
    }
  }
  
  deg_up <- deg[deg$log2FoldChange > 0,]
  deg_down <- deg[deg$log2FoldChange < 0,] %>% 
    dplyr::arrange(log2FoldChange)
  
  if(identical(user_design, reduced_user_design)) {
    
    out <- list(design = user_design,
                dds = dds,
                normalised_counts = normalised_counts,
                calculated_vals = calculated_vals,
                calculated_vals_description = calculated_vals_description,
                res = res,
                res_df = res_df,
                deg = deg,
                deg_up = deg_up,
                deg_down = deg_down)
    
  } else {
    
    out <- list(design = user_design,
                user_design = reduced_user_design,
                dds = dds,
                normalised_counts = normalised_counts,
                calculated_vals = calculated_vals,
                calculated_vals_description = calculated_vals_description,
                res = res,
                res_df = res_df,
                deg = deg,
                deg_up = deg_up,
                deg_down = deg_down)
    
  }
}

# recompute_padj_on_subset() ---------------------------------------------

recompute_padj_on_subset <- function(res,
                                     feature_id_filter_var,
                                     filter_var_pattern = ".*",
                                     alpha = 0.1,
                                     feature_id_var_name = "ensembl_gene_id_version",
                                     gene_metadata = NULL,
                                     ignore.case = F,
                                     fixed = F) {
  
  matches <- feature_id_filter_var[grepl(pattern = filter_var_pattern, 
                                         x = feature_id_filter_var[[2]],
                                         ignore.case = ignore.case,
                                         fixed = fixed),][[1]]
  
  res <- res[rownames(res) %in% matches,]
  
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  deg <- as.data.frame(res) 
  deg <- with(deg, deg[!is.na(padj) & padj < alpha,]) %>% 
    rownames_to_column(var = feature_id_var_name) %>% 
    dplyr::arrange(desc(log2FoldChange)) 
  
  if(!is.null(gene_metadata)) {
    if (!feature_id_var_name %in% colnames(gene_metadata)) {
      message("Column ", feature_id_var_name, " not present in gene_metadata: not joining deg to gene_metadata\n")
    } else {
      deg <- deg %>% 
        left_join(gene_metadata, by = feature_id_var_name) 
    }
  }
  
  deg_up <- deg[deg$log2FoldChange > 0,]
  deg_down <- deg[deg$log2FoldChange < 0,] %>% 
    dplyr::arrange(log2FoldChange)
  
  out <- list(filter_var_pattern = filter_var_pattern,
              res = res,
              deg = deg,
              deg_up = deg_up,
              deg_down = deg_down)
}

# plot_any_pc() -----------------------------------------------------------------------

plot_any_pc <- function(object,
                        intgroup,
                        ntop = 500,
                        first_pc = 1,
                        second_pc = 2,
                        returnData = F,
                        fill = NULL,
                        shape = NULL,
                        group = NULL,
                        fill_color_set = c("#00b4d8", "#0096c7", "#0077b6", "#ff4d6d", "#c9184a", "#a4133c", "#29a655"),
                        shape_set = c(21, 22)) {
  
  plotPCA.DESeqTransform <- function(object, 
                                     intgroup = "condition",
                                     ntop = 500, 
                                     returnData = FALSE, 
                                     pcsToUse = 1:2)
  {
    # calculate the variance for each gene
    rv <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(SummarizedExperiment::assay(object)[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(SummarizedExperiment::colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(SummarizedExperiment::colData(object)[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=":"))
    } else {
      SummarizedExperiment::colData(object)[[intgroup]]
    }
    
    # assembly the data for the plot
    pcs <- paste0("PC", pcsToUse)
    d <- data.frame(V1=pca$x[,pcsToUse[1]],
                    V2=pca$x[,pcsToUse[2]],
                    group=group, intgroup.df, name=colnames(object))
    colnames(d)[1:2] <- pcs
    
    if (returnData) {
      attr(d, "percentVar") <- percentVar[pcsToUse]
      return(d)
    }
    
    ggplot(data=d, aes_string(x=pcs[1], y=pcs[2], color="group")) +
      geom_point(size=3) + 
      xlab(paste0(pcs[1],": ",round(percentVar[pcsToUse[1]] * 100),"% variance")) +
      ylab(paste0(pcs[2],": ",round(percentVar[pcsToUse[2]] * 100),"% variance")) +
      coord_fixed()
  }
  
  pca_data <- plotPCA.DESeqTransform(object = object, 
                                     intgroup = intgroup,
                                     ntop = ntop,
                                     returnData = T,
                                     pcsToUse = c(first_pc, second_pc))
  
  pct_var <- attr(pca_data, "percentVar") * 100
  
  if (returnData) {
    pca_results <- list(pca_data = pca_data,
                        pct_var = pct_var)
  } else {
    pca_plot <- ggplot(data = pca_data,
                       mapping = aes(x = !!sym(paste0("PC", first_pc)), 
                                     y = !!sym(paste0("PC", second_pc)),
                                     fill = !!sym(fill),
                                     shape = !!sym(shape))) +
      geom_vline(xintercept = 0, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "orange", linetype = "dashed")
    
    if (!is.null(group)) {
      pca_plot <- pca_plot +
        geom_line(aes(group = !!sym(group)), color = "grey")
    }
    
    pca_plot <- pca_plot +
      geom_point(size = 5, color = "black") +
      theme_pubr(legend = "right") +
      scale_fill_manual(values = fill_color_set) +
      scale_shape_manual(values = shape_set) +
      guides(fill = guide_legend(override.aes = list(shape = shape_set[[1]]))) +
      labs(x = paste0("PC", first_pc, ": ", round(pct_var[1], digits = 2), "% variance"),
           y = paste0("PC", second_pc, ": ", round(pct_var[2], digits = 2), "% variance"))
    
    pca_plot
  }
}

# get_pc_influencing_features() -------------------------------------------

get_pc_influencing_features <- function(object,
                                        ntop = NULL,
                                        pc = 1,
                                        top_features_pct = 5,
                                        feature_id_var_name = "ensembl_gene_id_version",
                                        gene_metadata = NULL) {
  
  rv <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))
  
  ntop <- ifelse(is.null(ntop), length(rv), ntop)
  
  selection <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca <- prcomp(t(SummarizedExperiment::assay(object)[selection,]), 
                center = T,
                scale. = F)
  
  var_loadings <- as.data.frame(pca$rotation) %>% 
    rownames_to_column(var = feature_id_var_name) 
  
  pc <- paste0("PC", pc)
  
  var_loadings_thres <- quantile(x = abs(var_loadings[[pc]]), 
                                 probs = 1-(top_features_pct/100))
  
  influencing_features <- var_loadings[abs(var_loadings[[pc]]) >= var_loadings_thres, c(feature_id_var_name, pc)] 
  
  if(!is.null(gene_metadata)) {
    if (!feature_id_var_name %in% colnames(gene_metadata)) {
      message("Column ", feature_id_var_name, " not present in gene_metadata: not joining deg to gene_metadata\n")
    } else {
      influencing_features <- influencing_features %>% 
        left_join(gene_metadata, by = feature_id_var_name) 
    }
  }
  influencing_features
}

# import_salmon_meta_info() -----------------------------------------------

import_salmon_meta_info <- function(file) {
  
  meta_info <- rjson::fromJSON(file = file) 
  
  meta_info <- meta_info[c("library_types", 
                           "frag_dist_length", "frag_length_mean", "frag_length_sd", 
                           "seq_bias_correct", "gc_bias_correct", "keep_duplicates", 
                           "num_valid_targets", "num_decoy_targets", 
                           "num_processed", "num_mapped", 
                           "num_decoy_fragments", "num_dovetail_fragments", 
                           "num_fragments_filtered_vm", "num_alignments_below_threshold_for_mapped_fragments_vm", 
                           "percent_mapped")]
  
  meta_info <- as_tibble(meta_info) %>% 
    mutate(file = file) %>% 
    select(file, everything())
  
}

# summarise_average_transcript_length() -------------------------------------------------------------------------

summarise_average_transcript_length <- function(average_transcript_length_matrix,
                                                barcode_summary_variable_cell_type,
                                                cell_type,
                                                summary_variable) {
  
  # Subset the list of barcodes to only the cell type of interest
  filtered <- barcode_summary_variable_cell_type[barcode_summary_variable_cell_type$cell_type == cell_type,]

  # Summarize the average transcript length matrix to summary variable-level
  # => Compute the average transcript length per combination of gene and summary variable level
  out <- average_transcript_length_matrix %>% 
    select(filtered$barcode) %>% 
    rownames_to_column(var="ensembl_gene_id_version") %>% 
    pivot_longer(cols=-ensembl_gene_id_version, names_to = "barcode", values_to = "length") %>% 
    left_join(filtered, by = join_by(barcode)) %>% 
    nest(data = -c(all_of(summary_variable), ensembl_gene_id_version)) %>% 
    mutate(mean_length = map_dbl(.x = data,
                                 .f = ~ mean(.x$length))) %>% 
    select(-data) %>% 
    pivot_wider(names_from = all_of(summary_variable), values_from = mean_length)
  
}

# intersect_deg() ---------------------------------------------------------

intersect_deg <- function(deg_set1, deg_set2) {
  
  deg_set1_nm <- as.character(substitute(deg_set1))
  deg_set2_nm <- as.character(substitute(deg_set2))
  
  deg_set1 <- if(is.data.frame(deg_set1)) { deg_set1$gene_name } else { deg_set1 }
  deg_set2 <- if(is.data.frame(deg_set2)) { deg_set2$gene_name } else { deg_set2 }
  
  n_deg_set1 <- length(deg_set1)
  n_deg_set1_in_set2 <- sum(deg_set1 %in% deg_set2)
  
  n_deg_set2 <- length(deg_set2)
  n_deg_set2_in_set1 <- sum(deg_set2 %in% deg_set1)
  
  n_deg_set1_set2 <- length(union(deg_set1, deg_set2))
  intersecting_deg <- intersect(deg_set1, deg_set2)
  n_intersecting_deg <- length(intersecting_deg)
  
  message("DEG set 1: ", deg_set1_nm, "\n",
          "DEG set 2: ", deg_set2_nm, "\n",
          "n DEG set 1 in set 2: ", n_deg_set1_in_set2, "/", n_deg_set1, " (", round(n_deg_set1_in_set2/n_deg_set1*100, 1), "%)\n",
          "n DEG set 2 in set 1: ", n_deg_set2_in_set1, "/", n_deg_set2, " (", round(n_deg_set2_in_set1/n_deg_set2*100, 1), "%)\n",
          "n intersecting DEG: ", n_intersecting_deg, "/", n_deg_set1_set2, " (", round(n_intersecting_deg/n_deg_set1_set2*100, 1), "%)\n")
  
  return(intersecting_deg)
  
}
