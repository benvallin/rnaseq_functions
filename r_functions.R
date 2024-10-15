library(tidyverse)
library(ggpubr)
library(fgsea)
library(MAST)
library(limma)

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


# make_fgsea_pathways() ---------------------------------------------------

make_fgsea_pathways <- function(gene_metadata,
                                collections_names) {
  
  output <- lapply(X = collections_names,
                   FUN = function(x) {
                     
                     out_i <- as_tibble(gene_metadata) %>% 
                       select(ensembl_gene_id_version, all_of(x)) %>% 
                       unnest(all_of(x), keep_empty = F) %>% 
                       nest(data = ensembl_gene_id_version) %>% 
                       mutate(data = map(.x = data, .f = ~.x$ensembl_gene_id_version))
                     
                     out_i <- as.list(setNames(object = out_i[["data"]], nm = out_i[["set_name"]]))
                     
                   }) %>% setNames(nm = collections_names)
}

# batch_fgsea() -----------------------------------------------------------

batch_fgsea <- function(gene_metadata,
                        collections_names,
                        stats, 
                        min_set_size = 1,
                        max_set_size = length(stats) - 1,
                        padj_threshold = Inf,
                        ...) {
  
  # This does not make any difference: fgsea::fgsea() subsets pathways content to only the genes provided in stats
  # gene_metadata <- gene_metadata[gene_metadata$ensembl_gene_id_version %in% names(stats),] 
  
  collections <- make_fgsea_pathways(gene_metadata = gene_metadata,
                                     collections_names = collections_names)
  
  output <- lapply(X = names(collections),
                   FUN = function(x) {
                     
                     out_i <- fgsea(pathways = collections[[x]],
                                    stats = stats,
                                    minSize = min_set_size,
                                    maxSize = max_set_size,
                                    ...)
                     
                     out_i <- as_tibble(out_i) %>%
                       mutate(collection = x,
                              n_leadingEdge = map_int(.x = leadingEdge, 
                                                      .f = ~ length(.x)),
                              n_leadingEdge_over_total = paste0(n_leadingEdge, " / ", size),
                              pct_leadingEdge = (n_leadingEdge / size) * 100) %>% 
                       select(collection, set_name = pathway, everything()) %>% 
                       filter(padj < padj_threshold) %>% 
                       arrange(desc(NES), padj)
                     
                   }) %>% bind_rows()
}

# batch_ora() -------------------------------------------------------------

batch_ora <- function(gene_metadata,
                      collections_names,
                      genes,
                      universe,
                      min_set_size = 1,
                      max_set_size = length(stats) - 1,
                      padj_threshold = Inf,
                      ...) {
  
  collections <- make_fgsea_pathways(gene_metadata = gene_metadata,
                                     collections_names = collections_names)
  
  output <- lapply(X = names(collections),
                   FUN = function(x) {
                     
                     genes_collection_i <- gene_metadata %>% 
                       select(ensembl_gene_id_version, all_of(x)) %>% 
                       unnest(all_of(x), keep_empty = F)
                     
                     universe_i <- genes_collection_i %>% 
                       filter(ensembl_gene_id_version %in% universe) %>% 
                       pull(ensembl_gene_id_version) %>% 
                       unique()
                     
                     genes_i <- genes_collection_i %>% 
                       filter(ensembl_gene_id_version %in% genes) %>% 
                       pull(ensembl_gene_id_version) %>% 
                       unique()
                     
                     out_i <- fgsea::fora(pathways = collections[[x]],
                                          genes = genes_i,
                                          universe = universe_i,
                                          minSize = min_set_size,
                                          maxSize = max_set_size,
                                          ...)
                     
                     out_i <- as_tibble(out_i) %>%
                       mutate(collection = x,
                              n_overlapGenes = paste0(overlap, " / ", size),
                              pct_overlapGenes = (overlap / size) * 100) %>% 
                       select(collection, set_name = pathway, everything()) %>% 
                       filter(padj < padj_threshold) %>% 
                       arrange(padj)
                     
                   }) %>% bind_rows()
}

# batch_ids2indices() -----------------------------------------------------

batch_ids2indices <- function(gene_metadata,
                              collections_names,
                              identifiers) {
  
  output <- lapply(X = collections_names,
                   FUN = function(x) {
                     
                     out_i <- as_tibble(gene_metadata) %>% 
                       select(ensembl_gene_id_version, all_of(x)) %>% 
                       unnest(all_of(x), keep_empty = F) %>%
                       nest(data = -set_name)
                     
                     out_i <- setNames(object = map(.x = out_i$data,
                                                    .f = ~ .x$ensembl_gene_id_version),
                                       nm = out_i$set_name)
                     
                     out_i <- ids2indices(gene.sets = out_i, 
                                          identifiers = identifiers, 
                                          remove.empty = T)
                     
                   }) %>% setNames(nm = collections_names)
}

# batch_gseaAfterBoot() ---------------------------------------------------

batch_gseaAfterBoot <- function(sca,
                                gene_metadata,
                                collections_names,
                                zFit, 
                                boots, 
                                contrast, 
                                min_set_size = 1, 
                                max_set_size = dim(sca)[[1]] - 1, 
                                padj_threshold = Inf,
                                ...) {
  
  collections <- batch_ids2indices(gene_metadata = gene_metadata,
                                   collections_names = collections_names,
                                   identifiers = rownames(mcols(sca)))
  
  output <- lapply(X = names(collections),
                   FUN = function(x) {
                     
                     sets <- collections[[x]]
                     sets <- sets[map_lgl(.x = sets, .f = ~ length(.x) >= min_set_size)]
                     sets <- sets[map_lgl(.x = sets, .f = ~ length(.x) <= max_set_size)]
                     
                     if(length(sets) > 1) {
                       
                       out_i <- gseaAfterBoot(zFit = zFit, 
                                              boots = boots, 
                                              sets = sets, 
                                              hypothesis = CoefficientHypothesis(contrast),
                                              ...) 
                       
                       out_i <- summary(out_i, testType = "normal") %>% 
                         as_tibble() %>% 
                         mutate(collection = x) %>% 
                         select(collection, set_name = set, everything()) %>% 
                         filter(combined_adj < padj_threshold) %>% 
                         arrange(desc(combined_Z), combined_adj)
                       
                     } else {
                       
                       out_i <-  NULL
                       
                     }
                   }) %>% bind_rows()
}

# batch_limma_gsea() ------------------------------------------------------

batch_limma_gsea <- function(y,
                             gene_metadata,
                             collections_names,
                             design,
                             contrast,
                             min_set_size = 1, 
                             max_set_size = dim(y)[[1]] - 1, 
                             padj_threshold = Inf,
                             limma_test = "fry",
                             ...) {
  
  message("Performing ", limma_test, " test...")
  
  collections <- batch_ids2indices(gene_metadata = gene_metadata,
                                   collections_names = collections_names,
                                   identifiers = rownames(y))
  
  output <- lapply(X = names(collections),
                   FUN = function(x) {
                     
                     sets <- collections[[x]]
                     sets <- sets[map_lgl(.x = sets, .f = ~ length(.x) >= min_set_size)]
                     sets <- sets[map_lgl(.x = sets, .f = ~ length(.x) <= max_set_size)]
                     
                     if(length(sets) == 0) {
                       NULL
                     } else {
                       out_i <- do.call(what = limma_test,
                                        args = list(y = y,
                                                    index = sets,
                                                    design = design,
                                                    contrast = contrast,
                                                    ...))
                       
                       out_i <- rownames_to_column(out_i, var = "set_name") %>%
                         as_tibble() %>%
                         mutate(collection = x) %>%
                         select(collection, set_name, everything()) 
                       
                       if("FDR" %in% colnames(out_i)) {
                         out_i <- out_i %>% 
                           filter(FDR < padj_threshold) %>%
                           arrange(desc(Direction), FDR)
                       }
                     }
                   }) %>% bind_rows()
}

# get_genes_in_collection() -----------------------------------------------

get_genes_in_collection <- function(gene_metadata,
                                    collection,
                                    set_name) {
  
  gene_metadata %>% 
    select(ensembl_gene_id_version, gene_name, all_of(collection)) %>% 
    unnest(all_of(collection)) %>% 
    filter(!!set_name == set_name) %>% 
    mutate(collection = collection) %>% 
    select(collection, set_name, ensembl_gene_id_version, gene_name)
  
}

# get_leading_edge() ------------------------------------------------------

get_leading_edge <- function(fgsea_df,
                             set_name) {
  
  fgsea_df[fgsea_df$set_name == set_name, "leadingEdge"][[1]][[1]]
  
}

# plot_mast_lfc() ---------------------------------------------------------

plot_mast_lfc <- function(lfc_df,
                          genes,
                          key_genes = NULL,
                          title = NULL,
                          test_cond = "PFF",
                          ref_cond = "NT",
                          key_genes_name = "leading edge",
                          desc_log2fc = T,
                          xintercept = NULL) {
  
  lfc_df <- lfc_df %>% 
    filter(ensembl_gene_id_version %in% genes,
           !is.na(log2fc)) %>% 
    mutate(fdr_value = case_when(fdr < 0.01 ~ "FDR < 0.01",
                                 fdr < 0.05 ~ "FDR < 0.05", 
                                 T ~ paste0("FDR ", utf8::utf8_print("\u2265"), " 0.05")),
           fdr_color = case_when(fdr < 0.01 ~ "#e0007f",
                                 fdr < 0.05 ~ "#b892ff", 
                                 T ~ "#ff9100")) 
  
  if(!is.null(key_genes)) {
    lfc_df <- lfc_df %>% 
      mutate(in_key_genes = ifelse(ensembl_gene_id_version %in% key_genes, 
                                   paste0("in ", key_genes_name),
                                   paste0("not in ", key_genes_name)))
  }
  
  out <- ggplot(data = lfc_df,
                mapping = aes(x = log2fc,
                              y = fct_reorder(gene_name, log2fc, .desc = desc_log2fc),
                              fill = fdr_value)) +
    geom_hline(yintercept = lfc_df$gene_name,
               color = "grey", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#c9184a", linewidth = 0.5) +
    geom_vline(xintercept = xintercept, linetype = "dashed", color = "orange", linewidth = 0.5) +
    ggpubr::theme_pubr(legend = "bottom") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 14)) +
    scale_fill_manual(values = lfc_df %>% arrange(fdr) %>% pull(fdr_color) %>% unique()) +
    labs(title = title,
         x = paste0("log", utf8::utf8_print("\u2082"), " fold change (", test_cond, " vs ", ref_cond, ")"),
         y = NULL,
         fill = NULL,
         shape = NULL)
  
  if(is.null(key_genes)) {
    out +
      geom_point(shape = 21, color = "black", size = 5, stroke = 1)
  } else {
    if(length(unique(lfc_df$in_key_genes)) == 2) {
      shapes <- c(23, 21)
    } else if(unique(lfc_df$in_key_genes) == paste0("in ", key_genes_name)) {
      shapes <- 23
    } else {
      shapes <- 21
    }
    out +
      geom_point(aes(shape = in_key_genes), color = "black", size = 5, stroke = 1) +
      scale_shape_manual(values = shapes) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))
  }
}

# plot_mast_pct_exp() -----------------------------------------------------

plot_mast_pct_exp <- function(lfc_df, genes, ref_cond, test_cond, title, fill_colors = c("#ca3767", "#5d76cb")) {
  
  pct_exp_df <- lfc_df %>% 
    filter(ensembl_gene_id_version %in% genes, !is.na(log2fc)) 
  
  pct_exp_df <- pct_exp_df %>% 
    mutate(gene_name = fct_relevel(gene_name, 
                                   pct_exp_df %>% arrange(log2fc) %>% pull(gene_name)),
           reference = 100,
           test = (2^log2fc)*100) %>% 
    select(ensembl_gene_id_version, gene_name, reference, test) %>% 
    pivot_longer(cols = c(reference, test), 
                 names_to = "condition", 
                 values_to = "pct_expression") %>% 
    mutate(condition = case_when(condition == "reference" ~ ref_cond,
                                 condition == "test" ~ test_cond,
                                 T ~ NA_character_))
  
  if(mean(pct_exp_df$pct_expression) < 100) {
    out <- ggplot() +
      geom_col(data = pct_exp_df %>%
                 filter(condition == ref_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      geom_col(data = pct_exp_df %>%
                 filter(condition == test_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +  
      scale_y_continuous(limits = c(0L, 100L), n.breaks = 10L) 
  } else {
    out <- ggplot() +
      geom_col(data = pct_exp_df %>%
                 filter(condition == test_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      geom_col(data = pct_exp_df %>%
                 filter(condition == ref_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      scale_y_continuous(limits = c(0L, ceiling(max(pct_exp_df$pct_expression))), n.breaks = 12L) 
  }
  
  out <- out +
    ggpubr::theme_pubr(legend = "bottom") +
    scale_color_manual(values = fill_colors, aesthetics = "fill") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) +
    ggpubr::rotate_x_text(angle = 30) +
    labs(x = NULL,
         y = "% expression",
         title = title)
  
  out
}

# make_ggvenn_sets() ------------------------------------------------------

make_ggvenn_sets <- function(fgsea_df, set_names) {
  lapply(X = set_names,
         FUN = function(x) {
           get_leading_edge(fgsea_df = fgsea_df, 
                            set_name = x)
         }) %>% 
    setNames(nm = set_names)
}



# search_genes_in_leading_edges() -----------------------------------------

search_genes_in_leading_edges <- function(fgsea_df, 
                                          genes, 
                                          genes_group_name = "searched") {
  
  out <- fgsea_df %>% 
    mutate(leadingEdge_searched = map(.x = leadingEdge, .f = ~ intersect(.x, genes)),
           n_leadingEdge_searched = map_int(.x = leadingEdge_searched, .f = ~length(.x)),
           n_leadingEdge_searched_over_leadingEdge = paste0(n_leadingEdge_searched, " / ", n_leadingEdge),
           pct_leadingEdge_searched_in_leadingEdge = (n_leadingEdge_searched / n_leadingEdge)*100) %>% 
    arrange(desc(pct_leadingEdge_searched_in_leadingEdge), collection) %>% 
    mutate(set_number = seq_along(fgsea_df[[1]])) %>% 
    rename_at(.vars = c("leadingEdge_searched", 
                        "n_leadingEdge_searched", 
                        "n_leadingEdge_searched_over_leadingEdge", 
                        "pct_leadingEdge_searched_in_leadingEdge"),
              .funs = ~ str_replace(., "searched", genes_group_name)) %>% 
    select(set_number, everything())
  
}

# plot_genes_prop_in_leading_edges() --------------------------------------

plot_genes_prop_in_leading_edges <- function(fgsea_df, 
                                             genes, 
                                             x_label = "gene set",
                                             y_label_genes_group_name = "searched") {
  
  out_df <- search_genes_in_leading_edges(fgsea_df = fgsea_df,
                                          genes = genes,
                                          genes_group_name = "searched")
  
  out_df <- out_df %>%
    mutate(collection = case_when(collection == "msigdb_h" ~ "HALLMARK",
                                  collection == "msigdb_c2_cp_kegg" ~ "KEGG",
                                  collection == "msigdb_c2_cp_reactome" ~ "REACTOME",
                                  collection == "msigdb_c2_cp_wikipathways" ~ "Wikipathways",
                                  collection == "msigdb_c2_cp_biocarta" ~ "Biocarta",
                                  collection == "msigdb_c5_gobp" ~ "GO BP"),
           collection_color = case_when(collection == "HALLMARK" ~ "#ca3767",
                                        collection == "KEGG" ~ "#ff6d00",
                                        collection == "REACTOME" ~ "#5d76cb",
                                        collection == "Wikipathways" ~ "#29a655",
                                        collection == "Biocarta" ~ "#b370ff",
                                        collection == "GO BP" ~ "#ffc300"))
  
  collection_levels <- unique(out_df$collection)
  
  collection_levels <- collection_levels[match(x = c("HALLMARK", "KEGG", "REACTOME", 
                                                     "Wikipathways", "Biocarta", "GO BP"),
                                               table = collection_levels)]
  
  collection_levels <- na.omit(collection_levels)
  
  out_df <- out_df %>% 
    mutate(collection = fct_relevel(collection, collection_levels))
  
  collection_colors <- unique(out_df[, c("collection", "collection_color")])
  
  collection_colors <- collection_colors[match(x = collection_levels, 
                                               table = collection_colors$collection),
                                         "collection_color"][[1]]
  
  percentiles_x_pos <- c(0.25*length(out_df[[1]]),
                         0.5*length(out_df[[1]]),
                         0.75*length(out_df[[1]]))
  
  ggplot(data = out_df, 
         mapping = aes(x = set_number,
                       y = pct_leadingEdge_searched_in_leadingEdge)) +
    geom_vline(xintercept = out_df$set_number, color = "grey", alpha = 0.2) +
    geom_vline(xintercept = percentiles_x_pos, 
               linetype = "dashed", color = "orange", linewidth = 0.5) +
    geom_hline(yintercept = c(0, 25, 50, 75, 100), 
               linetype = "dotted", color = "#c9184a", linewidth = 0.5) +
    geom_label(data = tibble(x = percentiles_x_pos,
                             y = 100,
                             label = c("Q1", "Q2", "Q3")),
               aes(x = x, y = y, label = label)) + 
    geom_point(aes(fill = collection), 
               shape = 21, color = "black", size = 3, stroke = 1, alpha = 0.8) +
    ggpubr::theme_pubr(legend = "bottom") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16)) +
    scale_y_continuous(limits = c(0L, 100L), n.breaks = 10L) +
    scale_color_manual(values = collection_colors, aesthetics = "fill") +
    labs(x = x_label, 
         y = paste0("Proportion of leading edge corresponding to ", y_label_genes_group_name, " genes (%)"),
         fill = NULL)
}

# plot_lfc() --------------------------------------------------------------

plot_lfc <- function(lfc_df,
                     genes,
                     key_genes = NULL,
                     title = NULL,
                     test_cond = "test cond.",
                     ref_cond = "ref cond.",
                     key_genes_name = "leading edge",
                     desc_lfc = T,
                     xintercept = NULL,
                     padj_lab = T,
                     lab_nudge = 0.2,
                     input_pkg = "MAST",
                     lfcse_to_log2 = T) {
  
  if(input_pkg == "MAST") {
    lfc <- sym("log2fc")
    lfcse <- sym("lfcse")
    fdr <- sym("fdr")
  } else if(input_pkg == "DESeq2") {
    lfc <- sym("log2FoldChange")
    lfcse <- sym("lfcSE")
    fdr <- sym("padj")
  }
  
  lfc_df <- lfc_df %>% 
    filter(ensembl_gene_id_version %in% genes,
           !is.na(eval(lfc))) %>% 
    mutate(fdr_value = case_when(eval(fdr) < 0.01 ~ "adj. p < 0.01",
                                 eval(fdr) < 0.05 ~ "adj. p < 0.05", 
                                 eval(fdr) < 0.1 ~ "adj. p < 0.1", 
                                 eval(fdr) >= 0.1 ~ paste0("adj. p ", utf8::utf8_print("\u2265"), " 0.1"),
                                 is.na(eval(fdr)) ~ "adj. p not computed"),
           fdr_color = case_when(eval(fdr) < 0.01 ~ "#e0007f",
                                 eval(fdr) < 0.05 ~ "#b892ff", 
                                 eval(fdr) < 0.1 ~ "#29a655",
                                 eval(fdr) >= 0.1 ~ "#ff9100",
                                 is.na(eval(fdr)) ~ "grey"))
  
  if(input_pkg == "MAST") {
    
    if(lfcse_to_log2) {
      lfc_df <- lfc_df %>% 
        mutate(ci.hi = ci.hi / log(x = 2, base = exp(1)),
               ci.lo = ci.lo / log(x = 2, base = exp(1)))
    }
    
    lfc_df <- lfc_df %>% 
      mutate(lfcse = (ci.hi - ci.lo) / 3.92)
    
  }
  
  if(!is.null(key_genes)) {
    lfc_df <- lfc_df %>% 
      mutate(in_key_genes = ifelse(ensembl_gene_id_version %in% key_genes, 
                                   paste0("in ", key_genes_name),
                                   paste0("not in ", key_genes_name)))
  }
  
  out <- ggplot(data = lfc_df,
                mapping = aes(x = eval(lfc),
                              y = fct_reorder(gene_name, eval(lfc), .desc = desc_lfc),
                              fill = fdr_value)) +
    geom_hline(yintercept = lfc_df$gene_name,
               color = "grey", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#c9184a", linewidth = 0.5) +
    geom_vline(xintercept = xintercept, linetype = "dashed", color = "orange", linewidth = 0.5) +
    ggpubr::theme_pubr(legend = "bottom") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 14)) +
    scale_fill_manual(values = lfc_df %>% arrange(eval(fdr)) %>% pull(fdr_color) %>% unique()) +
    labs(title = title,
         x = paste0("log", utf8::utf8_print("\u2082"), " fold change (", test_cond, " vs ", ref_cond, ")"),
         y = NULL,
         fill = NULL,
         shape = NULL) +
    geom_linerange(mapping = aes(xmin = eval(lfc) - eval(lfcse),
                                 xmax = eval(lfc) + eval(lfcse))) 
  
  if(padj_lab) {
    out <- out +
      geom_text(mapping = aes(x = (eval(lfc) + eval(lfcse) + lab_nudge), 
                              label = ifelse(is.na(eval(fdr)), "", formatC(eval(fdr), format = "e", digits = 1))),
                hjust = 0)
  } 
  
  if(!is.null(xintercept)) {
    out <- out +
      scale_x_continuous(breaks = c(0, xintercept))
  }
  
  if(is.null(key_genes)) {
    out <- out +
      geom_point(shape = 21, color = "black", size = 5, stroke = 1)
  } else {
    if(length(unique(lfc_df$in_key_genes)) == 2) {
      shapes <- c(23, 21)
    } else if(unique(lfc_df$in_key_genes) == paste0("in ", key_genes_name)) {
      shapes <- 23
    } else {
      shapes <- 21
    }
    out <- out +
      geom_point(aes(shape = in_key_genes), color = "black", size = 5, stroke = 1) +
      scale_shape_manual(values = shapes) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))
  }
  out
}

# plot_pct_exp() ----------------------------------------------------------

plot_pct_exp <- function(lfc_df, 
                         genes, 
                         ref_cond, 
                         test_cond, 
                         title, 
                         fill_colors = c("#ca3767", "#5d76cb"), 
                         input_pkg = "MAST") {
  
  if(input_pkg == "MAST") {
    
    lfc <- sym("log2fc")
    
  } else if(input_pkg == "DESeq2") {
    
    lfc <- sym("log2FoldChange")
    
  }
  
  pct_exp_df <- lfc_df %>% 
    filter(ensembl_gene_id_version %in% genes, !is.na(eval(lfc))) 
  
  pct_exp_df <- pct_exp_df %>% 
    mutate(gene_name = fct_relevel(gene_name, 
                                   pct_exp_df %>% arrange(eval(lfc)) %>% pull(gene_name)),
           reference = 100,
           test = (2^eval(lfc))*100) %>% 
    select(ensembl_gene_id_version, gene_name, reference, test) %>% 
    pivot_longer(cols = c(reference, test), 
                 names_to = "condition", 
                 values_to = "pct_expression") %>% 
    mutate(condition = case_when(condition == "reference" ~ ref_cond,
                                 condition == "test" ~ test_cond,
                                 T ~ NA_character_) %>% 
             fct_relevel(c(ref_cond, test_cond)))
  
  if(mean(pct_exp_df$pct_expression) < 100) {
    out <- ggplot() +
      geom_col(data = pct_exp_df %>%
                 filter(condition == ref_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      geom_col(data = pct_exp_df %>%
                 filter(condition == test_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +  
      scale_y_continuous(limits = c(0L, 100L), n.breaks = 10L) 
  } else {
    out <- ggplot() +
      geom_col(data = pct_exp_df %>%
                 filter(condition == test_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      geom_col(data = pct_exp_df %>%
                 filter(condition == ref_cond),
               mapping = aes(x = gene_name,
                             y = pct_expression,
                             fill = condition),
               color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      scale_y_continuous(limits = c(0L, ceiling(max(pct_exp_df$pct_expression))), n.breaks = 12L) 
  }
  
  out <- out +
    ggpubr::theme_pubr(legend = "bottom") +
    scale_color_manual(values = fill_colors, aesthetics = "fill") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) +
    ggpubr::rotate_x_text(angle = 30) +
    labs(x = NULL,
         y = "% expression",
         title = title)
  
  out
}



# plot_norm_counts() ------------------------------------------------------

plot_norm_counts <- function(data, 
                             title, 
                             color = "precise_description", 
                             cont_colors = c("#fdab01", "#f85e00", "black"),
                             facet_var = NULL) {
  
  
  if(is.null(facet_var)) {
    split_counts <- F
  } else {
    split_counts <- T
    facet_var <- sym(facet_var)
  }
  
  if(split_counts) {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_norm_count,
                              y = gene_name,
                              fill = eval(facet_var),
                              color = eval(sym(color))))  +
      geom_col(width = 0.8, linewidth = 1) +
      geom_linerange(mapping = aes(xmin = mean_norm_count - sem_norm_count,
                                   xmax = mean_norm_count + sem_norm_count)) +
      geom_label(aes(x = 1.2, label = paste0(mean_norm_count_quartile, " - % ACTB: ", round(mean_norm_count_pct_actb, 2))), 
                 color = "black",
                 hjust = "left",
                 size = 3,
                 fill = "#e5e6e4") +
      scale_fill_manual(values = c("#0077b6", "#a4133c", "#d3d3d3")) +
      facet_wrap(~ eval(facet_var)) 
    
  } else {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_norm_count,
                              y = gene_name,
                              color = eval(sym(color))))  +
      geom_col(width = 0.8, linewidth = 1, fill = "#0077b6") +
      geom_linerange(mapping = aes(xmin = mean_norm_count - sem_norm_count,
                                   xmax = mean_norm_count + sem_norm_count)) +
      geom_label(aes(x = 1.2, label = paste0(mean_norm_count_quartile, " - % ACTB: ", round(mean_norm_count_pct_actb, 2))), 
                 color = "black",
                 hjust = "left",
                 size = 3,
                 fill = "#0077b6")
    
  }
  
  p +
    theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_x_log10(n.breaks = 12, labels = scales::comma) +
    scale_color_manual(values = cont_colors) +
    labs(x = "normalized count", 
         y = NULL,
         color = NULL,
         title = title) +
    guides(fill = "none",
           color = guide_legend(override.aes = list(fill = "white"))) +
    rotate_x_text(45) 
  
}

# sample_dist_heatmap() ---------------------------------------------------

sample_dist_heatmap <- function(data, labels, distance = "euclidian") {
  
  if(distance == "euclidian") {
    sampleDists <- dist(t(assay(data)))
  } else if(distance == "poisson") {
    sampleDists <- PoissonDistance(t(counts(data)))$dd
  }
  sampleDistMatrix <- as.matrix(sampleDists)
  
  if(!is.null(labels)) {
    if(length(labels) == 1) {
      rownames <- data[[labels]]
    } else {
      rownames <- data[[labels[[1]]]]
      for(i in labels[2:length(labels)]) {
        rownames <- paste0(rownames, " - ", data[[i]])
      }
    }
    rownames(sampleDistMatrix) <- rownames
  }
  colnames(sampleDistMatrix) <- NULL
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colorRampPalette( rev(brewer.pal(9, "Purples")) )(255))
  
}

# sample_dist_pca() -------------------------------------------------------

sample_dist_pca <- function(data, 
                            shape, 
                            fill,
                            shape_lab = NULL,
                            fill_lab = NULL,
                            shape_set = c(22, 21, 23), 
                            fill_color_set = NULL,
                            ntop = 500) {
  
  pca_data <- plotPCA(object = data,
                      intgroup = c(shape, fill),
                      ntop = ntop,
                      returnData = T,
                      pcsToUse = 1:2)
  
  pct_var <- attr(pca_data, "percentVar") * 100
  
  pca_plot <- ggplot(data = pca_data,
                     mapping = aes(x = PC1, 
                                   y = PC2,
                                   shape = eval(sym(shape)),
                                   fill = eval(sym(fill)))) +
    geom_point(size = 5, color = "black") +
    theme_pubr(legend = "right") +
    scale_shape_manual(values = shape_set) +
    guides(fill = guide_legend(override.aes = list(shape = shape_set[[1]]))) +
    labs(x = paste0("PC1: ", round(pct_var[[1]], digits = 1), "% variance"),
         y = paste0("PC2: ", round(pct_var[[2]], digits = 1), "% variance"),
         shape = shape_lab,
         fill = fill_lab)
  
  if(!is.null(fill_color_set)) {
    pca_plot <- pca_plot +
      scale_fill_manual(values = fill_color_set)
  }
  
  pca_plot
  
}

# run_pca() ---------------------------------------------------------------

run_pca <- function(cnt_mtx, ntop = 500, pcs = 1:2) {
  
  var <- rowVars(x = cnt_mtx)
  
  top_var <- order(var, decreasing = T)[seq_len(min(ntop, length(var)))]
  
  top_mtx <- t(cnt_mtx[top_var,])
  
  pca <- prcomp(x = top_mtx, center = T, scale. = T)
  
  pct_var <- (pca$sdev^2/sum(pca$sdev^2))*100
  
  pct_var <- setNames(object = pct_var[pcs], nm = paste0("PC", pcs))
  
  pcs <- names(pct_var)
  
  df <- as_tibble(x = pca$x[,pcs], rownames = "barcode")
  
  attr(df, "pct_var") <- pct_var
  
  return(df)
  
}

# plot_mean_tpms() --------------------------------------------------------

plot_mean_tpms <- function(data, 
                           title, 
                           x_label = "mean TPM",
                           add_geom_label = T,
                           color = "precise_protein_description", 
                           cont_colors = c("#5d76cb", "#ff9e00", "#29a655", "#ca3767", "#38a3a5"),
                           facet_var = NULL,
                           log_compatible = T) {
  
  if(log_compatible) {
    data <- data %>% 
      mutate(mean_tpm = mean_tpm + 1)
  }
  
  if(is.null(facet_var)) {
    split_counts <- F
  } else {
    split_counts <- T
    facet_var <- sym(facet_var)
  }
  
  if(split_counts) {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_tpm,
                              y = gene_name,
                              fill = eval(facet_var),
                              color = eval(sym(color))))  +
      geom_col(width = 0.9, linewidth = 1) +
      geom_linerange(mapping = aes(xmin = mean_tpm - sem_tpm,
                                   xmax = mean_tpm + sem_tpm))
    if(add_geom_label) {
      p <- p +
        geom_label(aes(x = 1.2, label = paste0(mean_tpm_quartile, " - % ACTB: ", round(mean_tpm_pct_actb, 3))), 
                   color = "black",
                   hjust = "left",
                   size = 4,
                   fill = "#e5e6e4")
    }
    p <- p +
      scale_fill_manual(values = c("#0077b6", "#a4133c", "#d3d3d3")) +
      facet_wrap(~ eval(facet_var)) 
    
  } else {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_tpm,
                              y = gene_name,
                              color = eval(sym(color))))  +
      geom_col(width = 0.9, linewidth = 1, fill = "#c0c0c0", alpha = 0.9) +
      geom_linerange(mapping = aes(xmin = mean_tpm - sem_tpm,
                                   xmax = mean_tpm + sem_tpm))
    if(add_geom_label) {
      p <- p +
        geom_label(aes(x = 1.2, label = paste0(mean_tpm_quartile, " - % ACTB: ", round(mean_tpm_pct_actb, 3))), 
                   color = "black",
                   hjust = "left",
                   size = 4,
                   fill = "#e5e6e4")
    }
    
  }
  
  p +
    theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          strip.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_x_log10(n.breaks = 8) +
    scale_color_manual(values = cont_colors) +
    labs(x = x_label, 
         y = NULL,
         color = NULL,
         title = title) +
    guides(fill = "none",
           color = guide_legend(override.aes = list(fill = "white")))
  
} 

# plot_mean_tpms2() -------------------------------------------------------

plot_mean_tpms2 <- function(data, 
                            title, 
                            x_label = "mean TPM",
                            add_geom_label = T,
                            include_pct_actb = F,
                            color = NULL, 
                            cont_colors = c("#5d76cb", "#ff9e00", "#29a655", "#ca3767", "#38a3a5"),
                            log_compatible = T) {
  
  if(log_compatible) {
    
    data <- data %>% 
      mutate(mean_tpm = mean_tpm + 1)
    
  }
  
  if(is.null(color)) {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_tpm,
                              y = gene_name))  +
      geom_col(width = 0.9, linewidth = 1, fill = "#5d76cb", alpha = 0.9, color = "black") +
      geom_linerange(mapping = aes(xmin = mean_tpm - sem_tpm,
                                   xmax = mean_tpm + sem_tpm))
    
  } else {
    
    p <- ggplot(data = data,
                mapping = aes(x = mean_tpm,
                              y = gene_name,
                              color = eval(sym(color))))  +
      geom_col(width = 0.9, linewidth = 1, fill = "#c0c0c0", alpha = 0.9) +
      geom_linerange(mapping = aes(xmin = mean_tpm - sem_tpm,
                                   xmax = mean_tpm + sem_tpm))
    
  }
  
  if(add_geom_label) {
    
    if(include_pct_actb) {
      
      p <- p +
        geom_label(aes(x = 1.2, 
                       label = paste0(mean_tpm_quartile, " - % ACTB: ", round(mean_tpm_pct_actb, 3))), 
                   color = "black",
                   hjust = "left",
                   size = 4,
                   fill = "#e5e6e4")
      
    } else {
      
      p <- p +
        geom_label(aes(x = 1.2, 
                       label = mean_tpm_quartile), 
                   color = "black",
                   hjust = "left",
                   size = 4,
                   fill = "#e5e6e4")
      
    }
    
  }
  
  p +
    theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          strip.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_x_log10(n.breaks = 8) +
    scale_color_manual(values = cont_colors) +
    labs(x = x_label, 
         y = NULL,
         color = NULL,
         title = title) +
    guides(fill = "none",
           color = guide_legend(override.aes = list(fill = "white")))
  
} 

# set_level_color() -------------------------------------------------------

set_level_color <- function(data, column, levels, colors) {
  
  factor_color <- tibble(factor = as_factor(x = levels),
                         color = as_factor(x = colors))
  
  filtered_factor_color <- factor_color %>% 
    dplyr::filter(factor %in% unique(data[[column]])) %>% 
    mutate(across(.cols = everything(), .fns = fct_drop))
  
  data <- data %>% 
    mutate(!!column := fct_relevel(eval(sym(column)), levels(filtered_factor_color$factor))) %>%
    left_join(filtered_factor_color, by = join_by(!!sym(column) == factor)) %>% 
    dplyr::rename(!!paste0(column, "_color") := color)
  
}

# split_splici_txi() ------------------------------------------------------

split_splici_txi <- function(txi, split_df) {
  
  out <- lapply(
    X = c("abundance", "counts", "length"),
    FUN = function(layer) {
      
      se <- SummarizedExperiment::SummarizedExperiment(assays = setNames(object = list(txi[[layer]]), nm = layer))
      
      txis <- tximeta::splitSE(se = se, splitDf = split_df, assayName = layer)
      
      setNames(object = list(SummarizedExperiment::assay(txis, colnames(split_df)[[1]]),
                             SummarizedExperiment::assay(txis, colnames(split_df)[[2]])),
               nm = c(colnames(split_df)[[1]], colnames(split_df)[[2]]))
      
    }) 
  
  out <- setNames(object = c(out, txi$countsFromAbundance), 
                  nm = c("abundance", "counts", "length", "countsFromAbundance"))
  
}