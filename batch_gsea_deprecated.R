library(tidyverse)
library(ggpubr)

# batch_ensembl_gene_id_annotations() -------------------------------------------------------------------------

batch_ensembl_gene_id_annotations <- function(res_df,
                                              gene_metadata,
                                              annotation_categories,
                                              keep_only_detected_genes = T) {
  
  lapply(X = annotation_categories,
         FUN = function(x) {
           
           # From gene_metadata, filter those genes that:
           # - are detected in the dataset res_df if keep_only_detected_genes = T
           # - have non-null annotation for the current gene set collection
           output <- as_tibble(gene_metadata) 
           if(keep_only_detected_genes) {
             output <- output %>% 
               dplyr::filter(ensembl_gene_id %in% res_df$ensembl_gene_id)
           }
           output <- output %>% 
             dplyr::select(ensembl_gene_id, all_of(x)) %>% 
             unnest(all_of(x), keep_empty = F)
         }) %>% 
    setNames(nm = annotation_categories)
  
}

# batch_ora_fgsea() -------------------------------------------

batch_ora_fgsea <- function(test_gene_set, 
                            ensembl_gene_id_annotations,
                            minSize = 10,
                            maxSize = 500,
                            p.adjust_threshold = Inf,
                            ...) {
  
  lapply(X = names(ensembl_gene_id_annotations),
         FUN = function(x) {
           
           # Define pathways
           # => All the set IDs in the current gene set collection 
           # that match one or more genes listed in ensembl_gene_id_annotations
           pathways_table <- ensembl_gene_id_annotations[[x]] %>% 
             nest(data = ensembl_gene_id) %>% 
             dplyr::mutate(data = map(.x = data, .f = ~.x$ensembl_gene_id))
           
           pathways <- as.list(setNames(object = pathways_table[["data"]],
                                        nm = pathways_table[["set_id"]]))
           
           # Define set of genes to test
           # => Genes in test_gene_set that have non-null annotation for the current gene set collection
           genes <- test_gene_set %>% 
             dplyr::filter(!map_lgl(.x = eval(sym(x)), .f = is.null)) %>% 
             pull(ensembl_gene_id)
           
           # Define universe
           # => All detected genes that have non-null annotation for the current gene set collection
           universe <- unique(ensembl_gene_id_annotations[[x]]$ensembl_gene_id)
           
           # Run fgsea::fora()
           fgsea::fora(pathways = pathways,
                       genes = genes,
                       universe = universe,
                       minSize = minSize,
                       maxSize = maxSize,
                       ...) %>% 
             # Tidy final output
             as_tibble() %>% 
             # Discard annotation terms that have no match in the set of tested genes and recompute adjusted pvalue
             dplyr::filter(overlap > 0) %>% 
             dplyr::mutate(padj = p.adjust(p = pval, method = "fdr")) %>% 
             # Bind annotation gene set name / description
             dplyr::rename(set_id = pathway) %>% 
             left_join(pathways_table[, c("set_id", "set_name")], 
                       by = join_by(set_id)) %>%
             dplyr::select(set_id, set_name, overlap, everything()) %>% 
             dplyr::filter(padj < p.adjust_threshold) %>% 
             mutate(n_gene_over_total = overlap / size,
                    n_gene_over_total_label = paste0(overlap, " / ", size))
           
         }) %>% 
    setNames(nm = names(ensembl_gene_id_annotations))
  
}

# batch_ora_clusterprofiler() -------------------------------------------

batch_ora_clusterprofiler <- function(test_gene_set, 
                                      ensembl_gene_id_annotations,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      pvalueCutoff = Inf,
                                      qvalueCutoff = Inf,
                                      p.adjust_threshold = Inf,
                                      ...) {
  
  out <- lapply(X = names(ensembl_gene_id_annotations),
                FUN = function(x) {
                  
                  # Define set of genes to test
                  # => Genes in test_gene_set that have non-null annotation for the current annotation category
                  gene <- test_gene_set %>% 
                    dplyr::filter(!map_lgl(.x = eval(sym(x)), .f = is.null)) %>% 
                    pull(ensembl_gene_id)
                  
                  # Define universe
                  # => All detected genes that have non-null annotation for the current annotation category
                  universe <- unique(ensembl_gene_id_annotations[[x]]$ensembl_gene_id)
                  
                  # Define TERM2GENE and TERM2NAME
                  # => Tables providing gene and annotation information only for the detected genes 
                  # which have non-null annotation for the current annotation category
                  # => TERM2GENE: column 1 = annotation term ID and column 2 = matching gene
                  TERM2GENE <- ensembl_gene_id_annotations[[x]][, c("set_id", "ensembl_gene_id")]
                  # => TERM2NAME: column 1 = annotation term ID and column 2 = annotation term name / description
                  TERM2NAME <- ensembl_gene_id_annotations[[x]][, c("set_id", "set_name")]
                  
                  # Run clusterProfiler::enricher()
                  res <- clusterProfiler::enricher(gene = gene,
                                                   universe = universe,
                                                   TERM2GENE = TERM2GENE,
                                                   TERM2NAME = TERM2NAME,
                                                   minGSSize = minGSSize,
                                                   maxGSSize = maxGSSize,
                                                   pvalueCutoff = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff,
                                                   ...) 
                  # Tidy final output if exists
                  if(!is.null(res)) {
                    res <- res %>% 
                      as_tibble() %>%
                      dplyr::rename(set_id = ID,
                                    set_name = Description) %>%
                      dplyr::select(set_id, set_name, Count, everything()) %>%
                      dplyr::filter(p.adjust < p.adjust_threshold) %>% 
                      mutate(overlap = Count,
                             size = as.integer(str_remove(BgRatio, "/.*$")),
                             n_gene_over_total = overlap / size,
                             n_gene_over_total_label = paste0(overlap, " / ", size))
                  }
                }) %>%
    setNames(nm = names(ensembl_gene_id_annotations))
  
  keep <- map_lgl(.x = out, .f = ~ !is.null(.x))
  out <- out[keep]
  
}

# batch_ora_goseq() -------------------------------------------

batch_ora_goseq <- function(test_gene_set,
                            ensembl_gene_id_annotations,
                            res_df, 
                            transcript_length,
                            p.adjust_threshold = Inf,
                            ...) {
  
  lapply(X = names(ensembl_gene_id_annotations),
         FUN = function(x) {
           
           # Define DEgenes
           # => Named vector of detected genes (0 = not DEG; 1 = DEG; names = ensembl_gene_id)
           # => Restricted to the detected genes that have non-null annotation for the current annotation category
           DEgenes <- res_df %>% 
             dplyr::filter(!map_lgl(.x = eval(sym(x)), .f = is.null)) %>%
             pull(ensembl_gene_id)
           
           DEgenes <- setNames(object = as.integer(DEgenes %in% test_gene_set$ensembl_gene_id),
                               nm = DEgenes)
           
           DEgenes <- DEgenes[order(names(DEgenes))]
           
           # Refine transcript_length
           # => Restrict to the detected genes that have non-null annotation for the current annotation category
           annotation_specific_transcript_length <- transcript_length[names(transcript_length) %in% names(DEgenes)]
           
           # Compute the probability weighting function for the set of detected genes based on the median transcript length
           pwf <- goseq::nullp(DEgenes = DEgenes,
                               bias.data = annotation_specific_transcript_length)
           
           # Define gene2cat
           # => List providing gene and annotation information only for the detected genes 
           # which have non-null annotation for the current annotation category
           # => Each entry is an annotation term ID with name = matching ensembl_gene_id
           gene2cat <- as.list(setNames(object = ensembl_gene_id_annotations[[x]][["set_id"]],
                                        nm = ensembl_gene_id_annotations[[x]][["ensembl_gene_id"]]))
           
           # Run goseq::goseq()
           goseq::goseq(pwf = pwf,
                        gene2cat = gene2cat,
                        ...) %>%
             # Tidy final output
             as_tibble() %>% 
             # Discard annotation terms that have no match in the set of tested genes and compute adjusted pvalues
             dplyr::filter(numDEInCat > 0) %>% 
             dplyr::mutate(over_represented_padj = p.adjust(over_represented_pvalue, method = "BH"),
                           under_represented_padj = p.adjust(under_represented_pvalue, method = "BH")) %>%
             # Bind annotation term name / description to the final output
             dplyr::rename(set_id = category) %>% 
             left_join(unique(ensembl_gene_id_annotations[[x]][, c("set_id", "set_name")]),
                       by = join_by(set_id)) %>% 
             dplyr::select(set_id, set_name, numDEInCat, everything(), -any_of("term")) %>%
             dplyr::arrange(over_represented_padj, desc(numDEInCat), set_id) %>% 
             dplyr::filter(over_represented_padj < p.adjust_threshold)
           
         }) %>% 
    setNames(nm = names(ensembl_gene_id_annotations))
  
}

# batch_ora() -------------------------------------------------------------

batch_ora <- function(test_gene_set,
                      ensembl_gene_id_annotations,
                      method = "fgsea::fora",
                      res_df = NULL, 
                      transcript_length = NULL,
                      p.adjust_threshold = Inf,
                      ...) {
  
  if (!method %in% c("fgsea::fora", "clusterProfiler::enricher", "goseq::goseq")) {
    stop("\nMethod must be \"fgsea::fora\", \"clusterProfiler::enricher\", or \"goseq::goseq\".\n")
  }
  
  if (method == "goseq::goseq" && (is.null(res_df) || is.null(transcript_length))) {
    stop("\nMethod \"goseq::goseq\" requires \"res_df\" and \"transcript_length\" to be provided.\n")
  }
  
  if (method != "goseq::goseq" && (!is.null(res_df) || !is.null(transcript_length))) {
    message("\nThe \"res_df\" and \"transcript_length\" arguments are ignored by method \"", method, "\".\n")
  }
  
  switch(EXPR = method,
         "fgsea::fora" = batch_ora_fgsea(test_gene_set = test_gene_set,
                                         ensembl_gene_id_annotations = ensembl_gene_id_annotations,
                                         p.adjust_threshold = p.adjust_threshold,
                                         ...),
         "clusterProfiler::enricher" = batch_ora_clusterprofiler(test_gene_set = test_gene_set,
                                                                 ensembl_gene_id_annotations = ensembl_gene_id_annotations,
                                                                 p.adjust_threshold = p.adjust_threshold,
                                                                 ...),
         "goseq::goseq" = batch_ora_goseq(test_gene_set = test_gene_set,
                                          ensembl_gene_id_annotations = ensembl_gene_id_annotations,
                                          res_df = res_df, 
                                          transcript_length = transcript_length,
                                          p.adjust_threshold = p.adjust_threshold,
                                          ...))
  
}

# batch_gsea_fgsea() -------------------------------------------

batch_gsea_fgsea <- function(ranked_genes, 
                             ensembl_gene_id_annotations,
                             minSize = 10,
                             maxSize = 500,
                             p.adjust_threshold = Inf,
                             ...) {
  
  lapply(X = names(ensembl_gene_id_annotations),
         FUN = function(x) {
           
           # Define pathways
           # => All the set IDs in the current gene set collection 
           # that match one or more genes listed in ensembl_gene_id_annotations
           pathways_table <- ensembl_gene_id_annotations[[x]] %>% 
             nest(data = ensembl_gene_id) %>% 
             dplyr::mutate(data = map(.x = data, .f = ~.x$ensembl_gene_id))
           
           pathways <- as.list(setNames(object = pathways_table[["data"]],
                                        nm = pathways_table[["set_id"]]))
           
           # Run fgsea::fgsea()
           output <- fgsea::fgsea(pathways = pathways,
                        stats = ranked_genes,
                        minSize = minSize,
                        maxSize = maxSize,
                        ...)
           
           # Tidy final output 
           output <- as_tibble(output) %>%
             dplyr::rename(set_id = pathway) %>% 
             # Bind annotation gene set name / description
             left_join(pathways_table[, c("set_id", "set_name")], 
                       by = join_by(set_id)) %>% 
             dplyr::filter(padj < p.adjust_threshold) %>% 
             dplyr::mutate(n_leading_edge = map2(.x = leadingEdge,
                                                 .y = size,
                                                 .f = ~ paste0(length(.x), " / ", .y))) 
           
         }) %>% 
    setNames(nm = names(ensembl_gene_id_annotations))
  
}

# batch_gsea_clusterprofiler() -------------------------------------------

batch_gsea_clusterprofiler <- function(ranked_genes, 
                                       ensembl_gene_id_annotations,
                                       minGSSize = 10,
                                       maxGSSize = 500,
                                       pvalueCutoff = Inf,
                                       p.adjust_threshold = Inf,
                                       ...) {
  
  out <- lapply(X = names(ensembl_gene_id_annotations),
                FUN = function(x) {
                  
                  # Define TERM2GENE and TERM2NAME
                  # => Tables providing gene and annotation information only for the detected genes 
                  # which have non-null annotation for the current annotation category
                  # => TERM2GENE: column 1 = annotation term ID and column 2 = matching gene
                  TERM2GENE <- ensembl_gene_id_annotations[[x]][, c("set_id", "ensembl_gene_id")]
                  # => TERM2NAME: column 1 = annotation term ID and column 2 = annotation term name / description
                  TERM2NAME <- ensembl_gene_id_annotations[[x]][, c("set_id", "set_name")]
                  
                  # Run clusterProfiler::GSEA()
                  res <- clusterProfiler::GSEA(geneList = ranked_genes,
                                               TERM2GENE = TERM2GENE,
                                               TERM2NAME = TERM2NAME,
                                               by = "fgsea",
                                               minGSSize = minGSSize,
                                               maxGSSize = maxGSSize,
                                               pvalueCutoff = pvalueCutoff,
                                               ...)  
                  # Tidy final output if exists
                  if(!is.null(res)) {
                    res <- res %>% 
                      as_tibble() %>%
                      dplyr::rename(set_id = ID,
                                    set_name = Description) %>% 
                      dplyr::filter(p.adjust < p.adjust_threshold) 
                    
                    if("core_enrichment" %in% colnames(res)) {
                      res <- res %>% 
                        dplyr::mutate(padj = p.adjust,
                                      size = setSize,
                                      leadingEdge = str_split(string = core_enrichment, pattern = "/"),
                                      n_leading_edge = map2(.x = leadingEdge,
                                                            .y = size,
                                                            .f = ~ paste0(length(.x), " / ", .y)))
                    }
                  }
                }) %>% 
    setNames(nm = names(ensembl_gene_id_annotations))
  
  keep <- map_lgl(.x = out, .f = ~ !is.null(.x))
  out <- out[keep]
  
}

# batch_gsea() -------------------------------------------

batch_gsea <- function(ranked_genes, 
                       ensembl_gene_id_annotations,
                       method = "fgsea::fgsea",
                       p.adjust_threshold = Inf,
                       ...) {
  
  if (!method %in% c("fgsea::fgsea", "clusterProfiler::GSEA")) {
    stop("\nMethod must be \"fgsea::fgsea\" or \"clusterProfiler::GSEA\".\n")
  }
  
  switch(EXPR = method,
         "fgsea::fgsea" = batch_gsea_fgsea(ranked_genes = ranked_genes, 
                                           ensembl_gene_id_annotations = ensembl_gene_id_annotations,
                                           p.adjust_threshold = p.adjust_threshold,
                                           ...),
         "clusterProfiler::GSEA" = batch_gsea_clusterprofiler(ranked_genes = ranked_genes, 
                                                              ensembl_gene_id_annotations = ensembl_gene_id_annotations,
                                                              p.adjust_threshold = p.adjust_threshold,
                                                              ...))
  
}

# summarise_batch_enrichment() ----------------------------------------------

summarise_batch_enrichment <- function(batch_enrichment_output) {
  
  out <- map2(.x = batch_enrichment_output,
              .y = names(batch_enrichment_output),
              .f = ~ mutate(.x, category = .y)) %>% 
    bind_rows() %>% 
    dplyr::select(category, set_id, set_name, everything()) %>% 
    arrange(category, padj)
}

# plot_gsea_results() -----------------------------------------------------

plot_gsea_results <- function(gsea_df, category, title) {
  
  ggplot(data = gsea_df[gsea_df$category == category,],
         mapping = aes(x = NES, y = fct_reorder(set_name, desc(NES)))) +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed", color = "orange", linewidth = 0.5) +
    geom_point(aes(fill = padj), size = 6, shape = 21) +
    geom_label(aes(label = n_leading_edge), size = 3, nudge_x = 0.3) +
    theme_pubr(legend = "right") +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 12)) +
    scale_fill_gradient2(guide = guide_colourbar(frame.colour = "black")) +
    labs(x = "NES", y = NULL, title = title, fill = "adj p")
  
}


# plot_ora_results() ------------------------------------------------------

plot_ora_results <- function(ora_df, category, title) {
  
  ggplot(data = ora_df[ora_df$category == category,]) +
    geom_col(mapping = aes(x = size, y = fct_reorder(set_name, size)),
             color = "black", fill = "grey", position = "dodge", width = 0.9, alpha = 0.5) +
    geom_col(mapping = aes(x = overlap, y = set_name, fill = p.adjust),
             color = "black", position = "dodge", width = 0.9, alpha = 0.5) +
    geom_label(aes(x = size, y = set_name, label = n_gene_over_total_label), size = 4, nudge_x = -50) +
    theme_pubr(legend = "right") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18)) +
    rotate_x_text(angle = 30) +
    scale_fill_gradient2(guide = guide_colourbar(frame.colour = "black")) +
    labs(x = NULL, y = NULL, title = title, fill = "adj p")
  
}

