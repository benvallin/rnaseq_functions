library(fgsea)
library(MAST)
library(limma)
library(tidyverse)

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
                              n_leadingEdge_over_total = paste0(n_leadingEdge, " / ", size)) %>% 
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
                              n_overlapGenes = paste0(overlap, " / ", size)) %>% 
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
                     test_cond = "PFF",
                     ref_cond = "NT",
                     key_genes_name = "leading edge",
                     desc_lfc = T,
                     xintercept = NULL,
                     input_pkg = "MAST") {
  
  
  if(input_pkg == "MAST") {
    
    lfc <- sym("log2fc")
    fdr <- sym("fdr")
    
  } else if(input_pkg == "DESeq2") {
    
    lfc <- sym("log2FoldChange")
    fdr <- sym("padj")
    
  }
  
  
  lfc_df <- lfc_df %>% 
    filter(ensembl_gene_id_version %in% genes,
           !is.na(eval(lfc))) %>% 
    mutate(fdr_value = case_when(eval(fdr) < 0.01 ~ "FDR < 0.01",
                                 eval(fdr) < 0.05 ~ "FDR < 0.05", 
                                 T ~ paste0("FDR ", utf8::utf8_print("\u2265"), " 0.05")),
           fdr_color = case_when(eval(fdr) < 0.01 ~ "#e0007f",
                                 eval(fdr) < 0.05 ~ "#b892ff", 
                                 T ~ "#ff9100")) 
  
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