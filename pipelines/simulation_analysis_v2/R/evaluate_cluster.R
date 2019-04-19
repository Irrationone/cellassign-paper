map_and_evaluate_clusters <- function(sce_sim, min_count = 0) {
  # Remove genes that are lowly expressed, and cells with no UMIs
  sce_sim <- sce_sim[rowSums(counts(sce_sim)) >= min_count,colSums(counts(sce_sim)) > 0]
  
  nclusts <- length(unique(sce_sim$cluster))
  
  if (nclusts > 1 & !all(is.na(sce_sim$cluster))) {
    sce_de <- map_clusters(sce_sim, method = "de2", FDR_cutoff = 0.05)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0)
  } else {
    max_group <- names(which.max(table(sce_sim$Group)))[1]
    
    sce_de <- sce %>% scater::mutate(inferred_group=max_group)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0)
  }
  
  de_evaluation_measures <- compute_evaluation_measures(sce_de, 
                                                        truth_labels = sce_de$Group,
                                                        inferred_labels = sce_de$inferred_group)
  corr_evaluation_measures <- compute_evaluation_measures(sce_correlation, 
                                                          truth_labels = sce_correlation$Group,
                                                          inferred_labels = sce_correlation$inferred_group)
  
  evaluation_measures <- rbind(data.frame(mapping_type = "de", de_evaluation_measures, stringsAsFactors = FALSE), data.frame(mapping_type = "correlation", corr_evaluation_measures, stringsAsFactors = FALSE))
  
  return(evaluation_measures)
}