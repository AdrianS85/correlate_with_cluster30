clust$python_cor_matrix_long_clusters <- purrr::pmap(
  .l = list(clust$python_cor_matrix_long, clust$python_medians$head, clust$python_medians$placenta),
  .f = function(sorted_, h_medians, p_medians){
    
    sorted__ <- purrr::map2(
      .x = sorted_$h_clusters,
      .y = sorted_$p_clusters,
      .f = function(h_clust, p_clust){
        
        df_ <- rlist::list.rbind( list(
          "head" = dplyr::filter(h_medians, h_medians[[colnames(h_medians)[grepl(".*sorted.*", colnames(h_medians))] ]] == h_clust),
          "placenta" = dplyr::filter(p_medians, p_medians[[colnames(p_medians)[grepl(".*sorted.*", colnames(p_medians))] ]] == p_clust ) ))
        
        df_$cluster_pair <- paste0("H", h_clust, "_P", p_clust)
        
        return(df_)
      })
  }) 


clust$python_cor_matrix_long_clusters_df <- purrr::map(
  .x = clust$python_cor_matrix_long_clusters,
  .f = function(list_){
    
    rlist::list.rbind(list_)
  })



clust$python_cor_matrix_long_clusters_df_for_ggplot <- purrr::map(
  .x = clust$python_cor_matrix_long_clusters_df,
  .f = function(df_){
    
    df_$tissue <- rownames(df_)
    
    df_$data <- NULL
    
    tidyr::pivot_longer(
      data = df_,
      cols = !dplyr::matches(".*sorted.*|tissue|cluster_pair"),)
  })



library(ggplot2)
clust$ploted_python_cor_matrix_long_clusters_df_for_ggplot <- purrr::map(
  .x = clust$python_cor_matrix_long_clusters_df_for_ggplot,
  .f = function(df_){
    
    
    purrr::map(
      .x = unique(df_$cluster_pair),
      .f = function(pair_){
        
        to_plot <- dplyr::filter(df_, cluster_pair == pair_)
        
        ggplot(data = to_plot, mapping = aes(x = name, y = value, group = tissue, color = tissue))+
          geom_line() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      })
  })



dir.create("corelated_clusters")
for (dfNb in seq_along(clust$ploted_python_cor_matrix_long_clusters_df_for_ggplot)) {
  
  dir_name <- paste0("corelated_clusters/", names(clust$ploted_python_cor_matrix_long_clusters_df_for_ggplot)[dfNb])
  
  dir.create(dir_name)
  
  for (clust_ in clust$ploted_python_cor_matrix_long_clusters_df_for_ggplot[[dfNb]]) {
    
    ggsave(
      filename = paste0(dir_name, "/", clust_$data$cluster_pair[[1]], ".png"),
      plot = clust_,
      device = "png",
      width = 297,
      height = 210,
      units = "mm")
  }
}
rm(clust_, dir_name)



clust$gpl21163_gemma <- readr::read_tsv(file = "GPL21163_noParents.an.txt", skip = 7)

clust$gpl21163_gemma$GeneSymbols <- stringr::str_remove(
  string = clust$gpl21163_gemma$GeneSymbols,
  pattern = "\\|.*")

clust$gpl21163_gemma$NCBIids <- stringr::str_remove(
  string = clust$gpl21163_gemma$NCBIids,
  pattern = "\\|.*")



dir.create("corelated_clusters")

clust$gene_data$wb <- purrr::map(
  .x = list("1000" = NA, "2500" = NA, "5000" = NA, "10000" = NA),
  .f = function(x) {openxlsx::createWorkbook()}) 
