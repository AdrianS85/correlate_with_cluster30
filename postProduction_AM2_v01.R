clust$python_lists <- purrr::map(
  .x = list("head" = ".*0\\.tsv", "placenta" = ".*0_p\\.tsv"),
  .f = function(tissue){
    
    purrr::map(
      .x = list.files(path = "python", pattern = tissue),
      .f = function(name_){
        
        df_ <- readr::read_tsv(file = paste0("python/", name_))
        
        df_$X1 <- NULL
        
        name__ <- stringr::str_remove(name_, "data_")
        name__ <- stringr::str_remove(name__, "_p\\.tsv")
        name__ <- stringr::str_remove(name__, "\\.tsv")
        
        colnames(df_)[1] <- name__
        
        df_ <- dplyr::select(df_, symbols_probes, dplyr::matches(match = "sorted.*"))
        
        df_
      })
  })



clust$python <- purrr::map(
  .x = clust$python_lists,
  .f = function(tissue){
    
    purrr::reduce(.x = tissue, .f = merge, by = "symbols_probes")
  })



clust$python_ <- purrr::map2(
  .x = clust$python,
  .y = list(clust$only_h, clust$only_p),
  .f = function(tissue_clust, tissue_input){
    
    merge(x = tissue_clust, y = tissue_input, by = "symbols_probes")
  })



### !!! Here we need to absolutize values. Reason: If you wont do it, then genes of reverse expression patterns will sum up to zero instead of providing actual median value. Also: we already work on clusters, which we know have similar expression or reverse. Therefore, we can assume that if we take absolute value of any other gene, it will have similar expression to any two genes.
clust$python_medians <- purrr::map(
  .x = clust$python_,
  .f = function(tissue){
    
    
    purrr::map(
      .x = list("sorted1000" = "sorted1000", "sorted2500" = "sorted2500", "sorted5000" = "sorted5000", "sorted10000" = "sorted10000"),
      .f = function(col_){
        
        df_ <- dplyr::select(.data = tissue, symbols_probes, rlang::sym(col_), dplyr::matches(".*ma.*"))
        
        df_nest <- tidyr::nest( dplyr::group_by(df_, !!rlang::sym(col_)) )
        
        
        
        medians_ <- purrr::map(
          .x = df_nest$data,
          .f = function(cluster_){
            
            cluster_ <- dplyr::select(cluster_, -symbols_probes)
            
            median_ <- numeric()
            for (colNb in seq(1, ncol(cluster_))) {
              
              median_[colNb] <- median(abs(cluster_[[colNb]])) ###!!! The input values are z-scores, so they are CENTERED. Here I order them also to be absolute
            }
            return( as.data.frame(t(as.data.frame(median_))))
          })
        
        medians_ <- rlist::list.rbind(medians_)
        
        colnames(medians_) <- colnames(df_)[3:17]
        
        rownames(medians_) <- df_nest[[1]]
        
        cbind(medians_, df_nest)
        
        return(cbind(medians_, df_nest))
      })
  })


### !!! Need better corelation method, as the current one tends to put too much enphasis on single large peaks
clust$python_cor_matrix <- purrr::map2(
  .x = clust$python_medians$head,
  .y = clust$python_medians$placenta,
  .f = function(ds_h, ds_p){
    
    ds_h <- ds_h[order(ds_h[[16]]),]
    ds_p <- ds_p[order(ds_p[[16]]),]
    
    assertthat::are_equal(colnames(ds_h)[1:15], colnames(ds_p)[1:15])
    
    # WRS2::winall(x = t(ds_h[,1:15]), y = t(ds_p[,1:15])) ### clust_wrs.save
    
    corrr::correlate(
      x = t(ds_h[,1:15]),
      y = t(ds_p[,1:15]),
      method = "spearman")
  })

# save(clust, file = "~/clust_wrs.save")
# load("~/clust_wrs.save")
# save(clust, file = "~/clust.save")
# load("~/clust.save")
# save(clust, file = "~/clust_spear.save")
# load("~/clust_spear.save")

### !!! FOR: corrr::correlate
clust$python_cor_matrix_long <- purrr::map(
  .x = clust$python_cor_matrix,
  .f = function(df_){

    temp1 <- df_

    temp1$h_clusters <- temp1$term

    temp1$term <- NULL

    temp1 <- tidyr::pivot_longer(data = temp1, cols = -h_clusters, names_to = "p_clusters", values_to = "corelation")

    dplyr::filter(temp1, abs(corelation) >= 0.9)
  })

clust$python_cor_matrix_long$sorted5000 <- dplyr::filter(clust$python_cor_matrix_long$sorted5000, abs(corelation) >= 0.92)
clust$python_cor_matrix_long$sorted10000 <- dplyr::filter(clust$python_cor_matrix_long$sorted10000, abs(corelation) >= 0.94)

### !!! FOR: WRS2::winall
# clust$python_cor_matrix_long <- purrr::map(
#   .x = clust$python_cor_matrix,
#   .f = function(df_){
# 
#     temp1 <- as.data.frame(df_$cor)
#     
#     temp1[temp1 == 1] <- NA ### !!! może to pomoże
# 
#     temp1$h_clusters <- rownames(temp1)
# 
#     temp1 <- tidyr::pivot_longer(data = temp1, cols = -h_clusters, names_to = "p_clusters", values_to = "corelation")
# 
#     dplyr::filter(temp1, abs(corelation) >= 0.95)
#   })



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
