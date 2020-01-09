#library(shiny)
#library(ggplot2)
#library(corrplot)
#library(reshape2)
#library(dplyr)
#source("Data_preparation.R")
#source("Clusterisation.R")
#source("Gene_expression.R")
#source("Correlation.R")

#' Creation of server part of Shiny App
#' @description the function forms the body of backend part of Shiny App.
#' It also shapes the architecture of the pipeline. there are pipeline
#' blocks here which reflect the steps of CyTOF analysis. All computational
#' functions are placed in other files. The blocks of this function request
#' computational functions from that files by suitable order.
#' @param input shiny input object with from shiny UI to server
#' @param output shiny output object from server to UI and back to server
#'
#' @return
#'
#' @import shiny shinyFiles visNetwork d3heatmap ggplot2 corrplot reshape2 dplyr
#' @examples
cytofCore_server <- function(input, output){

  ########################
  ### Data preparation ###
  ########################
  roots <- c(home = path.expand("~"))
  shinyFileChoose(input, 'fcs_files', roots=roots, filetypes=c('', 'fcs'))

  ##### Create "fcs_data" as reactive object to store the CyTOF data
  fcs_data <-reactiveValues()
  observeEvent(input$fcs_upload, {
    ## Get row data fcs files
    fcs_data$md <- get_fcs_metadata(parseFilePaths(roots, input$fcs_files)$datapath)
    fcs_data$fcs_raw <- get_fcs_raw(fcs_data$md)
    fcs_data$panel <- get_fcs_panel(fcs_data$fcs_raw)
    fcs_data$use_markers <- get_use_marker(fcs_data$panel)
    fcs_data$cell_number <- get_cell_number(fcs_data$fcs_raw)
    ## Transform row data to scaled data by set parameters
    if('asinh' %in% input$transformation_list){
      fcs_data$fcs_raw <- asinh_transformation(fcs_data$fcs_raw, input$cofactor)
    }
    if('outlier_by_quantile' %in% isolate(input$transformation_list)){
      fcs_data$fcs_raw <- outlier_by_quantile_transformation(fcs_data$fcs_raw, input$quantile)
    }
    ## Preparing data for heatmap
    sampling_size <- 0.5
    method <- "tSNE"
    if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
    if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
    fcs_data$tSNE <- sampled_tSNE(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method)
  })

  ##### Drawing the reactive tSNE plot
  output$tSNE_plot_data_preparation <- renderPlot({
    if(is.null(fcs_data$tSNE)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_data_preparation)){color_mk <- input$mk_data_preparation}
    ggplot(fcs_data$tSNE,  aes(x = tSNE1, y = tSNE2, color = eval(parse(text = color_mk)))) +
      geom_point(size = 0.2) +
      scale_color_gradient2(midpoint = 0.5, low = 'blue', mid = "white",  high = 'red') +
      labs(color = color_mk) +
      theme_bw()
  })

  ##### Create UI to choose target marker
  output$mk_target_data_preparation_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput("mk_data_preparation", label = h3("Ploted marker"),
                choices = names(fcs_data$use_markers),
                #choices = c(1,2,3),
                selected = 1)
  })

  ##### Create UI to choose excluded markers
  output$mk_subset_data_preparation_ui <- renderUI({
    color_mk <- names(fcs_data$use_markers)[1]
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput("exclude_mk_data_preparation", label = "exclude markers",
                choices = names(fcs_data$use_markers),
                multiple = TRUE)
  })

  ##### Redrawing plot after chanch number of draw cells ar methods
  observeEvent(input$redraw, {
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    sampling_size <- 0.5
    method <- "tSNE"
    if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
    if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
    fcs_data$tSNE <- sampled_tSNE(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method)
  })

  ##### Update reactive object "fcs_data" after excluding markers
  observeEvent(input$exclud_mk_button, {
    sampling_size <- 0.5
    method <- "tSNE"
    if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
    if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
    fcs_data$use_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_data_preparation)]
    fcs_data$tSNE <- sampled_tSNE(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method)
  })

  ##### Reactively show current use_markers from "fcs_data" object
  output$mk_rested_data_preparation <- renderPrint({
    fcs_data$use_markers
  })

  ##### Drawing the reactive histogram plot of marker expression
  output$mk_hist_data_preparation <- renderPlot({
    if(is.null(fcs_data$tSNE)){return(NULL)}
    color_mk <- "Marker was not choose"
    if(!is.null(input$mk_data_preparation)){color_mk <- input$mk_data_preparation}
    ggplot(fcs_data$tSNE, aes(x = eval(parse(text = color_mk)), y=..scaled..)) +
      geom_density(fill = 'black') +
      labs(x = color_mk)
  })

  ##### Drawing the reactive plot of cell number
  output$smpl_hist_preparation <- renderPlot({
    if(is.null(fcs_data$cell_number)){return(NULL)}
    ggplot(data = fcs_data$cell_number, aes(x = smpl, y = cell_number))+
      geom_bar(stat="identity") +
      coord_flip()
  })

  ########################
  ###  Clusterisation  ###
  ########################

  ##### Create "clusterisation" as reactive object to store the cluster information
  clusterisation <- reactiveValues()
  observeEvent(input$start_clusterization, {
    ## Clusterisation with set parameters
    clusterisation$clust_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_clusterisation)]
    som <- get_som(fcs_data$fcs_raw, clusterisation$clust_markers)
    if(input$mode_k_choice == 1){k <- input$maxK}
    if(input$mode_k_choice == 2){k <- input$k}
    mc <- get_consensusClust(som, maxK = k)
    if(input$mode_k_choice == 1){k <- get_optimal_clusters(mc, rate_var_expl = input$rate_var_explan)}
    clusterisation$cell_clustering_list <- get_cluster_annotation(fcs_data$fcs_raw, som, mc, k)
    clusterisation$cell_clustering <- get_cell_clustering_list(som, mc, k)
    ## Estimation of the distance between clusters
    clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
    ## Calculation of edges and nodes for graph
    clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
    clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
    ## Create a data frame to UMAP or tSNE plotting
    sampling_size <- 0.5
    method <- "UMAP"
    if(!is.null(input$n_cell_plot_clasterisation)){sampling_size <- as.numeric(input$n_cell_plot_clasterisation)}
    if(!is.null(input$method_plot_clasterisation)){method <- input$method_plot_clasterisation}
    tsne_inds <- get_inds_subset(fcs_data$fcs_raw, sampling_size = sampling_size)
    clusterisation$umap_df <- get_UMAP_dataframe(fcs_raw = fcs_data$fcs_raw, use_markers = fcs_data$use_markers,
                                                 clust_markers = clusterisation$clust_markers, tsne_inds = tsne_inds,
                                                 cell_clustering = clusterisation$cell_clustering, method = method)
    clusterisation$abundance_df <- get_abundance_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                           cell_clustering = clusterisation$cell_clustering)
  })

  ##### Redrawing plot after chanch number of draw cells ar methods
  observeEvent(input$redraw_clasterisation, {
    if(is.null(clusterisation$cell_clustering)){return(NULL)}
    sampling_size <- 0.5
    method <- "UMAP"
    if(!is.null(input$n_cell_plot_clasterisation)){sampling_size <- as.numeric(input$n_cell_plot_clasterisation)}
    if(!is.null(input$method_plot_clasterisation)){method <- input$method_plot_clasterisation}
    tsne_inds <- get_inds_subset(fcs_data$fcs_raw, sampling_size = sampling_size)
    clusterisation$umap_df <- get_UMAP_dataframe(fcs_raw = fcs_data$fcs_raw, use_markers = fcs_data$use_markers,
                                                 clust_markers = clusterisation$clust_markers, tsne_inds = tsne_inds,
                                                 cell_clustering = clusterisation$cell_clustering, method = method)
  })

  ##### Create UI to choose excluded markers from clusterisation
  output$mk_subset_clusterisation_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("exclude_mk_clusterisation", label = "exclude markers from clusterisation",
                choices = names(fcs_data$use_markers),
                multiple = TRUE)
  })

  ##### Create UI to choose clusters or marker to colour the UMAP plot
  output$mk_target_clusterisation_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("mk_target_clusterisation", label = h3("Ploted marker"),
                choices = c("cluster", names(fcs_data$use_markers)),
                selected = 1)
  })


  ##### Create UI to choose clusters to merge
  output$mergeing_clusterisation_ui <- renderUI({
    if(is.null(clusterisation$cell_clustering)){return(NULL)}
    wellPanel(
      h4("Merging clusters"),
      selectInput("cluster_to_merge_clusterisation", label = h5("Choose clusters"),
                  choices = unique(clusterisation$cell_clustering),
                  multiple = TRUE),
      actionButton("merge_clusterisation", label = "Merge")
    )
  })

  ##### Renew clusterisation reactive object after merging
  observeEvent(input$merge_clusterisation, {
    clusterisation$cell_clustering <- cluster_merging(clusterisation$cell_clustering, input$cluster_to_merge_clusterisation)
    clusterisation$cell_clustering_list <- lapply(clusterisation$cell_clustering_list, function(x)
                                                  cluster_merging(x,input$cluster_to_merge_clusterisation))
    ## Estimation of the distance between clusters
    clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
    ## Calculation of edges and nodes for graph
    clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
    clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
    clusterisation$umap_df$cluster <- cluster_merging(clusterisation$umap_df$cluster, input$cluster_to_merge_clusterisation)
    clusterisation$umap_df$cluster <- as.factor(clusterisation$umap_df$cluster)
  })

  ##### Create UI to rename clusters
  output$rename_clusterisation_ui <- renderUI({
    if(is.null(input$current_node_id)){return(NULL)}
    wellPanel(
      h4("Rename cluster"),
      textInput("new_cluster_name_clusterisation", label = h5("Write new name"),
                          value = as.character(input$current_node_id)),
      actionButton("rename_clusterisation", label = "Rename")
    )
  })

  ##### Renew clusterisation reactive object after cluster rename
  observeEvent(input$rename_clusterisation, {
    if(is.null(input$new_cluster_name_clusterisation)){return(NULL)}
    if(input$new_cluster_name_clusterisation == ""){return(NULL)}
    clusterisation$cell_clustering[clusterisation$cell_clustering==input$current_node_id] <- input$new_cluster_name_clusterisation
    clusterisation$cell_clustering_list <- lapply(clusterisation$cell_clustering_list, function(x)
      x[x==input$current_node_id] <- input$new_cluster_name_clusterisation)
    ## Estimation of the distance between clusters
    clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
    ## Calculation of edges and nodes for graph
    clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
    clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
    levels(clusterisation$umap_df$cluster)[grep(input$current_node_id, levels(clusterisation$umap_df$cluster))] <- input$new_cluster_name_clusterisation
    #clusterisation$umap_df$cluster[clusterisation$umap_df$cluster==input$current_node_id] <- input$new_cluster_name_clusterisation
    #clusterisation$umap_df$cluster <- as.factor(clusterisation$umap_df$cluster)
  })

  ##### Drawing the reactive and interactive UMAP plot
  output$tSNE_plot1_clust <- renderPlot({
    if(is.null(clusterisation$umap_df)){return(NULL)}
    focus_node <- input$current_node_id
    plt <- ggplot(clusterisation$umap_df,  aes(x = UMAP_1, y = UMAP_2, color = eval(parse(text = input$mk_target_clusterisation)))) +
      geom_point(size = 0.8)
    if(input$mk_target_clusterisation != 'cluster'){
      plt <- plt + scale_color_gradient2(midpoint=0.5, low='blue', mid='white', high='red')
    }
    plt + geom_point(data = clusterisation$umap_df[clusterisation$umap_df$cluster == focus_node,], colour = 'black', size = 1)+
      labs(color = input$mk_target_clusterisation)+
      theme_bw()
  })

  ##### Drawing the reactive abundance plot
  output$abundance_clust <- renderPlot({
    if(is.null(clusterisation$umap_df)){return(NULL)}
    plt <- ggplot(clusterisation$abundance_df, aes(x = sample_ids, y = abundance, fill = cluster)) +
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 90))
    return(plt)
  })

  ##### Drawing the reactive and interactive graph with clusters
  output$network <- renderVisNetwork({
    if(is.null(clusterisation$nodes)){return(NULL)}
    edges_threshold <- input$edges_threshold_clusterisation
    if(is.null(input$edges_threshold_clusterisation)){edges_threshold <- 0.5}
    gravity <- input$gravity_clusterisation
    if(is.null(input$gravity_clusterisation)){gravity <- -40}
    edges <- filter_edges(clusterisation$edges, edges_threshold)

    visNetwork(clusterisation$nodes, edges) %>%
      visInteraction(hover = TRUE) %>%
      visEvents(#hoverNode = "function(nodes) { Shiny.onInputChange('current_node_id', nodes);}",
                select = "function(nodes) { Shiny.onInputChange('current_node_id', nodes.nodes);}") %>%
      visPhysics(solver = "forceAtlas2Based",
                 forceAtlas2Based = list(gravitationalConstant = gravity))
  })


  ########################
  ### Gene expressions ###
  ########################

  ##### Create "gene_expression" as reactive object with grouped expression information
  gene_expression <- reactiveValues()
  observeEvent(input$start_clusterization, {
    summarize_mode <- input$method_summarize_expression
    if(is.null(input$method_summarize_expression)){summarize_mode <- 'median'}
    ## Preparing data for heatmap
    gene_expression$cluster_expr_median <- get_mk_clust_heatmap_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                          use_markers = fcs_data$use_markers,
                                                                          cell_clustering = clusterisation$cell_clustering,
                                                                          summarize_mode = summarize_mode)
    gene_expression$deconvol_expr_median <- get_deconvol_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                   use_markers = fcs_data$use_markers,
                                                                   cell_clustering = clusterisation$cell_clustering,
                                                                   summarize_mode = summarize_mode)
  })

  ##### Redrawing hm after chanch summarize mode
  observeEvent(input$redraw_expression, {
    summarize_mode <- input$method_summarize_expression
    if(is.null(input$method_summarize_expression)){summarize_mode <- 'median'}
    ## Preparing data for heatmap
    gene_expression$cluster_expr_median <- get_mk_clust_heatmap_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                          use_markers = fcs_data$use_markers,
                                                                          cell_clustering = clusterisation$cell_clustering,
                                                                          summarize_mode = summarize_mode)
    gene_expression$deconvol_expr_median <- get_deconvol_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                   use_markers = fcs_data$use_markers,
                                                                   cell_clustering = clusterisation$cell_clustering,
                                                                   summarize_mode = summarize_mode)
  })

  ##### Drawing expression heatmap of markers per clusters
  output$cluster_heatmap <- renderD3heatmap({
    if(is.null(gene_expression$cluster_expr_median)){return(NULL)}
    d3heatmap(gene_expression$cluster_expr_median,
              Rowv=NA,
              col=brewer.pal(9,"Reds")
              #scale="none"
              #RowSideColors=country_colours,
              #cellnote=clinical_data,
              #labRow=project_codes,
              #xaxis_font_size=10, yaxis_font_size=10, height=900
    )
  })

  ##### Create UI to choose marker to deconvolution
  output$mk_deconvol_gene_expression_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("mk_deconvol_gene_expression", label = h3("Marker for deconvolution"),
                choices = names(fcs_data$use_markers),
                selected = 1)
  })

  ##### Drawing deconvolution of marker to sample-cluster plot
  output$deconvol_expr <- renderPlot({
    if(is.null(gene_expression$deconvol_expr_median)){return(NULL)}
    if(is.null(input$mk_deconvol_gene_expression)){return(NULL)}
    ggplot(data = gene_expression$deconvol_expr_median,
           aes(x=sample_ids, y=cell_clustering, fill= eval(parse(text = input$mk_deconvol_gene_expression)))) +
      geom_tile() +
      theme_light() +
      labs(fill = input$mk_deconvol_gene_expression) +
      geom_text(aes(label = round(cell_rate, 3)), size = 5, color = 'white') +
      theme(axis.text.x = element_text(angle = 90))
  })

  ########################
  ####  Correlation   ####
  ########################

  ##### Create "correlation" as reactive object with correlation information
  correlation <- reactiveValues()
  observeEvent(input$start_clusterization, {
    ## Get the data of cluster abundance
    correlation$list_cell_ctDist <- get_list_cell_ctDist(clusterisation$cell_clustering_list)
    correlation$abundence_data <- get_abundance(correlation$list_cell_ctDist)
    ## Get correlations of clusters by its abundance
    correlation$abundance_correlation <- get_abundance_correlation(correlation$abundence_data)
  })

  ##### Drawing of triangle plot with correlations by cluster abundances
  output$abun_cor_plot <- renderPlot({
    if(is.null(correlation$abundance_correlation)){return(NULL)}
    corr_coef_matrix <- get_corr_coef_matrix(correlation$abundance_correlation)
    corr_pValue_matrix <- get_corr_pValue_matrix(correlation$abundance_correlation)
    #adj_corr_pValue_matrix <- get_adj_corr_pValue_matrix(correlation$abundance_correlation)
    corrplot(corr_coef_matrix,
             type = "lower",
             method = "circle",
             order = "hclust",
             #addrect = 4,
             tl.col = "black", tl.srt = 45,
             p.mat = corr_pValue_matrix,
             insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white"
    )
  })

  ##### UI for correlation analysis settings
  corr_settings <- reactiveValues(method = "spearman", pValue = 0.01, threshold = 0.1,
                                  bg_anova_alpha = 0.01, bg_ctCor_alpha = 0.01, gene_ctCor_alpha = 0.01)

  observeEvent(input$advanced_corr_settings, {
    output$corr_analysis_settings_ui <- renderUI({
      fluidRow(
        column(5, wellPanel(
          selectInput("method", "Select corr method:", c("spearman" = "spearman", "pearson" = "pearson")),
          numericInput("corr_pValue", "p-Value filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
          numericInput("corr_threshold", "Corr coefficient threshold:", min = -100, max = 100, value = 0.1, step = 0.1)
        )),
        column(5, wellPanel(
          numericInput("corr_bg_anova_alpha", "background anova filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
          numericInput("corr_bg_ctCor_alpha", "Background corr filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
          numericInput("corr_gene_ctCor_alpha", "Gene corr filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001)
        ))
      )
    })
  })
  observeEvent(input$simple_corr_settings,{
    corr_settings$method = "spearman"
    corr_settings$pValue = 0.01
    corr_settings$threshold = 0.1
    corr_settings$bg_anova_alpha = 0.01
    corr_settings$bg_ctCor_alpha = 0.01
    corr_settings$gene_ctCor_alpha = 0.01
  })

  observeEvent(input$advanced_corr_settings, {
    corr_settings$method = input$method
    corr_settings$pValue = input$corr_pValue
    corr_settings$threshold = input$corr_threshold
    corr_settings$bg_anova_alpha = input$corr_bg_anova_alpha
    corr_settings$bg_ctCor_alpha = input$corr_bg_ctCor_alpha
    corr_settings$gene_ctCor_alpha = input$corr_gene_ctCor_alpha
  })

  observeEvent(input$corr_analysis, {
    correlation$list_expData <- get_list_expData(fcs_data$fcs_raw)
    correlation$list_tt_expData <- tt_sample_aggregator(correlation$list_cell_ctDist, correlation$list_expData)
    ##### Background contrast computing
    correlation$bg_ctCor_data <- backgraund_correlation(list_cell_ctDist = correlation$list_cell_ctDist,
                                                        list_expData = correlation$list_expData,
                                                        method = corr_settings$method) ### Time-concuming (~30min)
    correlation$bg_anova <- background_anova(correlation$list_cell_ctDist, correlation$list_expData) ### Time-concuming (~40min)
    ##### Analys for each signalling call type
    correlation$signal_Stat <- signal_extractor(correlation$list_cell_ctDist, correlation$list_tt_expData,
                                                correlation$bg_ctCor_data, correlation$bg_anova,
                                                method = corr_settings$method,
                                                bg_anova_alpha = corr_settings$bg_anova_alpha,
                                                bg_ctCor_alpha = corr_settings$bg_ctCor_alpha,
                                                gene_ctCor_alpha = corr_settings$gene_ctCor_alpha) ### Time-concuming (~40min per cell type)
    correlation$signals <- signal_filter(correlation$signal_Stat, pValue = corr_settings$pValue,
                                         threshold = corr_settings$threshold)
    correlation$signals_between_clusters <- get_signals_between_clusters(correlation$signals)
    correlation$signals_in_cluster <- get_signals_in_cluster(correlation$signals)
  })

  ##### Drawing the reactive and interactive graph with clusters
  output$network2 <- renderVisNetwork({
    if(is.null(clusterisation$nodes)){return(NULL)}
    edges_threshold <- input$edges_threshold_abund_corr
    if(is.null(input$edges_threshold_abund_corr)){edges_threshold <- 0.5}
    gravity <- input$gravity_abund_corr
    if(is.null(input$gravity_abund_corr)){gravity <- -40}
    edges <- filter_edges(clusterisation$edges, edges_threshold)

    visNetwork(clusterisation$nodes, edges) %>%
      visInteraction(hover = T) %>%
      visEvents(select = "function(data) {
                            Shiny.onInputChange('current_node_id2', data.nodes)
                            Shiny.onInputChange('current_edges_id2', data.edges);
                          ;}") %>%
      visPhysics(solver = "forceAtlas2Based",
                 forceAtlas2Based = list(gravitationalConstant = gravity))
  })

  focus_corr <- reactive({
    list(input$current_node_id2, input$current_edges_id2)
  })

  output$gene_cor_table <- DT::renderDataTable(DT::datatable({
    if(is.null(correlation$signals_between_clusters)){return(NULL)}
    if(is.null(correlation$signals_in_cluster)){return(NULL)}
    focus_corr_data <- NULL
    if(is.null(focus_corr()[[1]]) & length(focus_corr()[[2]]==1)){
      target_clusters <- clusterisation$edges[clusterisation$edges$id == focus_corr()[[2]], c("from", "to")]
      signals_between_clusters_top <- get_signals_between_clusters_top(correlation$signals_between_clusters, fcs_data$use_markers)
      focus_corr_data <- signals_between_clusters_top[(signals_between_clusters_top$signaling_cluster %in% target_clusters) &
                                              (signals_between_clusters_top$targetet_cluster %in% target_clusters),]
      if(nrow(focus_corr_data) < 2){
        focus_corr_data <- signals_between_clusters[(signals_between_clusters$signaling_cluster %in% target_clusters) &
                                            (signals_between_clusters$targetet_cluster %in% target_clusters),]
      }
    }
    if(length(focus_corr()[[1]])==1){
      signals_in_cluster_top <- get_signals_in_cluster_top(correlation$signals_in_cluster, fcs_data$use_markers)
      focus_corr_data <- signals_in_cluster_top[signals_in_cluster_top$cluster == focus_corr()[[1]],]
      if (nrow(focus_corr_data) < 2) {focus_corr_data <- signals_in_cluster[signals_in_cluster$cluster == focus_corr()[[1]],]}
    }
    return(focus_corr_data)
  }))


}
