

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
#' @import shiny shinyFiles visNetwork d3heatmap ggplot2
#' @importFrom DT renderDataTable datatable
#' @importFrom magrittr "%>%"
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore write.flowSet
#' @importFrom RColorBrewer brewer.pal
#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @importFrom DT renderDataTable datatable
#' @importFrom neo4r neo4j_api
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.grab grid.draw gpar
#' @importFrom gridGraphics grid.echo
#' @importFrom fs path_home
#'
#' @examples
cytofBrowser_server <- function(input, output){

  ########################
  ###       iBox       ###
  ########################
  output$iBox_upload_dproc <- renderInfoBox({
    if(is.null(fcs_data$md)){
      return(infoBox("Data uploading", "0 files", icon = icon("arrow-alt-circle-up"), color = "red"))
    }
    if(!is.null(fcs_data$md)){
      return(infoBox("Data uploading", paste0(length(fcs_data$md$file_name), " files"), icon = icon("check"), color = "green", fill = T))
    }
  })

  output$iBox_preproc_dproc <- renderInfoBox({
    if(is.null(fcs_data$panel) | is.null(fcs_data$use_markers)){
      return(infoBox("Markers", "No panel", icon = icon("align-center"), color = "red"))
    }
    excluded_mk <- get_use_marker(fcs_data$panel)
    excluded_mk <- excluded_mk[!(excluded_mk %in% fcs_data$use_markers)]
    if(length(excluded_mk) == 0){
      return(infoBox("Markers", paste0(length(fcs_data$panel$name), " in panel"),
                     paste0(length(fcs_data$use_markers), " in use"),  icon = icon("align-right"), color = "yellow"))
    }
    if(length(excluded_mk) != 0){
      return(infoBox("Markers", paste0(length(fcs_data$panel$name), " in panel"),
                     paste0(length(fcs_data$use_markers), " in use, ", length(excluded_mk), " excluded"),
                     icon = icon("check"), color = "green", fill = T))
    }
  })

  output$iBox_clust_dproc <- renderInfoBox({
    if(is.null(clusterisation$cell_clustering)){
      return(infoBox("Clustering", "No clusters", icon = icon("spinner"), color = "red"))
    }
    if(!is.null(clusterisation$cell_clustering)){
      return(infoBox("Data uploading", paste0(length(unique(clusterisation$cell_clustering)), " clusters"),
                     icon = icon("check"), color = "green", fill = T))
    }
  })

  output$iBox_abund_corr <- renderInfoBox({
    if(is.null(correlation$signals_between_clusters) & is.null(correlation$signals_in_cluster)){
      return(infoBox("Abundance correlation analysis", "not done", icon = icon("window-close"), color = "red"))
    }
    if(!is.null(fcs_data$md)){
      return(infoBox("Abundance correlation analysis", paste0(ncol(correlation$signals_between_clusters)+
                                                                ncol(correlation$signals_in_cluster), " signals"), icon = icon("check"), color = "green", fill = T))
    }
  })

  output$iBox_mk_corr <- renderInfoBox({
    if(is.null(correlation$signal_between_cluster_mk) & is.null(correlation$signal_in_cluster_mk)){
      return(infoBox("Marker correlation analysis", "not done", icon = icon("window-close"), color = "red"))
    }
    if(!is.null(fcs_data$md)){
      return(infoBox("Marker correlation analysis", paste0(ncol(correlation$signal_between_cluster_mk)+
                                                                ncol(correlation$signal_in_cluster_mk), " signals"), icon = icon("check"), color = "green", fill = T))
    }
  })

  ########################
  ### Data preparation ###
  ########################
  #roots <- c(home = path.expand("~"))
  roots <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, 'choose_fcs_dp', roots=roots, filetypes=c('', 'fcs'))

  ##### Create "fcs_data" as reactive object to store the CyTOF data
  fcs_data <-reactiveValues()
  plots <-reactiveValues()
  data_prep_settings <- reactiveValues(perplexity = 30, theta = 0.5, max_iter = 1000, size_fuse = 5000)
  ctype <- reactiveValues()
  gates <- reactiveValues()

  observeEvent(input$butt_upload_dproc, {
    if(input$test_data_upload_dproc){fcs_data$md <- get_test_fcs_metadata(input$test_data_dproc)}
    if((length(input$choose_fcs_dp) <= 1) & !input$test_data_upload_dproc){return(NULL)}
    withProgress(message = "Extraction data", min =0, max = 7, value = 0,{
      ## Get row data fcs files
      if(!input$test_data_upload_dproc){
        fcs_data$md <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_dp)$datapath)}
      incProgress(1, detail = "Upload data" )
      fcs_data$fcs_raw <- get_fcs_raw(fcs_data$md)
      incProgress(1, detail = "Extraction panel")
      fcs_data$panel <- get_fcs_panel(fcs_data$fcs_raw)
      incProgress(1, detail = "Delete background markers" )
      fcs_data$use_markers <- get_use_marker(fcs_data$panel)
      incProgress(1, detail = "Cell number calculation" )
      fcs_data$cell_number <- get_cell_number(fcs_data$fcs_raw)
      incProgress(1, detail = "Plotting" )
      ## Preparing data for scatterplot
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter}
      sampling_size <- 0.5
      method <- "tSNE"
      if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
      if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
      fcs_data$tSNE <- scatter_plot_data_prep(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method,
                                    perplexity = data_prep_settings$perplexity, theta = data_prep_settings$theta,
                                    max_iter = data_prep_settings$max_iter, size_fuse = data_prep_settings$size_fuse)
      incProgress(1, detail = "Extraction fcs cluster info")
      if(input$extr_clust_dproc){
        withProgress(message = "Clusters from fcs files", min =0, max = 7, value = 0,{
          pattern <- "clust"
          if(!is.null(input$extr_clust_pattern_dproc)){pattern <- input$extr_clust_pattern_dproc}
          clusterisation$cell_clustering_list <- get_fcs_cluster_annotation(fcs_data$fcs_raw, pattern = pattern)
          incProgress(1, detail = "forming cluster lists")
          clusterisation$cell_clustering <- get_fcs_cell_clustering_vector(clusterisation$cell_clustering_list)
          incProgress(1, detail = "distance between clusters")
          clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
          incProgress(1, detail = "graph elements")
          ## Calculation of edges and nodes for graph
          clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
          clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
          incProgress(1, detail = "drawing scatter plot")
          ## Create a data frame to UMAP or tSNE plotting
          if(!is.null(input$cluster_perplexity)){cluster_settings$perplexity <- input$cluster_perplexity}
          if(!is.null(input$cluster_theta)){cluster_settings$theta <- input$cluster_theta}
          if(!is.null(input$cluster_max_iter)){cluster_settings$max_iter <- input$cluster_max_iter}
          sampling_size_clust <- 0.5
          method_clust <- "UMAP"
          if(!is.null(input$n_cell_plot_clasterisation)){sampling_size_clust <- as.numeric(input$n_cell_plot_clasterisation)}
          if(!is.null(input$method_plot_clasterisation)){method_clust <- input$method_plot_clasterisation}
          tsne_inds <- get_inds_subset(fcs_data$fcs_raw, sampling_size = sampling_size_clust, size_fuse = cluster_settings$size_fuse)
          clusterisation$umap_df <- get_UMAP_dataframe(fcs_raw = fcs_data$fcs_raw, use_markers = fcs_data$use_markers,
                                                       clust_markers = fcs_data$use_markers, tsne_inds = tsne_inds,
                                                       cell_clustering = clusterisation$cell_clustering, method = method_clust,
                                                       perplexity = cluster_settings$perplexity,
                                                       theta = cluster_settings$theta, max_iter = cluster_settings$max_iter)
          incProgress(1, detail = "drawing abundance plot")
          clusterisation$abundance_df <- get_abundance_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                 cell_clustering = clusterisation$cell_clustering)
        })
      }

      ## Add primer column to cell type data frame and gates data frame
      fcs_data$entire_panel <- get_entire_panel(fcs_raw = fcs_data$fcs_raw)
      gates$gates <- data.frame(all_cells = rep(TRUE, sum(fcs_data$cell_number$cell_nmbr)))
      ctype$ctype <- data.frame(all_cells = rep("cell", sum(fcs_data$cell_number$cell_nmbr)),
                                samples = rep(fcs_data$cell_number$smpl, fcs_data$cell_number$cell_nmbr))

      gates$antology <- data.frame(name = "all_cells", parent = NA)
      rownames(gates$antology) <- gates$antology$name
      if(!is.null(clusterisation$cell_clustering)){
        ctype$ctype$cell_clustering <- clusterisation$cell_clustering
      }
      incProgress(1)
    })
  })

  #### reaction to button "size_fuse" in "Data processing" section
  observeEvent(input$size_fuse_dproc, {if(!input$size_fuse_dproc){data_prep_settings$size_fuse <- NA}})

  #### reaction to button "Transform" in "Data processing" section: transformation fsc_data
  observeEvent(input$butt_trans_dproc, {
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    ## Transform row data to scaled data by set parameters
    withProgress(message = "Transformation", min =0, max = 4, value = 0,{
      incProgress(1, detail = "Transformation" )
      if('asinh' %in% input$transformation_list){
        fcs_data$fcs_raw <- asinh_transformation(fcs_data$fcs_raw, cofactor = input$cofactor, use_marker = fcs_data$use_markers)
      }
      incProgress(1, detail = "Outlier detection" )
      if('outlier_by_quantile' %in% isolate(input$transformation_list)){
        fcs_data$fcs_raw <- outlier_by_quantile_transformation(fcs_data$fcs_raw, quantile = input$quantile, use_marker = fcs_data$use_markers)
      }
      incProgress(1, detail = "Plotting" )
      ## Preparing data for scatterplot
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter}
      sampling_size <- 0.5
      method <- "tSNE"
      if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
      if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
      fcs_data$tSNE <- scatter_plot_data_prep(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method,
                                              perplexity = data_prep_settings$perplexity, theta = data_prep_settings$theta,
                                              max_iter = data_prep_settings$max_iter, size_fuse = data_prep_settings$size_fuse)
      incProgress(1)
    })
  })


  ##### Drawing the reactive tSNE plot
  output$scatter_plot_data_preparation <- renderPlot({
    if(is.null(fcs_data$tSNE)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_scatter_dp)){color_mk <- input$mk_scatter_dp}
    plots$scatter_dp <- ggplot(fcs_data$tSNE,  aes(x = tSNE1, y = tSNE2, color = fcs_data$tSNE[,color_mk])) +
      geom_point(size = 0.2) +
      scale_color_gradient2(midpoint = 0.5, low = 'blue', mid = "gray",  high = 'red') +
      labs(color = color_mk) +
      theme_bw()
    return(plots$scatter_dp)
  })

  ##### Download scatter plot data preparation
  output$dwn_scatter_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Scatter_plot_data_preparation", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$scatter_dp, device = ext)
    }
  )

  ##### Create UI to choose target marker
  output$mk_scatter_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput('mk_scatter_dp', label = h4("Plotted marker"),
                choices = names(fcs_data$use_markers),
                #choices = c(1,2,3),
                selected = 1)
  })

  ##### Create UI to choose excluded markers
  output$mk_subset_dp_ui <- renderUI({
    color_mk <- names(fcs_data$use_markers)[1]
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput("exclude_mk_data_preparation", label = "Exclude markers",
                         choices = names(fcs_data$use_markers),
                         multiple = TRUE),
             actionButton("exclud_mk_button", label = "Exclude markers")
      )
    )
  })

  ##### Redrawing plot after chanch number of draw cells ar methods
  observeEvent(input$redraw, {
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    withProgress(message = "Redrawing", min =0, max = 7, value = 0,{
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter}
      sampling_size <- 0.5
      method <- "tSNE"
      if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
      if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
      incProgress(3)
      fcs_data$tSNE <- scatter_plot_data_prep(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method,
                                    perplexity = data_prep_settings$perplexity, theta = data_prep_settings$theta,
                                    max_iter = data_prep_settings$max_iter, size_fuse = data_prep_settings$size_fuse)
      incProgress(4)
    })
  })

  ##### Update reactive object "fcs_data" after excluding markers
  observeEvent(input$exclud_mk_button, {
    withProgress(message = "Excluding", min =0, max = 7, value = 0,{
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter}
      sampling_size <- 0.5
      method <- "tSNE"
      if(!is.null(input$n_cell_plot_data_preparation)){sampling_size <- as.numeric(input$n_cell_plot_data_preparation)}
      if(!is.null(input$method_plot_data_preparation)){method <- input$method_plot_data_preparation}
      fcs_data$use_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_data_preparation)]
      incProgress(3, detail = "Excluding")
      fcs_data$tSNE <- scatter_plot_data_prep(fcs_data$fcs_raw, fcs_data$use_markers, sampling_size = sampling_size, method = method,
                                    perplexity = data_prep_settings$perplexity, theta = data_prep_settings$theta,
                                    max_iter = data_prep_settings$max_iter, size_fuse = data_prep_settings$size_fuse)
      incProgress(4, detail = "Redrawing")
    })
  })

  ##### Reactively show current use_markers from "fcs_data" object
  output$mk_rested_dp <- renderPrint({
    fcs_data$use_markers
  })
  output$mk_excluded_dp <- renderPrint({
    excluded_mk <- get_use_marker(fcs_data$panel)
    excluded_mk <- excluded_mk[!(excluded_mk %in% fcs_data$use_markers)]
    return(excluded_mk)
  })

  ##### Create UI to choose target marker for marker density plot
  output$mk_density_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput('mk_density_dp', label = h4("Plotted marker"),
                choices = names(fcs_data$use_markers),
                #choices = c(1,2,3),
                selected = 1)
  })

  ##### Drawing the reactive histogram plot of marker expression
  output$mk_density_plot_dp <- renderPlot({
    if(is.null(fcs_data$tSNE)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_density_dp)){color_mk <- input$mk_density_dp}
    plots$mk_hist <- ggplot(fcs_data$tSNE, aes(x = fcs_data$tSNE[,color_mk], y=..scaled..)) +
      geom_density(fill = 'black') +
      labs(x = color_mk)
    return(plots$mk_hist)
  })

  ##### Download plot of marker histogram data preparation
  output$dwn_mk_hist_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_mk_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      print(paste("Marker_histogram_plot_data_preparation", ext, sep = "."))
      return(paste("Marker_histogram_plot_data_preparation", ext, sep = ".")) },
    content = function(file) {
      ext <- input$dwn_mk_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$smpl_hist, device = ext)
    }
  )

  ##### Drawing the reactive plot of cell number
  output$smpl_hist_preparation <- renderPlot({
    if(is.null(fcs_data$cell_number)){return(NULL)}
    plots$smpl_hist <- ggplot(data = fcs_data$cell_number, aes(x = smpl, y = cell_nmbr))+
      geom_bar(stat="identity", fill = "black")
    return(plots$smpl_hist)
  })

  ##### Download plot of cell number data preparation
  output$dwn_smpl_hist_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_smpl_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Cell_number_plot_data_preparation", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_smpl_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$smpl_hist, device = ext)
    }
  )

  ########################
  ###      Gating      ###
  ########################

  ##### Create UI to set the gating plot
  output$gating_api_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(fcs_data$entire_panel)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput('gating_subset', label = h4("Which data use for gating"),
                         choices = colnames(gates$gates), selected = 1),
             selectInput('gating_mk1', label = h4("Marker for X"),
                         choices = names(fcs_data$entire_panel), selected = 1),
             selectInput('gating_mk2', label = h4("Marker for Y"),
                         choices = names(fcs_data$entire_panel), selected = 2),
             actionButton("butt_plot_for_gating", label = "Plot for gating")
      )
    )
  })

  ##### Create UI to draw the gating plot
  output$plot_for_gating_ui <- renderUI({
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             plotOutput('scatter_plot_gating', brush = "brush_gating")
      )
    )
  })

  ### Make subset of data to gatingplot
  observeEvent(input$butt_plot_for_gating, {
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gating_subset)){return(NULL)}
    if(is.null(input$gating_mk1)){return(NULL)}
    if(is.null(input$gating_mk2)){return(NULL)}
    withProgress(message = "Gate drawing", min =0, max = 2, value = 0,{
      incProgress(1)
      gates$gated_data_subset <- get_data_for_gating(fcs_raw = fcs_data$fcs_raw,
                                                     gating_subset = gates$gates[,input$gating_subset],
                                                     gating_mk1 = fcs_data$entire_panel[input$gating_mk1],
                                                     gating_mk2 = fcs_data$entire_panel[input$gating_mk2])
      incProgress(1)
    })
  })

  ### Drawing interactive dencity plot for gating
  output$scatter_plot_gating <- renderPlot({
    if(is.null(gates$gated_data_subset)){return(NULL)}
    colfunc <- grDevices::colorRampPalette(c("black", "red", "yellow"))
    min_mk1 <- min(gates$gated_data_subset$gating_mk1)
    max_mk1 <- max(gates$gated_data_subset$gating_mk1)
    min_mk2 <- min(gates$gated_data_subset$gating_mk2)
    max_mk2 <- max(gates$gated_data_subset$gating_mk2)
    plots$scatter_plot_gating <- ggplot(gates$gated_data_subset, aes(x = gating_mk1, y = gating_mk2)) +
      ylim((min_mk2 - 0.1*(max_mk2 - min_mk2)),
           (max_mk2 + 0.1*(max_mk2 - min_mk2))) +
      xlim((min_mk1 - 0.1*(max_mk1 - min_mk1)),
           (max_mk1 + 0.1*(max_mk1 - min_mk1))) +
      xlab(input$gating_mk1) +
      ylab(input$gating_mk2) +
      geom_point(size = 0.1) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon") +
      scale_fill_gradientn(colours=colfunc(400))+
      geom_density2d(colour="black") +
      theme_bw() +
      theme(legend.position = "none")
    return(plots$scatter_plot_gating)
  })

  output$gating_text <- renderPrint({
    if(is.null(input$brush_gating)){return(NULL)}
    list(colnames(gates$gates), colnames(ctype$ctype))
  })

  observeEvent(input$gete_chosen_cells, {
    if(is.null(input$brush_gating)){return(NULL)}
    if(is.null(gates$gated_data_subset)){return(NULL)}
    new_name <- paste0("Gate", as.character(ncol(gates$gates)+1), "_",
                      input$gating_mk1, "_", input$gating_mk2, "_from_", strsplit(input$gating_subset, "_")[[1]][1])

    if(!is.null(input$new_gate_name) & (input$new_gate_name != "marker1+/marker2+")){new_name <- input$new_gate_name}
    original_cell_coordinates <- brushedPoints(gates$gated_data_subset, input$brush_gating, xvar = "gating_mk1", yvar = "gating_mk2")
    original_cell_coordinates <- original_cell_coordinates$original_cell_coordinates
    gates$gates$new_gate <- FALSE
    gates$gates[original_cell_coordinates, "new_gate"] <- TRUE
    if(any(grepl(new_name, colnames(gates$gates)))){
      new_name <- paste0(new_name,"_",sum(grepl(new_name, colnames(gates$gates)))+1)}
    colnames(gates$gates)[which(colnames(gates$gates) ==  "new_gate")] <- new_name

    gates$antology <- rbind(gates$antology, data.frame(name = new_name, parent = input$gating_subset))
    rownames(gates$antology) <- gates$antology$name
  })

  output$gate_antology_graph <- renderVisNetwork({
    if(is.null(gates$antology)){return(NULL)}
    nodes <- data.frame(id = gates$antology$name, label = gates$antology$name)
    edges <- data.frame(from = gates$antology$parent, to = gates$antology$name)
    edges <- edges[!is.na(edges$from),]

    visNetwork(nodes, edges) %>%
      visInteraction(hover = TRUE) %>%
      visEvents(select = "function(nodes) { Shiny.onInputChange('gated_node_id', nodes.nodes);}") %>%
      visHierarchicalLayout()
  })

  ##### Create UI to convert gates to cell types
  output$mergeing_gates_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Convert gates to cell annotation"),
             selectInput("gate_converting_method", label = h5("converting method"),
                         choices = list("Pure selection" = 'pure', "Gates squeezing" = 'squeeze'),
                         selected = 1),
             selectInput("gates_to_convert_gating", label = h5("Choose gates"),
                         choices = gates$antology$name[-1], multiple = TRUE),
             actionButton("convert_gates", label = "Convert")
      )
    )
  })

  ##### Renew clusterisation reactive object after merging
  observeEvent(input$convert_gates, {
    ctype$ctype <- get_cell_type_from_gates(gates$gates, input$gates_to_convert_gating,
                                            ctype$ctype, method = input$gate_converting_method)

  })

  ##### Create UI to rename gates
  output$rename_gates_ui <- renderUI({
    if(is.null(input$gated_node_id)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Rename gates"),
             textInput('new_gate_name_gating', label = h5("Write new name"),
                       value = as.character(input$gated_node_id)),
             actionButton("rename_gates", label = "Rename")
      )
    )
  })

  ##### Renew gate reactive object after gate rename
  observeEvent(input$rename_gates, {
    if(is.null(input$new_gate_name_gating)){return(NULL)}
    if(input$new_gate_name_gating == ""){return(NULL)}
    levels(gates$antology$name)[levels(gates$antology$name) == input$gated_node_id] <- input$new_gate_name_gating
    levels(gates$antology$parent)[levels(gates$antology$parent) == input$gated_node_id] <- input$new_gate_name_gating
    rownames(gates$antology) <- gates$antology$name
    colnames(gates$gates)[colnames(gates$gates) == input$gated_node_id] <- input$new_gate_name_gating
  })

  ##### Create UI to convert cell type to clusters
  output$convert_cells_annotation_ui <- renderUI({
    if(is.null(ctype$ctype)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Convert cell annotation to clusters"),
             selectInput("ctype_to_convert_gates", label = h5("Choose cell annotation to clusters"),
                         choices = colnames(ctype$ctype), multiple = FALSE),
             actionButton("convert_ctype_to_clusters", label = "Convert")
      )
    )
  })

  ###### Set choosed cell type to clusters and update objects
  #observeEvent(input$convert_ctype_to_clusters, {
  #  withProgress(message = "Convert", min =0, max = 7, value = 0,{
  #    clusterisation$cell_clustering <- ctype$ctype[, input$ctype_to_convert_gates]
  #    incProgress(1)
  #    clusterisation$cell_clustering_list <- lapply(unique(ctype$ctype$samples), function(x)
  #      ctype$ctype[ctype$ctype$samples == x, input$ctype_to_convert_gates])
  #    incProgress(1)
  #    ## Estimation of the distance between clusters
  #    clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
  #    incProgress(1)
  #    ## Calculation of edges and nodes for graph
  #    clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
  #    incProgress(1)
  #    clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
  #    incProgress(1)
  #    clusterisation$umap_df$cluster <- ctype$ctype[, input$ctype_to_convert_gates] ### Mistake
  #    incProgress(1)
  #    clusterisation$umap_df$cluster <- as.factor(clusterisation$umap_df$cluster)
  #    ## Renew data in cell type data frame
  #    ctype$ctype$cell_clustering <- clusterisation$cell_clustering
  #    incProgress(1)
  #  })
  #})

  ##### UI for extraction data from fcs to add to cell annotation
  output$extract_cell_annotation_ui <- renderUI({
    if(is.null(ctype$ctype)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Extract cell annotations"),
             selectInput("col_to_ctype_gates", label = h5("Choose data to add to cell annotation"),
                         choices = names(fcs_data$entire_panel), multiple = TRUE),
             actionButton("btm_col_to_ctype_gates", label = "Extract")
      )
    )
  })

  ##### Extraction data from fcs to add to cell annotation
  observeEvent(input$btm_col_to_ctype_gates, {
    if(is.null(ctype$ctype)){return(NULL)}
    ctype$ctype <- get_add_cell_annotation_from_data(entire_panel = fcs_data$entire_panel,
                                                     fcs_raw = fcs_data$fcs_raw, cell_annotation = ctype$ctype)
  })

  ##### UI for saving fcs files with cell annotations
  output$save_cell_annaotation <- renderUI({
    if(is.null(ctype$ctype)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Cell annotations to save in fcs"),
             selectInput("save_cell_ann_col_names", label = h5("Choose cell annotation"),
                         choices = colnames(ctype$ctype)[-1], multiple = TRUE),
             h5("Saving the data for samples as FCS files with cell annotation"),
             shinyDirButton('choose_panel_gate', "Folder choose", "Select a folder to save panel samples"),
             hr(),
             h5("Saved set of files can be used as panel data in a cross-panel analysiso"),
             textInput("panel_name_gate", label = h4("Panel name"), value = "Panel1"),
             hr(),
             actionButton('dwn_panel_gate', label = "Save panel")
      )
    )
  })

  ##### Save sample data with cluster info as panal files
  shinyDirChoose(input, 'choose_panel_gate', roots = roots)
  observeEvent(input$dwn_panel_gate, {
    if(is.null(input$choose_panel_gate)){return(NULL)}
    if(is.null(ctype$ctype)){return(NULL)}
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    panel_folder_path <- parseDirPath(roots, input$choose_panel_gate)
    panel_name <- input$panel_name_gate
    if(is.null(panel_name)){panel_name <- "UNNAMED_PANEL"}
    if(!is.null(panel_folder_path) | length(panel_folder_path) !=0 ){
      clustered_fcs <- get_cell_annotation_fcs_files(fcs_raw = fcs_data$fcs_raw,
                                                     cell_annotation = ctype$ctype,
                                                     column_name = input$save_cell_ann_col_names)
      filenames <- paste0("CyBr_", as.character(panel_name),"_", sampleNames(clustered_fcs), ".fcs")
      flowCore::write.flowSet(clustered_fcs, outdir = panel_folder_path, filename = filenames)
    }
  })

  ########################
  ###  Clusterisation  ###
  ########################

  ##### Create "clusterisation" as reactive object to store the cluster information
  clusterisation <- reactiveValues()
  cluster_settings <- reactiveValues(perplexity = 30, theta = 0.5, max_iter = 1000, size_fuse = 5000)

  observeEvent(input$start_clusterization, {
    withProgress(message = "Clustering", min =0, max = 7, value = 0,{
      ## Clusterisation with set parameters
      clusterisation$clust_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_clusterisation)]
      som <- get_som(fcs_data$fcs_raw, clusterisation$clust_markers)
      if(input$mode_k_choice == 1){k <- input$maxK}
      if(input$mode_k_choice == 2){k <- input$k}
      incProgress(1, detail = "clustering")
      mc <- get_consensusClust(som, maxK = k)
      incProgress(1, detail = "forming cluster lists")
      if(input$mode_k_choice == 1){k <- get_optimal_clusters(mc, rate_var_expl = input$rate_var_explan)}
      clusterisation$cell_clustering_list <- get_cluster_annotation(fcs_data$fcs_raw, som, mc, k)
      clusterisation$cell_clustering <- get_cell_clustering_vector(som, mc, k)
      incProgress(1, detail = "distance between clusters")
      ## Estimation of the distance between clusters
      clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
      incProgress(1, detail = "graph elements")
      ## Calculation of edges and nodes for graph
      clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
      clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
      incProgress(1, detail = "drawing scatter plot")
      ## Create a data frame to UMAP or tSNE plotting
      if(!is.null(input$cluster_perplexity)){cluster_settings$perplexity <- input$cluster_perplexity}
      if(!is.null(input$cluster_theta)){cluster_settings$theta <- input$cluster_theta}
      if(!is.null(input$cluster_max_iter)){cluster_settings$max_iter <- input$cluster_max_iter}
      sampling_size <- 0.5
      method <- "UMAP"
      if(!is.null(input$n_cell_plot_clasterisation)){sampling_size <- as.numeric(input$n_cell_plot_clasterisation)}
      if(!is.null(input$method_plot_clasterisation)){method <- input$method_plot_clasterisation}
      tsne_inds <- get_inds_subset(fcs_data$fcs_raw, sampling_size = sampling_size, size_fuse = cluster_settings$size_fuse)
      clusterisation$umap_df <- get_UMAP_dataframe(fcs_raw = fcs_data$fcs_raw, use_markers = fcs_data$use_markers,
                                                   clust_markers = clusterisation$clust_markers, tsne_inds = tsne_inds,
                                                   cell_clustering = clusterisation$cell_clustering, method = method,
                                                   perplexity = cluster_settings$perplexity,
                                                   theta = cluster_settings$theta, max_iter = cluster_settings$max_iter)
      clusterisation$abundance_df <- get_abundance_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                             cell_clustering = clusterisation$cell_clustering)
      ## Add cluster data to cell type data frame
      ctype$ctype$cell_clustering <- clusterisation$cell_clustering
      incProgress(1)
    })
  })

  #### reaction to button "size_fuse" in "Clustering" section
  observeEvent(input$size_fuse_clust, {if(!input$size_fuse_clust){cluster_settings$size_fuse <- NA}})

  ##### Redrawing plot after chanch number of draw cells ar methods
  observeEvent(input$redraw_clasterisation, {
    withProgress(message = "Redrawing", min =0, max = 2, value = 0,{
      if(is.null(clusterisation$cell_clustering)){return(NULL)}
      if(!is.null(input$cluster_perplexity)){cluster_settings$perplexity <- input$cluster_perplexity}
      if(!is.null(input$cluster_theta)){cluster_settings$theta <- input$cluster_theta}
      if(!is.null(input$cluster_max_iter)){cluster_settings$max_iter <- input$cluster_max_iter}
      sampling_size <- 0.5
      method <- "UMAP"
      if(!is.null(input$n_cell_plot_clasterisation)){sampling_size <- as.numeric(input$n_cell_plot_clasterisation)}
      if(!is.null(input$method_plot_clasterisation)){method <- input$method_plot_clasterisation}
      tsne_inds <- get_inds_subset(fcs_data$fcs_raw, sampling_size = sampling_size, size_fuse = cluster_settings$size_fuse)
      incProgress(1)
      clusterisation$umap_df <- get_UMAP_dataframe(fcs_raw = fcs_data$fcs_raw, use_markers = fcs_data$use_markers,
                                                   clust_markers = clusterisation$clust_markers, tsne_inds = tsne_inds,
                                                   cell_clustering = clusterisation$cell_clustering, method = method,
                                                   perplexity = cluster_settings$perplexity,
                                                   theta = cluster_settings$theta, max_iter = cluster_settings$max_iter)
      incProgress(1)
    })
  })

  ##### Create UI to choose excluded markers from clusterisation
  output$mk_subset_clusterisation_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("exclude_mk_clusterisation", label = "Exclude markers from clustering",
                choices = names(fcs_data$use_markers),
                multiple = TRUE)
  })

  ##### Create UI to choose clusters to merge
  output$mergeing_clusterisation_ui <- renderUI({
    if(is.null(clusterisation$cell_clustering)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Merging clusters"),
             selectInput("cluster_to_merge_clusterisation", label = h5("Choose clusters"),
                         choices = unique(clusterisation$cell_clustering), multiple = TRUE),
             actionButton("merge_clusterisation", label = "Merge")
      )
    )
  })

  ##### Renew clusterisation reactive object after merging
  observeEvent(input$merge_clusterisation, {
    withProgress(message = "Merging", min =0, max = 7, value = 0,{
      clusterisation$cell_clustering <- cluster_merging(clusterisation$cell_clustering, input$cluster_to_merge_clusterisation)
      incProgress(1)
      clusterisation$cell_clustering_list <- lapply(clusterisation$cell_clustering_list, function(x)
        cluster_merging(x,input$cluster_to_merge_clusterisation))
      incProgress(1)
      ## Estimation of the distance between clusters
      clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
      incProgress(1)
      ## Calculation of edges and nodes for graph
      clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
      incProgress(1)
      clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
      incProgress(1)
      clusterisation$umap_df$cluster <- cluster_merging(clusterisation$umap_df$cluster, input$cluster_to_merge_clusterisation)
      incProgress(1)
      clusterisation$umap_df$cluster <- as.factor(clusterisation$umap_df$cluster)
      ## Renew data in cell type data frame
      ctype$ctype$cell_clustering <- clusterisation$cell_clustering
      incProgress(1)
    })
  })

  ##### Create UI to rename clusters
  output$rename_clusterisation_ui <- renderUI({
    if(is.null(input$current_node_id)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Rename cluster"),
             textInput("new_cluster_name_clusterisation", label = h5("Write new name"),
                       value = as.character(input$current_node_id)),
             actionButton("rename_clusterisation", label = "Rename")
      )
    )
  })

  ##### Renew clusterisation reactive object after cluster rename
  observeEvent(input$rename_clusterisation, {
    if(is.null(input$new_cluster_name_clusterisation)){return(NULL)}
    if(input$new_cluster_name_clusterisation == ""){return(NULL)}
    withProgress(message = "Renaming", min =0, max = 5, value = 0,{
      clusterisation$cell_clustering[clusterisation$cell_clustering==input$current_node_id] <- input$new_cluster_name_clusterisation
      clusterisation$cell_clustering_list <- lapply(clusterisation$cell_clustering_list, function(x)
        x[x==input$current_node_id] <- input$new_cluster_name_clusterisation)
      incProgress(1)
      ## Estimation of the distance between clusters
      clusterisation$clus_euclid_dist <- get_euclid_dist(fcs_data$fcs_raw, fcs_data$use_markers, clusterisation$cell_clustering)
      incProgress(1)
      ## Calculation of edges and nodes for graph
      clusterisation$edges <- get_edges(clusterisation$clus_euclid_dist)
      incProgress(1)
      clusterisation$nodes <- get_nodes(clusterisation$edges, clusterisation$cell_clustering)
      incProgress(1)
      levels(clusterisation$umap_df$cluster)[grep(input$current_node_id, levels(clusterisation$umap_df$cluster))] <- input$new_cluster_name_clusterisation
      incProgress(1)
      #clusterisation$umap_df$cluster[clusterisation$umap_df$cluster==input$current_node_id] <- input$new_cluster_name_clusterisation
      #clusterisation$umap_df$cluster <- as.factor(clusterisation$umap_df$cluster)
    })
  })

  ##### Drawing the reactive and interactive UMAP plot
  output$scatter_plot_clust <- renderPlot({
    if(is.null(clusterisation$umap_df)){return(NULL)}

    focus_node <- input$current_node_id
    plt <- ggplot(clusterisation$umap_df,  aes(x = UMAP_1, y = UMAP_2, color = clusterisation$umap_df[,input$mk_target_clusterisation])) +
      geom_point(size = 0.8)
    if(input$mk_target_clusterisation == 'cluster'){
      plt <- plt + scale_color_manual(values = as.character(clusterisation$nodes$color))
    }
    if(input$mk_target_clusterisation != 'cluster'){
      plt <- plt + scale_color_gradient2(midpoint=0.5, low='blue', mid='gray', high='red')
    }
    plt <- plt + geom_point(data = clusterisation$umap_df[clusterisation$umap_df$cluster == focus_node,], colour = 'black', size = 1)+
      labs(color = input$mk_target_clusterisation)+
      theme_bw()
    plots$scatter_clust <- plt
    return(plt)
  })

  ##### Create Scatter plot and mk choice ui for dpsection
  output$scatter_plot_dp_ui <- renderUI({
    if(is.null(fcs_data$tSNE)){return(NULL)}
    if(is.null(clusterisation$umap_df)){
      return(
        box(
          fluidRow(
            column(2, actionBttn(inputId = "redraw", style = "material-circle", color = "default" ,icon = icon("redo"))),
            column(2,
                   dropdownButton(
                     tags$h4("Options of plotting"),
                     numericInput("n_cell_plot_data_preparation",
                                  label = h5("Cell fraction to display"), value = 0.5, step = 0.1),
                     materialSwitch(inputId = 'size_fuse_dproc', label = h4("Size fuse"), value = TRUE),
                     selectInput("method_plot_data_preparation", label = h5("Visualisation method"),
                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                 selected = "tSNE"),
                     icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                   )
            ),
            column(2,
                   dropdownButton(
                     tags$h4("Advanced options"),
                     uiOutput('advanced_opt_dp_ui'),
                     icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                   )

            ),
            column(2,
                   dropdownButton(
                     selectInput('dwn_scatter_dp_ext', label = NULL,
                                 choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                     downloadButton('dwn_scatter_dp', ""),
                     icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                   )

            )
          ),
          plotOutput('scatter_plot_data_preparation'),
          selectInput('mk_scatter_dp', label = h4("Plotted marker"),
                      choices = names(fcs_data$use_markers), selected = 1)
        )
      )
    }
    if(!is.null(clusterisation$umap_df)){
      return(
        box(
          fluidRow(
            column(2, actionBttn(inputId = "redraw_clasterisation", style = "material-circle", color = "default" ,icon = icon("redo"))),
            column(2,
                   dropdownButton(
                     tags$h4("Options of plotting"),
                     numericInput("n_cell_plot_clasterisation",
                                  label = h4("Cell fraction to display"), value = 0.5, step = 0.1),
                     materialSwitch(inputId = 'size_fuse_clust', label = h4("Size fuse"), value = TRUE),
                     selectInput("method_plot_clasterisation", label = h4("Visualisation method"),
                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                 selected = "UMAP"),
                     icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                   )
            ),
            column(2,
                   dropdownButton(
                     tags$h4("Advanced options"),
                     uiOutput('advanced_opt_clust_ui'),
                     icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                   )
            ),
            column(2,
                   dropdownButton(
                     selectInput('dwn_scatter_clust_ext', label = NULL,
                                 choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                     downloadButton('dwn_scatter_clust', ""),
                     icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                   )

            )
          ),
          plotOutput('scatter_plot_clust', click = "plot_click"),
          selectInput("mk_target_clusterisation", label = h4("Plotted marker"),
                      choices = c("cluster", names(fcs_data$use_markers)),
                      selected = 1)
        )
      )
    }
  })

  ##### UI for advanced ortions data preparation
  output$advanced_opt_dp_ui <- renderUI({
    if(input$method_plot_data_preparation == 'tSNE'){return(
      fluidRow(
        column(1),
        column(10,
               numericInput("data_prep_perplexity", "tSNE Perplexity", min = 0, max = 200, value = 30, step = 5),
               numericInput("data_prep_theta", "tSNE Theta", min = 0, max = 1, value = 0.5, step = 0.1),
               numericInput("data_prep_max_iter", "tSNE Iterations", value = 1000, step = 500)
        )
      )
    )}
  })

  ##### UI for advanced ortions data preparation
  output$advanced_opt_clust_ui <- renderUI({
    if(input$method_plot_clasterisation == 'tSNE'){return(
      fluidRow(
        column(1),
        column(10,
               numericInput("cluster_perplexity", "Perplexity", min = 0, max = 200, value = 30, step = 5),
               numericInput("cluster_theta", "Theta", min = 0, max = 1, value = 0.5, step = 0.1),
               numericInput("cluster_max_iter", "Iterations", value = 1000, step = 500)
        )
      )
    )}
  })

  ##### Download scatter plot for clustering
  output$dwn_scatter_clust <- downloadHandler(
    filename = function() {
      ext <- input$dwn_scatter_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Scatter_plot_clustering", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_scatter_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$scatter_clust, device = ext)
    }
  )

  ##### Drawing the reactive abundance plot
  output$abundance_clust <- renderPlot({
    if(is.null(clusterisation$umap_df)){return(NULL)}
    plots$abundance_clust <- ggplot(clusterisation$abundance_df, aes(x = sample_ids, y = abundance, fill = cluster)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = as.character(clusterisation$nodes$color)) +
      theme(axis.text.x = element_text(angle = 90))
    return(plots$abundance_clust)
  })

  ##### Download abundance plot for clustering
  output$dwn_abundance_clust <- downloadHandler(
    filename = function() {
      ext <- input$dwn_abundance_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Abundance_plot_clustering", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_abundance_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$abundance_clust, device = ext)
    }
  )


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

  ##### Create network of clusters ui for navigation
  output$network_clust_ui <- renderUI({
    if (is.null(clusterisation$edges)){return(NULL)}
    box(
      visNetworkOutput("network"),
      sliderInput('edges_threshold_clusterisation', "Edge weight threshold for graph",
                  min =0, max = 1, value = 0.5, step = 0.01),
      sliderInput('gravity_clusterisation', "Gravity for graph",
                  min = -100, max = 0, value = -40, step = 1)
    )
  })

  ##### Save sample data with cluster info as panal files
  shinyDirChoose(input, 'choose_panel_clust', roots = roots)
  observeEvent(input$dwn_panel_clust, {
    if(is.null(clusterisation$cell_clustering_list)){return(NULL)}
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    panel_folder_path <- parseDirPath(roots, input$choose_panel_clust)
    panel_name <- input$panel_name_clust
    if(is.null(panel_name)){panel_name <- "UNNAMED_PANEL"}
    print(length(panel_folder_path))
    if(!is.null(panel_folder_path) | length(panel_folder_path) !=0 ){
      clustered_fcs <- get_clustered_fcs_files(fcs_data$fcs_raw, clusterisation$cell_clustering_list, column_name = "cluster")
      filenames <- paste0("CyBr_", as.character(panel_name),"_", sampleNames(clustered_fcs), ".fcs")
      flowCore::write.flowSet(clustered_fcs, outdir = panel_folder_path, filename = filenames)
    }
  })


  ########################
  ### Gene expressions ###
  ########################

  ##### Create "gene_expression" as reactive object with grouped expression information
  gene_expression <- reactiveValues()
  observeEvent(input$start_clusterization, {
    withProgress(message = "Heatmap drawing", min =0, max = 3, value = 0,{
      summarize_mode <- input$method_summarize_expression
      if(is.null(input$method_summarize_expression)){summarize_mode <- 'median'}
      incProgress(1)
      ## Preparing data for heatmap
      gene_expression$cluster_expr_median <- get_mk_clust_heatmap_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                            use_markers = fcs_data$use_markers,
                                                                            cell_clustering = clusterisation$cell_clustering,
                                                                            summarize_mode = summarize_mode)
      incProgress(1)
      gene_expression$deconvol_expr_median <- get_deconvol_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                     use_markers = fcs_data$use_markers,
                                                                     cell_clustering = clusterisation$cell_clustering,
                                                                     summarize_mode = summarize_mode)
      incProgress(1)
    })

  })

  ##### Redrawing hm after chanch summarize mode
  observeEvent(input$redraw_expression, {
    withProgress(message = "Redrawing", min =0, max = 3, value = 0,{
      summarize_mode <- input$method_summarize_expression
      if(is.null(input$method_summarize_expression)){summarize_mode <- 'median'}
      incProgress(1)
      ## Preparing data for heatmap
      gene_expression$cluster_expr_median <- get_mk_clust_heatmap_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                            use_markers = fcs_data$use_markers,
                                                                            cell_clustering = clusterisation$cell_clustering,
                                                                            summarize_mode = summarize_mode)
      incProgress(1)
      gene_expression$deconvol_expr_median <- get_deconvol_dataframe(fcs_raw = fcs_data$fcs_raw,
                                                                     use_markers = fcs_data$use_markers,
                                                                     cell_clustering = clusterisation$cell_clustering,
                                                                     summarize_mode = summarize_mode)
      incProgress(1)
    })
  })

  ##### Drawing expression heatmap of markers per clusters
  output$cluster_heatmap <- renderD3heatmap({
    if(is.null(gene_expression$cluster_expr_median)){return(NULL)}
    d3heatmap(gene_expression$cluster_expr_median,
              #Rowv = NA,
              col = RColorBrewer::brewer.pal(9,"Reds")
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
    selectInput("mk_deconvol_gene_expression", label = h4("Marker for deconvolution"),
                choices = names(fcs_data$use_markers),
                selected = 1)
  })

  ##### Drawing deconvolution of marker to sample-cluster plot
  output$deconvol_expr <- renderPlot({
    if(is.null(gene_expression$deconvol_expr_median)){return(NULL)}
    if(is.null(input$mk_deconvol_gene_expression)){return(NULL)}
    plots$deconvol_expr <- ggplot(data = gene_expression$deconvol_expr_median,
           aes(x=sample_ids, y=cell_clustering, fill= gene_expression$deconvol_expr_median[, input$mk_deconvol_gene_expression])) +
      geom_tile() +
      theme_light() +
      labs(fill = input$mk_deconvol_gene_expression) +
      geom_text(aes(label = round(cell_rate, 3)), size = 5, color = 'white') +
      theme(axis.text.x = element_text(angle = 90))
    return(plots$deconvol_expr)
  })

  ##### Download deconvolution plot of marker for marker expression
  output$dwn_deconvol_expr <- downloadHandler(
    filename = function() {
      ext <- input$dwn_deconvol_expr_ext
      mk <- input$mk_deconvol_gene_expression
      if(is.null(ext)){ext <- "pdf"}
      if(is.null(mk)){mk <- ""}
      paste0("Deconvolution_plot_clustering_", mk, ".", ext) },
    content = function(file) {
      ext <- input$dwn_deconvol_expr_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$deconvol_expr, device = ext)
    }
  )

  observeEvent(input$drawn_cluster_hm_expr, {
    if(is.null(gene_expression$cluster_expr_median)){return(NULL)}
    plots$cluster_hm_expr <- ComplexHeatmap::Heatmap(as.matrix(gene_expression$cluster_expr_median),
                                                     col = RColorBrewer::brewer.pal(9,"Reds"),
                                                     heatmap_legend_param = list(title = "expression"))
  })

  ##### Drawing expression heatmap of markers per clusters
  output$cluster_hm_expr <- renderPlot({plots$cluster_hm_expr})

  ##### Download heatmap expression plot marker expressions
  output$dwn_drawn_cluster_hm_expr <- downloadHandler(
    filename = function() {
      ext <- input$dwn_drawn_cluster_hm_expr_ext
      if(is.null(ext)){ext <- "pdf"}
      return(paste("Clusters_hm_expression", ext, sep = ".")) },
    content = function(file) {
      ext <- input$dwn_drawn_cluster_hm_expr_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        ComplexHeatmap::draw(plots$cluster_hm_expr)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        ComplexHeatmap::draw(plots$cluster_hm_expr)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        ComplexHeatmap::draw(plots$cluster_hm_expr)
        dev.off()
      }
      #ggsave(file, plot = plots$cluster_hm_expr, device = ext)
    }
  )

  ########################
  ####  Correlation   ####
  ########################

  ##### Create "correlation" as reactive object with correlation information
  correlation <- reactiveValues()
  observeEvent(input$start_clusterization, {
    withProgress(message = "Cluster abundances", min =0, max = 3, value = 0,{
      ## Get the data of cluster abundance
      correlation$list_cell_ctDist <- get_list_cell_ctDist(clusterisation$cell_clustering_list)
      incProgress(1)
      correlation$abundence_data <- get_abundance(correlation$list_cell_ctDist)
      incProgress(1)
      ## Get correlations of clusters by its abundance
      correlation$abundance_correlation <- get_abundance_correlation(correlation$abundence_data)
      incProgress(1)
    })
  })

  ##### Drawing of triangle plot with correlations by cluster abundances
  output$abun_cor_plot <- renderPlot({
    if(is.null(correlation$abundance_correlation)){return(NULL)}
    corr_coef_matrix <- get_corr_coef_matrix(correlation$abundance_correlation)
    corr_pValue_matrix <- get_corr_pValue_matrix(correlation$abundance_correlation)
    #adj_corr_pValue_matrix <- get_adj_corr_pValue_matrix(correlation$abundance_correlation)
    corrplot::corrplot(corr_coef_matrix,
             type = "lower",
             method = "circle",
             order = "hclust",
             #addrect = 4,
             tl.col = "black", tl.srt = 45,
             p.mat = corr_pValue_matrix,
             insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = 1.2, pch.col = "white"
    )
    gridGraphics::grid.echo()
    plots$abun_corr <- grid::grid.grab()
    return(plots$abun_corr)
  })

  ##### Download abundance correlation plot
  output$dwn_abund_corr <- downloadHandler(
    filename = function() {
      ext <- input$dwn_abund_corr_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Abunadnce_correlations", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_abund_corr_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        grid::grid.draw(plots$abun_corr)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        grid::grid.draw(plots$abun_corr)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        grid::grid.draw(plots$abun_corr)
        dev.off()
      }
      #ggsave(file, plot = plots$abun_corr, device = ext)
    }
  )

  ##### UI for correlation analysis settings
  corr_settings <- reactiveValues(method = "spearman", pValue = 0.01, threshold = 0.1,
                                  bg_anova_alpha = 0.01, bg_ctCor_alpha = 0.01, gene_ctCor_alpha = 0.01)

  observeEvent(input$simple_corr_settings,{
    output$corr_analysis_settings_ui <- NULL
    corr_settings$method <- "spearman"
    corr_settings$pValue <- 0.01
    corr_settings$threshold <- 0.1
    corr_settings$bg_anova_alpha <- 0.01
    corr_settings$bg_ctCor_alpha <- 0.01
    corr_settings$gene_ctCor_alpha <- 0.01
  })

  observeEvent(input$advanced_corr_settings, {
    if(!is.null(input$method)){corr_settings$method <- input$method}
    if(!is.null(input$corr_pValue)){corr_settings$pValue <- input$corr_pValue}
    if(!is.null(input$corr_threshold)){corr_settings$threshold <- input$corr_threshold}
    if(!is.null(input$corr_bg_anova_alpha)){corr_settings$bg_anova_alpha <- input$corr_bg_anova_alpha}
    if(!is.null(input$corr_bg_ctCor_alpha)){corr_settings$bg_ctCor_alpha <- input$corr_bg_ctCor_alpha}
    if(!is.null(input$corr_gene_ctCor_alpha)){corr_settings$gene_ctCor_alpha <- input$corr_gene_ctCor_alpha}
  })

  observeEvent(input$corr_analysis, {
    options(warn=-1)
    withProgress(message = "Abundances correlation", min =0, max = 6, value = 0,{
      correlation$list_expData <- get_list_expData(fcs_data$fcs_raw)
      correlation$list_tt_expData <- tt_sample_aggregator(correlation$list_cell_ctDist, correlation$list_expData)
      incProgress(1, detail = "background correlations")
      ##### Background contrast computing
      correlation$bg_ctCor_data <- backgraund_correlation(list_cell_ctDist = correlation$list_cell_ctDist,
                                                          list_expData = correlation$list_expData,
                                                          method = corr_settings$method) ### Time-concuming (~30min)
      incProgress(1, detail = "background anova")
      correlation$bg_anova <- background_anova(correlation$list_cell_ctDist, correlation$list_expData) ### Time-concuming (~40min)
      incProgress(1, detail = "signals extraction")
      ##### Analys for each signalling call type
      correlation$signal_Stat <- signal_extractor(correlation$list_cell_ctDist, correlation$list_tt_expData,
                                                  correlation$bg_ctCor_data, correlation$bg_anova,
                                                  method = corr_settings$method,
                                                  bg_anova_alpha = corr_settings$bg_anova_alpha,
                                                  bg_ctCor_alpha = corr_settings$bg_ctCor_alpha,
                                                  gene_ctCor_alpha = corr_settings$gene_ctCor_alpha) ### Time-concuming (~40min per cell type)
      incProgress(1, detail = "signals filtration")
      correlation$signals <- signal_filter(correlation$signal_Stat, pValue = corr_settings$pValue,
                                           threshold = corr_settings$threshold)
      incProgress(1, detail = "connection to graph")
      correlation$signals_between_clusters <- get_signals_between_clusters(correlation$signals)
      correlation$signals_in_cluster <- get_signals_in_cluster(correlation$signals)
      incProgress(1)

    })

  })

  ##### Drawing the reactive and interactive graph network2 with clusters
  output$network_corr <- renderVisNetwork({
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

  #### Choice a focused node from the network2
  focus_corr <- reactive({
    list(input$current_node_id2, input$current_edges_id2)
  })

  #### Create the table of correlation signal  by a focused node from the network2
  output$gene_cor_table <- DT::renderDataTable(DT::datatable({
    if(is.null(correlation$signals_between_clusters)){return(NULL)}
    if(is.null(correlation$signals_in_cluster)){return(NULL)}
    focus_corr_data <- NULL
    signals_between_clusters_top <- get_signals_between_clusters_top(correlation$signals_between_clusters, fcs_data$use_markers, adj_pValue = 0.01)
    signals_in_cluster_top <- get_signals_in_cluster_top(correlation$signals_in_cluster, fcs_data$use_markers, adj_pValue = 0.01)
    if(is.null(focus_corr()[[1]]) & length(focus_corr()[[2]]==1)){
      target_clusters <- clusterisation$edges[clusterisation$edges$id == focus_corr()[[2]], c("from", "to")]
      focus_corr_data <- signals_between_clusters_top[(signals_between_clusters_top$signaling_cluster %in% target_clusters) &
                                              (signals_between_clusters_top$targetet_cluster %in% target_clusters),]
      if(nrow(focus_corr_data) < 2){
        focus_corr_data <- signals_between_clusters[(signals_between_clusters$signaling_cluster %in% target_clusters) &
                                            (signals_between_clusters$targetet_cluster %in% target_clusters),]
      }
    }
    if(length(focus_corr()[[1]])==1){
      focus_corr_data <- signals_in_cluster_top[signals_in_cluster_top$cluster == focus_corr()[[1]],]
      if (nrow(focus_corr_data) < 2) {focus_corr_data <- signals_in_cluster[signals_in_cluster$cluster == focus_corr()[[1]],]}
    }
    if(is.null(focus_corr()[[2]]) & is.null(focus_corr()[[1]])){
      focus_corr_data <- rbind(signals_between_clusters_top, data.frame(signaling_cluster = signals_in_cluster_top$cluster,
                                                                  targetet_cluster = signals_in_cluster_top$cluster,
                                                                  gene_in_target_cluster = signals_in_cluster_top$gene,
                                                                  signals_in_cluster_top[,-c(1, 2)]))
    }
    correlation$focus_corr_data <- focus_corr_data
    signals_in_cluster_top
    return(focus_corr_data)
  }))

  ##### Download abundance corelation table
  output$dwn_gene_cor_acorr <- downloadHandler(
    filename = function() {paste0("abund_", input$dwn_gene_cor_acorr_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_gene_cor_acorr_ext
      if(is.null(mode)){mode <- 'signals'}
      if(mode == 'signals_in_cluster'){write.csv(correlation$signals_in_cluster, file)}
      if(mode == 'signals_between_clusters'){write.csv(correlation$signals_between_clusters, file)}
      if(mode == 'signals'){write.csv(correlation$signals, file)}
      if(mode == 'focus_corr_data'){
        if(is.null(correlation$focus_corr_data)){return(NULL)}
        write.csv(correlation$focus_corr_data, file)}
    }
  )

  ########################
  ## Marker correlation ##
  ########################

  mk_corr_settings <- reactiveValues(method = "spearman", pValue = 0.01, threshold = 0.1,
                                  bg_anova_alpha = 0.01, bg_ctCor_alpha = 0.01, gene_ctCor_alpha = 0.01)

  observeEvent(input$simple_mk_corr_settings,{
    output$corr_analysis_settings_ui <- NULL
    mk_corr_settings$method <- "spearman"
    mk_corr_settings$pValue <- 0.1
    mk_corr_settings$threshold <- 0
    mk_corr_settings$bg_anova_alpha <- 0.01
    mk_corr_settings$bg_ctCor_alpha <- 0.01
    mk_corr_settings$gene_ctCor_alpha <- 0.01
  })

  observeEvent(input$advanced_mk_corr_settings, {
    if(!is.null(input$mk_corr_method)){mk_corr_settings$method <- input$mk_corr_method}
    if(!is.null(input$mk_corr_pValue)){mk_corr_settings$pValue <- input$mk_corr_pValue}
    if(!is.null(input$mk_corr_threshold)){mk_corr_settings$threshold <- input$mk_corr_threshold}
    if(!is.null(input$mk_corr_bg_anova_alpha)){mk_corr_settings$bg_anova_alpha <- input$mk_corr_bg_anova_alpha}
    if(!is.null(input$mk_corr_bg_ctCor_alpha)){mk_corr_settings$bg_ctCor_alpha <- input$mk_corr_bg_ctCor_alpha}
    if(!is.null(input$mk_corr_gene_ctCor_alpha)){mk_corr_settings$gene_ctCor_alpha <- input$mk_corr_gene_ctCor_alpha}
  })

  observeEvent(input$mk_corr_analysis, {
    options(warn=-1)
    print("Correlation analysis started...")
    withProgress(message = "Marker correlations", min =0, max = 7, value = 0,{
      correlation$list_expData <- get_list_expData(fcs_data$fcs_raw)
      correlation$list_tt_expData <- tt_sample_aggregator(correlation$list_cell_ctDist, correlation$list_expData)
      incProgress(1, detail = "background correlations")
      ##### Background contrast computing
      if(length(fcs_data$fcs_raw) < 4){
        correlation$bg_corr_mk <- get_bg_naive_corr_land_mk(list_expData = correlation$list_expData,
                                                            use_markers = fcs_data$use_markers, method = mk_corr_settings$method)
        correlation$contrast_corr_land_mk <- get_contrast_naive_corr_land_mk(list_cell_ctDist = correlation$list_cell_ctDist,
                                                                             list_expData = correlation$list_expData,
                                                                             use_markers = fcs_data$use_markers,
                                                                             method = mk_corr_settings$method)
      }
      if(length(fcs_data$fcs_raw) >= 4){
        correlation$bg_corr_mk <- get_bg_corr_land_mk(list_expData = correlation$list_expData,
                                                      use_markers = fcs_data$use_markers, method = mk_corr_settings$method)
        correlation$contrast_corr_land_mk <- get_contrast_corr_land_mk(list_cell_ctDist = correlation$list_cell_ctDist,
                                                                       list_expData = correlation$list_expData,
                                                                       use_markers = fcs_data$use_markers,
                                                                       method = mk_corr_settings$method)
      }
      incProgress(1, detail = "anova analysis")
      correlation$anova_mk <- get_anova_mk(list_cell_ctDist = correlation$list_cell_ctDist,
                                           list_expData = correlation$list_expData, use_markers = fcs_data$use_markers)
      incProgress(1, detail = "complex correlations")
      ##### Analys for each signalling call type
      correlation$mk_to_mk_corr <- htest_data_extractor_mk(list_tt_expData = correlation$list_tt_expData,
                                                           use_markers = fcs_data$use_markers,
                                                           method = mk_corr_settings$method)
      incProgress(1, detail = "background comparison")
      sig_summary <- comparator_mk(gene_to_gene_cor = correlation$mk_to_mk_corr, anova_mk = correlation$anova_mk,
                                   bg_corr_mk = correlation$bg_corr_mk, contrast_corr_land_mk = correlation$contrast_corr_land_mk,
                                   anova_mk_alpha = mk_corr_settings$bg_anova_alpha,
                                   bg_corr_mk_alpha = mk_corr_settings$bg_ctCor_alpha,
                                   contrast_corr_land_mk_alpha = mk_corr_settings$bg_ctCor_alpha)
      incProgress(1, detail = "signal filtering")
      correlation$signals_mk <- filter_mk(gene_to_gene_cor = correlation$mk_to_mk_corr, sig_summary = sig_summary,
                                          threshold = mk_corr_settings$threshold, pValue = mk_corr_settings$pValue)
      incProgress(1, detail = "signal graph association")
      ##### Dividing signals on signals between and in clusters
      correlation$signal_in_cluster_mk <- get_signal_in_cluster_mk(correlation$signals_mk)
      correlation$signal_in_cluster_mk$marker_1 <- names(fcs_data$use_markers[
        match(correlation$signal_in_cluster_mk$marker_1, fcs_data$use_markers)])
      correlation$signal_in_cluster_mk$marker_2 <- names(fcs_data$use_markers[
        match(correlation$signal_in_cluster_mk$marker_2, fcs_data$use_markers)])

      correlation$signal_between_cluster_mk <- get_signal_between_cluster_mk(correlation$signals_mk)
      correlation$signal_between_cluster_mk$signaling_marker <- names(fcs_data$use_markers[
        match(correlation$signal_between_cluster_mk$signaling_marker, fcs_data$use_markers)])
      correlation$signal_between_cluster_mk$target_marker <- names(fcs_data$use_markers[
        match(correlation$signal_between_cluster_mk$target_marker, fcs_data$use_markers)])
      incProgress(1)
    })
    print("... Correlation analysis finished")
  })

  #### Choice a focused node from the network_mk
  focus_mk_corr <- reactive({
    return(list(input$current_node_id2, input$current_edges_id2))
  })

  #### Create the table of correlation signal by a focused node from the network_mk
  output$marker_mk_corr_table <- DT::renderDataTable(DT::datatable({
    if(is.null(correlation$signal_in_cluster_mk) & is.null(correlation$signal_between_cluster_mk)){return(NULL)}
    focus_mk_corr_data <- NULL
    if(is.null(focus_mk_corr()[[1]]) & length(focus_mk_corr()[[2]]==1)){
      target_clusters <- clusterisation$edges[clusterisation$edges$id == focus_mk_corr()[[2]], c("from", "to")]
      focus_mk_corr_data <- correlation$signal_between_cluster_mk[(
        correlation$signal_between_cluster_mk$signaling_cluster %in% target_clusters) &
          (correlation$signal_between_cluster_mk$target_cluster %in% target_clusters),]

    }
    if(length(focus_mk_corr()[[1]])==1){
      focus_mk_corr_data <- correlation$signal_in_cluster_mk[correlation$signal_in_cluster_mk$cluster == focus_mk_corr()[[1]],]
    }
    correlation$focus_mk_corr_data <- focus_mk_corr_data
    return(focus_mk_corr_data)
  }))

  ##### Download marker corelation table
  output$dwn_marker_table_mk_corr <- downloadHandler(
    filename = function() {paste0("mk_",input$dwn_marker_table_mk_corr_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_marker_table_mk_corr_ext
      if(is.null(mode)){mode <- 'signals'}
      if(mode == 'signals_in_cluster'){write.csv(correlation$signal_in_cluster_mk, file)}
      if(mode == 'signals_between_clusters'){write.csv(correlation$signal_between_cluster_mk, file)}
      if(mode == 'signals'){write.csv(correlation$signals_mk, file)}
      if(mode == 'focus_corr_data'){
        if(is.null(correlation$focus_mk_corr_data)){return(NULL)}
        write.csv(correlation$focus_mk_corr_data, file)}
    }
  )

  ########################
  ####  Cross-panel   ####
  ########################

  ##### Creating reactive Values with panel data
  crosspanel <- reactiveValues(md_set = list(), fcs_raw_set = list(), panel_set = list(), use_marker_set = list(),
                               cell_clustering_list_set = list())

  shinyFileChoose(input, 'choose_fcs_cross_p', roots=roots, filetypes=c('', 'fcs'))
  observeEvent(input$butt_upload_cross_p, {
    withProgress(message = "Clusters from fcs files", min =0, max = 4, value = 0,{
      panel_name <- input$panel_name_cross_p
      if(is.null(panel_name) | panel_name == ""){
        f_name <- basename(parseFilePaths(roots, input$choose_fcs_cross_p)$datapath[1])
        panel_name <- gsub(".fcs", "", f_name)
        if(grepl("CyBr_", f_name)){panel_name <- strsplit(f_name, split = "_")[[1]][2]}
      }
      crosspanel$md_set[[panel_name]] <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_cross_p)$datapath)
      incProgress(1, detail = "Upload data" )
      crosspanel$fcs_raw_set[[panel_name]] <- get_fcs_raw(crosspanel$md_set[[panel_name]])
      incProgress(1, detail = "Extraction panel")
      crosspanel$panel_set[[panel_name]] <- get_fcs_panel(crosspanel$fcs_raw_set[[panel_name]])
      incProgress(1, detail = "Delete background markers" )
      crosspanel$use_marker_set[[panel_name]] <- get_use_marker(crosspanel$panel_set[[panel_name]])
      incProgress(1, detail = "Extract cluster info" )
      pattern <- "clust"
      if(!is.null(input$extr_clust_pattern_cross_p)){pattern <- input$extr_clust_pattern_cross_p}
      crosspanel$cell_clustering_list_set[[panel_name]] <- get_fcs_cluster_annotation(crosspanel$fcs_raw_set[[panel_name]], pattern = pattern)
    })
  })

  ##### Show the overlapped samples of the panels
  observe({
    if(length(crosspanel$md_set) == 0){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    crosspanel$use_samples <- rownames(panel_content)[apply(panel_content, 1, function(x) {all(x == "presented")})]
  })

  ##### Show the overlapped markers of the panels
  observe({
    if(length(crosspanel$use_marker_set) == 0){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    crosspanel$use_common_mk <- rownames(mk_panel_content)[apply(mk_panel_content, 1, function(x) {all(x == "presented")})]
  })

  ##### Create ui for removing panels
  output$remove_panel_cross_p_ui <- renderUI({
    if(length(crosspanel$fcs_raw_set) < 1){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Remove panels"),
             selectInput('panels_to_remove_cross_p', label = h5("Choose panels"),
                         choices = names(crosspanel$fcs_raw_set), multiple = TRUE),
             actionButton('remove_cross_p', label = "Remove")
      )
    )
  })

  ##### Removing panels
  observeEvent(input$remove_cross_p, {
    if(length(crosspanel$fcs_raw_set) < 1){return(NULL)}
    if(is.null(input$panels_to_remove_cross_p)){return(NULL)}
    stay_panel <- !(names(crosspanel$fcs_raw_set) %in% input$panels_to_remove_cross_p)
    crosspanel$md_set <- crosspanel$md_set[stay_panel]
    crosspanel$fcs_raw_set <- crosspanel$fcs_raw_set[stay_panel]
    crosspanel$panel_set <- crosspanel$panel_set[stay_panel]
    crosspanel$use_marker_set <- crosspanel$use_marker_set[stay_panel]
    crosspanel$cell_clustering_list_set <- crosspanel$cell_clustering_list_set[stay_panel]
  })

  ##### Show the sample content of the panels
  output$content_cross_p <- renderPlot({
    if(length(crosspanel$md_set) < 1){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    ComplexHeatmap::Heatmap(panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top",
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  ##### Show the marker content of the panels
  output$mk_content_cross_p <- renderPlot({
    if(length(crosspanel$use_marker_set) < 1){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    ComplexHeatmap::Heatmap(mk_panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top", row_names_gp = grid::gpar(fontsize = 8),
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  ##### UI for choosing p-valueadjusting method for correlation cross-panel analysis
  output$abund_corr_padj_method_cross_p_ui <- renderUI({
    if(input$abund_corr_p_mode_cross_p == 'padj'){return((
      selectInput('abund_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                    "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  ##### Make abundance cross panel correlation analysis
  observeEvent(input$abund_corr_cross_p, {
    crosspanel$abund_cross_p <- get_abund_cross_panel(cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                           use_samples = crosspanel$use_samples)
    crosspanel$abund_corr_info <- get_matrices_corr_info(crosspanel$abund_cross_p, method_cortest = input$abund_corr_method_cross_p,
                                                         method_padj = input$abund_corr_padj_method_cross_p)
  })

  #### Drawing abundance correlation plot for cross panel analysis
  output$plot_abund_corr_cross_p <- renderPlot({
    if(is.null(crosspanel$abund_corr_info))return(NULL)
    signif_matrix <- crosspanel$abund_corr_info$padj
    if(input$abund_corr_p_mode_cross_p == 'pval'){signif_matrix <- crosspanel$abund_corr_info$p_value}
    corrplot::corrplot(crosspanel$abund_corr_info$corr_coef, method = "square", outline = T, order="hclust", type = 'lower',
             p.mat = signif_matrix,
             insig = 'label_sig', sig.level = c(.001, .01, .05), pch.cex = 0.7, pch.col = 'white',
             addgrid.col = "darkgray", tl.col = "black", tl.cex = 1, cl.cex = 1,
             col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
    gridGraphics::grid.echo()
    plots$abun_corr_cross_p <- grid::grid.grab()
    return(plots$abun_corr_cross_p)
  })

  ##### Download abundance correlation plot for cross panel analysis
  output$dwn_abund_corr_cross_p <- downloadHandler(
    filename = function() {
      ext <- input$dwn_abund_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Abunadnce_correlations_cross_panel", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_abund_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      #ggsave(file, plot = plots$abun_corr_cross_p, device = ext)
    }
  )

  ##### Download abundance corelation tables  for cross-panel analusis
  output$dwn_table_abund_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_", input$dwn_table_abund_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_abund_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'abund_data'){write.csv(crosspanel$abund_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$abund_corr_info), file)}
    }
  )

  ##### UI to choose marker for expressing cell fraction correlations
  output$mk_exp_cell_f_cross_p_ui <- renderUI({
    if(is.null(crosspanel$use_common_mk)){return(NULL)}
    selectInput('mk_exp_cell_f_cross_p', label = h5("Marker for cross-panel analysis"),
                choices = crosspanel$use_common_mk,
                selected = 1)
  })

  ##### UI to choose the threshold to expressing cell fraction allocation
  output$threshold_exp_cell_f_cross_p_ui <- renderUI({
    if(input$exp_cell_f_corr_method != 'threshold'){return(NULL)}
    if(is.null(crosspanel$fcs_raw_set)){return(NULL)}
    if(is.null(input$mk_exp_cell_f_cross_p)){return(NULL)}
    mk_exp_data_cross_p <- get_mk_exp_data_cross_p(fcs_raw_set = crosspanel$fcs_raw_set , use_samples = crosspanel$use_samples,
                                                   use_marker_set = crosspanel$use_marker_set, target_mk = input$mk_exp_cell_f_cross_p)
    mk_exp_data_cross_p <<- mk_exp_data_cross_p
    numericInput('threshold_exp_cell_f_cross_p', "threshold",min = min(mk_exp_data_cross_p), max = max(mk_exp_data_cross_p),
                 value = stats::median(mk_exp_data_cross_p), step = 0.1)
  })

  ##### UI for choosing p-value adjusting method for exp cell fraction cross-panel analysis
  output$exp_cell_f_corr_padj_method_cross_p_ui <- renderUI({
    if(input$exp_cell_f_corr_p_mode_cross_p == 'padj'){return((
      selectInput('exp_cell_f_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                              "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  ##### Make expressing cell fraction cross panel correlation analysis
  observeEvent(input$exp_cell_f_corr_cross_p, {
    exp_cell_f_corr_method <- input$exp_cell_f_corr_method
    if(is.null(exp_cell_f_corr_method)){exp_cell_f_corr_method <- 'clustering'}
    if(exp_cell_f_corr_method == 'clustering'){
      crosspanel$exp_cell_f_cross_p <- get_clustering_exp_cell_f_cross_panel(fcs_raw_set = crosspanel$fcs_raw_set,
                                                                  cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                                  use_samples = crosspanel$use_samples,
                                                                  use_marker_set = crosspanel$use_marker_set,
                                                                  target_mk = input$mk_exp_cell_f_cross_p,
                                                                  min_cell_number = input$exp_cell_f_corr_min_cell_cross_p)
    }
    if(exp_cell_f_corr_method == 'threshold'){
      crosspanel$exp_cell_f_cross_p <- get_threshold_exp_cell_f_cross_panel(fcs_raw_set = crosspanel$fcs_raw_set,
                                                                             cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                                             use_samples = crosspanel$use_samples,
                                                                             use_marker_set = crosspanel$use_marker_set,
                                                                             target_mk = input$mk_exp_cell_f_cross_p,
                                                                             threshold = input$threshold_exp_cell_f_cross_p,
                                                                            min_cell_number = input$exp_cell_f_corr_min_cell_cross_p)
    }
    crosspanel$exp_cell_f_corr_info <- get_matrices_corr_info(crosspanel$exp_cell_f_cross_p,
                                                              method_cortest = input$exp_cell_f_corr_method_cross_p,
                                                              method_padj = input$exp_cell_f_corr_padj_method_cross_p)
  })

  #### Drawing abundance correlation plot for cross panel analysis
  output$plot_exp_cell_f_corr_cross_p <- renderPlot({
    if(is.null(crosspanel$exp_cell_f_corr_info))return(NULL)
    signif_matrix <- crosspanel$exp_cell_f_corr_info$padj
    if(input$exp_cell_f_corr_p_mode_cross_p == 'pval'){signif_matrix <- crosspanel$exp_cell_f_corr_info$p_value}
    corrplot::corrplot(crosspanel$exp_cell_f_corr_info$corr_coef, method = "square", outline = T, order="hclust", type = 'lower',
                       p.mat = signif_matrix,
                       insig = 'label_sig', sig.level = c(.001, .01, .05), pch.cex = 0.7, pch.col = 'white',
                       addgrid.col = "darkgray", tl.col = "black", tl.cex = 1, cl.cex = 1,
                       col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
    gridGraphics::grid.echo()
    plots$exp_cell_f_corr_cross_p <- grid::grid.grab()
    return(plots$exp_cell_f_corr_cross_p)
  })

  ##### Drawing expression density plot for a marker in cross-panel analysis
  output$mk_density_plot_cross_p <- renderPlot({
    if(is.null(crosspanel$fcs_raw_set)){return(NULL)}
    if(is.null(input$mk_exp_cell_f_cross_p)){return(NULL)}
    mk_exp_data_cross_p <- get_mk_exp_data_cross_p(fcs_raw_set = crosspanel$fcs_raw_set , use_samples = crosspanel$use_samples,
                                                   use_marker_set = crosspanel$use_marker_set, target_mk = input$mk_exp_cell_f_cross_p)
    ggplot(data.frame(expr = mk_exp_data_cross_p), aes(x = expr)) +
      geom_density(fill = 'black')
  })


  ##### Download abundance correlation plot for cross panel analysis
  output$dwn_exp_cell_f_corr_cross_p <- downloadHandler(
    filename = function() {
      ext <- input$dwn_exp_cell_f_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      paste0("Expr_cell_f_correlations_", input$mk_exp_cell_f_cross_p, "_cross_panel.", ext) },
    content = function(file) {
      ext <- input$dwn_exp_cell_f_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      #ggsave(file, plot = plots$exp_cell_f_corr_cross_p, device = ext)
    }
  )

  ##### Download expressing cell fraction corelation tables for cross-panel analusis
  output$dwn_table_exp_cell_f_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_",input$mk_exp_cell_f_cross_p,"_", input$dwn_table_exp_cell_f_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_exp_cell_f_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'exp_cell_f_data'){write.csv(crosspanel$exp_cell_f_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$exp_cell_f_corr_info), file)}
    }
  )



  ########################
  ####      GDB       ####
  ########################

  observeEvent(input$neo4j_export,{
    user <- "neo4j"
    password <- "password"
    if(!is.null(input$user_neo4j)){user <- input$user_neo4j}
    if(!is.null(input$password_neo4j)){password <- input$password_neo4j}
    gdb <- get_neo_api(user = input$user_neo4j, password = input$password_neo4j)
    ping_answer <- neo_api_ping(gdb)
    if(!is.null(ping_answer)){
      showModal(modalDialog(title = "Error with Neo4j", ping_answer, easyClose = TRUE))
      return(NULL)
    }
    if(!is.null(fcs_data$fcs_raw)){add_sample_GDB(fcs_data$fcs_raw, gdb)}
    if(!is.null(fcs_data$use_markers)){add_marker_GDB(fcs_data$use_markers, gdb)}
    if(!is.null(clusterisation$cell_clustering)){
      add_cluster_GDB(clusterisation$cell_clustering, gdb)
      add_populatio_GDB(clusterisation$cell_clustering_list, gdb)
      add_observation_GDB(gdb)
      add_phenounite_GDB(gdb)
    }
    if(!is.null(correlation$signals_between_clusters)){add_signals_between_clusters_GDB(correlation$signals_between_clusters, gdb)}
    if(!is.null(correlation$signals_in_cluster)){add_signals_in_cluster_GDB(correlation$signals_in_cluster, gdb)}
  })


}

