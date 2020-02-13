
#' The main function to run cutofCore graphic interface
#' @description Run graphical user interface for cytofanalyzer.
#' The function runs the Shiny App in browser. The current package
#' was created with assumption that it will be run with GUI mode
#' as Shiny App, That is why the function is the most appropriate
#' way to use the cytofanalyzer.
#'
#' @return The function runs the Shiny App in browser.
#' @import shiny shinyFiles visNetwork d3heatmap
#' @export
#'
#' @examples
#' \dontrun{
#' cytofCoreGUI()
#' }
cytofCoreGUI <-function(){
  cytofCore_ui <- fluidPage(
    navbarPage("CyTOF pipeline",
               ################################################################### Tab 1
               tabPanel("Data preparation",
                        sidebarLayout(
                          sidebarPanel(
                            shinyFilesButton('fcs_files', label='Select fcs files', title='Please select fcs files', multiple=TRUE),
                            checkboxGroupInput("transformation_list", label = h4("Transformations"),
                                               choices = list("asinh" = "asinh", "outlier by quartile" = "outlier_by_quantile"),
                                               selected = c("asinh", "outlier_by_quantile")),
                            conditionalPanel(
                              condition = "input.transformation_list.includes('asinh')",
                              numericInput("cofactor", label = h4("Cofactor for dividing"), value = 5)
                            ),
                            conditionalPanel(
                              condition = "input.transformation_list.includes('outlier_by_quantile')",
                              numericInput("quantile", label = h4("Quantile for outlier removing"), value = 0.01)
                            ),
                            hr(),
                            actionButton("fcs_upload", label = "Upload")
                          ),
                          mainPanel(
                            fluidRow(
                              column(4,
                                     numericInput("n_cell_plot_data_preparation",
                                                  label = h4("Cell fraction to display"), value = 0.5, step = 0.1)
                                     ),
                              column(4,
                                     selectInput("method_plot_data_preparation", label = h4("Visualisation method"),
                                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                                 selected = "tSNE")
                                     ),
                              column(2, actionButton("redraw", label = "Redraw")),
                              column(2,
                                     actionButton("draw_simple_data_prep", label = "Simple"),
                                     actionButton("draw_advance_data_prep", label = "Advance")
                                     )
                            ),
                            fluidRow(uiOutput('dvance_data_prep_ui')),
                            fluidRow(
                              column(8,
                                     plotOutput('scatter_plot_data_preparation'),
                                     fluidRow(
                                       column(2, downloadButton('dwn_scatter_dp', "")),
                                       column(4,
                                              selectInput('dwn_scatter_dp_ext', label = NULL,
                                                          choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                         'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                       )
                                     )
                                     ),
                              column(4,
                                     plotOutput("smpl_hist_preparation"),
                                     fluidRow(
                                       column(4, downloadButton('dwn_smpl_hist_dp', "")),
                                       column(7,
                                              selectInput('dwn_smpl_hist_dp_ext', label = NULL,
                                                          choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                         'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                       )
                                     )
                                     )
                            ),
                            fluidRow(
                              column(6,
                                     uiOutput("mk_target_data_preparation_ui"),
                                     uiOutput("mk_subset_data_preparation_ui")
                                     ),
                              column(6,
                                     plotOutput('mk_hist_data_preparation'),
                                     fluidRow(
                                       column(3, downloadButton('dwn_mk_hist_dp', "")),
                                       column(5,
                                              selectInput('dwn_mk_hist_dp_ext', label = NULL,
                                                          choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                         'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                       )
                                     )
                                     )
                            ),
                            fluidRow(
                              column(6,
                                     h4("Analysed markers"),
                                     verbatimTextOutput('mk_rested_data_preparation')
                                     ),
                              column(6,
                                     h4("Excluded markers"),
                                     verbatimTextOutput('mk_excluded_data_preparation')
                                     )
                            )
                          )
                        )
               ),
               ################################################################### Tab 2
               tabPanel("Clustering",
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("mode_k_choice", label = h4("Choose of number of clusters"),
                                         choices = list("Automatically detect optimum" = 1, "Manually choose" = 2),
                                         selected = 1),
                            conditionalPanel(
                              condition = "input.mode_k_choice == 1",
                              numericInput("rate_var_explan", label = h4("Rate of explained variance"), value = 0.9),
                              numericInput("maxK", label = h4("Max number of clusters"), value = 20)
                            ),
                            conditionalPanel(
                              condition = "input.mode_k_choice == 2",
                              numericInput("k", label = h4("Choose number of clusters"), value = 8)
                            ),
                            uiOutput("mk_subset_clusterisation_ui"),
                            actionButton("start_clusterization", label = "Clustering")
                          ),
                          mainPanel(
                            fluidRow(
                              column(3,
                                     numericInput("n_cell_plot_clasterisation",
                                                  label = h4("Cell fraction to display"), value = 0.5, step = 0.1)
                              ),
                              column(3,
                                     selectInput("method_plot_clasterisation", label = h4("Visualisation method"),
                                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                                 selected = "UMAP")
                              ),
                              column(2, actionButton("redraw_clasterisation", label = "Redraw")),
                              column(2,
                                     actionButton("draw_simple_claster", label = "Simple"),
                                     actionButton("draw_advance_claster", label = "Advance")
                              )
                            ),
                            fluidRow(uiOutput('dvance_cluster_ui')),
                            fluidRow(
                              column(8,
                                     plotOutput('scatter_plot_clust', click = "plot_click"),
                                     fluidRow(
                                       column(2, downloadButton('dwn_scatter_clust', "")),
                                       column(4,
                                              selectInput('dwn_scatter_clust_ext', label = NULL,
                                                          choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                         'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                       )
                                     )
                                     ),
                              column(4,
                                     plotOutput('abundance_clust'),
                                     fluidRow(
                                       column(4, downloadButton('dwn_abundance_clust', "")),
                                       column(7,
                                              selectInput('dwn_abundance_clust_ext', label = NULL,
                                                          choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                         'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                       )
                                     )
                                     )
                            ),
                            fluidRow(
                              column(6,
                                     uiOutput("mk_target_clusterisation_ui"),
                                     uiOutput("mergeing_clusterisation_ui"),
                                     uiOutput("rename_clusterisation_ui")
                              ),
                              column(6,
                                     sliderInput('edges_threshold_clusterisation', "Edge weight threshold for graph",
                                                 min =0, max = 1, value = 0.5, step = 0.01),
                                     sliderInput('gravity_clusterisation', "Gravity for graph",
                                                 min = -100, max = 0, value = -40, step = 1),
                                     visNetworkOutput("network")
                              )
                            )
                          )
                        )
               ),
               ################################################################### Tab 3
               tabPanel("Marker level",
                        fluidRow(
                          column(3,
                                 selectInput("method_summarize_expression", label = h4("summarise method"),
                                             choices = list('median' = "median", 'mean' = "mean"),
                                             selected = 'median')
                          ),
                          column(2, actionButton("redraw_expression", label = "Redraw"))
                        ),
                        fluidRow(
                          d3heatmapOutput('cluster_heatmap')
                        ),
                        fluidRow(
                          column(6,
                                 uiOutput('mk_deconvol_gene_expression_ui')
                          ),
                          column(6,
                                 plotOutput("deconvol_expr"),
                                 fluidRow(
                                   column(2, downloadButton('dwn_deconvol_expr', "")),
                                   column(4,
                                          selectInput('dwn_deconvol_expr_ext', label = NULL,
                                                      choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                     'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp"))
                                   )
                                 )
                          )
                        ),
                        fluidRow(
                          column(3, actionButton('drawn_cluster_hm_expr', label = "Drawing heatmap")),
                          column(3, downloadButton('dwn_drawn_cluster_hm_expr', "")),
                          column(3,
                                 selectInput('dwn_drawn_cluster_hm_expr_ext', label = NULL,
                                             choices = list("pdf" = 'pdf', "jpeg" = 'jpeg', 'png' = "png"))
                          )
                        ),
                        fluidRow(plotOutput('cluster_hm_expr'))
               ),
               ################################################################### Tab 4
               tabPanel("Abundance correlation",
                        fluidRow(
                          column(6,
                                 plotOutput("abun_cor_plot", click = "abun_cor_click"),
                                 fluidRow(
                                   column(3, downloadButton('dwn_abund_corr', "")),
                                   column(5,
                                          selectInput('dwn_abund_corr_ext', label = NULL,
                                                      choices = list("pdf" = 'pdf', "jpeg" = 'jpeg', 'png' = "png"))
                                   )
                                 )
                          ),
                          column(6,
                                 sliderInput('edges_threshold_abund_corr', "Edge weight threshold for graph",
                                             min =0, max = 1, value = 0.5, step = 0.01),
                                 sliderInput('gravity_abund_corr', "Gravity for graph",
                                             min = -100, max = 0, value = -40, step = 1),
                                 visNetworkOutput("network2")
                          )
                        ),
                        fluidRow(
                          column(3,
                                 h4("Correlation settings:")
                          ),
                          column(3,
                                 actionButton("simple_corr_settings", "Simple"),
                                 actionButton("advanced_corr_settings", "Advanced")
                          ),
                          column(3,
                                 actionButton('corr_analysis', label = "Start correlation analysis")
                          )
                        ),
                        fluidRow(
                          uiOutput("corr_analysis_settings_ui")
                        ),
                        fluidRow(
                          fluidRow(
                            column(2, downloadButton('dwn_gene_cor_acorr', "")),
                            column(4,
                                   selectInput('dwn_gene_cor_acorr_ext', label = NULL,
                                               choices = list("chosen correlations" = 'focus_corr_data',
                                                              "correlations within clusters" = 'signals_in_cluster',
                                                              "correlations between clusters" = 'signals_between_clusters',
                                                              "all correlations" = 'signals'))
                            )
                          ),
                          DT::dataTableOutput("gene_cor_table")
                        )
               ),
               ################################################################### Tab 5
               tabPanel("Marker correlation",
                        fluidRow(
                          column(6,
                                 sliderInput('edges_threshold_mk_corr', "Edge weight threshold for graph",
                                             min =0, max = 1, value = 0.5, step = 0.01),
                                 sliderInput('gravity_mk_corr', "Gravity for graph",
                                             min = -100, max = 0, value = -40, step = 1),
                                 visNetworkOutput("network_mk_corr")
                                 ),
                          column(6,
                                 visNetworkOutput('network_mk')
                                 )
                        ),
                        fluidRow(
                          column(3,
                                 h4("Correlation settings:")
                          ),
                          column(3,
                                 actionButton("simple_mk_corr_settings", "Simple"),
                                 actionButton("advanced_mk_corr_settings", "Advanced")
                          ),
                          column(3,
                                 actionButton('mk_corr_analysis', label = "Start correlation analysis")
                          )
                        ),
                        fluidRow(
                          uiOutput('mk_corr_analysis_settings_ui')
                        ),
                        fluidRow(
                          fluidRow(
                            column(2, downloadButton('dwn_marker_table_mk_corr', "")),
                            column(4,
                                   selectInput('dwn_marker_table_mk_corr_ext', label = NULL,
                                               choices = list("chosen correlations" = 'focus_corr_data',
                                                              "correlations within clusters" = 'signals_in_cluster',
                                                              "correlations between clusters" = 'signals_between_clusters',
                                                              "all correlations" = 'signals'))
                            )
                          ),
                          DT::dataTableOutput('marker_mk_corr_table')
                        )

               )

    )
  )

  shinyApp(ui = cytofCore_ui, server = cytofCore_server)
}
