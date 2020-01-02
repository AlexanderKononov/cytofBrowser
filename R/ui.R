#library(shiny)

#library(visNetwork)
#library(d3heatmap)
#library(shinyFiles)

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
                            shinyFilesButton('fcs_files', label='File fcs select', title='Please select fcs file', multiple=TRUE),
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
                                                  label = h4("Number of cell to print"), value = 2000)
                                     ),
                              column(4,
                                     selectInput("method_plot_data_preparation", label = h4("Visualisation method"),
                                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                                 selected = "tSNE")
                                     ),
                              column(2, actionButton("redraw", label = "Redraw"))
                            ),
                            fluidRow(
                              column(8,
                                     plotOutput("tSNE_plot_data_preparation")
                                     ),
                              column(4,
                                     plotOutput("smpl_hist_preparation")
                                     )
                            ),
                            fluidRow(
                              column(6,
                                     uiOutput("mk_target_data_preparation_ui"),
                                     uiOutput("mk_subset_data_preparation_ui"),
                                     actionButton("exclud_mk_button", label = "Exclude markers"),
                                     hr(),
                                     h5("Analyzed markers"),
                                     verbatimTextOutput("mk_rested_data_preparation")
                                     ),
                              column(6,
                                     plotOutput("mk_hist_data_preparation")
                                     )
                            )
                          )
                        )
               ),
               ################################################################### Tab 2
               tabPanel("Clusterisation",
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("mode_k_choice", label = h3("Choose of number of clusters"),
                                         choices = list("Automatically detect optimum" = 1, "Manually choose" = 2),
                                         selected = 1),
                            conditionalPanel(
                              condition = "input.mode_k_choice == 1",
                              numericInput("rate_var_explan", label = h3("Rate of explained variance"), value = 0.9),
                              numericInput("maxK", label = h3("Max number of clusters"), value = 20)
                            ),
                            conditionalPanel(
                              condition = "input.mode_k_choice == 2",
                              numericInput("k", label = h3("Choose number of clusters"), value = 8)
                            ),
                            uiOutput("mk_subset_clusterisation_ui"),
                            #selectInput("exclude_mk_clusterisation", label = "exclude markers from clusterisation",
                            #            choices = names(fcs_data$use_markers),
                            #            multiple = TRUE),
                            actionButton("start_clusterization", label = "Clusterization")
                          ),
                          mainPanel(
                            fluidRow(
                              column(3,
                                     numericInput("n_cell_plot_clasterisation",
                                                  label = h4("Number of cell to print"), value = 2000)
                              ),
                              column(3,
                                     selectInput("method_plot_clasterisation", label = h4("Visualisation method"),
                                                 choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),
                                                 selected = "UMAP")
                              ),
                              column(2, actionButton("redraw_clasterisation", label = "Redraw"))
                            ),
                            fluidRow(
                              plotOutput("tSNE_plot1_clust", click = "plot_click")
                            ),
                            fluidRow(
                              column(6,
                                     uiOutput("mk_target_clusterisation_ui"),
                                     uiOutput("mergeing_clusterisation_ui"),
                                     uiOutput("rename_clusterisation_ui")
                              ),
                              column(6,
                                     visNetworkOutput("network")
                              )
                            )
                          )
                        )
               ),
               ################################################################### Tab 3
               tabPanel("Gene expression",
                        fluidRow(
                          d3heatmapOutput('cluster_heatmap')
                        ),
                        fluidRow(
                          column(6,
                                 uiOutput('mk_deconvol_gene_expression_ui')
                          ),
                          column(6,
                                 plotOutput("deconvol_expr")
                          )
                        )
               ),
               ################################################################### Tab 4
               tabPanel("Abundance correlation",
                        fluidRow(
                          column(6,
                                 plotOutput("abun_cor_plot", click = "abun_cor_click")
                          ),
                          column(6,
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
                          DT::dataTableOutput("gene_cor_table")
                        )
               )

    )
  )

  shinyApp(ui = cytofCore_ui, server = cytofCore_server)
}

