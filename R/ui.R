
#' The main function to run cutofCore graphic interface
#' @description Run graphical user interface for cytofBrowser.
#' The function runs the Shiny App in browser. The current package
#' was created with assumption that it will be run with GUI mode
#' as Shiny App, That is why the function is the most appropriate
#' way to use the cytofBrowser.
#'
#' @return The function runs the Shiny App in browser.
#' @import shiny shinyFiles visNetwork d3heatmap shinydashboard shinyWidgets
#' @export
#'
#' @examples
#' \dontrun{
#' cytofBrowserGUI()
#' }
cytofBrowserGUI <-function(){
  cytofBrowser_ui <- dashboardPage(
    dashboardHeader(title = "cytofBrowser"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Data", tabName = 'data_processing', icon = icon("bars")),
        menuItem("Exploration", tabName = 'data_exploration', icon = icon("chart-area")),
        menuItem("Correlation", tabName = 'data_correlation', icon = icon("braille")),
        menuItem("Cross-panel", tabName = 'data_crosspanel', icon = icon("clone"))
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = 'data_processing',
                fluidRow(
                  infoBoxOutput('iBox_upload_dproc'),
                  infoBoxOutput('iBox_preproc_dproc'),
                  infoBoxOutput('iBox_clust_dproc')
                ),
                fluidRow(
                  tabBox(
                    tabPanel("Uploading",
                             shinyFilesButton('choose_fcs_dp', label='Select FCS files', title='Please select FCS files', multiple=TRUE),
                             hr(),
                             materialSwitch(inputId = 'extr_clust_dproc', label = "extract cluster info"),
                             conditionalPanel(
                               condition = "input.extr_clust_dproc == true",
                               textInput("extr_clust_pattern_dproc",
                                         label = h5("full or part column name with clusters info (for cytofBrowser and cytofkit : <cluster>)"), value = "cluster")
                             ),
                             hr(),
                             actionButton('butt_upload_dproc', label = "Upload")
                    ),
                    tabPanel("Transforming",
                             checkboxGroupInput("transformation_list", label = h4("Transformations"),
                                                choices = list("asinh" = 'asinh', "outlier by quartile" = 'outlier_by_quantile',
                                                               "extract cluster info" = 'extract_cluster_info'),
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
                             actionButton('butt_trans_dproc', label = "Transform")
                    ),
                    tabPanel("Markers",
                             uiOutput('mk_subset_dp_ui')
                    ),
                    tabPanel("Clustering",
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
                             actionButton('start_clusterization', label = "Clustering")
                    ),
                    tabPanel("Cluster management",
                             uiOutput("mergeing_clusterisation_ui"),
                             hr(),
                             uiOutput("rename_clusterisation_ui")
                    ),
                    tabPanel("Save",
                             h5("Saving the data for samples as FCS files with cluster-info"),
                             shinyDirButton('choose_panel_clust', "Folder choose", "Select a folder to save panel samples"),
                             hr(),
                             h5("Saved set of files can be used as panel data in a cross-panel analysiso"),
                             textInput("panel_name_clust", label = h4("Panel name"), value = "Panel1"),
                             hr(),
                             actionButton('dwn_panel_clust', label = "Save panel")

                    )
                  ),
                  uiOutput('scatter_plot_dp_ui'),
                  tabBox(
                    tabPanel("Cell number",
                             fluidRow(
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_smpl_hist_dp_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_smpl_hist_dp', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                              ),
                             plotOutput("smpl_hist_preparation")
                    ),
                    tabPanel("Expression",
                             fluidRow(
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_mk_hist_dp_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_mk_hist_dp', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                                )
                             ),
                             plotOutput('mk_density_plot_dp'),
                             uiOutput('mk_density_dp_ui')
                    ),
                    tabPanel("Markers",
                             h4("Analysed markers"),
                             verbatimTextOutput('mk_rested_dp'),
                             h4("Excluded markers"),
                             verbatimTextOutput('mk_excluded_dp')
                    ),
                    tabPanel("Abundance",
                             fluidRow(
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_abundance_clust_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_abundance_clust', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             plotOutput('abundance_clust')
                    )
                  ),
                 uiOutput('network_clust_ui')
                )
        ),

        # Second tab content
        tabItem(tabName = 'data_exploration',
                fluidRow(
                  box(
                    fluidRow(
                      column(1, actionBttn(inputId = "redraw_expression", style = "material-circle", color = "default" ,icon = icon("redo"))),
                      column(1,
                             dropdownButton(
                               tags$h4("Options of plotting"),
                               selectInput("method_summarize_expression", label = h4("summarise method"),
                                           choices = list('median' = "median", 'mean' = "mean"),
                                           selected = 'median'),
                               icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                      ),
                      column(1,
                             dropdownButton(
                               selectInput('dwn_drawn_cluster_hm_expr_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                               downloadButton('dwn_drawn_cluster_hm_expr', ""),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
                    d3heatmapOutput('cluster_heatmap'),
                    width = 12
                  ),
                  box(
                    fluidRow(
                      column(2,
                             dropdownButton(
                               selectInput('dwn_deconvol_expr_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                          'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                               downloadButton('dwn_deconvol_expr', ""),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
                    plotOutput("deconvol_expr"),
                    uiOutput('mk_deconvol_gene_expression_ui')
                  )
                )
        ),

        # Third tab content
        tabItem(tabName = 'data_correlation',
                fluidRow(
                  infoBoxOutput('iBox_abund_corr'),
                  infoBoxOutput('iBox_mk_corr')
                ),
                fluidRow(
                  box(
                    sliderInput('edges_threshold_abund_corr', "Edge weight threshold for graph",
                                min =0, max = 1, value = 0.5, step = 0.01),
                    sliderInput('gravity_abund_corr', "Gravity for graph",
                                min = -100, max = 0, value = -40, step = 1),
                    visNetworkOutput('network_corr')
                  ),
                  box(
                    fluidRow(
                      column(2,
                             dropdownButton(
                               selectInput('dwn_abund_corr_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                               downloadButton('dwn_abund_corr', ""),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
                    plotOutput("abun_cor_plot", click = "abun_cor_click")
                  ),
                  tabBox(
                    tabPanel("Abundance correlations",
                             fluidRow(
                               column(1, actionBttn(inputId = "corr_analysis", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(1,
                                      dropdownButton(
                                        tags$h4("Advanced settings"),
                                        selectInput("method", "Select correlation method:", c("spearman" = "Spearman", "pearson" = "Pearson")),
                                        numericInput("corr_pValue", "p-Value filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        numericInput("corr_threshold", "Correlation coefficient threshold:", min = -100, max = 100, value = 0.1, step = 0.1),
                                        numericInput("corr_bg_anova_alpha", "Background anova filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        numericInput("corr_bg_ctCor_alpha", "Background correlation filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        numericInput("corr_gene_ctCor_alpha", "Marker correlation filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation settings")
                                      )
                               ),
                               column(1,
                                      dropdownButton(
                                        tags$h4("Neo4j export"),
                                        h5("Check that neo4j database is created and activated"),
                                        textInput('user_neo4j', label = h4("User name"), value = "neo4j"),
                                        textInput('password_neo4j', label = h4("Password"), value = "password"),
                                        actionButton('neo4j_export', "Export"),
                                        icon = icon("share-alt"), status = "primary", tooltip = tooltipOptions(title = "export database to neo4j")
                                      )
                               ),
                               column(1,
                                      dropdownButton(
                                        selectInput('dwn_gene_cor_acorr_ext', label = NULL,
                                                    choices = list("chosen correlations" = 'focus_corr_data',
                                                                   "correlations within clusters" = 'signals_in_cluster',
                                                                   "correlations between clusters" = 'signals_between_clusters',
                                                                   "all correlations" = 'signals')),
                                        downloadButton('dwn_gene_cor_acorr', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             DT::dataTableOutput("gene_cor_table")
                    ),
                    tabPanel("Marker correlations",
                             fluidRow(
                               column(1, actionBttn(inputId = "mk_corr_analysis", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(1,
                                      dropdownButton(
                                        tags$h4("Advanced settings"),
                                        selectInput("mk_corr_method", "Select correlation method:", c("spearman" = "Spearman", "pearson" = "Pearson")),
                                        numericInput("mk_corr_pValue", "p-Value filter (alpha):", min = 0, max = 1, value = 0.1, step = 0.001),
                                        numericInput("mk_corr_threshold", "Correlation coefficient threshold:", min = 0, max = 100, value = 0, step = 0.1),
                                        numericInput("mk_corr_bg_anova_alpha", "Background anova filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        numericInput("mk_corr_bg_ctCor_alpha", "Background correlation filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        numericInput("mk_corr_gene_ctCor_alpha", "Marker correlation filter (alpha):", min = 0, max = 1, value = 0.01, step = 0.001),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation settings")
                                      )
                               ),
                               column(1,
                                      dropdownButton(
                                        selectInput('dwn_marker_table_mk_corr_ext', label = NULL,
                                                    choices = list("chosen correlations" = 'focus_corr_data',
                                                                   "correlations within clusters" = 'signals_in_cluster',
                                                                   "correlations between clusters" = 'signals_between_clusters',
                                                                   "all correlations" = 'signals')),
                                        downloadButton('dwn_marker_table_mk_corr', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             DT::dataTableOutput('marker_mk_corr_table')
                    ),
                    width = 12
                  )
                )
        ),

        # Fourth tab content
        tabItem(tabName = 'data_crosspanel',
                fluidRow(
                  tabBox(
                    tabPanel("Upload panel",
                             h4("FCS files"),
                             h5("Choose FSC files to upload them as one panel"),
                             shinyFilesButton('choose_fcs_cross_p', label= h5("Select FCS files"),
                                              title='Please select clustered FCS files to upload as panel', multiple=TRUE),
                             hr(),
                             h4("Name the panel"),
                             textInput("panel_name_cross_p", label = h5("You can name the new panel (try to use a short name or even one or few letters)")),
                             hr(),
                             textInput("extr_clust_pattern_cross_p",
                                       label = h5("full or part column name with clusters info (for cytofBrowser and cytofkit : <cluster>)"), value = "cluster"),
                             hr(),
                             actionButton('butt_upload_cross_p', label = "Upload")

                    ),
                    tabPanel("Remove panel",
                             uiOutput('remove_panel_cross_p_ui')
                    )
                  ),
                  tabBox(
                    tabPanel("Samples content",
                             plotOutput('content_cross_p')
                    ),
                    tabPanel("Markers content",
                             plotOutput('mk_content_cross_p')
                    )
                  ),
                  tabBox(
                    tabPanel("Abundance cross-panel correlation",
                             fluidRow(
                               column(2, actionBttn(inputId = "abund_corr_cross_p", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(2,
                                       dropdownButton(
                                         tags$h4("Advanced options"),
                                         selectInput('abund_corr_method_cross_p', "Select correlation method:", c("spearman" = 'spearman', "pearson" = 'pearson')),
                                         selectInput('abund_corr_p_mode_cross_p', "p-value", c("p-value" = "pval", "adjusted p-value" = "padj"), selected = 'padj'),
                                         uiOutput('abund_corr_padj_method_cross_p_ui'),
                                         icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation test options")
                                       )
                               ),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_abund_corr_cross_p_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                                        downloadButton('dwn_abund_corr_cross_p', ""),
                                        hr(),
                                        selectInput('dwn_table_abund_corr_cross_p_ext', label = NULL,
                                                    choices = list("abundance data" = 'abund_data', "correlation data" = 'corr_data'),
                                                    selected = 'corr_data'),
                                        downloadButton('dwn_table_abund_corr_cross_p', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             plotOutput('plot_abund_corr_cross_p')
                    ),
                    tabPanel("Expressing cell fraction Correlation",
                             fluidRow(
                               column(2, actionBttn(inputId = "exp_cell_f_corr_cross_p", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(4, uiOutput('mk_exp_cell_f_cross_p_ui')),

                               column(2,
                                       dropdownButton(
                                         tags$h4("Settings"),
                                         selectInput("exp_cell_f_corr_method", label =  "Choose the method to get expression cells fraction",
                                                     choices =c("clustering" = 'clustering', "threshold" = 'threshold'), selected = 'clustering'),
                                         uiOutput('threshold_exp_cell_f_cross_p_ui'),
                                         icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "Cell fraction allocation methods")
                                       )
                               ),
                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        selectInput('exp_cell_f_corr_method_cross_p', "Select correlation method:", c("spearman" = 'spearman', "pearson" = 'pearson')),
                                        selectInput('exp_cell_f_corr_p_mode_cross_p', "p-value", c("p-value" = "pval", "adjusted p-value" = "padj"), selected = 'padj'),
                                        uiOutput('exp_cell_f_corr_padj_method_cross_p_ui'),
                                        numericInput('exp_cell_f_corr_min_cell_cross_p', "Min expressing cell number ", min = 0, value = 5, step = 1),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation test options")
                                      )

                               ),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_exp_cell_f_corr_cross_p_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                                        downloadButton('dwn_exp_cell_f_corr_cross_p', ""),
                                        hr(),
                                        selectInput('dwn_table_exp_cell_f_corr_cross_p_ext', label = NULL,
                                                    choices = list("expressing cell fractions" = 'exp_cell_f_data', "correlation data" = 'corr_data'),
                                                    selected = 'corr_data'),
                                        downloadButton('dwn_table_exp_cell_f_corr_cross_p', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             fluidRow(
                               column(7, plotOutput('plot_exp_cell_f_corr_cross_p')),
                               column(5, plotOutput('mk_density_plot_cross_p'))
                             )

                    ),
                    width = 12
                  )

                )
        )
      )
    )
  )

  shinyApp(ui = cytofBrowser_ui, server = cytofBrowser_server)
}
