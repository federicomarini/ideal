# ideal.R

#' ideal: Interactive Differential Expression Analysis
#'
#' ideal makes differential expression analysis interactive, easy and reproducible.
#' This function launches the main application included in the package.
#'
#' @param dds_obj A \code{\link{DESeqDataSet}} object. If not provided, then a
#' \code{countmatrix} and a \code{expdesign} need to be provided. If none of
#' the above is provided, it is possible to upload the data during the
#' execution of the Shiny App
#' @param res_obj  A \code{\link{DESeqResults}} object. If not provided, it can
#' be computed during the execution of the application
#' @param annotation_obj A \code{data.frame} object, with row.names as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. If not provided, it can be constructed during
#' the execution via the org.eg.XX.db packages - these need to be installed
#' @param countmatrix A count matrix, with genes as rows and samples as columns.
#' If not provided, it is possible to upload the data during the execution of
#' the Shiny App
#' @param expdesign A \code{data.frame} containing the info on the covariates
#' of each sample. If not provided, it is possible to upload the data during the
#' execution of the Shiny App
#' @param gene_signatures A list of vectors, one for each pathway/signature. This 
#' is for example the output of the \code{\link{read_gmt}} function. The provided
#' object can also be replaced during runtime in the dedicated upload widget.
#'
#' @return A Shiny App is launched for interactive data exploration and
#' differential expression analysis
#'
#' @export
#'
#' @examples
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
#' cm <- counts(dds)
#' cd <- colData(dds)
#'
#' # with the well known airway package...
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' \dontrun{
#'
#' ideal()
#' ideal(dds)
#' ideal(dds_airway)
#'
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#' ideal(dds_airway, res_airway)
#' }
#'
ideal<- function(dds_obj = NULL,
                 res_obj = NULL,
                 annotation_obj = NULL,
                 countmatrix = NULL,
                 expdesign = NULL,
                 gene_signatures = NULL){

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }


  # create environment for storing inputs and values
  ## i need the assignment like this to export it up one level - i.e. "globally"
  ideal_env <<- new.env(parent = emptyenv())

  ## upload max 300mb files - can be changed if necessary
  options(shiny.maxRequestSize=300*1024^2)
  options(shiny.launch.browser = TRUE)

  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  #   # components defined in separated .R files
  #   shinyApp(ui = ideal_ui, server = ideal_server)

  # ui definition -----------------------------------------------------------
  ideal_ui <- shinydashboard::dashboardPage(
    title = "ideal - Interactive Differential Expression AnaLysis",
    # header definition -----------------------------------------------------------
    shinydashboard::dashboardHeader(
      title = tags$span(
        img(src = "ideal/ideal_logo_v2.png", height = "50px"),
        paste0("ideal - Interactive Differential Expression AnaLysis ",
                     packageVersion("ideal"))),
      titleWidth = 600,

      # TODO:
      # http://stackoverflow.com/questions/31440564/adding-a-company-logo-to-shinydashboard-header
      # replace text with image
      # ideal_header$children[[2]]$children <- tags$a(href='https://github.com/federicomarini/ideal',
      # tags$img(src='ideal_logo_v2.png',height='50',width='200'))
      # title = tags$a(href='https://github.com/federicomarini/ideal',
      #                tags$img(src='ideal_logo_v2.png',height='50',width='200')),

      # task menu for saving state to environment or binary data
      shinydashboard::dropdownMenu(
        type = "tasks",icon = icon("cog"),
        badgeStatus = NULL, 
        headerText = "ideal Tasks menu",
        notificationItem(
          text = actionButton("task_exit_and_save","Exit ideal & save",
                              class = "btn_no_border",
                              onclick = "setTimeout(function(){window.close();}, 100); "),
          icon = icon("sign-out"),status = "primary"),
        menuItem(
          text = downloadButton("task_state_save","Save State as .RData"))
      )
    ), # end of dashboardHeader

    # sidebar definition -----------------------------------------------------------
    dashboardSidebar(
      width = 280,
      menuItem("App settings",
               icon = icon("cogs"),
               startExpanded = TRUE,
               uiOutput("color_by"),
               shinyBS::bsTooltip(
                 "color_by", 
                 paste0("Select the group(s) of samples to stratify the analysis, and ideally match the contrast of interest. Can also assume multiple values, in this case the interaction of the factors is used."),
                 "right", options = list(container = "body")),
               uiOutput("available_genes"),
               shinyBS::bsTooltip(
                 "available_genes", 
                 paste0("Select one or more features (genes) from the list to inspect. Autocompletion is provided, so you can easily find your genes of interest by started typing their names. Defaults to the row names if no annotation object is provided."),
                 "right", options = list(container = "body")),
               numericInput("FDR","False Discovery Rate",value = 0.05, min = 0, max = 1, step = 0.01),
               shinyBS::bsTooltip(
                 "FDR", 
                 paste0("Select the alpha level at which you would like to control the FDR (False Discovery Rate) for the set of multiple tests in your dataset. The sensible choice of 0.05 is provided as default, 0.1 is more liberal, while 0.01 is more stringent - keep in mind this does not tell anything on the effect size for the expression change."),
                 "right", options = list(container = "body"))
               
      ),
      menuItem("Plot export settings", 
               icon = icon("paint-brush"),
               startExpanded = TRUE,
               numericInput("export_width",label = "Width of exported figures (cm)",value = 16,min = 2),
               shinyBS::bsTooltip(
                 "export_width", paste0("Width of the figures to export, expressed in cm"),
                 "right", options = list(container = "body")),
               numericInput("export_height",label = "Height of exported figures (cm)",value = 10,min = 2),
               shinyBS::bsTooltip(
                 "export_height", paste0("Height of the figures to export, expressed in cm"),
                 "right", options = list(container = "body"))
               ),
      menuItem("Quick viewer", 
               icon = icon("flash"), 
               startExpanded = TRUE,
               id = "qvmenu",
               fluidRow(
                 fluidRow(column(6,p("Count matrix")), column(6,uiOutput("ok_cm"))),
                 fluidRow(column(6,p("Experimental design")), column(6,uiOutput("ok_ed"))),
                 fluidRow(column(6,p("DESeqDataset")), column(6,uiOutput("ok_dds"))),
                 fluidRow(column(6,p("Annotation")), column(6,uiOutput("ok_anno"))),
                 fluidRow(column(6,p("Results")), column(6,uiOutput("ok_resu")))
               )),
      menuItem("First steps help", 
               icon = icon("question-circle"),
               startExpanded = TRUE,
               actionButton("btn", "Click me for a quick tour", icon("info"),
                            style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4")
      )
    ), # end of dashboardSidebar

    # body definition -----------------------------------------------------------
    dashboardBody(
      introjsUI(),

      ## Define output size and style of error messages, and also the style of the icons e.g. check
      ## plus, define the myscrollbox div to prevent y overflow when page fills up
      tags$head(
        tags$style(HTML("
                        .shiny-output-error-validation {
                        font-size: 15px;
                        color: forestgreen;
                        text-align: center;
                        }
                        .icon-done {
                        color: green;
                        }
                        #myScrollBox{
                        overflow-y: scroll;

                        .dataTables_wrapper{
                        overflow-x: scroll;
                        }
                        }
                        #myAnchorBox{}
                        "))
        ),

      # value boxes to always have an overview on the available data
      fluidRow(
        valueBoxOutput("box_ddsobj"),
        valueBoxOutput("box_annobj"),
        valueBoxOutput("box_resobj")
      ),

      ## main structure of the body for the dashboard
      div(
        id = "myScrollBox", # trick to have the y direction scrollable
        tabBox(
          width=12,

          # ui panel welcome -----------------------------------------------------------
          tabPanel(
            title = "Welcome!",  icon = icon("home"), value="tab-welcome",
            
            fluidRow(
              column(
                width = 8,
                includeMarkdown(system.file("extdata", "welcome.md",package = "ideal")),
                br(),br(),
                p("If you see a grey box like this one open below..."),

                shinyBS::bsCollapse(
                  id = "help_welcome",open = "Help", 
                  shinyBS::bsCollapsePanel(
                    "Help", 
                    includeMarkdown(system.file("extdata", "help_welcome.md",package = "ideal"))
                  )
                ),

                actionButton("introexample", "If you see a button like this...", icon("info"),
                             style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
                p("... you can click on that to start a tour based on introJS"),
                br(),br(),
                
                uiOutput("ui_instructions")
              )
            )
          ), # end of Welcome panel

          # ui panel data setup -----------------------------------------------------------
          tabPanel(
            "Data Setup",icon = icon("upload"), # value="tab-ds",
            value = "tab-datasetup",
            headerPanel("Setup your data for the analysis"),
            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(
                  id = "help_datasetup",open = NULL, 
                  shinyBS::bsCollapsePanel(
                    "Help",
                    includeMarkdown(system.file("extdata", "help_datasetup.md",package = "ideal"))
                  )
                )
              )
            ),

            actionButton("tour_datasetup", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

            box(
              width = 12, 
              title = "Step 1", status = "danger", solidHeader = TRUE,
              h2("Upload your count matrix and the info on the experimental design"),

              fluidRow(
                column(
                  width = 4,
                  uiOutput("upload_count_matrix"),
                  uiOutput("upload_metadata"),
                  br(),
                  "... or you can also ",
                  actionButton("btn_loaddemo", "Load the demo airway data", icon = icon("play-circle"),
                               class = "btn btn-info"),br(), p()
                ),
                column(
                  width = 4,
                  br(),
                  actionButton("help_format",label = "",icon = icon("question-circle"),
                               style="color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"),
                  shinyBS::bsTooltip(
                    "help_format", 
                    "How to provide your input data to ideal",
                    "bottom", options = list(container = "body")
                  )
                )
              ),
              
              fluidRow(
                column(
                  width = 6,
                  box(width = NULL, title = "Count matrix preview",status = "primary",
                      solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
                      fluidRow(
                        column(
                          width = 12,
                          offset = 0.5,
                          DT::dataTableOutput("dt_cm"))
                      )
                  )
                ),
                column(
                  width = 6,
                  box(width = NULL, title = "Experimental design preview",status = "primary",
                      solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
                      fluidRow(
                        column(
                          width = 12,
                          offset = 0.5,
                          DT::dataTableOutput("dt_ed"))
                      )
                  )
                )
              )
            ),
            uiOutput("ui_step2"),
            fluidRow(
              column(
                width = 6,
                uiOutput("ui_stepanno")
                ## this ideally populates also the list of genes of interest to choose among
              ),
              column(
                width = 6,
                uiOutput("ui_stepoutlier")
              )
            ),
            uiOutput("ui_step3")
          ), # end of Data Setup panel
          # ui panel counts overview -----------------------------------------------------------
          tabPanel(
            "Counts Overview",
            icon = icon("eye"),
            conditionalPanel(
              condition="!output.checkdds",
              headerPanel("Get an overview on your data"),
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_countsoverview",open = NULL, 
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_overview.md",package = "ideal")))
                  )
                  )
                ),

              actionButton("tour_countsoverview", "Click me for a quick tour of the section", 
                           icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
              br(),
              selectInput("countstable_unit", label = "Data scale in the table",
                          choices = list("Counts (raw)" = "raw_counts",
                                         "Counts (normalized)" = "normalized_counts",
                                         "Regularized logarithm transformed" = "rlog_counts",
                                         "Log10 (pseudocount of 1 added)" = "log10_counts",
                                         "TPM (Transcripts Per Million)" = "tpm_counts")),

              DT::dataTableOutput("showcountmat"),
              downloadButton("downloadData","Download", class = "btn btn-success"),
              hr(),
              fluidRow(
                column(
                  width = 8,
                  h3("Basic summary for the counts"),
                  p("Number of uniquely aligned reads assigned to each sample"),
                  # verbatimTextOutput("reads_summary"),
                  wellPanel(
                    fluidRow(
                      column(
                        width = 6,
                        numericInput("threshold_rowsums","Threshold on the row sums of the counts",value = 0, min = 0)),
                      column(
                        width = 6,
                        numericInput("threshold_rowmeans","Threshold on the row means of the normalized counts",value = 0, min = 0))
                    )),
                  p("According to the selected filtering criteria, this is an overview on the provided count data"),
                  verbatimTextOutput("detected_genes"),

                  selectInput("filter_crit",label = "Choose the filtering criterium",
                              choices = c("row means", "row sums"), selected = "row means"),

                  actionButton("featfilt_dds", "Filter the DDS object",class = "btn btn-primary")
                )
              ),

              h3("Sample to sample scatter plots"),
              selectInput("corr_method","Correlation method",choices = list("pearson","spearman", "kendall")),
              checkboxInput(inputId = "corr_uselogs",
                            label = "Use log2 values for plot axes and values",
                            value = TRUE),
              checkboxInput(inputId = "corr_usesubset",
                            label = "Use a subset of max 1000 genes (quicker to plot)",
                            value = TRUE),
              p("Compute sample to sample correlations on the normalized counts - warning, it can take a while to plot all points (depending mostly on the number of samples you provided)."),
              actionButton("compute_pairwisecorr", "Run", class = "btn btn-primary"),
              uiOutput("pairwise_plotUI"),
              uiOutput("heatcorr_plotUI")
            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )
          ), # end of Counts Overview panel
          # ui panel extract results -----------------------------------------------------------
          tabPanel(
            "Extract Results", icon = icon("table"),
            # see: http://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value?noredirect=1&lq=1
            conditionalPanel(
              condition="!output.checkdds",
              headerPanel("Extract and inspect the DE results"),
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_extractresults",open = NULL,
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_results.md",package = "ideal")))
                  )
                )
              ),

              actionButton("tour_results", "Click me for a quick tour of the section", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
              br(),
              fluidRow(
                column(
                  width = 6,
                  uiOutput("choose_fac")
                )
              ),
              fluidRow(
                column(
                  width = 4,
                  # factor as covariate
                  wellPanel(
                    width = 4, id = "factor_opts",
                    uiOutput("fac1"),
                    uiOutput("fac2"),
                    # continuous covariate
                    uiOutput("facnum")
                  )

                ),
                column(
                  width = 4,
                  # factor with > 2 levels
                  wellPanel(
                    width = 4,
                    uiOutput("lrtavailable"),
                    uiOutput("lrtfull"),
                    uiOutput("lrtreduced")
                  ),

                  uiOutput("runlrt")
                )
              ),

              ## general options for result function
              # alpha is set via FDR on the left side
              fluidRow(
                column(
                  width = 4,
                  wellPanel(id = "resu_opts",
                    selectInput("resu_indfil",label = "Apply independent filtering automatically",
                                choices = c(TRUE,FALSE), selected = TRUE),
                    selectInput("resu_lfcshrink",label = "Shrink the log fold change for the contrast of interest",
                                choices = c(TRUE,FALSE), selected = TRUE),
                    selectInput("resu_ihw", "Use Independent Hypothesis Weighting (IHW) as a filtering function",
                                choices = c(TRUE, FALSE), selected = FALSE)
                  )
                )
              ),
              #, evtl also the *filter* parameter of the function, i.e. baseMean if not specified
              fluidRow(
                column(
                  width = 6,
                  uiOutput("runresults"),
                  uiOutput("store_result"),
                  verbatimTextOutput("diyres_summary")
                )
              ),

              DT::dataTableOutput("table_res"),
              downloadButton("downloadTblResu","Download", class = "btn btn-success"),
              fluidRow(
                h3("Diagnostic plots"),
                column(
                  width = 6,
                  plotOutput("pvals_hist"),
                  div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                      downloadButton("download_plot_pvals_hist", "Download Plot"),
                      textInput("filename_plot_pvals_hist",label = "Save as...",value = "plot_pvals_hist.pdf"))
                ),
                column(
                  width = 6,
                  plotOutput("pvals_hist_strat"),
                  div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                      downloadButton("download_plot_pvals_hist_strat", "Download Plot"),
                      textInput("filename_plot_pvals_hist_strat",label = "Save as...",value = "plot_pvals_hist_strat.pdf"))
                ),
                column(
                  width = 6,
                  plotOutput("pvals_ss"),
                  div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                      downloadButton("download_plot_pvals_ss", "Download Plot"),
                      textInput("filename_plot_pvals_ss",label = "Save as...",value = "plot_pvals_ss.pdf"))
                ),
                column(
                  width = 6,
                  plotOutput("logfc_hist"),
                  div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                      downloadButton("download_plot_logfc_hist", "Download Plot"),
                      textInput("filename_plot_logfc_hist",label = "Save as...",value = "plot_logfc_hist.pdf"))
                )
              )
            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )

          ), # end of Extract Results panel
          # ui panel summary plots -----------------------------------------------------------
          tabPanel(
            "Summary Plots", icon = icon("photo"),
            conditionalPanel(
              condition="!output.checkresu",
              headerPanel("Interactive graphical exploration of the results"),
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_summaryplots",open = NULL, 
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_plots.md",package = "ideal")))
                  )
                )
              ),

              actionButton("tour_plots", "Click me for a quick tour of the section", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
              
              br(),
              fluidRow(
                column(6,
                       h4("MA plot - Interactive!"),
                       plotOutput('plotma', brush = 'ma_brush'),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_ma", "Download Plot"),
                           textInput("filename_plot_ma",label = "Save as...",value = "plot_ma.pdf"))),
                column(6,
                       h4("Zoomed section"),
                       plotOutput("mazoom",click= 'mazoom_click'),
                       numericInput('size_genelabels', label = 'Labels size: ', value = 4,min = 1,max = 8),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_mazoom", "Download Plot"),
                           textInput("filename_plot_mazoom",label = "Save as...",value = "plot_mazoom.pdf")))
              ),
              fluidRow(
                column(6,
                       h4("Selected gene"),
                       checkboxInput("ylimZero_genes","Set y axis limit to 0",value=TRUE),
                       plotOutput("genefinder_plot"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_genefinder", "Download Plot"),
                           textInput("filename_plot_genefinder",label = "Save as...",value = "plot_genefinder.pdf"))
                ),
                column(6,
                       h4("Gene infobox"),
                       htmlOutput("rentrez_infobox"))
              ),

              fluidRow(
                column(6,
                       h4("volcano plot"),
                       plotOutput("volcanoplot"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_volcanoplot", "Download Plot"),
                           textInput("filename_plot_volcanoplot",label = "Save as...",value = "plot_volcanoplot.pdf"))
              )),

              fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
              fluidRow(
                column(4,
                       checkboxInput("rowscale",label = "Scale by rows",value = TRUE)),
                column(4,
                       checkboxInput("pseudocounts","use log2(1+counts)",value = TRUE))
              ),
              fluidRow(
                column(6,
                       plotOutput("heatbrush"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_heatbrush", "Download Plot"),
                           textInput("filename_plot_heatbrush",label = "Save as...",value = "plot_heatbrush.pdf"))
                ),
                column(6,
                       d3heatmapOutput("heatbrushD3"))
              ),

              box(
                title = "Brushed table", status = "primary", solidHeader = TRUE,
                id = "box_brushedtbl",
                collapsible = TRUE, collapsed = TRUE, width = 12,
                fluidRow(DT::dataTableOutput("ma_brush_out"),
                         downloadButton("downloadTblMabrush","Download", class = "btn btn-success")))
            ),
            conditionalPanel(
              condition="output.checkresu",
              h2("You did not create the result object yet. Please go the dedicated tab and generate it")
            )
          ), # end of Summary Plots panel
          # ui panel gene finder -----------------------------------------------------------
          tabPanel(
            "Gene Finder", icon = icon("crosshairs"),
            conditionalPanel(
              condition="!output.checkdds",
              headerPanel("Find your gene(s) of interest"),
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_genefinder",open = NULL,
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_genefinder.md",package = "ideal")))
                  )
                )
              ),
              actionButton("tour_genefinder", "Click me for a quick tour of the section", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
              br(),
              fluidRow(
                column(6,checkboxInput("ylimZero_genefinder","Set y axis limit to 0",value=TRUE))),
              fluidRow(
                column(6,
                       plotOutput("bp1"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plotbp1", "Download Plot"),
                           textInput("filename_plotbp1",label = "Save as...",value = "plotbp1.pdf"))
                ),
                column(6,
                       plotOutput("bp2"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plotbp2", "Download Plot"),
                           textInput("filename_plotbp2",label = "Save as...",value = "plotbp2.pdf")))
              ),
              fluidRow(
                column(6,
                       plotOutput("bp3"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plotbp3", "Download Plot"),
                           textInput("filename_plotbp3",label = "Save as...",value = "plotbp3.pdf"))
                ),
                column(6,
                       plotOutput("bp4"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plotbp4", "Download Plot"),
                           textInput("filename_plotbp4",label = "Save as...",value = "plotbp4.pdf")))
              ),

              plotOutput("ma_highlight"),
              div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                  downloadButton("download_plot_mahighlight", "Download Plot"),
                  textInput("filename_plot_mahighlight",label = "Save as...",value = "plot_mahighlight.pdf")),
              DT::dataTableOutput("table_combi"),
              downloadButton("downloadTblCombi","Download", class = "btn btn-success"),

              fileInput(inputId = "gl_ma",
                        label = "Upload a gene list file",
                        accept = c("text/csv", "text/comma-separated-values",
                                   "text/tab-separated-values", "text/plain",
                                   ".csv", ".tsv"), multiple = FALSE),
              plotOutput("ma_hl_list"),
              div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                  downloadButton("download_plot_mahllist", "Download Plot"),
                  textInput("filename_plot_mahllist",label = "Save as...",value = "plot_mahllist.pdf")),
              DT::dataTableOutput("table_combi_list"),
              downloadButton("downloadTblCombiList","Download", class = "btn btn-success")
            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )
          ), # end of Gene Finder panel
          # ui panel functional analysis ----------------------------------------------------------
          tabPanel(
            "Functional Analysis", icon = icon("list-alt"),
            conditionalPanel(
              condition="!output.checkresu",
              headerPanel("Find functions enriched in gene sets"),
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_functionalanalysis",open = NULL,
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_funcanalysis.md",package = "ideal")))
                  )
                )
              ),
              actionButton("tour_funcanalysis", "Click me for a quick tour of the section", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

              selectInput("go_cats",label = "Select the GO category(ies) of interest",
                          choices = list("GO Biological Process" = "BP", "GO Molecular Function" = "MF", "GO Cellular Component" = "CC"),
                          selected = "BP",multiple = TRUE
              ),

              div(
                id="myAnchorBox",
                tabBox(
                  width = NULL,
                  id="gse_tabbox",
                  tabPanel("UPregu", icon = icon("arrow-circle-up"),
                           fluidRow(column(width = 6,actionButton("button_enrUP", "Perform gene set enrichment analysis on the upregulated genes",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrUP_goseq", "Perform gene set enrichment analysis on the upregulated genes - goseq",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrUP_topgo", "Perform gene set enrichment analysis on the upregulated genes - topGO",class = "btn btn-primary"))),
                           DT::dataTableOutput("DT_gse_up"),
                           DT::dataTableOutput("DT_gse_up_goseq"),
                           fluidRow(
                             column(width = 9, DT::dataTableOutput("DT_gse_up_topgo"),
                                    downloadButton("downloadGOTbl_up","Download", class = "btn btn-success")),
                             column(width = 3, plotOutput("goterm_heatmap_up_topgo"))
                           )
                  ),
                  tabPanel("DOWNregu", icon = icon("arrow-circle-down"),
                           fluidRow(column(width = 6,actionButton("button_enrDOWN", "Perform gene set enrichment analysis on the downregulated genes",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrDOWN_goseq", "Perform gene set enrichment analysis on the downregulated genes - goseq",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrDOWN_topgo", "Perform gene set enrichment analysis on the downregulated genes - topGO",class = "btn btn-primary"))),
                           DT::dataTableOutput("DT_gse_down"),
                           DT::dataTableOutput("DT_gse_down_goseq"),
                           fluidRow(
                             column(width = 9, DT::dataTableOutput("DT_gse_down_topgo"),
                                    downloadButton("downloadGOTbl_down","Download", class = "btn btn-success")),
                             column(width = 3, plotOutput("goterm_heatmap_down_topgo"))
                           )
                  ),
                  tabPanel("UPDOWN", icon = icon("arrows-v"),
                           fluidRow(column(width = 6,actionButton("button_enrUPDOWN", "Perform gene set enrichment analysis on the up- and downregulated genes",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrUPDOWN_goseq", "Perform gene set enrichment analysis on the up- and downregulated genes - goseq",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrUPDOWN_topgo", "Perform gene set enrichment analysis on the up- and downregulated genes - topGO",class = "btn btn-primary"))),
                           DT::dataTableOutput("DT_gse_updown"),
                           DT::dataTableOutput("DT_gse_updown_goseq"),
                           fluidRow(
                             column(width = 9, DT::dataTableOutput("DT_gse_updown_topgo"),
                                    downloadButton("downloadGOTbl_updown","Download", class = "btn btn-success")),
                             column(width = 3, plotOutput("goterm_heatmap_updown_topgo"))
                           )
                  ),
                  tabPanel("List1", icon = icon("list"),
                           fileInput(inputId = "gl1",
                                     label = "Upload a gene list file",
                                     accept = c("text/csv", "text/comma-separated-values",
                                                "text/tab-separated-values", "text/plain",
                                                ".csv", ".tsv"), multiple = FALSE),
                           fluidRow(column(width = 6,actionButton("button_enrLIST1", "Perform gene set enrichment analysis on the genes in list1",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrLIST1_goseq", "Perform gene set enrichment analysis on the list1 genes - goseq",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrLIST1_topgo", "Perform gene set enrichment analysis on the list1 genes - topGO",class = "btn btn-primary"))),
                           DT::dataTableOutput("DT_gse_list1"),
                           DT::dataTableOutput("DT_gse_list1_goseq"),
                           fluidRow(
                             column(width = 9, DT::dataTableOutput("DT_gse_list1_topgo"),
                                    downloadButton("downloadGOTbl_l1","Download", class = "btn btn-success")),
                             column(width = 3, plotOutput("goterm_heatmap_l1_topgo"))
                           )
                           
                  ),
                  tabPanel("List2", icon = icon("list-alt"),
                           fileInput(inputId = "gl2",
                                     label = "Upload a gene list file",
                                     accept = c("text/csv", "text/comma-separated-values",
                                                "text/tab-separated-values", "text/plain",
                                                ".csv", ".tsv"), multiple = FALSE),
                           fluidRow(column(width = 6,actionButton("button_enrLIST2", "Perform gene set enrichment analysis on the genes in list2",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrLIST2_goseq", "Perform gene set enrichment analysis on the list2 genes - goseq",class = "btn btn-primary"))),
                           fluidRow(column(width = 6,actionButton("button_enrLIST2_topgo", "Perform gene set enrichment analysis on the list2 genes - topGO",class = "btn btn-primary"))),
                           DT::dataTableOutput("DT_gse_list2"),
                           DT::dataTableOutput("DT_gse_list2_goseq"),
                           fluidRow(
                             column(width = 9, DT::dataTableOutput("DT_gse_list2_topgo"),
                                    downloadButton("downloadGOTbl_l2","Download", class = "btn btn-success")),
                             column(width = 3, plotOutput("goterm_heatmap_l2_topgo"))
                           )
                  )
                )
              ),

              ## will put collapsible list elements? or multi tab panel? or something to select on the left, and operate output-wise on the right e.g. venn diagrams or table for gene set enrichment
              # h3("custom list 3 - handpicked") # use the select input from the left column?
              # ,verbatimTextOutput("debuggls"),

              # verbatimTextOutput("printUPgenes"),
              # verbatimTextOutput("debuglists"),

              h2("Intersection of gene sets"),

              fluidRow(
                column(width = 4,
                       checkboxInput("toggle_updown","Use up and down regulated genes", TRUE),
                       checkboxInput("toggle_up","Use up regulated genes", FALSE),
                       checkboxInput("toggle_down","Use down regulated genes", FALSE)
                ),
                column(width = 4,
                       checkboxInput("toggle_list1","Use list1 genes", TRUE),
                       checkboxInput("toggle_list2","Use list2 genes", FALSE),
                       checkboxInput("toggle_list3","Use list3 genes", FALSE)
                )
              ),

              fluidRow(
                column(width = 6,plotOutput("vennlists"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_vennlists", "Download Plot"),
                           textInput("filename_plot_vennlists",label = "Save as...",value = "plot_vennlists.pdf")),
                       offset = 3)),
              fluidRow(
                column(width = 6,plotOutput("upsetLists"),
                       div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                           downloadButton("download_plot_upsetlists", "Download Plot"),
                           textInput("filename_plot_upsetlists",label = "Save as...",value = "plot_upsetlists.pdf")),
                       offset = 3))

            ),
            conditionalPanel(
              condition="output.checkresu",
              h2("You did not create the result object yet. Please go the dedicated tab and generate it")
            )
          ), # end of Functional Analysis panel
          # ui panel signatures explorer ---------------------------------------------------------
          tabPanel(
            "Signatures Explorer",
            icon = icon("map"),
            conditionalPanel(
              condition="!output.checkdds",
              fluidRow(
                column(
                  width = 8,
                  shinyBS::bsCollapse(
                    id = "help_signatureexplorer",open = NULL,
                    shinyBS::bsCollapsePanel(
                      "Help",
                      includeMarkdown(system.file("extdata", "help_signatureexplorer.md",package = "ideal")))
                  )
                )
              ),
              actionButton("tour_signatureexplorer", "Click me for a quick tour of the section", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),
              
              fluidRow(
                column(
                  width = 6,
                  h4("Setup options"),
                  wellPanel(
                    uiOutput("sig_ui_gmtin"),
                    uiOutput("sig_ui_nrsigs"),
                    actionButton("sig_button_computevst",
                                 label = "Compute the variance stabilized transformed data", 
                                 icon = icon("spinner"), class = "btn btn-success")
                  )
                ),
                column(
                  width = 6,
                  h4("Conversion options"),
                  wellPanel(
                    uiOutput("sig_ui_id_data"),
                    uiOutput("sig_ui_id_sigs"),
                    uiOutput("sig_ui_orgdbpkg"),
                    actionButton("sig_convert_setup",
                                 label = "Apply id conversion between data and signatures")
                  ),
                  verbatimTextOutput("sig_convcheck")
                  
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  wellPanel(
                    uiOutput("sig_ui_selectsig"),
                    uiOutput("sig_ui_annocoldata"),
                    checkboxInput("sig_useDEonly",
                                  label = "Use only DE genes in the signature",value = FALSE)
                  )
                  # ,
                  # verbatimTextOutput("sig_sigmembers")
                ),
                column(
                  width = 6,
                  wellPanel(
                    checkboxInput("sig_clusterrows",label = "Cluster rows", value = TRUE),
                    checkboxInput("sig_clustercols", label = "Cluster columns"),
                    checkboxInput("sig_centermean", label = "Center mean",value = TRUE),
                    checkboxInput("sig_scalerow", label = "Standardize by row")
                  )
                )
              ),
              fluidRow(
                column(
                  width = 8, offset = 2,
                  plotOutput("sig_heat")
                )
              )
            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )
          ), # end of Signatures Explorer panel
          
          # ui panel report editor -----------------------------------------------------------
          tabPanel(
            "Report Editor",
            icon = icon("pencil"),
            headerPanel("Create, view and export a report of your analysis"),
            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(
                  id = "help_reporteditor",open = NULL, 
                  shinyBS::bsCollapsePanel(
                    "Help",
                    includeMarkdown(system.file("extdata", "help_report.md",package = "ideal")))
                )
              )
            ),

            actionButton("tour_report", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

            fluidRow(
              column(
                width = 6,
                box(
                  title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
                  id = "md_opts",
                  radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = TRUE),
                  textInput("report_title", "Title: "),
                  textInput("report_author", "Author: "),
                  radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false")),
                  radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
                  selectInput("report_theme", "Theme", choices = list("Default" = "default", "Cerulean" = "cerulean",
                                                                      "Journal" = "journal", "Flatly" = "flatly",
                                                                      "Readable" = "readable", "Spacelab" = "spacelab",
                                                                      "United" = "united", "Cosmo" = "cosmo")),
                  radioButtons("report_echo", "Echo the commands in the output", choices = list("Yes" = "TRUE", "No" = "FALSE")))),
              column(
                width = 6,
                box(
                  title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
                  id = "editor_opts",
                  checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
                  conditionalPanel(
                    "input.enableAutocomplete",
                    wellPanel(
                      checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
                      checkboxInput("enableRCompletion", "R code completion", TRUE)
                    )
                  ),

                  selectInput("mode", "Mode: ", choices=shinyAce::getAceModes(), selected="markdown"),
                  selectInput("theme", "Theme: ", choices=shinyAce::getAceThemes(), selected="solarized_light"))
              )
            ),
            fluidRow(
              column(3,
                     actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
              ),
              column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success")),
              column(3, uiOutput("ui_iSEEexport"))
            ),

            tabBox(
              width = NULL,
              id="report_tabbox",
              tabPanel("Report preview",
                       icon = icon("file-text"),
                       htmlOutput("knitDoc")
              ),

              tabPanel("Edit report",
                       icon = icon("pencil-square-o"),
                       aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
                                 value=readLines(system.file("extdata", "irt.Rmd",package = "ideal")),
                                 height="800px"))
            )
          ), # end of Report Editor panel
          # ui panel about -----------------------------------------------------------
          tabPanel(
            "About", icon = icon("institution"),

            # headerPanel("Information on ideal/session"),

            fluidRow(
              column(
                width = 8,
                includeMarkdown(system.file("extdata", "about.md",package = "ideal")),

                verbatimTextOutput("sessioninfo")
              )
            )
          ) # end of About panel
        ) # end of box
      ) # end of myScrollBox
      ,footer()
    ), # end of dashboardBody
    skin="black"
  ) # end of dashboardPage


  # server definition -----------------------------------------------------------
  ideal_server <- shinyServer(function(input, output, session) {

    # server tours setup -----------------------------------------------------------
    
    # here will go the coded - i.e. not explicitly wrapped in introBox - steps
    intro_firsttour <- read.delim(system.file("extdata", "intro_firsttour.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_datasetup <- read.delim(system.file("extdata", "intro_datasetup.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_countsoverview <- read.delim(system.file("extdata", "intro_countsoverview.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_results <- read.delim(system.file("extdata", "intro_results.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_plots <- read.delim(system.file("extdata", "intro_plots.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_genefinder <- read.delim(system.file("extdata", "intro_genefinder.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_funcanalysis <- read.delim(system.file("extdata", "intro_funcanalysis.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_signatureexplorer <- read.delim(system.file("extdata", "intro_signatureexplorer.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)
    intro_report <- read.delim(system.file("extdata", "intro_report.txt",package = "ideal"), sep=";", stringsAsFactors = FALSE)

    observeEvent(input$btn, {
      introjs(session,
              options = list(steps = intro_firsttour)
      )
    })
    
    observeEvent(input$introexample, {
      intro_example <- data.frame(element=c("#introexample","#introexample"),
                                  intro=c("Tour elements can be anchored to elements of the UI that are intended to be highlighted. You can proceed to the next step by using the button, or also pushing the right arrow key.",
                                          "Well done. This is how a tour can look like. Click outside of this window to close the tour, or on the corresponding button."))
      introjs(session,
              options = list(steps = intro_example)
      )
    })
    
    observeEvent(input$tour_datasetup, {
      introjs(session,
              options = list(steps = intro_datasetup)
      )
    })

    observeEvent(input$tour_countsoverview, {
      introjs(session,
              options = list(steps = intro_countsoverview)
      )
    })

    observeEvent(input$tour_results, {
      introjs(session,
              options = list(steps = intro_results)
      )
    })

    observeEvent(input$tour_plots, {
      introjs(session,
              options = list(steps = intro_plots)
      )
    })

    observeEvent(input$tour_genefinder, {
      introjs(session,
              options = list(steps = intro_genefinder)
      )
    })

    observeEvent(input$tour_funcanalysis, {
      introjs(session,
              options = list(steps = intro_funcanalysis)
      )
    })
    
    observeEvent(input$tour_signatureexplorer, {
      introjs(session,
              options = list(steps = intro_signatureexplorer)
      )
    })
    
    observeEvent(input$tour_report, {
      introjs(session,
              options = list(steps = intro_report)
      )
    })



    ## Update directory
    userdir <- tempfile()
    dir.create(userdir, recursive = TRUE)
    # sapply(file.path(newuserdir, dir(newuserdir)[grep("code_", dir(newuserdir))]), file.remove)
    # file.copy(file.path(userdir, "code_All.R"), newuserdir)
    # userdir <- newuserdir
    # dir.create(file.path(userdir, "data"))

    # server setup reactivevalues -----------------------------------------------------------
    ## placeholder for the figures to export
    exportPlots <- reactiveValues()
    # expfig_fig1 <- NULL
    # )

    # will store all the reactive values relevant to the app
    values <- reactiveValues()

    values$countmatrix <- countmatrix
    values$expdesign <- expdesign

    values$dds_obj <- dds_obj
    values$res_obj <- res_obj
    values$annotation_obj <- annotation_obj
    values$gene_signatures <- gene_signatures


    # this part sets the "matching" objects if something is provided that is depending on these
    if(!is.null(dds_obj)){
      values$countmatrix <- counts(dds_obj, normalized = FALSE)
      values$expdesign <- as.data.frame(colData(dds_obj))
    }
    
    # server welcome home ---------------------------------------------------------
    output$ui_instructions <- renderUI({
      box(width = 12, 
          title = "Instructions", status = "info", solidHeader = TRUE, 
          collapsible = TRUE, collapsed = TRUE,
          includeMarkdown(system.file("extdata", "instructions.md",package = "ideal"))
      )
    })


    # server info boxes -----------------------------------------------------------
    output$box_ddsobj <- renderUI({
      if(!is.null(values$dds_obj))
        return(valueBox("dds object",
                        paste0(nrow(values$dds_obj), " genes - ",ncol(values$dds_obj)," samples"),
                        icon = icon("list"),
                        color = "green",width = NULL))
      else
        return(valueBox("dds object",
                        "yet to create",
                        icon = icon("list"),
                        color = "red",width = NULL))

    })

    output$box_annobj <- renderUI({
      if(!is.null(values$annotation_obj))
        return(valueBox("Annotation",
                        paste0(nrow(values$annotation_obj), " genes - ",ncol(values$annotation_obj)," ID types"),
                        icon = icon("book"),
                        color = "green",width = NULL))
      else
        return(valueBox("Annotation",
                        "yet to create",
                        icon = icon("book"),
                        color = "red",width = NULL))
    })

    output$box_resobj <- renderUI({
      if(!is.null(values$res_obj)){
        DEregu <- sum(values$res_obj$padj < input$FDR & values$res_obj$log2FoldChange != 0, na.rm = TRUE)
        return(valueBox("DE genes",
                        paste0(DEregu, " DE genes - out of ",nrow(values$res_obj),""),
                        icon = icon("list-alt"),
                        color = "green",width = NULL))
      } else
        return(valueBox("DE genes",
                        "yet to create",
                        icon = icon("list-alt"),
                        color = "red",width = NULL))
    })


    # if i want to focus a little more on the ihw object
    values$ihwres <- NULL

    # server uploading data -----------------------------------------------------------
    ## count matrix
    output$upload_count_matrix <- renderUI({
      if (!is.null(dds_obj) | !is.null(countmatrix)) {
        return(fluidRow(column(
          width = 12,
          tags$li("You already provided a count matrix or a DESeqDataSet object as input. You can check your input data in the collapsible box here below."), offset = 2)))
      } else {
        return(fileInput(inputId = "uploadcmfile",
                         label = "Upload a count matrix file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readCountmatrix <- reactive({
      if (is.null(input$uploadcmfile))
        return(NULL)
      guessed_sep <- sepguesser(input$uploadcmfile$datapath)
      cm <- utils::read.delim(input$uploadcmfile$datapath, header = TRUE,
                              as.is = TRUE, sep = guessed_sep, quote = "",
                              row.names = 1, # https://github.com/federicomarini/pcaExplorer/issues/1
                              ## TODO: tell the user to use tsv, or use heuristics
                              ## to check what is most frequently occurring separation character? -> see sepGuesser.R
                              check.names = FALSE)

      return(cm)
    })

    ## exp design
    output$upload_metadata <- renderUI({
      if (!is.null(dds_obj) | !is.null(expdesign)) {
        return(fluidRow(column(
          width = 12,
          tags$li("You already provided a matrix/data.frame with the experimental covariates or a DESeqDataSet object as input. You can check your input data in the collapsible box here below."), offset = 2)))

      } else {
        return(fileInput(inputId = "uploadmetadatafile",
                         label = "Upload a sample metadata matrix file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readMetadata <- reactive({
      if (is.null(input$uploadmetadatafile))
        return(NULL)
      guessed_sep <- sepguesser(input$uploadmetadatafile$datapath)
      expdesign <- utils::read.delim(input$uploadmetadatafile$datapath, header = TRUE,
                                     as.is = TRUE, sep = guessed_sep, quote = "",
                                     check.names = FALSE)

      return(expdesign)
    })

    # load the demo data
    observeEvent(input$btn_loaddemo,withProgress(
                 message = "Loading demo data",
                 detail = "Loading airway count and metadata information", value = 0,
                 {
                   aw <- requireNamespace("airway",quietly = TRUE)
                   incProgress(0.2,detail = "`airway` package loaded")
                   if(aw) {
                     data(airway,package="airway",envir = environment())

                   cm_airway <- assay(airway)
                   incProgress(0.7, detail = "Count matrix loaded")
                   ed_airway <- as.data.frame(colData(airway))

                   values$countmatrix <- cm_airway
                   values$expdesign <- ed_airway
                   incProgress(0.3, detail = "Experimental metadata loaded")
                   # just to be sure, erase the annotation and the rest
                   values$dds_obj <- NULL
                   values$annotation_obj <- NULL
                   values$res_obj <- NULL
                   showNotification("All components for generating the DESeqDataset object have been loaded, proceed to Step 2!",
                                    type = "message")
                   } else {
                     showNotification("The 'airway' package is currently not installed. Please do so by executing BiocManager::install('airway') before launching ideal()",type = "warning")
                   }
                 })
    )
    
    observeEvent(input$help_format, {
      showModal(modalDialog(
        title = "Format specifications for ideal",
        includeMarkdown(system.file("extdata", "datainput.md",package = "ideal")),
        h4("Example:"),
        tags$img(
          src = base64enc::dataURI(file = system.file("www", "help_dataformats.png",package = "pcaExplorer"), mime = "image/png"),
          width = 750
        ),
        easyClose = TRUE,
        footer = NULL,
        size = "l"
      ))
    })

    output$ddsdesign <- renderUI({
      if(is.null(values$expdesign))
        return(NULL)
      poss_covars <- colnames(values$expdesign)
      selectInput('dds_design', label = 'Select the design for your experiment: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })

    # server ui steps -----------------------------------------------------------
    output$ui_step2 <- renderUI({
      if (is.null(values$expdesign) | is.null(values$countmatrix))
        return(NULL)
      box(width = 12, title = "Step 2", status = "warning", solidHeader = TRUE,
          tagList(
            # as in https://groups.google.com/forum/#!topic/shiny-discuss/qQ8yICfvDu0
            h2("Select the DE design and create the DESeqDataSet object"),
            fluidRow(
              column(
                width = 6,
                uiOutput("ddsdesign"),
                uiOutput("ui_diydds"),
                hr(),
                # uiOutput("ok_dds"),
                verbatimTextOutput("debugdiy")
              )
            )
          ))
    })

    output$ui_stepanno <- renderUI({
      if (is.null(values$dds_obj)) ### and not provided already with sep annotation?
        return(NULL)

      box(width = 12, title = "Optional Step", status = "info", solidHeader = TRUE,
          tagList(
            h2("Create the annotation data frame for your dataset"),

            fluidRow(
              column(
                width = 8,
                uiOutput("ui_selectspecies"),
                verbatimTextOutput("speciespkg"),
                uiOutput("ui_idtype"),
                verbatimTextOutput("printDIYanno")

              )
            )
            ,
            uiOutput("ui_getanno")
          )
      )
    })

    output$ui_stepoutlier <- renderUI({
      if (is.null(values$dds_obj)) ### and not provided already with sep annotation?
        return(NULL)

      box(
        width = 12, title = "Optional Step", status = "info", solidHeader = TRUE,
        tagList(
          h2("Remove sample(s) from the current dataset - suspected outliers!"),
          
          fluidRow(
            column(
              width = 8,
              uiOutput("ui_selectoutliers"),
              uiOutput("outliersout"),
              verbatimTextOutput("printremoved")
            )
          )
        )
      )
    })

    output$ui_diydds <- renderUI({
      if (is.null(values$expdesign) | is.null(values$countmatrix) | is.null(input$dds_design))
        return(NULL)
      actionButton("button_diydds","Generate the dds object", class = "btn btn-success")
    })

    output$ui_getanno <- renderUI({
      if (is.null(values$dds_obj) ) ### and not provided already with sep annotation?
        return(NULL)
      shiny::validate(
        need(input$speciesSelect != "",
             "Select a species first in the panel")
      )
      actionButton("button_getanno","Retrieve the gene symbol annotation for the uploaded data", class = "btn btn-primary")
    })

    output$ui_nrcores <- renderUI({
      mincores <- 1
      maxcores <- BiocParallel::multicoreWorkers()
      sliderInput("nrcores",label = "Choose how many cores to use for computing:",
                  min = mincores, max = maxcores,value = 1,step = 1)
    })

    output$ui_step3 <- renderUI({
      if (is.null(values$dds_obj)) #
        return(NULL)
      box(width = 12, title = "Step 3", status = "success", solidHeader = TRUE,
          tagList(
            h2("Run DESeq!"),

            fluidRow(
              column(
                width = 4,
                uiOutput("ui_nrcores")
              )
            ),

            uiOutput("rundeseq"),
            verbatimTextOutput("printDIYresults"),
            uiOutput("ui_stepend")
          )
      )
    })

    output$ui_stepend <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      if (!"results" %in% mcols(mcols(values$dds_obj))$type) #
        return(NULL)
      tagList(
        h2("Good to go!"),
        box(width = 6, title = "Diagnostic plot", status = "info", solidHeader = TRUE,
            collapsible = TRUE, collapsed = TRUE,
            plotOutput("diagno_dispests"))
      )
    })

    output$diagno_dispests <- renderPlot({
      plotDispEsts(values$dds_obj)
    })

    # server ok objects -----------------------------------------------------------
    output$ok_cm <- renderUI({
      if (is.null(values$countmatrix))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_ed <- renderUI({
      if (is.null(values$expdesign))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_dds <- renderUI({
      if (is.null(values$dds_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_anno <- renderUI({
      if (is.null(values$annotation_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_resu <- renderUI({
      if (is.null(values$res_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })

    output$checkdds <- reactive({
      is.null(values$dds_obj)
    })
    output$checkresu<-reactive({
      is.null(values$res_obj)
    })

    outputOptions(output, 'checkresu', suspendWhenHidden=FALSE)
    outputOptions(output, 'checkdds', suspendWhenHidden=FALSE)

    output$dt_cm <- DT::renderDataTable({
      if(is.null(values$countmatrix))
        return(NULL)
      datatable(values$countmatrix,options = list(scrollX = TRUE,scrollY = "400px"))
    })

    output$dt_ed <- DT::renderDataTable({
      if(is.null(values$expdesign))
        return(NULL)
      datatable(values$expdesign,options = list(scrollX = TRUE))
    })

    # http://stackoverflow.com/questions/17024685/how-to-use-a-character-string-in-formula
    # http://stats.stackexchange.com/questions/29477/how-to-write-a-linear-model-formula-with-100-variables-in-r
    # http://stackoverflow.com/questions/7666807/anova-test-fails-on-lme-fits-created-with-pasted-formula/7668846#7668846
    diyDDS <- reactive({
      if(is.null(values$countmatrix) | is.null(values$expdesign) | is.null(input$dds_design))
        return(NULL)

      dds <- DESeqDataSetFromMatrix(countData = values$countmatrix,
                                    colData = values$expdesign,
                                    design=as.formula(paste0("~",paste(input$dds_design, collapse=" + "))))
      dds <- estimateSizeFactors(dds)
      return(dds)
    })

    observeEvent(input$button_diydds,
                 {
                   if(!is.null(values$countmatrix) & !is.null(values$expdesign))
                     values$dds_obj <- diyDDS()
                 })

    output$debugdiy <- renderPrint({
      if(!is.null(values$dds_obj)){
        print(values$dds_obj)
        print(design(values$dds_obj))
      }
    })

    # as in http://stackoverflow.com/questions/29716868/r-shiny-how-to-get-an-reactive-data-frame-updated-each-time-pressing-an-actionb
    observeEvent(input$uploadcmfile,
                 {
                   values$countmatrix <- readCountmatrix()
                 })

    observeEvent(input$uploadmetadatafile,
                 {
                   values$expdesign <- readMetadata()
                 })


    # server retrieving anno --------------------------------------------------
    annoSpecies_df <- 
      data.frame(species=c("","Anopheles","Arabidopsis","Bovine","Worm",
                           "Canine","Fly","Zebrafish","E coli strain K12",
                           "E coli strain Sakai","Chicken","Human","Mouse",
                           "Rhesus","Malaria","Chimp","Rat",
                           "Yeast","Streptomyces coelicolor", "Pig","Toxoplasma gondii",
                           "Xenopus"),
                 pkg=c("","org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
                       "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
                       "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
                       "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
                       "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
                       "org.Xl.eg.db"),
                 stringsAsFactors = FALSE)

    annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species),]
    # this one is relevant for creating links to the genes
    annoSpecies_df$ensembl_db <- c("","","","Bos_taurus","Canis_familiaris","Gallus_gallus","Pan_troglodytes",
                                   "","","Drosophila_melanogaster","Homo_sapiens","","Mus_musculus",
                                   "Sus_scrofa","Rattus_norvegicus","Macaca_mulatta","","","Caenorhabditis_elegans",
                                   "Xenopus_tropicalis","Saccharomyces_cerevisiae","Danio_rerio"
    )
    # this one is the shortcut for the limma::goana function
    annoSpecies_df$species_short[grep(pattern = "eg.db",annoSpecies_df$pkg)] <- gsub(".eg.db","",gsub("org.","",annoSpecies_df$pkg))[grep(pattern = "eg.db",annoSpecies_df$pkg) ]
    # to match to the goseq genome setting
    annoSpecies_df$goseq_shortcut <- c("","anoGam1","Arabidopsis","bosTau8","canFam3","galGal4","panTro4","E. coli K12","E. coli Sakai",
                                       "dm6","hg38","Malaria","mm10","susScr3","rn6","rheMac","","","ce11","xenTro","sacCer3","danRer10")
    rownames(annoSpecies_df) <- annoSpecies_df$species # easier to access afterwards
    # annoSpecies_df <- annoSpecies_df[annoSpecies_df$species %in% c("","Human", "Mouse", "Rat", "Fly", "Chimp"),]

    output$ui_selectspecies <- renderUI({
      if (is.null(values$dds_obj)) #
        return(NULL)
      selectInput("speciesSelect",label = "Select the species of your samples - it will also be used for enhancing result tables",
                  choices = annoSpecies_df$species,selected="")
    })

    output$ui_idtype <- renderUI({
      if (is.null(values$dds_obj)) #
        return(NULL)
      
      std_choices <- c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")
      if (input$speciesSelect!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
        require(annopkg,character.only=TRUE)
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput("idtype", "select the id type in your data", choices=std_choices)
    })

    output$speciespkg <- renderText({
      if (is.null(values$dds_obj)) #
        return(NULL)
      shiny::validate(
        need(input$speciesSelect!="",
             "Select a species - requires the corresponding annotation package"
        )
      )
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      shiny::validate(
        need(require(annopkg,character.only=TRUE),
             paste0("The package ",annopkg, " is not installed/available. Try installing it with BiocManager::install('",annopkg,"')"))
      )
      retmsg <- paste0(annopkg," - package available and loaded")
      # if (!require(annopkg,character.only=TRUE)) {
      # stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
      # }
      retmsg <- paste0(retmsg," - ",gsub(".eg.db","",gsub("org.","",annopkg)))
      retmsg
    })

    # server outliers --------------------------------------------------------
    output$ui_selectoutliers <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      else
        selectInput("selectoutliers","Select the samples to remove - candidate outliers",
                    choices = colnames(values$dds_obj), selected = NULL,multiple = TRUE
        )
    })

    output$outliersout <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      else
        actionButton("button_outliersout","Recompute the dds without some samples",class = "btn btn-primary")
    })

    observeEvent(input$button_outliersout,{
      withProgress({
        allsamples <- colnames(values$dds_obj)
        outliersamples <- input$selectoutliers

        keptsamples <- setdiff(allsamples,outliersamples)
        dds <- DESeqDataSetFromMatrix(countData = values$countmatrix[,keptsamples],
                                      colData = values$expdesign[keptsamples,],
                                      design= design(values$dds_obj)
                                      # design=as.formula(paste0("~",paste(input$dds_design, collapse=" + ")))
        )
        dds <- estimateSizeFactors(dds)

        # return(dds)
        # re-create the dds and keep track of which samples were removed
        values$removedsamples <- input$selectoutliers

        curr_species <- input$speciesSelect
        values$dds_obj <- dds
        updateSelectInput(session, inputId = "speciesSelect", selected = curr_species)
        # accordingly, reset the results
        values$res_obj <- NULL},
        message = "Removing selected samples from the current dataset")
    })

    output$printremoved <- renderPrint({
      print(values$removedsamples)
    })

    # server run deseq --------------------------------------------------------
    output$rundeseq <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      else
        actionButton("button_rundeseq","Run DESeq!", icon = icon("spinner"), class = "btn btn-success")
    })

    observeEvent(input$button_rundeseq,
                 {
                   withProgress(message="Running DESeq on your data...",
                                detail = "This step might take a while", value = 0,{
                                  # trick to keep species info while still changing the dds_obj
                                  curr_species <- input$speciesSelect
                                  incProgress(0.1)

                                  if(input$nrcores == 1)
                                    values$dds_obj <- DESeq(values$dds_obj)
                                  else
                                    # leave open option for computing in parallel?
                                    values$dds_obj <- DESeq(values$dds_obj,
                                                            parallel = TRUE,
                                                            BPPARAM = MulticoreParam(workers = input$nrcores))
                                  incProgress(0.89)
                                  updateSelectInput(session, inputId = "speciesSelect", selected = curr_species)
                                })
                 })
    
    observeEvent(input$speciesSelect,
                 {
                   curr_idtype <- values$cur_type
                   updateSelectInput(session, inputId = "idtype", selected = curr_idtype)
                 }
                 )

    output$printDIYresults <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Provide or construct a dds object")
      )
      shiny::validate(
        need("results" %in% mcols(mcols(values$dds_obj))$type ,
             "dds object provided, but couldn't find results. you should first run DESeq() with the button up here"
        )
      )
      summary(results(values$dds_obj), alpha = input$FDR)
    })

    # server counts overview --------------------------------------------------------
    current_countmat <- reactive({
      if(input$countstable_unit=="raw_counts")
        return(counts(values$dds_obj,normalized=FALSE))
      if(input$countstable_unit=="normalized_counts")
        return(counts(values$dds_obj,normalized=TRUE))
      if(input$countstable_unit=="rlog_counts")
        return(NULL) ## see if it is worth to keep in here or explore possibility with fast vst
      if(input$countstable_unit=="log10_counts")
        return(log10(1 + counts(values$dds_obj,normalized=TRUE)))
      if(input$countstable_unit=="tpm_counts")
        return(NULL) ## TODO!: assumes length of genes/exons as known, and is currently not required in the dds

    })

    output$showcountmat <- DT::renderDataTable({
      datatable(current_countmat())
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste0(input$countstable_unit,"table.csv")
      },
      content = function(file) {
        write.csv(current_countmat(), file)
      }
    )

    output$corrplot <- renderPlot({
      if(input$compute_pairwisecorr)
        withProgress(
          pair_corr(current_countmat(),
                    method=input$corr_method,
                    log = input$corr_uselogs,
                    use_subset = input$corr_usesubset),
          message = "Preparing the plot",
          detail = "this can take a while..."
        )
    })

    output$heatcorr <- renderPlot({
      if(input$compute_pairwisecorr)
        pheatmap(cor(current_countmat()))
    })

    output$pairwise_plotUI <- renderUI({
      if(!input$compute_pairwisecorr) return()

      plotOutput("corrplot", height = "1000px")
      # )
    })

    output$heatcorr_plotUI <- renderUI({
      if(!input$compute_pairwisecorr) return()

      plotOutput("heatcorr")
    })

    # overview on number of detected genes on different threshold types
    output$detected_genes <- renderPrint({
      t1 <- rowSums(counts(values$dds_obj))
      t2 <- rowMeans(counts(values$dds_obj,normalized=TRUE))

      thresh_rowsums <- input$threshold_rowsums
      thresh_rowmeans <- input$threshold_rowmeans
      abs_t1 <- sum(t1 > thresh_rowsums)
      rel_t1 <- 100 * mean(t1 > thresh_rowsums)
      abs_t2 <- sum(t2 > thresh_rowmeans)
      rel_t2 <- 100 * mean(t2 > thresh_rowmeans)

      cat("Number of detected genes:\n")
      cat(abs_t1,"genes have at least a sample with more than",thresh_rowsums,"counts\n")
      cat(paste0(round(rel_t1,3),"%"), "of the",nrow(values$dds_obj),
          "genes have at least a sample with more than",thresh_rowsums,"counts\n")
      cat(abs_t2,"genes have more than",thresh_rowmeans,"counts (normalized) on average\n")
      cat(paste0(round(rel_t2,3),"%"), "of the",nrow(values$dds_obj),
          "genes have more than",thresh_rowsums,"counts (normalized) on average\n")
      cat("Counts are ranging from", min(counts(values$dds_obj)),"to",max(counts(values$dds_obj)))
    })

    observeEvent(input$featfilt_dds,
                 {
                   t1 <- rowSums(counts(values$dds_obj))
                   t2 <- rowMeans(counts(values$dds_obj,normalized=TRUE))

                   thresh_rowsums <- input$threshold_rowsums
                   thresh_rowmeans <- input$threshold_rowmeans

                   if(input$filter_crit == "row sums") {
                     filt_dds <- values$dds_obj[t1 > thresh_rowsums, ]
                   } else {
                     filt_dds <- values$dds_obj[t2 > thresh_rowmeans, ]
                   }

                   # TODO: see if re-estimation of size factors is required
                   filt_dds <- estimateSizeFactors(filt_dds)

                   curr_species <- input$speciesSelect
                   values$dds_obj <- filt_dds
                   updateSelectInput(session, inputId = "speciesSelect", selected = curr_species)
                 })

    # server managing gene lists --------------------------------------------------------
    ## gene lists upload

    observeEvent(input$gl1,
                 {
                   mydf <- as.data.frame(gl1(),stringsAsFactors=FALSE)
                   names(mydf) <- "Gene Symbol"
                   values$genelist1 <- mydf
                 })

    gl1 <- reactive({
      if (is.null(input$gl1)) {
        # User has not uploaded a file yet
        return(data.frame())
      } else {
        gl1 <- readLines(input$gl1$datapath)
        return(gl1)
      }
    })

    observeEvent(input$gl2,
                 {
                   mydf <- as.data.frame(gl2(),stringsAsFactors=FALSE)
                   names(mydf) <- "Gene Symbol"
                   values$genelist2 <- mydf
                 })

    gl2 <- reactive({
      if (is.null(input$gl2)) {
        # User has not uploaded a file yet
        return(data.frame())
      } else {
        gl2 <- readLines(input$gl2$datapath)
        return(gl2)
      }
    })

    observeEvent(input$gl_ma,
                 {
                   mydf <- as.data.frame(gl_ma(),stringsAsFactors=FALSE)
                   names(mydf) <- "Gene Symbol"
                   values$genelist_ma <- mydf
                 })

    gl_ma <- reactive({
      if (is.null(input$gl_ma)) {
        # User has not uploaded a file yet
        return(data.frame())
      } else {
        gl_ma <- readLines(input$gl_ma$datapath)
        return(gl_ma)
      }
    })

    output$debuggls <- renderPrint({
      values$genelist1
      # values$genelist2
    })

    values$genelistUP <- reactive({

      res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = input$FDR)
      res_tbl_UP <- res_tbl[res_tbl$log2FoldChange > 0 & !is.na(res_tbl$padj),]
      # res_tbl_DOWN <- res_tbl[res_tbl$log2FoldChange < 0 & !is.na(res_tbl$padj),]

      # this will have to be modified!
      # res_tbl_UP$symbol <- anno_df$gene_name[match(res_tbl_UP$id,anno_df$gene_id)]

      listUP <- res_tbl_UP$symbol
      return(listUP)
    })

    values$genelistDOWN <- reactive({

      res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = input$FDR)
      # res_tbl_UP <- res_tbl[res_tbl$log2FoldChange > 0 & !is.na(res_tbl$padj),]
      res_tbl_DOWN <- res_tbl[res_tbl$log2FoldChange < 0 & !is.na(res_tbl$padj),]

      # this will have to be modified!
      # res_tbl_DOWN$symbol <- anno_df$gene_name[match(res_tbl_DOWN$id,anno_df$gene_id)]

      listDOWN <- res_tbl_DOWN$symbol
      return(listDOWN)
    })

    values$genelistUPDOWN <- reactive({

      res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = input$FDR)

      # this will have to be modified!
      # res_tbl$symbol <- anno_df$gene_name[match(res_tbl$id,anno_df$gene_id)]

      listUPDOWN <- res_tbl$symbol
      return(listUPDOWN)
    })

    ## list of gene lists
    gll <- reactive({
      mylist <- list(listUP = values$genelistUP(),
                     listDOWN = values$genelistDOWN(),
                     listUPDOWN = values$genelistUPDOWN(),
                     list1 = as.character(values$genelist1$`Gene Symbol`),
                     list2 = as.character(values$genelist2$`Gene Symbol`),
                     list3 = NULL) # will be changed to be the ones selected by the user


      gll_nonempty <- mylist[!sapply(mylist,is.null)]

      # plus, add toggles to selectively keep only some lists?

      lists_tokeep <- names(mylist)[which(c(input$toggle_up,
                                            input$toggle_down,
                                            input$toggle_updown,
                                            input$toggle_list1,
                                            input$toggle_list2,
                                            input$toggle_list3))]
      gll_final <- gll_nonempty[match(lists_tokeep,names(gll_nonempty))]
    })

    output$debuglists <- renderPrint({
      # length(gll_nonempty)
      # length(gll())
      # lapply(gll(),length)
      print(gll())
    })

    output$vennlists <- renderPlot({
      shiny::validate(
        need(all(sapply(gll(),function(arg) !is.null(arg))),
             message = "Some lists are empty - make sure you extracted the results using the annotation object")

      )

      gplots::venn(gll())
    })

    output$upsetLists <- renderPlot({
      shiny::validate(
        need(sum(sapply(gll(),function(arg) length(arg)>0)) > 1,
             message = "Make sure you provide at least two sets")
      )
      UpSetR::upset(fromList(gll()))
    })

    observeEvent(input$button_getanno,
                 {
                   withProgress(message="Retrieving the annotation...",
                                detail = "Locating package", value = 0,{

                                  annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
                                  incProgress(0.1,detail = "Matching identifiers")
                                  annotation_obj <- get_annotation_orgdb(values$dds_obj,orgdb_species = annopkg, idtype = input$idtype)
                                  values$annotation_obj <- annotation_obj
                                  # and also, set the species in the reactiveValues
                                  values$cur_species <- input$speciesSelect
                                  values$cur_type <- input$idtype
                                })
                 })

    output$printDIYanno <- renderPrint({
      print(head(values$annotation_obj))
    })

    output$printUPgenes <- renderPrint({
      print(head(values$genelistUP()))
      print(str(values$genelistUP()))

      organism <- annoSpecies_df[values$cur_species,]$species_short
      backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
      inputType <- "SYMBOL" # will be replaced by input$...
      # annopkg <- paste0("org.",organism,".eg.db")
      annopkg <- annoSpecies_df[values$cur_species,]$pkg
      listGenesEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUP(),
                                                            column="ENTREZID", keytype=inputType))
      listBackgroundEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                 column="ENTREZID", keytype=input$idtype))

      # print(values$genelistUP())
      print(str(listGenesEntrez))
      print(class(listGenesEntrez))
      print(str(listBackgroundEntrez))
      print(class(listBackgroundEntrez))
      print(head(listGenesEntrez))
      print(head(listBackgroundEntrez))

      # values$gse_up <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
      # ontology="BP", # could be ideally replaced by input$
      # number=200)
    })


    ### UP
    observeEvent(input$button_enrUP,
                 {
                   withProgress(message="Performing Gene Set Enrichment on upregulated genes...",value = 0,{
                     organism <- annoSpecies_df[values$cur_species,]$species_short
                     backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
                     inputType <- "SYMBOL" # will be replaced by input$...
                     # annopkg <- paste0("org.",organism,".eg.db")
                     annopkg <- annoSpecies_df[values$cur_species,]$pkg
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUP(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype=input$idtype))
                     incProgress(0.1, detail = "IDs mapped")
                     values$gse_up <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
                                                   ontology= input$go_cats[1],
                                                   number=200)

                     incProgress(0.7, detail = "adding gene names to GO terms") # good indicator for showing it has progressed
                     go_ids <- rownames(values$gse_up)
                     allegs_list <- lapply(go_ids, function(arg) AnnotationDbi::get(arg, get(paste0("org.",organism,".egGO2ALLEGS"))))
                     genes_list <- lapply(allegs_list, function(arg) unlist(AnnotationDbi::mget(arg,get(paste0("org.",organism,".egSYMBOL")))))
                     degenes <- values$genelistUP()
                     DEgenes_list <- lapply(genes_list, function(arg) intersect(arg,degenes))

                     values$gse_up$genes <- unlist(lapply(DEgenes_list,function(arg) paste(arg,collapse=",")))
                   })
                 })

    observeEvent(input$button_enrUP_goseq,
                 {
                   withProgress(message="GOSEQ - Performing Gene Set Enrichment on upregulated genes...",value = 0,{

                     de.genes <- values$genelistUP() # assumed to be in symbols
                     assayed.genes.ids <- rownames(values$dds_obj) # as IDs, but then to be converted back
                     assayed.genes <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                             keys=assayed.genes.ids,
                                             column="SYMBOL",
                                             keytype=input$idtype,
                                             multiVals="first")
                     de.genes.ids <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                            keys=de.genes,
                                            column="ENSEMBL",
                                            keytype="SYMBOL",
                                            multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")

                     values$gse_up_goseq <- goseqTable(de.genes.ids,
                                                       assayed.genes.ids,
                                                       genome = annoSpecies_df[values$cur_species,]$goseq_short,
                                                       id= "ensGene",
                                                       testCats=paste0("GO:",input$go_cats),
                                                       FDR_GO_cutoff = 1,
                                                       nTop = 200,
                                                       addGeneToTerms=TRUE,
                                                       orgDbPkg = annoSpecies_df[values$cur_species,]$pkg # ,
                     )

                     incProgress(0.89)

                   })
                 })

    observeEvent(input$button_enrUP_topgo,
                 {
                   withProgress(message="TOPGO - Performing Gene Set Enrichment on upregulated genes...",value = 0,{

                     de_symbols <- values$genelistUP() # assumed to be in symbols
                     bg_ids <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
                     bg_symbols <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                          keys=bg_ids,
                                          column="SYMBOL",
                                          keytype=input$idtype,
                                          multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")
                     # library(topGO)
                     # requireNamespace("topGO")
                     values$topgo_up <- pcaExplorer::topGOtable(de_symbols, bg_symbols,
                                                                ontology = input$go_cats[1],
                                                                mapping = annoSpecies_df[values$cur_species,]$pkg,
                                                                geneID = "symbol",addGeneToTerms = TRUE)
                     incProgress(0.89)
                   })
                 })

    ### DOWN
    observeEvent(input$button_enrDOWN,
                 {
                   withProgress(message="Performing Gene Set Enrichment on downregulated genes...",value = 0,{
                     organism <- annoSpecies_df[values$cur_species,]$species_short
                     backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
                     inputType <- "SYMBOL" # will be replaced by input$...
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistDOWN(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype=input$idtype))
                     incProgress(0.1, detail = "IDs mapped")
                     values$gse_down <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
                                                     ontology=input$go_cats[1],
                                                     number=200)


                     incProgress(0.7, detail = "adding gene names to GO terms") # good indicator for showing it has progressed
                     go_ids <- rownames(values$gse_down)
                     allegs_list <- lapply(go_ids, function(arg) AnnotationDbi::get(arg, get(paste0("org.",organism,".egGO2ALLEGS"))))
                     genes_list <- lapply(allegs_list, function(arg) unlist(AnnotationDbi::mget(arg,get(paste0("org.",organism,".egSYMBOL")))))
                     degenes <- values$genelistDOWN()
                     DEgenes_list <- lapply(genes_list, function(arg) intersect(arg,degenes))

                     values$gse_down$genes <- unlist(lapply(DEgenes_list,function(arg) paste(arg,collapse=",")))
                   })
                 })

    observeEvent(input$button_enrDOWN_goseq,
                 {
                   withProgress(message="GOSEQ - Performing Gene Set Enrichment on downregulated genes...",value = 0,{

                     de.genes <- values$genelistDOWN() # assumed to be in symbols
                     assayed.genes.ids <- rownames(values$dds_obj) # as IDs, but then to be converted back
                     assayed.genes <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                             keys=assayed.genes.ids,
                                             column="SYMBOL",
                                             keytype=input$idtype,
                                             multiVals="first")
                     de.genes.ids <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                            keys=de.genes,
                                            column="ENSEMBL",
                                            keytype="SYMBOL",
                                            multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")

                     values$gse_down_goseq <- goseqTable(de.genes.ids,
                                                         assayed.genes.ids,
                                                         genome = annoSpecies_df[values$cur_species,]$goseq_short,
                                                         id= "ensGene",
                                                         testCats=paste0("GO:",input$go_cats),
                                                         FDR_GO_cutoff = 1,
                                                         nTop = 200,
                                                         addGeneToTerms=TRUE,
                                                         orgDbPkg = annoSpecies_df[values$cur_species,]$pkg # ,
                     )

                     incProgress(0.89)

                   })
                 })

    observeEvent(input$button_enrDOWN_topgo,
                 {
                   withProgress(message="TOPGO - Performing Gene Set Enrichment on downregulated genes...",value = 0,{

                     de_symbols <- values$genelistDOWN() # assumed to be in symbols
                     bg_ids <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
                     bg_symbols <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                          keys=bg_ids,
                                          column="SYMBOL",
                                          keytype=input$idtype,
                                          multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")
                     # library(topGO)
                     # requireNamespace("topGO")
                     values$topgo_down <- pcaExplorer::topGOtable(de_symbols, bg_symbols,
                                                                  ontology = input$go_cats[1], # will take the first ontology
                                                                  mapping = annoSpecies_df[values$cur_species,]$pkg,
                                                                  geneID = "symbol",addGeneToTerms = TRUE)
                     incProgress(0.89)


                   })
                 })

    ### UPDOWN
    observeEvent(input$button_enrUPDOWN,
                 {
                   withProgress(message="Performing Gene Set Enrichment on up- and downregulated genes...",value = 0,{
                     organism <- annoSpecies_df[values$cur_species,]$species_short
                     backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
                     inputType <- "SYMBOL" # will be replaced by input$...
                     # annopkg <- paste0("org.",organism,".eg.db")
                     annopkg <- annoSpecies_df[values$cur_species,]$pkg
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUPDOWN(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype=input$idtype))
                     incProgress(0.1, detail = "IDs mapped")
                     values$gse_updown <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
                                                       ontology=input$go_cats[1],
                                                       number=200)

                     incProgress(0.7, detail = "adding gene names to GO terms") # good indicator for showing it has progressed
                     go_ids <- rownames(values$gse_updown)
                     allegs_list <- lapply(go_ids, function(arg) AnnotationDbi::get(arg, get(paste0("org.",organism,".egGO2ALLEGS"))))
                     genes_list <- lapply(allegs_list, function(arg) unlist(AnnotationDbi::mget(arg,get(paste0("org.",organism,".egSYMBOL")))))
                     degenes <- values$genelistDOWN()
                     DEgenes_list <- lapply(genes_list, function(arg) intersect(arg,degenes))

                     # values$gse_down$genes[1:20] <- DEgenes_list
                     # lapply(values$gse_down,class)
                     values$gse_updown$genes <- unlist(lapply(DEgenes_list,function(arg) paste(arg,collapse=",")))
                   })
                 })

    observeEvent(input$button_enrUPDOWN_goseq,
                 {
                   withProgress(message="GOSEQ - Performing Gene Set Enrichment on up and downregulated genes...",value = 0,{

                     de.genes <- values$genelistUPDOWN() # assumed to be in symbols
                     assayed.genes.ids <- rownames(values$dds_obj) # as IDs, but then to be converted back
                     assayed.genes <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                             keys=assayed.genes.ids,
                                             column="SYMBOL",
                                             keytype=input$idtype,
                                             multiVals="first")
                     de.genes.ids <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                            keys=de.genes,
                                            column="ENSEMBL",
                                            keytype="SYMBOL",
                                            multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")

                     values$gse_updown_goseq <- goseqTable(de.genes.ids,
                                                           assayed.genes.ids,
                                                           genome = annoSpecies_df[values$cur_species,]$goseq_short,
                                                           id= "ensGene",
                                                           testCats=paste0("GO:",input$go_cats),
                                                           FDR_GO_cutoff = 1,
                                                           nTop = 200,
                                                           addGeneToTerms=TRUE,
                                                           orgDbPkg = annoSpecies_df[values$cur_species,]$pkg # ,
                     )

                     incProgress(0.89)

                   })
                 })

    observeEvent(input$button_enrUPDOWN_topgo,
                 {
                   withProgress(message="TOPGO - Performing Gene Set Enrichment on up and downregulated genes...",value = 0,{

                     de_symbols <- values$genelistUPDOWN() # assumed to be in symbols
                     bg_ids <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
                     bg_symbols <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                          keys=bg_ids,
                                          column="SYMBOL",
                                          keytype=input$idtype,
                                          multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")
                     # library(topGO)
                     # requireNamespace("topGO")
                     values$topgo_updown <- pcaExplorer::topGOtable(de_symbols, bg_symbols,
                                                                    ontology = input$go_cats[1],
                                                                    mapping = annoSpecies_df[values$cur_species,]$pkg,
                                                                    geneID = "symbol",addGeneToTerms = TRUE)
                     incProgress(0.89)
                   })
                 })

    ### LIST1
    observeEvent(input$button_enrLIST1,
                 {
                   withProgress(message="Performing Gene Set Enrichment on upregulated genes...",value = 0,{
                     organism <- annoSpecies_df[values$cur_species,]$species_short
                     backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
                     inputType <- "SYMBOL" # will be replaced by input$...
                     # annopkg <- paste0("org.",organism,".eg.db")
                     annopkg <- annoSpecies_df[values$cur_species,]$pkg
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
                     }
                     listGenesEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = as.character(values$genelist1$`Gene Symbol`),
                                                              column="ENTREZID", keytype=inputType)
                     listBackgroundEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                   column="ENTREZID", keytype=input$idtype)
                     incProgress(0.1, detail = "IDs mapped")
                     values$gse_list1 <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
                                                      ontology=input$go_cats[1],
                                                      number=200)

                     incProgress(0.7, detail = "adding gene names to GO terms") # good indicator for showing it has progressed
                     go_ids <- rownames(values$gse_list1)
                     allegs_list <- lapply(go_ids, function(arg) AnnotationDbi::get(arg, get(paste0("org.",organism,".egGO2ALLEGS"))))
                     genes_list <- lapply(allegs_list, function(arg) unlist(AnnotationDbi::mget(arg,get(paste0("org.",organism,".egSYMBOL")))))
                     degenes <- values$genelistDOWN()
                     DEgenes_list <- lapply(genes_list, function(arg) intersect(arg,degenes))

                     values$gse_list1$genes <- unlist(lapply(DEgenes_list,function(arg) paste(arg,collapse=",")))
                   })
                 })

    observeEvent(input$button_enrLIST1_goseq,
                 {
                   withProgress(message="GOSEQ - Performing Gene Set Enrichment on list 1 genes...",value = 0,{

                     de.genes <- values$genelist1$`Gene Symbol` # assumed to be in symbols
                     assayed.genes.ids <- rownames(values$dds_obj) # as IDs, but then to be converted back
                     assayed.genes <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                             keys=assayed.genes.ids,
                                             column="SYMBOL",
                                             keytype=input$idtype,
                                             multiVals="first")
                     de.genes.ids <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                            keys=de.genes,
                                            column="ENSEMBL",
                                            keytype="SYMBOL",
                                            multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")

                     values$gse_list1_goseq <- goseqTable(de.genes.ids,
                                                          assayed.genes.ids,
                                                          genome = annoSpecies_df[values$cur_species,]$goseq_short,
                                                          id= "ensGene",
                                                          testCats=paste0("GO:",input$go_cats),
                                                          FDR_GO_cutoff = 1,
                                                          nTop = 200,
                                                          addGeneToTerms=TRUE,
                                                          orgDbPkg = annoSpecies_df[values$cur_species,]$pkg # ,
                     )

                     incProgress(0.89)

                   })
                 })

    observeEvent(input$button_enrLIST1_topgo,
                 {
                   withProgress(message="TOPGO - Performing Gene Set Enrichment on list1 genes...",value = 0,{

                     de_symbols <- values$genelist1$`Gene Symbol` # assumed to be in symbols
                     bg_ids <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
                     bg_symbols <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                          keys=bg_ids,
                                          column="SYMBOL",
                                          keytype=input$idtype,
                                          multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")
                     # library(topGO)
                     # requireNamespace("topGO")
                     values$topgo_list1 <- pcaExplorer::topGOtable(de_symbols, bg_symbols,
                                                                   ontology = input$go_cats[1],
                                                                   mapping = annoSpecies_df[values$cur_species,]$pkg,
                                                                   geneID = "symbol",addGeneToTerms = TRUE)
                     incProgress(0.89)
                   })
                 })
    ### LIST2
    observeEvent(input$button_enrLIST2,
                 {
                   withProgress(message="Performing Gene Set Enrichment on upregulated genes...",value = 0,{
                     organism <- annoSpecies_df[values$cur_species,]$species_short
                     backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
                     inputType <- "SYMBOL" # will be replaced by input$...
                     # annopkg <- paste0("org.",organism,".eg.db")
                     annopkg <- annoSpecies_df[values$cur_species,]$pkg
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
                     }
                     listGenesEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = as.character(values$genelist2$`Gene Symbol`),
                                                              column="ENTREZID", keytype=inputType)
                     listBackgroundEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                   column="ENTREZID", keytype=input$idtype)
                     incProgress(0.1, detail = "IDs mapped")
                     values$gse_list2 <- limma::topGO(limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
                                                      ontology=input$go_cats[1],
                                                      number=200)
                     incProgress(0.7, detail = "adding gene names to GO terms") # good indicator for showing it has progressed
                     go_ids <- rownames(values$gse_list2)
                     allegs_list <- lapply(go_ids, function(arg) AnnotationDbi::get(arg, get(paste0("org.",organism,".egGO2ALLEGS"))))
                     genes_list <- lapply(allegs_list, function(arg) unlist(AnnotationDbi::mget(arg,get(paste0("org.",organism,".egSYMBOL")))))
                     degenes <- values$genelistDOWN()
                     DEgenes_list <- lapply(genes_list, function(arg) intersect(arg,degenes))

                     values$gse_list2$genes <- unlist(lapply(DEgenes_list,function(arg) paste(arg,collapse=",")))
                   })
                 })

    observeEvent(input$button_enrLIST2_goseq,
                 {
                   withProgress(message="GOSEQ - Performing Gene Set Enrichment on list 2 genes...",value = 0,{

                     de.genes <- values$genelist2$`Gene Symbol` # assumed to be in symbols
                     assayed.genes.ids <- rownames(values$dds_obj) # as IDs, but then to be converted back
                     assayed.genes <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                             keys=assayed.genes.ids,
                                             column="SYMBOL",
                                             keytype=input$idtype,
                                             multiVals="first")
                     de.genes.ids <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                            keys=de.genes,
                                            column="ENSEMBL",
                                            keytype="SYMBOL",
                                            multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")

                     values$gse_list2_goseq <- goseqTable(de.genes.ids,
                                                          assayed.genes.ids,
                                                          genome = annoSpecies_df[values$cur_species,]$goseq_short,
                                                          id= "ensGene",
                                                          testCats=paste0("GO:",input$go_cats),
                                                          FDR_GO_cutoff = 1,
                                                          nTop = 200,
                                                          addGeneToTerms=TRUE,
                                                          orgDbPkg = annoSpecies_df[values$cur_species,]$pkg # ,
                     )

                     incProgress(0.89)

                   })
                 })

    observeEvent(input$button_enrLIST2_topgo,
                 {
                   withProgress(message="TOPGO - Performing Gene Set Enrichment on list2 genes...",value = 0,{

                     de_symbols <- values$genelist2$`Gene Symbol` # assumed to be in symbols
                     bg_ids <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
                     bg_symbols <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                                          keys=bg_ids,
                                          column="SYMBOL",
                                          keytype=input$idtype,
                                          multiVals="first")
                     incProgress(0.1, detail = "IDs mapped")
                     # library(topGO)
                     # requireNamespace("topGO")
                     values$topgo_list2 <- pcaExplorer::topGOtable(de_symbols, bg_symbols,
                                                                   ontology = input$go_cats[1],
                                                                   mapping = annoSpecies_df[values$cur_species,]$pkg,
                                                                   geneID = "symbol",addGeneToTerms = TRUE)
                     incProgress(0.89)
                   })
                 })




    # server gse datatables --------------------------------------------------------
    output$DT_gse_up <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_up))
        return(NULL)
      mytbl <- values$gse_up
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_down <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_down))
        return(NULL)
      mytbl <- values$gse_down
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_updown <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_updown))
        return(NULL)
      mytbl <- values$gse_updown
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list1 <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_list1))
        return(NULL)
      mytbl <- values$gse_list1
      # mytbl$GOid <- rownames(mytbl)
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list2 <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_list2))
        return(NULL)
      mytbl <- values$gse_list2
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
      datatable(mytbl,escape=FALSE)
    })


    output$DT_gse_up_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_up))
        return(NULL)
      mytbl <- values$topgo_up
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      DT::datatable(mytbl,escape=FALSE, selection=list(mode="single"))
    })
    output$DT_gse_down_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_down))
        return(NULL)
      mytbl <- values$topgo_down
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      DT::datatable(mytbl,escape=FALSE, selection=list(mode="single"))
    })
    output$DT_gse_updown_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_updown))
        return(NULL)
      mytbl <- values$topgo_updown
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      DT::datatable(mytbl,escape=FALSE, selection=list(mode="single"))
    })
    output$DT_gse_list1_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_list1))
        return(NULL)
      mytbl <- values$topgo_list1
      # mytbl$GOid <- rownames(mytbl)
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      DT::datatable(mytbl,escape=FALSE, selection=list(mode="single"))
    })
    output$DT_gse_list2_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_list2))
        return(NULL)
      mytbl <- values$topgo_list2
      # mytbl$GOid <- rownames(mytbl)
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      DT::datatable(mytbl,escape=FALSE, selection=list(mode="single"))
    })


    output$DT_gse_up_goseq <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_up_goseq))
        return(NULL)
      mytbl <- values$gse_up_goseq
      mytbl$category <- createLinkGO(mytbl$category)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_down_goseq <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_down_goseq))
        return(NULL)
      mytbl <- values$gse_down_goseq
      mytbl$category <- createLinkGO(mytbl$category)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_updown_goseq <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_updown_goseq))
        return(NULL)
      mytbl <- values$gse_updown_goseq
      mytbl$category <- createLinkGO(mytbl$category)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list1_goseq <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_list1_goseq))
        return(NULL)
      mytbl <- values$gse_list1_goseq
      # mytbl$GOid <- rownames(mytbl)
      mytbl$category <- createLinkGO(mytbl$category)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list2_goseq <- DT::renderDataTable({
      # if not null...
      if(is.null(values$gse_list2_goseq))
        return(NULL)
      mytbl <- values$gse_list2_goseq
      # mytbl$GOid <- rownames(mytbl)
      mytbl$category <- createLinkGO(mytbl$category)
      datatable(mytbl,escape=FALSE)
    })


    # server gse heatmaps --------------------------------------------------------
    output$goterm_heatmap_up_topgo <- renderPlot({

      s <- input$DT_gse_up_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      # allow only one selected line
      mygenes <- values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,input$idtype,"SYMBOL",multiVals="first")
      log2things <- assay(normTransform(values$dds_obj))
      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...
      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)

    })

    output$goterm_heatmap_down_topgo <- renderPlot({

      s <- input$DT_gse_down_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      # allow only one selected line
      mygenes <- values$topgo_down[input$DT_gse_down_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_down[input$DT_gse_down_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_down[input$DT_gse_down_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,input$idtype,"SYMBOL",multiVals="first")
      log2things <- assay(normTransform(values$dds_obj))
      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...
      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)

    })

    output$goterm_heatmap_updown_topgo <- renderPlot({

      s <- input$DT_gse_updown_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      values$topgo_updown[input$DT_gse_updown_topgo_rows_selected,]$genes

      # allow only one selected line
      mygenes <- values$topgo_updown[input$DT_gse_updown_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_updown[input$DT_gse_updown_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_updown[input$DT_gse_updown_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,input$idtype,"SYMBOL",multiVals="first")
      log2things <- assay(normTransform(values$dds_obj))
      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...
      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)

    })

    output$goterm_heatmap_l1_topgo <- renderPlot({

      s <- input$DT_gse_list1_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      # allow only one selected line
      mygenes <- values$topgo_list1[input$DT_gse_list1_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_list1[input$DT_gse_list1_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_list1[input$DT_gse_list1_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,input$idtype,"SYMBOL",multiVals="first")
      log2things <- assay(normTransform(values$dds_obj))
      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...
      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)

    })

    output$goterm_heatmap_l2_topgo <- renderPlot({

      s <- input$DT_gse_list2_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      # allow only one selected line
      mygenes <- values$topgo_list2[input$DT_gse_list2_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_list2[input$DT_gse_list2_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_list2[input$DT_gse_list2_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,input$idtype,"SYMBOL",multiVals="first")
      log2things <- assay(normTransform(values$dds_obj))
      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...
      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)

    })

    # server signature explorer ------------------------------------------------------
    output$sig_ui_gmtin <- renderUI({
      fileInput("sig_gmtin","gmt input file")
    })
    
    loaded_gmt <- reactive({
      if (is.null(input$sig_gmtin))
        return(NULL)
      mysigs <- read_gmt(input$sig_gmtin$datapath)
      return(mysigs)
    })
    
    observeEvent(input$sig_gmtin,
                 {
                   values$gene_signatures <- loaded_gmt()
                 })
    
    output$sig_ui_nrsigs <- renderUI({
      if(!is.null(values$gene_signatures))
        return(valueBox("Gene signatures",
                        paste0(length(values$gene_signatures), " gene signatures"),
                        icon = icon("list"),
                        color = "green",width = NULL))
      else
        return(valueBox("Gene signatures",
                        "yet to be loaded",
                        icon = icon("list"),
                        color = "red",width = NULL))
    })
    
    observeEvent(input$sig_button_computevst,
                 {
                   withProgress(message="Computing the variance stabilized transformed data...",
                                detail = "This step can take a little while",
                                value = 0,{
                                  values$vst_obj <- vst(values$dds_obj)
                                })
                 })

    output$sig_ui_selectsig <- renderUI({
      if(!is.null(values$gene_signatures))
        return(selectizeInput("sig_selectsig", label = "Select the gene signature",
                              choices = NULL, selected = NULL, multiple = FALSE))
      else
        return(NULL)
    })

    observe({
      updateSelectizeInput(session = session, inputId = 'sig_selectsig', choices = c(Choose = '', names(values$gene_signatures)), server = TRUE)
    })
    
    output$sig_sigmembers <- renderPrint({
      values$gene_signatures[[input$sig_selectsig]]
    })
    
    output$sig_ui_annocoldata <- renderUI({
      if(!is.null(values$dds_obj))
        return(selectizeInput("sig_annocoldata", label = "Select the colData to decorate",
                              choices = names(colData(values$dds_obj)),
                              selected = NULL, multiple = TRUE))
      else
        return(NULL)
    })


    output$sig_ui_id_data <- renderUI({
      if (is.null(values$dds_obj)) #
        return(NULL)
      validate(
        need(!is.null(input$speciesSelect), message = "Please specify the species in the Data Setup panel")
      )
      
      std_choices <- c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")
      if (input$speciesSelect!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput("sig_id_data", "select the id type in your dds data", choices=std_choices)
    })
    
    output$sig_ui_id_sigs <- renderUI({
      if (is.null(values$gene_signatures)) #
        return(NULL)
      validate(
        need(!is.null(input$speciesSelect), message = "Please specify the species in the Data Setup panel")
      )
      
      std_choices <- c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")
      if (input$speciesSelect!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput("sig_id_sigs", "select the id type in your signatures", choices=std_choices)
    })
    
    available_orgdb <- rownames(installed.packages())[
      grep(pattern = "^org.*db$",rownames(installed.packages()))]
    
    output$sig_ui_orgdbpkg <- renderUI({
      selectInput("sig_orgdbpkg", "Select the organism package for matching", 
                  choices=c("",available_orgdb),selected = "")
    })
    
    observeEvent(input$speciesSelect,
                 {
                   suggested_orgdb <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
                   if(suggested_orgdb %in% available_orgdb)
                     updateSelectInput(session, inputId = "sig_orgdbpkg", selected = suggested_orgdb)
                 })
    
    observeEvent(input$sig_convert_setup,
                 {
                   withProgress(message="Matching the identifiers",
                                detail = "Locating package", value = 0,{
                                  require(input$sig_orgdbpkg,character.only=TRUE)
                                  incProgress(0.1,detail = "Matching identifiers")
                                  
                                  x <- get(input$sig_orgdbpkg)
                                  values$anno_vec <- mapIds(x, rownames(values$dds_obj),
                                                            column = input$sig_id_sigs,
                                                            keytype = input$sig_id_data)
                                  
                                })
                 })
    
    output$sig_convcheck <- renderPrint({
      head(values$anno_vec)
    })
    
    output$sig_heat <- renderPlot({
      validate(
        need(!is.null(values$gene_signatures), message = "Please provide some gene signatures in gmt format"),
        need(!is.null(values$vst_obj), message = "Compute the vst transformed data"),
        need(!is.null(values$anno_vec), message = "Setup the conversion between data ids and signature ids"),
        need((!is.null(values$res_obj) | !input$sig_useDEonly),
             message = "Please compute the results first if you want to subset to DE genes only"),
        need(input$sig_selectsig!="", message = "Select a signature")
      )
      
      print(
        sig_heatmap(
          values$vst_obj,
          my_signature = values$gene_signatures[[input$sig_selectsig]],
          res_data = values$res_obj,
          FDR = input$FDR,
          de_only = input$sig_useDEonly,
          annovec = values$anno_vec,
          # anno_colData = colData(values$vst_obj)[,input$sig_annocoldata, drop = FALSE],
          title = names(values$gene_signatures)[match(input$sig_selectsig,names(values$gene_signatures))],
          cluster_rows = input$sig_clusterrows,
          cluster_cols = input$sig_clustercols,
          center_mean = input$sig_centermean,
          scale_row = input$sig_scalerow
        ))
              
    })


    # server ui update/observers --------------------------------------------------------
    output$color_by <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      poss_covars <- names(colData(values$dds_obj))
      selectInput('color_by', label = 'Group/color by: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })

    # this trick speeds up the populating of the select(ize) input widgets,
    # see http://stackoverflow.com/questions/38438920/shiny-selectinput-very-slow-on-larger-data-15-000-entries-in-browser
    observe({
      updateSelectizeInput(
        session = session,
        inputId = 'avail_ids',
        choices = c(Choose = '', rownames(values$dds_obj)),
        server = TRUE)
    })

    observe({
      updateSelectizeInput(
        session = session,
        inputId = 'avail_symbols',
        choices = c(Choose = '', values$annotation_obj$gene_name[match(rownames(values$dds_obj), values$annotation_obj$gene_id)]),
        server = TRUE)
    })

    output$available_genes <- renderUI({
      if(!is.null(values$annotation_obj)) {
        selectizeInput("avail_symbols", label = "Select the gene(s) of interest",
                       choices = NULL, selected = NULL, multiple = TRUE)
      } else { # else use the rownames as identifiers
        selectizeInput("avail_ids", label = "Select the gene(s) of interest - ids",
                       choices = NULL, selected = NULL, multiple = TRUE)
      }
    })

    design_factors <- reactive({
      rev(attributes(terms.formula(design(values$dds_obj)))$term.labels)
    })

    output$choose_fac <- renderUI({
      selectInput("choose_expfac",label = "Choose the experimental factor to build the contrast upon (must be in the design formula)",
                  choices = c("",design_factors()), selected = "")
    })

    observe({
      updateSelectizeInput(session = session, inputId = 'color_by', selected = input$choose_expfac)
    })


    # server DE results --------------------------------------------------------
    # nrl <- reactive
    output$lrtavailable <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      nrl <- length(levels(colData(values$dds_obj)[,fac1]))

      if(nrl > 2)
        p("I can perform a LRT test on the chosen factor, select the full and the reduced model")

    })

    output$lrtfull <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      nrl <- length(levels(colData(values$dds_obj)[,fac1]))

      if(nrl > 2)
        selectInput("choose_lrt_full",label = "Choose the factors for the full model",
                    choices = c("",design_factors()), selected = "", multiple = TRUE)

    })

    output$lrtreduced <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      nrl <- length(levels(colData(values$dds_obj)[,fac1]))

      if(nrl > 2)
        selectInput("choose_lrt_reduced",label = "Choose the factor(s) for the reduced model",
                    choices = c("",design_factors()), selected = "", multiple = TRUE)
    })


    output$runlrt <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      nrl <- length(levels(colData(values$dds_obj)[,fac1]))

      if(nrl > 2)
        actionButton("button_runlrt",label = "(re)Run LRT for the dataset",class = "btn btn-primary")
    })

    observeEvent(input$button_runlrt,{
      withProgress(message="Computing the LRT results...",
                   detail = "This step can take a little while",
                   value = 0,{

                     values$ddslrt <- DESeq(values$dds_obj,test = "LRT",
                                            full = as.formula(paste0("~",paste(input$choose_lrt_full, collapse=" + "))),
                                            reduced = as.formula(paste0("~",paste(input$choose_lrt_reduced, collapse=" + "))))

                     values$reslrt <- results(values$ddslrt)


                     if(!is.null(values$annotation_obj))
                       values$reslrt$symbol <- values$annotation_obj$gene_name[match(rownames(values$reslrt),
                                                                                     rownames(values$annotation_obj))]
                   })

    })

    # copy this in the report for debugging purposes or so
    # # section title
    #
    # ```{r setup, include=FALSE}
    # knitr::opts_chunk$set(echo = TRUE)
    # ```
    #
    # ## R Markdown
    #
    # This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.
    #
    # When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
    #
    #   ```{r cars}
    # values$reslrt
    # summary(values$reslrt)
    #
    # deseqresult2DEgenes(values$reslrt)
    # plotCounts(dds_airway_lrt,intgroup="cell",gene="ENSG00000262902")
    # plotCounts(dds_airway_lrt,intgroup="cell",gene="ENSG00000123243")
    # resultsNames(dds_airway_lrt)
    # ```
    #
    #
    # ```{r}
    # footertemplate()
    # ```

    ## TODO; think if we want to allow for a continuous factor in the results, if so then do something like building
    # the ui elements accordingly

    output$fac1 <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             "Please select an experimental factor to generate the results"
        )
      )
      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]

      fac1_levels <- levels(fac1_vals)
      if(is.factor(colData(values$dds_obj)[,fac1]))
        selectInput("fac1_c1","Select the name of the numerator level for the fold change",choices = c("",fac1_levels), selected = "")
      # selectInput("fac1_c2","c2",choices = fac1_levels)
    })

    output$fac2 <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]
      fac1_levels <- levels(fac1_vals)
      if(is.factor(colData(values$dds_obj)[,fac1]))
        # selectInput("fac1_c1","c1",choices = fac1_levels)
        selectInput("fac1_c2","Select the name of the denominator level for the fold change (must be different from the numerator)",choices = c("",fac1_levels), selected = "")
    })

    output$facnum <- renderPrint({
      shiny::validate(
        need(input$choose_expfac!="",
             ""
        )
      )
      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]

      # fac1_levels <- levels(fac1_vals)
      if(class(colData(values$dds_obj)[,fac1]) %in% c("integer","numeric"))
        print("numeric/integer factor provided")

      # selectInput("fac1_num","num/int",choices = c("",fac1_levels), selected = "")
      # selectInput("fac1_c2","c2",choices = fac1_levels)
    })

    output$runresults <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             "Select a factor for the contrast first")
      )

      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]

      if(!(class(colData(values$dds_obj)[,fac1]) %in% c("integer","numeric"))){

        shiny::validate(
          need(input$fac1_c1 != "" & input$fac1_c2 != "" & input$fac1_c1 != input$fac1_c2,
               "Select two different levels of the factor for the contrast")
        )
      }

      # if((class(colData(values$dds_obj)[,fac1]) %in% c("integer","numeric"))){
      # 
      #   shiny::validate(
      #     need(input$resu_lfcshrink==FALSE,
      #          "Set the Add the unshrunken MLE to FALSE")
      #   )
      # }
      shiny::validate(
        need("results" %in% mcols(mcols(values$dds_obj))$type ,
             "I couldn't find results. you should first run DESeq() with the button up here"
        )
      )

      # if(input$choose_expfac=="" | input$fac1_c1 == "" | input$fac1_c2 == "" | input$fac1_c1 == input$fac1_c2)
      #   return(NULL)
      # else
      actionButton("button_runresults","Extract the results!", icon = icon("spinner"), class = "btn btn-success")
    })

    observeEvent(input$button_runresults, {
      withProgress(message="Computing the results...",
                   detail = "DE table on its way!",
                   value = 0,{
                     # handling the experimental covariate correctly to extract the results...
                     if(is.factor(colData(values$dds_obj)[,input$choose_expfac])) {
                       if(input$resu_ihw) {
                         values$res_obj <- results(values$dds_obj,
                                                   contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                   independentFiltering = input$resu_indfil, 
                                                   alpha = input$FDR,
                                                   filterFun = ihw)
                         
                         if(input$resu_lfcshrink) {
                           incProgress(amount = 0.15,detail = "Results extracted. Shrinking the logFC now...")
                           values$res_obj <- lfcShrink(values$dds_obj,
                                                       contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                       res = values$res_obj)
                           incProgress(amount = 0.8,detail = "logFC shrunken, adding annotation info...")
                         } else {
                           incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                         }
                       } else {
                         values$res_obj <- results(values$dds_obj,
                                                   contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                   independentFiltering = input$resu_indfil, 
                                                   alpha = input$FDR)
                         if(input$resu_lfcshrink) {
                           incProgress(amount = 0.15,detail = "Results extracted. Shrinking the logFC now...")
                           values$res_obj <- lfcShrink(values$dds_obj,
                                                       contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                       res = values$res_obj)
                           incProgress(amount = 0.8,detail = "logFC shrunken, adding annotation info...")
                         } else {
                           incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                         }
                       }
                     }
                     
                     if(class(colData(values$dds_obj)[,input$choose_expfac]) %in% c("integer","numeric"))
                       values$res_obj <- results(values$dds_obj,name = input$choose_expfac,
                                                 independentFiltering = input$resu_indfil, 
                                                 alpha = input$FDR
                                                 # , addMLE = input$resu_lfcshrink
                                                 )
                     
                     # adding info from the annotation
                     if(!is.null(values$annotation_obj))
                       values$res_obj$symbol <- values$annotation_obj$gene_name[
                         match(rownames(values$res_obj),
                               rownames(values$annotation_obj))]
                   })
    })

    output$diyres_summary <- renderPrint({
      shiny::validate(
        need(input$choose_expfac!="" & input$fac1_c1 != "" & input$fac1_c2 != "" & input$fac1_c1 != input$fac1_c2 ,
             "Please select the factor to build the contrast upon, and two different levels to build the contrast"
        )
      )
      shiny::validate(
        need(!is.null(values$res_obj), "Parameters selected, please compute the results first")
      )
      # summary(results(values$dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2)))
      summary(values$res_obj,alpha = input$FDR)
    })

    output$printdds <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please provide a count matrix/dds object"
        )
      )

      values$dds_obj
      design(values$dds_obj)
    })

    output$printres <- renderPrint({
      shiny::validate(
        need(!is.null(values$res_obj),
             "Please provide a DESeqResults object"
        )
      )
      print(sub(".*p-value: (.*)","\\1",mcols(values$res_obj, use.names=TRUE)["pvalue","description"]))
      summary(values$res_obj,alpha = input$FDR) # use fdr shiny widget
    })


    output$store_result <- renderUI({
      if(is.null(values$res_obj))
        return(NULL)
      actionButton("button_store_result", "Store current results",class = "btn btn-primary")
    })

    observeEvent(input$button_store_result,
                 {
                   values$stored_res <- values$res_obj
                   # this is in such a way to store & compare later if some parameters are edited
                 })

    output$table_res <- DT::renderDataTable({
      if(is.null(values$res_obj))
        return(NULL)
      mydf <- as.data.frame(values$res_obj[order(values$res_obj$padj),])#[1:500,]
      rownames(mydf) <- createLinkENS(rownames(mydf),species = annoSpecies_df$ensembl_db[match(input$speciesSelect,annoSpecies_df$species)]) ## TODO: check what are the species from ensembl and
      ## TODO: add a check to see if wanted?
      mydf$symbol <- createLinkGeneSymbol(mydf$symbol)
      datatable(mydf, escape = FALSE)
    })

    # server resu diagnostics --------------------------------------------------------
    output$pvals_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      p <- ggplot(res_df, aes_string("pvalue")) +
        geom_histogram(binwidth = 0.01, boundary = 0) + theme_bw()
      
      # for visual estimation of the false discovery proportion in the first bin
      alpha <- binw <- input$FDR
      pi0 <- 2*mean(res_df$pvalue > 0.5)
      p <- p + geom_hline(yintercept = pi0 * binw * nrow(res_df), col = "steelblue") + 
        geom_vline(xintercept = alpha, col = "red")
      
      p <- p + ggtitle(
        label = "p-value histogram",
        subtitle = paste0(
          "Expected nulls = ", pi0 * binw * nrow(res_df), 
          " - #elements in the selected bins = ", sum(res_df$pvalue < alpha)
        ))
      
      exportPlots$plot_pvals_hist <- p
      p
      
    })
    
    output$pvals_hist_strat <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      res_df <- mutate(
        res_df, 
        stratum = cut(baseMean, include.lowest = TRUE, 
                      breaks = signif(quantile(baseMean, probs = seq(0,1, length.out = 10)),2)))
  
      p <- ggplot(res_df, aes_string("pvalue")) +
        geom_histogram(binwidth = 0.01, boundary = 0) + 
        facet_wrap(~stratum) + 
        theme_bw()
      
      p <- p + ggtitle(
        label = "p-value histogram",
        subtitle = "stratified on the different value classes of mean expression values")
      
      exportPlots$plot_pvals_hist_strat <- p
      p
    })
    
    output$pvals_ss <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      phi <- input$FDR
      res_df <- mutate(res_df, rank = rank(pvalue))
      m <- nrow(res_df)
      
      p <- ggplot(filter(res_df, rank <= 6000), 
                  aes_string(x = "rank", y = "pvalue")) + 
        geom_line() + 
        geom_abline(slope = phi/m, col = "red") + 
        theme_bw()
      
      p <- p + ggtitle(
        label = "Schweder-Spjotvoll plot",
        subtitle = paste0(
          "Intersection point at rank ", with(arrange(res_df,rank), last(which(pvalue <= phi * rank / m))))
        )
      exportPlots$plot_pvals_ss <- p
      p
    })
    

    output$logfc_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      p <- ggplot(res_df, aes_string("log2FoldChange")) +
        geom_histogram(binwidth = 0.1) + theme_bw()
      
      p <- p + ggtitle(
        "Histogram of the log2 fold changes"
      )

      exportPlots$plot_logfc_hist <- p
      p
    })


    output$dds_design <- renderPrint({
      design(values$dds_obj)
    })

    output$res_names <- renderPrint({
      resultsNames(values$dds_obj)
    })

    output$explore_res <- renderPrint({
      expfac <- attributes(terms.formula(design(values$dds_obj)))$term.labels
      expfac # plus, support up to four factors that are either there or not according to the length
    })

    output$plotma <- renderPlot({
      p <- plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = input$FDR)
      exportPlots$plot_ma <- p
      p
    })

    output$mazoom <- renderPlot({
      if(is.null(input$ma_brush)) return(ggplot() + annotate("text",label="click and drag to zoom in",0,0) + theme_bw())

      if(!is.null(values$annotation_obj))
        p <- plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = input$FDR) +
        coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                        ylim = c(input$ma_brush$ymin,input$ma_brush$ymax)) +
        geom_text(aes_string(label="genename"),size=input$size_genelabels,hjust=0.25, vjust=-0.75)
      else
        p <-  plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = input$FDR) +
        coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                        ylim = c(input$ma_brush$ymin,input$ma_brush$ymax))
      exportPlots$plot_mazoom <- p
      p
    })

    output$ma_highlight <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "Please generate the results object to display the plot and show the combined tables")
      )

      if("symbol" %in% names(values$res_obj)) {
        p <- plot_ma(values$res_obj,
                intgenes = input$avail_symbols,annotation_obj = values$annotation_obj,FDR = input$FDR)
      } else {
        p <- plot_ma(values$res_obj,
                intgenes = input$avail_ids,annotation_obj = values$annotation_obj,FDR = input$FDR)
      }

      exportPlots$plot_mahighlight <- p
      p
    })

    output$ma_hl_list <- renderPlot({
      if(is.null(values$genelist_ma))
        return(NULL)
      if("symbol" %in% names(values$res_obj)) {
        p <- plot_ma(values$res_obj,
                intgenes = values$genelist_ma$`Gene Symbol`,annotation_obj = values$annotation_obj,FDR = input$FDR)
      } else {
        # plot_ma(values$res_obj,
        # intgenes = values$genelist_ma,annotation_obj = values$annotation_obj)
        return(NULL)
      }
      exportPlots$plot_mahllist <- p
      p
    })


    curData <- reactive({
      mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      res <- brushedPoints(mama, input$ma_brush,xvar="logmean",yvar="lfc")
      res
    })


    curDataClick <- reactive({
      mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      res <- nearPoints(mama, input$mazoom_click,threshold = 20, maxpoints = 1,
                        addDist = TRUE)
      res
    })




    output$ma_brush_out <- DT::renderDataTable({
      if(nrow(curData())==0)
        return(NULL)
      datatable(curData(),options=list(pageLength=100))
    })

    output$heatbrush <- renderPlot({
      if((is.null(input$ma_brush))|is.null(values$dds_obj)) return(NULL)

      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]

      if(input$pseudocounts) toplot <- log2(1+toplot)
      mat_rowscale <- function(x)
      {
        m <- apply(x, 1, mean, na.rm = TRUE)
        s <- apply(x, 1, sd, na.rm = TRUE)
        return((x - m)/s)
      }
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
    })


    output$heatbrushD3 <- renderD3heatmap({
      if((is.null(input$ma_brush))|is.null(values$dds_obj)) return(NULL)
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding
      if(input$pseudocounts) toplot <- log2(1+toplot)
      mat_rowscale <- function (x)
      {
        m = apply(x, 1, mean, na.rm = TRUE)
        s = apply(x, 1, sd, na.rm = TRUE)
        return((x - m)/s)
      }
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
    })


    output$deb <- renderPrint({
      # curDataClick()
      selectedGene <- curDataClick()$ID
      #         selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
      #         # plotCounts(dds_cleaner,)
      #         genedata <- plotCounts(dds_cleaner,gene=selectedGene,intgroup = "condition",returnData = T)
      #         genedata
      # str(as.character(selectedGene))
      selectedGene
    })

    output$volcanoplot <- renderPlot({
      p <- plot_volcano(values$res_obj, FDR = input$FDR)
      exportPlots$plot_volcanoplot <- p
      p
    })

    # server genefinder --------------------------------------------------------
    output$genefinder_plot <- renderPlot({

      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )

      if(is.null(input$ma_brush)) return(NULL)

      if(is.null(input$mazoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

      selectedGene <- as.character(curDataClick()$ID)
      selectedGeneSymbol <- values$annotation_obj$gene_name[match(selectedGene,values$annotation_obj$gene_id)]

      p <- ggplotCounts(values$dds_obj, selectedGene, intgroup = input$color_by,annotation_obj=values$annotation_obj)

      if(input$ylimZero_genes)
        p <- p + ylim(0.1, NA)

      exportPlots$plot_genefinder <- p
      p
    })

    output$rentrez_infobox <- renderUI({
      shiny::validate(
        need(
          (nrow(curDataClick()) > 0),
          "Select a gene first to display additional info (retrieved from the NCBI/ENTREZ db website)"
        )
      )
      shiny::validate(
        need(
          (!is.null(values$cur_species)),
          "Select a species first in the Data Setup panel"
        )
      )

      selectedGene <- as.character(curDataClick()$ID)
      selgene_entrez <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                               selectedGene, "ENTREZID", input$idtype)
      fullinfo <- geneinfo(selgene_entrez)

      ## TODO: build up link manually to paste under the info!
      #
      link_pubmed <- paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',
                            selgene_entrez,
                            '" target="_blank" >Click here to see more at NCBI</a>')

      if(fullinfo$summary == "")
        return(HTML(paste0("<b>",fullinfo$name, "</b><br/><br/>",
                           fullinfo$description,"<br/><br/>",
                           link_pubmed
        )))
      else
        return(HTML(paste0("<b>",fullinfo$name, "</b><br/><br/>",
                           fullinfo$description, "<br/><br/>",
                           fullinfo$summary, "<br/><br/>",
                           link_pubmed
        )))
    })


    cur_combires <- reactive({

      if(is.null(values$res_obj))
        return(NULL)

      normCounts <- as.data.frame(counts(estimateSizeFactors(values$dds_obj),normalized=TRUE))
      normCounts$id <- rownames(normCounts)
      res_df <- deseqresult2tbl(values$res_obj)

      combi_obj <- dplyr::inner_join(res_df,normCounts,by="id")
      combi_obj$symbol <- values$annotation_obj$gene_name[match(combi_obj$id,values$annotation_obj$gene_id)]


      if("symbol" %in% names(values$res_obj)) {
        sel_genes <- input$avail_symbols
        sel_genes_ids <- values$annotation_obj$gene_id[match(sel_genes,values$annotation_obj$gene_name)]
      } else {
        sel_genes_ids <- input$avail_ids
      }

      if(length(sel_genes_ids) > 0) {
        combi_obj[match(sel_genes_ids,combi_obj$id),]
      } else {
        combi_obj
      }
    })

    output$table_combi <- DT::renderDataTable({
      datatable(cur_combires(),options = list(scrollX=TRUE))
    })

    cur_combires_list <- reactive({
      if(is.null(values$res_obj))
        return(NULL)

      normCounts <- as.data.frame(counts(estimateSizeFactors(values$dds_obj),normalized=TRUE))
      normCounts$id <- rownames(normCounts)
      res_df <- deseqresult2tbl(values$res_obj)

      combi_obj <- dplyr::inner_join(res_df,normCounts,by="id")
      combi_obj$symbol <- values$annotation_obj$gene_name[match(combi_obj$id,values$annotation_obj$gene_id)]


      if("symbol" %in% names(values$res_obj)) {
        sel_genes <- values$genelist_ma$`Gene Symbol`
        sel_genes_ids <- values$annotation_obj$gene_id[match(sel_genes,values$annotation_obj$gene_name)]
      } else {
        # sel_genes_ids <- values$genelist_ma$`Gene Symbol`
      }

      if(length(sel_genes_ids) > 0) {
        combi_obj[match(sel_genes_ids,combi_obj$id),]
      } else {
        combi_obj
      }
    })

    output$table_combi_list <- DT::renderDataTable({
      if(is.null(values$genelist_ma))
        return(NULL)
      datatable(cur_combires_list(),options = list(scrollX=TRUE))
    })

    output$bp1 <- renderPlot({
      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          (length(input$avail_symbols)>0 | length(input$avail_ids)>0),
          "Select at least a gene to plot"
        )
      )
      if(length(input$avail_symbols)>0) {
        # got the symbol, look for the id
        mysym <- input$avail_symbols[1]
        myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
      } else {
        myid <- input$avail_ids[1]
        # make it optional if annot is available
        if(!is.null(values$annotation_obj)) {
          mysim <- values$annotation_obj$gene_name[match(myid, values$annotation_obj$gene_id)]
        } else {
          mysim <- ""
        }
      }
      p <- ggplotCounts(values$dds_obj, myid, intgroup = input$color_by,annotation_obj=values$annotation_obj)
      if(input$ylimZero_genefinder)
        p <- p + ylim(0.1, NA)
      exportPlots$plotbp1 <- p
      p
    })

    output$bp2 <- renderPlot({
      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          (length(input$avail_symbols)>1 | length(input$avail_ids)>1),
          "Select at least a second gene to plot"
        )
      )
      if(length(input$avail_symbols)>0) {
        # got the symbol, look for the id
        mysym <- input$avail_symbols[2]
        myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
      } else {
        myid <- input$avail_ids[2]
        # make it optional if annot is available
        if(!is.null(values$annotation_obj)) {
          mysim <- values$annotation_obj$gene_name[match(myid, values$annotation_obj$gene_id)]
        } else {
          mysim <- ""
        }
      }
      p <- ggplotCounts(values$dds_obj, myid, intgroup = input$color_by,annotation_obj=values$annotation_obj)
      if(input$ylimZero_genefinder)
        p <- p + ylim(0.1, NA)
      exportPlots$plotbp2 <- p
      p
    })

    output$bp3 <- renderPlot({
      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          (length(input$avail_symbols)>2 | length(input$avail_ids)>2),
          "Select at least a third gene to plot"
        )
      )
      if(length(input$avail_symbols)>0) {
        # got the symbol, look for the id
        mysym <- input$avail_symbols[3]
        myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
      } else {
        myid <- input$avail_ids[3]
        # make it optional if annot is available
        if(!is.null(values$annotation_obj)) {
          mysim <- values$annotation_obj$gene_name[match(myid, values$annotation_obj$gene_id)]
        } else {
          mysim <- ""
        }
      }
      p <- ggplotCounts(values$dds_obj, myid, intgroup = input$color_by,annotation_obj=values$annotation_obj)
      if(input$ylimZero_genefinder)
        p <- p + ylim(0.1, NA)
      exportPlots$plotbp3 <- p
      p
    })

    output$bp4 <- renderPlot({
      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          (length(input$avail_symbols)>3 | length(input$avail_ids)>3),
          "Select at least a fourth gene to plot"
        )
      )
      if(length(input$avail_symbols)>0) {
        # got the symbol, look for the id
        mysym <- input$avail_symbols[4]
        myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
      } else {
        myid <- input$avail_ids[4]
        # make it optional if annot is available
        if(!is.null(values$annotation_obj)) {
          mysim <- values$annotation_obj$gene_name[match(myid, values$annotation_obj$gene_id)]
        } else {
          mysim <- ""
        }
      }
      p <- ggplotCounts(values$dds_obj, myid, intgroup = input$color_by,annotation_obj=values$annotation_obj)
      if(input$ylimZero_genefinder)
        p <- p + ylim(0.1, NA)
      exportPlots$plotbp4 <- p
      p
    })


    # server report editor --------------------------------------------------------
    ### yaml generation
    rmd_yaml <- reactive({
      paste0("---",
             "\ntitle: '", input$report_title,
             "'\nauthor: '", input$report_author,
             "'\ndate: '", Sys.Date(),
             "'\noutput:\n  html_document:\n    toc: ", input$report_toc, "\n    number_sections: ", input$report_ns, "\n    theme: ", input$report_theme, "\n---\n\n",collapse = "\n")
    })


    # rmd_full <- reactive({
    #   paste0(rmd_yaml(),"\n",
    #          readLines("reportTemplate.Rmd"))
    # })
    # output$loadedRmd <- renderPrint({
    #   # rmd_yaml() # or rmd_full()
    #   paste0(
    #     # rmd_yaml(),
    #     paste0(readLines("reportTemplate.Rmd"),collapse = "\n"))
    #   # head(paste0(rmd_yaml(),
    #   # readLines("reportTemplate.Rmd")),collapse="\n")
    # })

    ### loading report template
    # update aceEditor module
    observe({
      # loading rmd report from disk
      inFile <- system.file("extdata", "irt.Rmd",package = "ideal")

      isolate({
        if(!is.null(inFile) && !is.na(inFile)) {

          rmdfilecontent <- paste0(readLines(inFile),collapse="\n")

          shinyAce::updateAceEditor(session, "acereport_rmd", value = rmdfilecontent)
        }
      })
    })


    ### ace editor options
    observe({
      autoComplete <- if(input$enableAutocomplete) {
        if(input$enableLiveCompletion) "live" else "enabled"
      } else {
        "disabled"
      }

      updateAceEditor(session, "acereport_rmd", autoComplete = autoComplete,theme=input$theme, mode=input$mode)
      # updateAceEditor(session, "plot", autoComplete = autoComplete)
    })

    #Enable/Disable R code completion
    rmdOb <- aceAutocomplete("acereport_rmd")
    observe({
      if(input$enableRCompletion) {
        rmdOb$resume()
      } else {
        rmdOb$suspend()
      }
    })

    ## currently not working as I want with rmarkdown::render, but can leave it like this - the yaml will be taken in the final version only
    output$knitDoc <- renderUI({
      input$updatepreview_button


      ## TODO: this does what it should do but messes up with CSS and so
      #
      #     # error_I <- 0
      #     withProgress(message = 'Processing', value = 0, {
      #       isolate({
      #         fileConn<-file("www/tmp.Rmd")
      #         tmp_content <-
      #           paste0(rmd_yaml(),
      #                  input$acereport_rmd,collapse = "\n")
      #         writeLines(tmp_content, fileConn)
      #         close(fileConn)
      #         incProgress(0.5, detail = "Synthesizing report...")
      #         # tryCatch({
      #           rmarkdown::render(input = "www/tmp.Rmd", output_format = "html_document", output_file = "../www/Rmd_preview.html",quiet = TRUE) #},
      #           # error = function(e) {
      #           #   # error_I <<- 1
      #           # }
      #         # )
      #       })
      #       setProgress(1)
      #     })
      #
      #     return(isolate(includeHTML("www/Rmd_preview.html")))
      #     # return(isolate(includeHTML("<iframe src='www/Rmd_preview.html', width='100%', height='800'></iframe>")))
      #     # return(isolate(HTML("<iframe src='www/Rmd_preview.html', width='100%', height='800'></iframe>")))


      return(
        withProgress({
          # temporarily switch to the temp dir, in case you do not have write
          # permission to the current working directory
          owd <- setwd(tempdir())
          on.exit(setwd(owd))
          tmp_content <- paste0(rmd_yaml(),input$acereport_rmd,collapse = "\n")
          isolate(HTML(knit2html(text = tmp_content, fragment.only = TRUE, quiet = TRUE)))
        },
        message = "Updating the report in the app body",
        detail = "This can take some time"
        )
      )
    })

    # Generate and Download module
    output$saveRmd <- downloadHandler(
      filename = function() {
        if(input$rmd_dl_format == "rmd") {
          "report.Rmd"
        } else {
          "report.html" # TODO: maybe add Sys.time() to the filename to improve traceability?
        }
      },
      content = function(file) {
        
        # knit2html(text = input$rmd, fragment.only = TRUE, quiet = TRUE))
        
        tmp_content <-
          paste0(rmd_yaml(),
                 input$acereport_rmd,collapse = "\n")
        if(input$rmd_dl_format == "rmd") {
          cat(tmp_content,file=file,sep="\n")
        } else {
          if(input$rmd_dl_format == "html") {
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            cat(tmp_content,file="ideal_tempreport.Rmd",sep="\n")
            withProgress(rmarkdown::render(input = "ideal_tempreport.Rmd",
                                           output_file = file,
                                           # fragment.only = TRUE,
                                           quiet = TRUE),
                         message = "Generating the html report",
                         detail = "This can take some time")
          }
        }
      }
    )
    
    output$ui_iSEEexport <- renderUI({
      validate(
        need(((!is.null(values$dds_obj)) & (!is.null(values$res_obj))),
             message = "Please build and compute the dds and res object to export as 
             SummarizedExperiment for use in iSEE")
      )
      return(
        tagList(
          textInput(
            "se_export_name",label = "Choose a filename for the serialized .rds object",
            value = "se_ideal_toiSEE.rds"),
          downloadButton(
            "button_iSEEexport",
            label = "Export as serialized SummarizedExperiment"
          )
        )
      )
    })

    output$button_iSEEexport <- downloadHandler(
      filename = function() {
        # paste0("se_ideal_toiSEE_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".rds")
        input$se_export_name
      }, content = function(file) {
        se <- wrapup_for_iSEE(values$dds_obj, values$res_obj)
        saveRDS(se, file = file)
      }
    )

    # server state saving --------------------------------------------------------
    ### to environment
    observe({
      if(is.null(input$task_exit_and_save) || input$task_exit_and_save ==0 ) return()

      # quit R, unless you are running an interactive session
      if(interactive()) {
        # flush input and values to the environment in two distinct objects (to be reused later?)
        isolate({

          # ideal_env <<- new.env(parent = emptyenv())
          cur_inputs <- reactiveValuesToList(input)
          cur_values <- reactiveValuesToList(values)
          tstamp <- gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))

          # myvar <- "frfr"
          # assign("test", myvar, ideal_env)

          # better practice rather than assigning to global env - notify users of this
          assign(paste0("ideal_inputs_", tstamp),cur_inputs, envir = ideal_env)
          assign(paste0("ideal_values_", tstamp),cur_values, envir = ideal_env)
          stopApp("ideal closed, state successfully saved to global R environment.")

          # assign(paste0("ideal_inputs_",
          #               gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
          #        reactiveValuesToList(input), envir = .GlobalEnv)
          # assign(paste0("ideal_values_",
          #               gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
          #        reactiveValuesToList(values), envir = .GlobalEnv)
          # stopApp("ideal closed, state successfully saved to global R environment.")
        })
      } else {
        stopApp("ideal closed")
        q("no")
      }
    })

    ### to binary data
    saveState <- function(filename) {
      isolate({
        LiveInputs <- reactiveValuesToList(input)
        # values[names(LiveInputs)] <- LiveInputs
        r_data <- reactiveValuesToList(values)
        save(LiveInputs, r_data , file = filename)
      })
    }

    output$task_state_save <- downloadHandler(
      filename = function() {
        paste0("idealState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".RData")
      },
      content = function(file) {
        saveState(file)
      }
    )

    output$sessioninfo <- renderPrint({
      sessionInfo()
    })

    # server export plots and tables --------------------------------------------------------
    
    ## here, all export of plots and tables
    output$download_plot_pvals_hist <- downloadHandler(filename = function() {
      input$filename_plot_pvals_hist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_pvals_hist, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_logfc_hist <- downloadHandler(filename = function() {
      input$filename_plot_logfc_hist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_logfc_hist, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_ma <- downloadHandler(filename = function() {
      input$filename_plot_ma
    }, content = function(file) {
      ggsave(file, exportPlots$plot_ma, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_mazoom <- downloadHandler(filename = function() {
      input$filename_plot_mazoom
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mazoom, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_mahighlight <- downloadHandler(filename = function() {
      input$filename_plot_mahighlight
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mahighlight, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_mahllist <- downloadHandler(filename = function() {
      input$filename_plot_mahllist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mahllist, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_volcanoplot <- downloadHandler(filename = function() {
      input$filename_plot_volcanoplot
    }, content = function(file) {
      ggsave(file, exportPlots$plot_volcanoplot, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plot_genefinder <- downloadHandler(filename = function() {
      input$filename_plot_genefinder
    }, content = function(file) {
      ggsave(file, exportPlots$plot_genefinder, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plotbp1 <- downloadHandler(filename = function() {
      input$filename_plotbp1
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp1, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plotbp2 <- downloadHandler(filename = function() {
      input$filename_plotbp2
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp2, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plotbp3 <- downloadHandler(filename = function() {
      input$filename_plotbp3
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp3, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    output$download_plotbp4 <- downloadHandler(filename = function() {
      input$filename_plotbp4
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp4, width = input$export_width,
             height = input$export_height, units = "cm")
    })

    # tbls
    output$downloadTblResu <- downloadHandler(
      filename = function() {
        "table_results.csv"
      },
      content = function(file) {
        mydf <- as.data.frame(values$res_obj[order(values$res_obj$padj),])
        write.csv(mydf, file)
      }
    )
    output$downloadTblMabrush <- downloadHandler(
      filename = function() {
        "table_mabrush.csv"
      },
      content = function(file) {
        write.csv(curData(), file)
      }
    )
    output$downloadTblCombi <- downloadHandler(
      filename = function() {
        "table_combi.csv"
      },
      content = function(file) {
        write.csv(cur_combires(), file)
      }
    )
    output$downloadTblCombiList <- downloadHandler(
      filename = function() {
        "table_combilist.csv"
      },
      content = function(file) {
        write.csv(cur_combires_list(), file)
      }
    )

    # base graphics plots
    output$download_plot_heatbrush <- downloadHandler(filename = function() {
      input$filename_plot_heatbrush
    }, content = function(file) {
      pdf(file)
      brushedObject <- curData()

      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]

      if(input$pseudocounts) toplot <- log2(1+toplot)

      mat_rowscale <- function(x)
      {
        m <- apply(x, 1, mean, na.rm = TRUE)
        s <- apply(x, 1, sd, na.rm = TRUE)
        return((x - m)/s)
      }

      if(input$rowscale) toplot <- mat_rowscale(toplot)

      pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
      dev.off()
    })

    output$download_plot_vennlists <- downloadHandler(filename = function() {
      input$filename_plot_vennlists
    }, content = function(file) {
      pdf(file)
      gplots::venn(gll())
      dev.off()
    })

    output$download_plot_upsetlists <- downloadHandler(filename = function() {
      input$filename_plot_upsetlists
    }, content = function(file) {
      pdf(file)
      UpSetR::upset(fromList(gll()))
      dev.off()
    })

    ## GO tbls topGO
    output$downloadGOTbl_up <- downloadHandler(
      filename = function() {
        "table_GOresults_up.csv"
      },
      content = function(file) {
        write.csv(values$topgo_up, file)
      }
    )
    output$downloadGOTbl_down <- downloadHandler(
      filename = function() {
        "table_GOresults_down.csv"
      },
      content = function(file) {
        write.csv(values$topgo_down, file)
      }
    )
    output$downloadGOTbl_updown <- downloadHandler(
      filename = function() {
        "table_GOresults_updown.csv"
      },
      content = function(file) {
        write.csv(values$topgo_updown, file)
      }
    )
    output$downloadGOTbl_l1 <- downloadHandler(
      filename = function() {
        "table_GOresults_list1.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list1, file)
      }
    )
    output$downloadGOTbl_l2 <- downloadHandler(
      filename = function() {
        "table_GOresults_list2.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list2, file)
      }
    )
  }) # end of server function definition

  # launch the app!
  shinyApp(ui = ideal_ui, server = ideal_server)
}
