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
#'
#' ideal()
#' ideal(dds)
#' ideal(dds_airway)
#'
#' dds_airway <- DESeq(dds_airway)
#' res_airway <- results(dds_airway)
#' ideal(dds_airway, res_airway)
#' }
#'
ideal<- function(dds_obj = NULL,
                 res_obj = NULL,
                 annotation_obj = NULL,
                 countmatrix = NULL,
                 expdesign = NULL){

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }


  # create environment for storing inputs and values
  ## i need the assignment like this to export it up one level - i.e. "globally"
  ideal_env <<- new.env(parent = emptyenv())

  ## upload max 300mb files - can be changed if necessary
  options(shiny.maxRequestSize=300*1024^2)
  #
  options(shiny.launch.browser = TRUE)

  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  # source(file.path("R/ideal_ui.R"),local = TRUE)
  # source(file.path("R/ideal_server.R"), local= TRUE)



#
#   # components defined in separated .R files
#   shinyApp(ui = ideal_ui, server = ideal_server)
#



  # runApp(appDir = system.file("shiny", package = "ideal"))



  ideal_ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = paste0("ideal - Interactive Differential Expression AnaLysis ",
                     packageVersion("ideal")),
      titleWidth = 900,

      # TODO:
      # http://stackoverflow.com/questions/31440564/adding-a-company-logo-to-shinydashboard-header
      # replace text with image
      # ideal_header$children[[2]]$children <- tags$a(href='https://github.com/federicomarini/ideal',
      # tags$img(src='ideal_logo_v2.png',height='50',width='200'))
      # title = tags$a(href='https://github.com/federicomarini/ideal',
      #                tags$img(src='ideal_logo_v2.png',height='50',width='200')),

      # task menu for saving state to environment or binary data
      shinydashboard::dropdownMenu(type = "tasks",icon = icon("cog"),badgeStatus = NULL, # something to change the top message? maybe file an issue @shinydashboard development
                                   notificationItem(
                                     text = actionButton("task_exit_and_save","Exit ideal & save",
                                                         class = "btn_no_border",
                                                         onclick = "setTimeout(function(){window.close();}, 100); "),
                                     icon = icon("sign-out"),status = "primary"),
                                   menuItem(
                                     text = downloadButton("task_state_save","Save State as .RData"))
      )
    ),


    dashboardSidebar(
      width = 280,
      menuItem("App settings",icon = icon("cogs"),
               uiOutput("color_by"),
               uiOutput("available_genes"),

               #
               #              selectizeInput(
               #                inputId = 'available_genes', label = 'Select Something',
               #                choices = NULL,
               #                multiple = TRUE,
               #                selected = 1
               #              ),

               numericInput("FDR","False Discovery Rate",value = 0.05, min = 0, max = 1, step = 0.01)

      ),
      menuItem("Plot export settings", icon = icon("paint-brush")),
      menuItem("Quick viewer", icon = icon("flash"),
               fluidRow(
                 fluidRow(column(6,p("Count matrix")), column(6,uiOutput("ok_cm"))),
                 fluidRow(column(6,p("Experimental design")), column(6,uiOutput("ok_ed"))),
                 fluidRow(column(6,p("DESeqDataset")), column(6,uiOutput("ok_dds"))),
                 fluidRow(column(6,p("Annotation")), column(6,uiOutput("ok_anno"))),
                 fluidRow(column(6,p("Results")), column(6,uiOutput("ok_resu")))
               ))
    ),



    dashboardBody(
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
}
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


          tabPanel(
            # "Welcome!",  icon = icon("info-circle"),
            "Welcome!",  icon = icon("home"),
            # carouselPanel(
            #   img(src = "www/ideal_logo_v2.png"),
            #   img(src = "ideal_logo_v2.png"),
            #   img(src = "ideal_logo_v2.png")
            # ),
            includeMarkdown(system.file("extdata", "welcome.md",package = "ideal")),

            ## TODO: explore possibility to put a carousel of images: https://github.com/dcurrier/carouselPanel/
            # carouselPanel(
            #   plotOutput("distPlot1"),
            #   plotOutput("distPlot2")
            # ),
            # img(src = "ideal_logo_v2.png"),
            includeMarkdown(system.file("extdata", "instructions.md",package = "ideal"))
          ),


          tabPanel(
            "Data Setup",icon = icon("upload"),

            box(width = 12, title = "Step 1", status = "danger", solidHeader = TRUE,
                h2("Upload your count matrix and the info on the experimental design"),

                fluidRow(
                  column(
                    width = 4,
                    uiOutput("upload_count_matrix"),
                    uiOutput("upload_metadata")
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
            # h2("Step 1: Upload your count matrix and the info on the experimental design"),


            uiOutput("ui_step2"),

            # hr(),
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
            # hr(),

            uiOutput("ui_step3")
          ),




          tabPanel(
            "Counts Overview",
            icon = icon("eye"),
            conditionalPanel(
              condition="!output.checkdds",
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
              selectInput("corr_method","Correlation method palette",choices = list("pearson","spearman")),
              p("Compute sample to sample correlations on the normalized counts - warning, it can take a while to plot all points (depending mostly on the number of samples you provided)."),
              actionButton("compute_pairwisecorr", "Run", class = "btn btn-primary"),
              uiOutput("pairwise_plotUI"),
              uiOutput("heatcorr_plotUI")

            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )
          ),



          tabPanel(
            "Extract Results", icon = icon("table"),

            # see: http://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value?noredirect=1&lq=1
            conditionalPanel(
              condition="!output.checkdds",

              # conditionalPanel(
              #   condition="output.checkresu==0",
              #   h2('RESU not provided')
              # ),

              # "Data Overview", icon = icon("eye"),

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
                    width = 4,
                    uiOutput("fac1"),
                    uiOutput("fac2")
                  ),
                  # continuous covariate
                  uiOutput("facnum")
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
                  wellPanel(
                    selectInput("resu_indfil",label = "Apply independent filtering automatically",
                                choices = c(TRUE,FALSE), selected = TRUE),
                    selectInput("resu_addmle",label = "Add the unshrunken MLE of log2 fold change",
                                choices = c(TRUE,FALSE), selected = TRUE),
                    selectInput("resu_ihw", "Use Independent Hypothesis Weighting as a filtering function",
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
              fluidRow(
                column(
                  width = 6,
                  plotOutput("pvals_hist")
                ),
                column(
                  width = 6,
                  plotOutput("logfc_hist")
                )
              )
            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )

          ),

          tabPanel(
            "Gene Lists", icon = icon("list-alt"),
            conditionalPanel(
              condition="!output.checkresu",
              h2("Gene Set Enrichment on the lists"),

              selectInput("go_cats",label = "Select the GO category(ies) of interest",
                          choices = list("GO Biological Process" = "BP", "GO Molecular Function" = "MF", "GO Cellular Component" = "CC"),
                          selected = "BP",multiple = TRUE
              ),


              tabBox(
                width = NULL,
                id="gse_tabbox",
                tabPanel("UPregu", icon = icon("arrow-circle-up"),
                         fluidRow(column(width = 6,actionButton("button_enrUP", "Perform gene set enrichment analysis on the upregulated genes",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrUP_goseq", "Perform gene set enrichment analysis on the upregulated genes - goseq",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrUP_topgo", "Perform gene set enrichment analysis on the upregulated genes - topGO",class = "btn btn-primary"))),
                         DT::dataTableOutput("DT_gse_up"),
                         DT::dataTableOutput("DT_gse_up_goseq"),
                         DT::dataTableOutput("DT_gse_up_topgo")
                ),
                tabPanel("DOWNregu", icon = icon("arrow-circle-down"),
                         fluidRow(column(width = 6,actionButton("button_enrDOWN", "Perform gene set enrichment analysis on the downregulated genes",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrDOWN_goseq", "Perform gene set enrichment analysis on the downregulated genes - goseq",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrDOWN_topgo", "Perform gene set enrichment analysis on the downregulated genes - topGO",class = "btn btn-primary"))),
                         DT::dataTableOutput("DT_gse_down"),
                         DT::dataTableOutput("DT_gse_down_goseq"),
                         DT::dataTableOutput("DT_gse_down_topgo")
                ),
                tabPanel("UPDOWN", icon = icon("arrows-v"),
                         fluidRow(column(width = 6,actionButton("button_enrUPDOWN", "Perform gene set enrichment analysis on the up- and downregulated genes",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrUPDOWN_goseq", "Perform gene set enrichment analysis on the up- and downregulated genes - goseq",class = "btn btn-primary"))),
                         fluidRow(column(width = 6,actionButton("button_enrUPDOWN_topgo", "Perform gene set enrichment analysis on the up- and downregulated genes - topGO",class = "btn btn-primary"))),
                         DT::dataTableOutput("DT_gse_updown"),
                         DT::dataTableOutput("DT_gse_updown_goseq"),
                         DT::dataTableOutput("DT_gse_updown_topgo")
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
                         DT::dataTableOutput("DT_gse_list1_topgo")
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
                         DT::dataTableOutput("DT_gse_list2_topgo")
                )
              ),

              verbatimTextOutput("goterm_heatmap_up_topgo"),





              ## will put collapsible list elements? or multi tab panel? or something to select on the left, and operate output-wise on the right e.g. venn diagrams or table for gene set enrichment
              # h3("custom list 3 - handpicked") # use the select input from the left column?
              # ,verbatimTextOutput("debuggls"),

              # verbatimTextOutput("printUPgenes"),
              # verbatimTextOutput("debuglists"),

              h2("Intersection of gene sets"),
              checkboxInput("toggle_updown","Use up and down regulated genes", TRUE),
              checkboxInput("toggle_up","Use up regulated genes", FALSE),
              checkboxInput("toggle_down","Use down regulated genes", FALSE),
              checkboxInput("toggle_list1","Use list1 genes", TRUE),
              checkboxInput("toggle_list2","Use list2 genes", FALSE),
              checkboxInput("toggle_list3","Use list3 genes", FALSE),


              plotOutput("vennlists"),
              plotOutput("upsetLists")
            ),
            conditionalPanel(
              condition="output.checkresu",
              h2("You did not create the result object yet. Please go the dedicated tab and generate it")
            )
          ),

          tabPanel(
            "MA Plot", icon = icon("photo"),
            conditionalPanel(
              condition="!output.checkresu",

              headerPanel("MA plot interactive exploration"),

              fluidRow(column(6,
                              h4("MA plot - Interactive!"),
                              plotOutput('plotma', brush = 'ma_brush')),
                       column(6,
                              h4("Zoomed section"),
                              plotOutput("mazoom",click= 'mazoom_click'))
                       # ,
                       # column(4,
                       #        h4("Boxplot for the selected gene"),
                       #        plotOutput("geneplot")
                       # )
              ),
              fluidRow(column(6,
                              h4("Selected gene"),
                              checkboxInput("ylimZero_genes","Set y axis limit to 0",value=TRUE),
                              plotOutput("genefinder_plot")
              ),
              column(6,
                     h4("Gene infobox"),
                     htmlOutput("rentrez_infobox"))

              ),


              fluidRow(column(6,
                              h4("volcano plot"),
                              plotOutput("volcanoplot")
              )),


              # plotOutput("volcanoplot"),
              fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
              fluidRow(
                column(4,
                       checkboxInput("rowscale",label = "Scale by rows",value = TRUE)),
                column(4,
                       checkboxInput("pseudocounts","use log2(1+counts)",value = TRUE))
              ),
              fluidRow(
                column(6,
                       plotOutput("heatbrush")
                ),
                column(6,
                       d3heatmapOutput("heatbrushD3"))
              ),


              box(
                title = "Brushed table", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, collapsed = TRUE, width = 12,
                fluidRow(DT::dataTableOutput("ma_brush_out")))
            ),
            conditionalPanel(
              condition="output.checkresu",
              h2("You did not create the result object yet. Please go the dedicated tab and generate it")
            )
          ),
          tabPanel(
            "Gene Finder", icon = icon("crosshairs"),
            conditionalPanel(
              condition="!output.checkdds",
              fluidRow(
                column(6,checkboxInput("ylimZero_genefinder","Set y axis limit to 0",value=TRUE))),
              fluidRow(
                column(6,
                       plotOutput("bp1")
                ),
                column(6,
                       plotOutput("bp2"))
              ),
              fluidRow(
                column(6,
                       plotOutput("bp3")
                ),
                column(6,
                       plotOutput("bp4"))
              ),

              plotOutput("ma_highlight"),
              DT::dataTableOutput("table_combi")

            ),
            conditionalPanel(
              condition="output.checkdds",
              h2("You did not create the dds object yet. Please go the main tab and generate it")
            )
          ),

          tabPanel(
            "Report Editor",
            icon = icon("pencil"),


            fluidRow(
              column(
                width = 6,
                box(
                  title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
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
              # ,
              # column( # kept for debugging purposes!
              #   width = 6,
              #   verbatimTextOutput("loadedRmd")
              # )
            ),
            fluidRow(
              column(3,
                     actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
              ),
              column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success"))
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
          ),
          tabPanel(
            "About", icon = icon("institution"),
            includeMarkdown(system.file("extdata", "about.md",package = "ideal")),
            hr(),

            h4("Session Info"),
            verbatimTextOutput("sessioninfo")
          )
          ,tabPanel(
            "devel", icon = icon("github")
            # ,
            # verbatimTextOutput("debugihw"),
            #
            # plotOutput("ihwp1"),
            # plotOutput("ihwp2"),
            # plotOutput("ihwp3"),
            # plotOutput("ihwp4"),
            # plotOutput("ihwp5"),
            # plotOutput("ihwp6"),
            # plotOutput("ihwp7"),
            # plotOutput("ihwp8")

          )

        )
      )
      ,footer() # TODO: decide where to place the footer
        ),


    skin="blue"
  )



  ideal_server <- shinyServer(function(input, output, session) {

    ## Update directory
    userdir <- tempfile()
    dir.create(userdir, recursive = TRUE)
    # sapply(file.path(newuserdir, dir(newuserdir)[grep("code_", dir(newuserdir))]), file.remove)
    # file.copy(file.path(userdir, "code_All.R"), newuserdir)
    # userdir <- newuserdir
    # dir.create(file.path(userdir, "data"))



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


    # this part sets the "matching" objects if something is provided that is depending on these
    if(!is.null(dds_obj)){
      values$countmatrix <- counts(dds_obj, normalized = FALSE)
      values$expdesign <- as.data.frame(colData(dds_obj))
    }

    # info boxes, to keep on top of the page  on the left side?

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

      # "", paste0(25 + input$count, "%"), icon = icon("list"),
      # color = "purple"
      # )
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



    ###### uploading data
    ## count matrix
    output$upload_count_matrix <- renderUI({
      if (!is.null(dds_obj) | !is.null(countmatrix)) {
        NULL
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
        NULL
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





    output$ddsdesign <- renderUI({
      if(is.null(values$expdesign))
        return(NULL)
      poss_covars <- colnames(values$expdesign)
      selectInput('dds_design', label = 'Select the design for your experiment: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })



    output$ui_step2 <- renderUI({
      if (is.null(values$expdesign) | is.null(values$countmatrix))
        return(NULL)
      box(width = 12, title = "Step 2", status = "warning", solidHeader = TRUE,
          tagList(
            # as in https://groups.google.com/forum/#!topic/shiny-discuss/qQ8yICfvDu0
            h2("Step 2: Select the DE design and create the DESeqDataSet object"),
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
            h2("Optional Step: Create the annotation data frame for your dataset"),

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

      box(width = 12, title = "Optional Step", status = "info", solidHeader = TRUE,
          tagList(
            h2("Optional Step: Remove sample(s) from the current dataset - suspected outliers!"),

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
            h2("Step 3: Run DESeq!"),

            fluidRow(
              column(
                width = 4,
                uiOutput("ui_nrcores")
              )
            ),

            uiOutput("rundeseq"),

            verbatimTextOutput("printDIYresults"),

            # verbatimTextOutput("printdds"),
            # verbatimTextOutput("printres"),
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
    # > fff <- c("cell","dex")
    # > fff
    # [1] "cell" "dex"
    # > DESeqDataSetFromMatrix(countData = ccmm,
    #                          +                        colData = eedd,
    #                          +                        design=~1) -> mydds
    # > DESeqDataSetFromMatrix(countData = ccmm,
    #                          +                        colData = eedd,
    #                          +                        design=~1) -> mydds ; design(mydds)
    # ~1
    # > DESeqDataSetFromMatrix(countData = ccmm,
    #                          +                        colData = eedd,
    #                          +                        design= paste0("~",paste(fff, collapse=" + ")) ) -> mydds ; design(mydds)
    # Error: $ operator is invalid for atomic vectors
    # > paste0("~",paste(fff, collapse=" + "))
    # [1] "~cell + dex"
    # > as.formula(paste0("~",paste(fff, collapse=" + ")))
    # ~cell + dex
    # > DESeqDataSetFromMatrix(countData = ccmm,
    #                          +                        colData = eedd,
    #                          +                        design= as.formula(paste0("~",paste(fff, collapse=" + "))) ) -> mydds ; design(mydds)
    # ~cell + dex
    # >



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
                   # if(!is.null(values$expdesign)){
                   #   withProgress(message="Computing the objects...",value = 0,{
                   #
                   #     values$dds_object <- DESeqDataSetFromMatrix(countData = values$countmatrix,
                   #                                                 colData = values$expdesign,
                   #                                                 design=~1)})
                   # }
                 })

    observeEvent(input$uploadmetadatafile,
                 {
                   values$expdesign <- readMetadata()
                   # if(!is.null(values$countmatrix)){
                   #   withProgress(message="Computing the objects...",value = 0,{
                   #
                   #     values$dds_object <- DESeqDataSetFromMatrix(countData = values$countmatrix,
                   #                                                 colData = values$expdesign,
                   #                                                 design=~1)})
                   # }
                 })



    # for retrieving the annotation
    annoSpecies_df <- data.frame(species=c("","Anopheles","Arabidopsis","Bovine","Worm",
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
      selectInput("idtype", "select the id type in your data", choices=c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL"))
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
             paste0("The package ",annopkg, " is not installed/available. Try installing it with biocLite('",annopkg,"')"))
      )

      retmsg <- paste0(annopkg," - package available and loaded")
      # if (!require(annopkg,character.only=TRUE)) {
      # stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
      # }
      retmsg <- paste0(retmsg," - ",gsub(".eg.db","",gsub("org.","",annopkg)))
      retmsg

    })


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


                                  if(input$nrcores == 1)
                                    values$dds_obj <- DESeq(values$dds_obj)
                                  else
                                    # leave open option for computing in parallel?
                                    values$dds_obj <- DESeq(values$dds_obj,
                                                            parallel = TRUE,
                                                            BPPARAM = MulticoreParam(workers = input$nrcores))


                                  updateSelectInput(session, inputId = "speciesSelect", selected = curr_species)
                                })
                 })

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



    ###### counts overview

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
          pair_corr(current_countmat(),method=input$corr_method),
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












    #### MANAGING THE GENE LISTS
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
      annopkg <- paste0("org.",organism,".eg.db")
      listGenesEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUP(),
                                                            column="ENTREZID", keytype=inputType))
      listBackgroundEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                 column="ENTREZID", keytype="ENSEMBL"))

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
                     annopkg <- paste0("org.",organism,".eg.db")
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUP(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype="ENSEMBL"))
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
                                             keytype="ENSEMBL",
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
                                          keytype="ENSEMBL",
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
                     annopkg <- paste0("org.",organism,".eg.db")
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistDOWN(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype="ENSEMBL"))
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
                                             keytype="ENSEMBL",
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
                                          keytype="ENSEMBL",
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
                     annopkg <- paste0("org.",organism,".eg.db")
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
                     }
                     listGenesEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUPDOWN(),
                                                                            column="ENTREZID", keytype=inputType))
                     listBackgroundEntrez <-  as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                                 column="ENTREZID", keytype="ENSEMBL"))
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
                                             keytype="ENSEMBL",
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
                                          keytype="ENSEMBL",
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
                     annopkg <- paste0("org.",organism,".eg.db")
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
                     }
                     listGenesEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = as.character(values$genelist1$`Gene Symbol`),
                                                              column="ENTREZID", keytype=inputType)
                     listBackgroundEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                   column="ENTREZID", keytype="ENSEMBL")
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
                                             keytype="ENSEMBL",
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
                                          keytype="ENSEMBL",
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
                     annopkg <- paste0("org.",organism,".eg.db")
                     if (!require(annopkg,character.only=TRUE)) {
                       stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
                     }
                     listGenesEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = as.character(values$genelist2$`Gene Symbol`),
                                                              column="ENTREZID", keytype=inputType)
                     listBackgroundEntrez <- AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                   column="ENTREZID", keytype="ENSEMBL")
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
                                             keytype="ENSEMBL",
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
                                          keytype="ENSEMBL",
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
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_down_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_down))
        return(NULL)
      mytbl <- values$topgo_down
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_updown_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_updown))
        return(NULL)
      mytbl <- values$topgo_updown
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list1_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_list1))
        return(NULL)
      mytbl <- values$topgo_list1
      # mytbl$GOid <- rownames(mytbl)
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      datatable(mytbl,escape=FALSE)
    })
    output$DT_gse_list2_topgo <- DT::renderDataTable({
      # if not null...
      if(is.null(values$topgo_list2))
        return(NULL)
      mytbl <- values$topgo_list2
      # mytbl$GOid <- rownames(mytbl)
      mytbl$GO.ID <- createLinkGO(mytbl$GO.ID)
      datatable(mytbl,escape=FALSE)
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




    output$goterm_heatmap_up_topgo <- renderPrint({

      s <- input$DT_gse_up_topgo_rows_selected
      if(length(s) == 0)
        return(NULL)

      print(s)

      print(values$topgo_up[s,])


      values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$genes

      # allow only one selected line
      mygenes <- values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$genes[1]
      myterm <- paste0(
        values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$`GO.ID`, " - ",
        values$topgo_up[input$DT_gse_up_topgo_rows_selected,]$Term)

      genevec <- unlist(strsplit(mygenes,split=","))
      genevec

      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]

      # genevec_ids <- mapIds(org.Hs.eg.db,genevec,"ENSEMBL","SYMBOL",multiVals="first")
      genevec_ids <- mapIds(eval(parse(text=annopkg)),genevec,"ENSEMBL","SYMBOL",multiVals="first")

      genevec_ids



      log2things <- assay(normTransform(values$dds_obj))
      #log2things

      selectedLogvalues <- log2things[genevec_ids,]

      # check that I do not have nas or similar...

      if(length(genevec_ids)==length(genevec)){
        rowlabs <- genevec
      } else {
        rowlabs <- genevec_ids
        # rowlabs <- ifelse(, genevec, genevec_ids)
      }
      pheatmap(selectedLogvalues,scale="row",labels_row=rowlabs,main = myterm)







      # values$topgo_down



    })














    ## some dev section on usage of ihw
    # ihwres <- reactive({
    #   de_res <- as.data.frame(res_obj)
    #   ihw_res <- ihw(pvalue ~ baseMean,  data = de_res, alpha = 0.1)
    #   ihw_res
    # })
    #
    # output$debugihw <- renderPrint({
    #   ihw_res # replace then here with ihwres()
    # })
    #
    # output$ihwp1 <- renderPlot({
    #   plot(ihw_res)
    # })
    #
    # output$ihwp2 <- renderPlot({
    #   plot(ihw_res, what = "decisionboundary")
    # })
    #
    # output$ihwp3 <- renderPlot({
    #   gg <- ggplot(as.data.frame(ihw_res), aes(x = pvalue, y = adj_pvalue, col = group)) +
    #     geom_point(size = 0.25) +
    #     scale_colour_hue(l = 70, c = 150, drop = FALSE)
    #
    #   gg
    # })
    #
    # output$ihwp4 <- renderPlot({
    #   de_res <- res_airway ## CHANGE!!
    #   de_res <- na.omit(de_res)
    #   de_res$geneid <- as.numeric(gsub("ENSG[+]*", "", rownames(de_res)))
    #
    #   # set up data frame for plotting
    #   df <- rbind(data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$baseMean)/nrow(de_res),
    #                          covariate_type="base mean"),
    #               data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$geneid)/nrow(de_res),
    #                          covariate_type="gene id"))
    #
    #   ggplot(df, aes(x=covariate, y = -log10(pvalue))) +
    #     geom_hex(bins = 100) +
    #     facet_grid( . ~ covariate_type)
    #
    # })
    #
    #
    #
    # output$ihwp5 <- renderPlot({
    #   de_res <- as.data.frame(res_airway) ## CHANGE!!
    #   de_res <- na.omit(de_res)
    #   ggplot(de_res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
    # })
    #
    # output$ihwp6 <- renderPlot({
    #   de_res <- as.data.frame(res_airway) ## CHANGE!!
    #   de_res <- na.omit(de_res)
    #   de_res$baseMean_group <- groups_by_filter(de_res$baseMean, 8)
    #
    #   ggplot(de_res, aes(x=pvalue)) +
    #     geom_histogram(binwidth = 0.025, boundary = 0) +
    #     facet_wrap( ~ baseMean_group, nrow = 2)
    # })
    #
    # output$ihwp7 <- renderPlot({
    #   de_res <- as.data.frame(res_airway) ## CHANGE!!
    #   de_res <- na.omit(de_res)
    #   de_res$baseMean_group <- groups_by_filter(de_res$baseMean, 8)
    #   ggplot(de_res, aes(x = pvalue, col = baseMean_group)) + stat_ecdf(geom = "step")
    # })
    #
    # output$ihwp8 <- renderPlot({
    #   de_res <- as.data.frame(res_airway) ## CHANGE!!
    #   de_res <- na.omit(de_res)
    #
    #   de_res$lfc_group <- groups_by_filter(abs(de_res$log2FoldChange),8)
    #
    #   ggplot(de_res, aes(x = pvalue)) +
    #     geom_histogram(binwidth = 0.025, boundary = 0) +
    #     facet_wrap( ~ lfc_group, nrow=2)
    # })




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
      updateSelectizeInput(session = session, inputId = 'avail_ids', choices = c(Choose = '', rownames(values$res_obj)), server = TRUE)
    })

    observe({
      updateSelectizeInput(session = session, inputId = 'avail_symbols', choices = c(Choose = '', unname(values$res_obj$symbol)), server = TRUE) # using unname to prevent issues if res is modified with mapIds
    })



    output$available_genes <- renderUI({
      if("symbol" %in% names(values$res_obj)) {
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
      selectInput("choose_expfac",label = "choose the experimental factor to build the contrast upon",
                  choices = c("",design_factors()), selected = "")
    })


    ## LRT test...
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
             "Please select one level of the factor to build the contrast upon - contrast1"
        )
      )
      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]

      fac1_levels <- levels(fac1_vals)
      if(class(colData(values$dds_obj)[,fac1]) == "factor")
        selectInput("fac1_c1","c1",choices = c("",fac1_levels), selected = "")
      # selectInput("fac1_c2","c2",choices = fac1_levels)
    })

    output$fac2 <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             "Please select the other level of the factor to build the contrast upon - contrast2"
        )
      )
      fac1 <- input$choose_expfac
      fac1_vals <- colData(values$dds_obj)[,fac1]
      fac1_levels <- levels(fac1_vals)
      if(class(colData(values$dds_obj)[,fac1]) == "factor")
        # selectInput("fac1_c1","c1",choices = fac1_levels)
        selectInput("fac1_c2","c2",choices = c("",fac1_levels), selected = "")
    })

    output$facnum <- renderPrint({
      shiny::validate(
        need(input$choose_expfac!="",
             "Please select one level of the factor to build the contrast upon - contrast1_NUM!"
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

      if((class(colData(values$dds_obj)[,fac1]) %in% c("integer","numeric"))){

        shiny::validate(
          need(input$resu_addmle==FALSE,
               "Set the Add the unshrunken MLE to FALSE")
        )
      }
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
                     if(class(colData(values$dds_obj)[,input$choose_expfac]) == "factor") {
                       if(input$resu_ihw)
                         values$res_obj <- results(values$dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                   independentFiltering = input$resu_indfil, alpha = input$FDR, addMLE = input$resu_addmle,
                                                   filterFun = ihw)
                       else
                         values$res_obj <- results(values$dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                   independentFiltering = input$resu_indfil, alpha = input$FDR, addMLE = input$resu_addmle)
                     }
                     ## think more whether to include or not the IHW and ihw as filterFun...
                     if(class(colData(values$dds_obj)[,input$choose_expfac]) %in% c("integer","numeric"))
                       values$res_obj <- results(values$dds_obj,name = input$choose_expfac,
                                                 independentFiltering = input$resu_indfil, alpha = input$FDR, addMLE = input$resu_addmle)
                     if(!is.null(values$annotation_obj))
                       values$res_obj$symbol <- values$annotation_obj$gene_name[match(rownames(values$res_obj),
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
      summary(values$res_obj,alpha = 0.05) # use fdr shiny widget

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


    output$pvals_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )

      res_df <- as.data.frame(values$res_obj)
      p <- ggplot(res_df, aes_string("pvalue")) +
        geom_histogram(binwidth = 0.01) + theme_bw()

      p

    })

    output$logfc_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )

      res_df <- as.data.frame(values$res_obj)
      p <- ggplot(res_df, aes_string("log2FoldChange")) +
        geom_histogram(binwidth = 0.1) + theme_bw()

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
      plot_ma(values$res_obj,annotation_obj = values$annotation_obj)
    })

    output$mazoom <- renderPlot({
      if(is.null(input$ma_brush)) return(ggplot() + annotate("text",label="click and drag to zoom in",0,0) + theme_bw())

      if(!is.null(values$annotation_obj))
        plot_ma(values$res_obj,annotation_obj = values$annotation_obj) +
        coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                        ylim = c(input$ma_brush$ymin,input$ma_brush$ymax)) +
        geom_text(aes_string(label="genename"),size=3,hjust=0.25, vjust=-0.75)
      else
        plot_ma(values$res_obj,annotation_obj = values$annotation_obj) +
        coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                        ylim = c(input$ma_brush$ymin,input$ma_brush$ymax))
    })


    output$ma_highlight <- renderPlot({
      if("symbol" %in% names(values$res_obj)) {
        plot_ma_highlight(values$res_obj,
                          intgenes = input$avail_symbols,annotation_obj = values$annotation_obj)
      } else {
        plot_ma_highlight(values$res_obj,
                          intgenes = input$avail_ids,annotation_obj = values$annotation_obj)
      }
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
      plot_volcano(values$res_obj, FDR = input$FDR)
    })




    ## TODO: same thing but with the volcano plot?


    output$geneplot <- renderPlot({

      # if(length(input$color_by_G)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
      if(is.null(input$ma_brush)) return(NULL)

      if(is.null(input$mazoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

      selectedGene <- as.character(curDataClick()$ID)
      selectedGeneSymbol <- values$annotation_obj$gene_name[match(selectedGene,rownames(values$annotation_obj))]
      # plotCounts(dds_cleaner,)
      genedata <- plotCounts(values$dds_obj,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

      # onlyfactors <- genedata[match(input$color_by_G,colnames(genedata))]
      # onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]
      onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]

      ## intgroup can be a vector of factors. then i need interactions of the two factors
      # plotCounts(ddsmf_global,gene="ENSMUSG00000026080",intgroup = c("tissue","condition"),returnData = T) -> dh
      # dh$f1f2 <- interaction(dh$tissue,dh$condition)
      # dh  %>% ggplot(aes(x=f1f2,y=count,fill=f1f2)) + geom_boxplot()



      # genedata$plotby <- lapply(1:ncol(onlyfactors),function(arg) onlyfactors[,arg]) %>% interaction()
      genedata$plotby <- interaction(onlyfactors)



      ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) + geom_boxplot() + scale_y_log10(name="Normalized counts") + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes_string(x="plotby",y="count"),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
      # exportPlots$genesZoom <- res
      # res
    })


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

      p



    })



    output$rentrez_infobox <- renderUI({
      shiny::validate(
        need(
          (nrow(curDataClick()) > 0),
          "Select a gene first to display additional info (ENTREZ db)"
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




      # return(fullinfo$summary)




    })






    #
    #   if(!is.null(values$res_obj) & !is.null(values$dds_obj)) {
    #
    #
    #   }
    #

    cur_combires <- reactive({

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
      p
    })






    ## REPORT EDITOR
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
          tmp_content <- paste0(rmd_yaml(),input$acereport_rmd,collapse = "\n")
          isolate(HTML(knit2html(text = tmp_content, fragment.only = TRUE, quiet = TRUE)))
        },
        ## TODO: I should just "render on the fly" the Rmd without messing with the folders...
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
        # input$acereport_rmd
        if(input$rmd_dl_format == "rmd") {
          cat(tmp_content,file=file,sep="\n")
        } else {
          # write it somewhere too keeping the source
          # tmpfile <- tempfile()
          # file.create(tmpfile)
          # fileConn<- file(tempfile())
          # writeLines(tmp_content, fileConn)
          # close(fileConn)
          if(input$rmd_dl_format == "html") {
            cat(tmp_content,file="ideal_tempreport.Rmd",sep="\n")
            ## TODO: use writeLines maybe instead of this?


            # attempt?
            ## seems to work - otherwise revert to previous
            ## TODO: see how this can work for the live rendering!?
            src <- normalizePath("ideal_tempreport.Rmd")
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            cwd <- getwd()
            owd <- setwd(tempdir())
            # on.exit(setwd(owd))
            file.copy(src, "ideal_tempreport.Rmd",overwrite=TRUE)


            withProgress(rmarkdown::render(input = "ideal_tempreport.Rmd",
                                           output_file = file,


                                           # fragment.only = TRUE,
                                           quiet = TRUE),
                         message = "Generating the html report",
                         detail = "This can take some time")
            setwd(cwd)
          }
        }
      })



    ## STATE SAVING
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

    #
    #   # try out the carousel
    #   nBins = reactive({
    #     13
    #   })
    #
    #   output$distPlot1 <- renderPlot({
    #
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = nBins() + 1)
    #
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white')
    #
    #   })
    #
    #   output$distPlot2 <- renderPlot({
    #
    #     #browser()
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = nBins() + 1)
    #
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'steelblue', border = 'white')
    #
    #   })







  })





  shinyApp(ui = ideal_ui, server = ideal_server)


}






