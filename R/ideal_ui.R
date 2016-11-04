ideal_ui <- shinydashboard::dashboardPage(
  dashboardHeader(
    title = paste0("ideal - Interactive Differential Expression AnaLysis ",
                   packageVersion("ideal")),
    titleWidth = 900,

    # task menu for saving state to environment or binary data
    shinydashboard::dropdownMenu(type = "tasks",icon = icon("cog"),badgeStatus = "success",
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
                      "))
      ),

    ## main structure of the body for the dashboard
    fluidRow(
      valueBoxOutput("box_ddsobj"),
      valueBoxOutput("box_annobj"),
      valueBoxOutput("box_resobj")
    ),

    tabBox(
      width=12,

      tabPanel(
        "Data Setup",icon = icon("upload"),

        hr(),

        h2("Step 1: Upload your count matrix and the info on the experimental design"),

        uiOutput("upload_count_matrix"),
        uiOutput("upload_metadata"),

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
        ),
        # verbatimTextOutput("eddesign"),
        hr(),

        uiOutput("ui_step2"),
        uiOutput("ddsdesign"),
        verbatimTextOutput("debugdesign"),


        hr(),
        uiOutput("ui_step3"),
        uiOutput("ui_diydds"),
        # uiOutput("ok_dds"),
        verbatimTextOutput("debugdiy"),

        hr(),
        uiOutput("ui_stepanno"),
        uiOutput("ui_selectspecies"),
        verbatimTextOutput("speciespkg"),
        selectInput("idtype", "select the id type in your data", choices=c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")),
        # selectInput("orgdbspecies","select the species for your data",
        #             choices = c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db"),
        #             selected=""),
        verbatimTextOutput("printDIYanno"),
        # descr = c("Anopheles","Arabidopsis","Bovine","Worm","Canine","Fly","Zebrafish","E coli strain K12","E coli strain Sakai","Chicken","Human","org.Hs.ipi.db","Mouse","Rhesus","Malaria","Chimp","Rat","Yeast","Streptomyces coelicolor","Pig","Toxoplasma gondii","Xenopus"))
        uiOutput("ui_getanno"),
        ## this ideally populates also the list of genes of interest to choose among

        hr(),
        uiOutput("ui_stepoutlier"),
        uiOutput("ui_selectoutliers"),
        uiOutput("outliersout"),
        verbatimTextOutput("printremoved"),



        hr(),
        uiOutput("ui_step4"),
        uiOutput("rundeseq"),

        verbatimTextOutput("printDIYresults"),



        # verbatimTextOutput("printdds"),
        # verbatimTextOutput("printres"),
        uiOutput("ui_stepend")
      ),



      tabPanel(
        "Instructions",  icon = icon("info-circle"),
        includeMarkdown(system.file("extdata", "instructions.md",package = "ideal")),
        footer()
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
          h3("Basic summary for the counts"),
          p("Number of uniquely aligned reads assigned to each sample"),
          # verbatimTextOutput("reads_summary"),
          wellPanel(
            fluidRow(
              column(
                width = 4,
                numericInput("threshold_rowsums","Threshold on the row sums of the counts",value = 0, min = 0)),
              column(
                width = 4,
                numericInput("threshold_rowmeans","Threshold on the row means of the normalized counts",value = 0, min = 0))
            )),
          p("According to the selected filtering criteria, this is an overview on the provided count data"),
          verbatimTextOutput("detected_genes"),

          ## TODO: section to filter out features manually?

          selectInput("filter_crit",label = "Choose the filtering criterium",
                      choices = c("row means", "row sums"), selected = "row means"),

          actionButton("featfilt_dds", "Filter the DDS object"),




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
        "View Results", icon = icon("table"),

        # verbatimTextOutput("id"),
        # see: http://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value?noredirect=1&lq=1
        conditionalPanel(
          condition="!output.checkdds",

          # conditionalPanel(
          #   condition="output.checkresu==0",
          #   h2('RESU not provided')
          # ),

          ## TODO: exploe if possible to use conditional panels?

          # "Data Overview", icon = icon("eye"),

          uiOutput("choose_fac"),
          uiOutput("fac1"),
          uiOutput("fac2"),
          uiOutput("runresults"),

          uiOutput("facnum"),

          uiOutput("lrtavailable"),
          uiOutput("lrtfull"),
          uiOutput("lrtreduced"),

          uiOutput("runlrt"),

          uiOutput("store_result"),

          ## general options for result function
          # alpha is set via FDR on the left side
          selectInput("resu_indfil",label = "Apply independent filtering automatically",
                      choices = c(TRUE,FALSE), selected = TRUE),
          selectInput("resu_addmle",label = "Add the unshrunken MLE of log2 fold change",
                      choices = c(TRUE,FALSE), selected = TRUE),
          selectInput("resu_ihw", "Use Independent Hypothesis Weighting as a filtering function",
                      choices = c(TRUE, FALSE), selected = FALSE),
          #, evtl also the *filter* parameter of the function, i.e. baseMean if not specified

          verbatimTextOutput("diyres_summary"),

          DT::dataTableOutput("table_res"),
          plotOutput("pvals_hist"),
          plotOutput("logfc_hist")#,
          ,h4('Visible')
        ),
        conditionalPanel(
          condition="output.checkdds",
          h2("You did not create the dds object yet. Please go the main tab and generate it")
        )

        # verbatimTextOutput("diyres_summary"),
        # verbatimTextOutput("diyres")
      ),

      tabPanel(
        "Gene Lists", icon = icon("list-alt"),
        conditionalPanel(
          condition="!output.checkresu",
          tabBox(
            width = NULL,
            id="gse_tabbox",
            tabPanel("UPregu", icon = icon("arrow-circle-up"),
                     actionButton("button_enrUP", "Perform gene set enrichment analysis on the upregulated genes"),
                     DT::dataTableOutput("DT_gse_up")
            ),
            tabPanel("DOWNregu", icon = icon("arrow-circle-down"),
                     actionButton("button_enrDOWN", "Perform gene set enrichment analysis on the downregulated genes"),
                     actionButton("button_enrDOWN_goseq", "Perform gene set enrichment analysis on the downregulated genes"),
                     actionButton("button_enrDOWN_topgo", "Perform gene set enrichment analysis on the downregulated genes"),
                     DT::dataTableOutput("DT_gse_down"),
                     DT::dataTableOutput("DT_gse_down_goseq"),
                     DT::dataTableOutput("DT_gse_down_topgo")
            ),
            tabPanel("UPDOWN", icon = icon("arrows-v"),
                     actionButton("button_enrUPDOWN", "Perform gene set enrichment analysis on the up- and downregulated genes"),
                     DT::dataTableOutput("DT_gse_updown")
            ),
            tabPanel("List1", icon = icon("list"),
                     fileInput(inputId = "gl1",
                               label = "Upload a gene list file",
                               accept = c("text/csv", "text/comma-separated-values",
                                          "text/tab-separated-values", "text/plain",
                                          ".csv", ".tsv"), multiple = FALSE),
                     actionButton("button_enrLIST1", "Perform gene set enrichment analysis on the genes in list1"),
                     actionButton("button_enrLIST1_topgo", "Perform gene set enrichment analysis on the list1 genes"),
                     DT::dataTableOutput("DT_gse_list1"),
                     DT::dataTableOutput("DT_gse_list1_topgo")
            ),
            tabPanel("List2", icon = icon("list-alt"),
                     fileInput(inputId = "gl2",
                               label = "Upload a gene list file",
                               accept = c("text/csv", "text/comma-separated-values",
                                          "text/tab-separated-values", "text/plain",
                                          ".csv", ".tsv"), multiple = FALSE),
                     actionButton("button_enrLIST2", "Perform gene set enrichment analysis on the genes in list2"),
                     DT::dataTableOutput("DT_gse_list2")
            )
          ),
          ## will put collapsible list elements? or multi tab panel? or something to select on the left, and operate output-wise on the right e.g. venn diagrams or table for gene set enrichment
          h3("custom list 3 - handpicked") # use the select input from the left column?
          ,verbatimTextOutput("debuggls"),
          h2("Enrichment on the lists"),

          verbatimTextOutput("printUPgenes"),
          verbatimTextOutput("debuglists"),
          checkboxInput("toggle_updown","Use up and down regulated genes", TRUE),
          checkboxInput("toggle_up","Use up regulated genes", FALSE),
          checkboxInput("toggle_down","Use down regulated genes", FALSE),
          checkboxInput("toggle_list1","Use list1 genes", TRUE),
          checkboxInput("toggle_list2","Use list2 genes", FALSE),
          checkboxInput("toggle_list3","Use list3 genes", FALSE),


          plotOutput("vennlists")
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
          fluidRow(verbatimTextOutput("deb")),
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
          plotOutput("genefinder_plot"),
          plotOutput("volcanoplot"),
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
          )
          ,
          # fluidRow(dataTableOutput('ma_brush_out')),

          box(
            title = "Brushed table", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, collapsed = TRUE, width = 12,
            fluidRow(dataTableOutput("ma_brush_out")))
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
          # verbatimTextOutput("d1"),
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
              title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
              radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = T),
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
              title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
              checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
              conditionalPanel(
                "input.enableAutocomplete",
                wellPanel(
                  checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
                  checkboxInput("enableRCompletion", "R code completion", TRUE)
                )
              ),

              selectInput("mode", "Mode: ", choices=modes, selected="markdown"),
              selectInput("theme", "Theme: ", choices=themes, selected="solarized_light"))
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
        #             shiny::verbatimTextOutput("showuploaded1"),
        #             shiny::verbatimTextOutput("showuploaded2"),
        #             shiny::verbatimTextOutput("showuploaded3"),
        #             shiny::verbatimTextOutput("showuploaded4"),

        h4("Session Info"),
        verbatimTextOutput("sessioninfo"),
        footer()
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
      ),
  skin="blue"
      )
