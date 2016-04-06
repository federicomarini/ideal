# ideal.R

ideal <- function(dds=NULL,
                        rlt=NULL,
                        countmatrix=NULL,
                        coldata=NULL,
                        pca2go=NULL,
                        annotation=NULL){
  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  # library("shinyURL")

  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  myui <-
    shinydashboard::dashboardPage(
      dashboardHeader(
        title = paste0("ideal - Interactive Differential Expression AnaLysis",
                       " - version ",packageVersion("ideal")),
        titleWidth = 900),

      dashboardSidebar(
        width = 280,
        menuItem("Data upload",icon = icon("upload"),
                 uiOutput("upload_count_matrix"),
                 shinyBS::bsTooltip(
                   "upload_count_matrix", paste0("Select file containing the count matrix"),
                   "right", options = list(container = "body")),
                 uiOutput("upload_metadata"),
                 shinyBS::bsTooltip(
                   "upload_metadata", paste0("Select file containing the samples metadata"),
                   "right", options = list(container = "body")),
                 uiOutput("upload_annotation"),
                 shinyBS::bsTooltip(
                   "upload_annotation", paste0("Select file containing the annotation data"),
                   "right", options = list(container = "body"))),
        menuItem("App settings",icon = icon("cogs"),
                 selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8, selected = 1),
                 shinyBS::bsTooltip(
                   "pc_x", paste0("Select the principal component to display on the x axis"),
                   "right", options = list(container = "body")),
                 selectInput('pc_y', label = 'y-axis PC: ', choices = 1:8, selected = 2),
                 shinyBS::bsTooltip(
                   "pc_y", paste0("Select the principal component to display on the y axis"),
                   "right", options = list(container = "body")),
                 uiOutput("color_by"),
                 shinyBS::bsTooltip(
                   "color_by", paste0("Select the group of samples to stratify the analysis. Can also assume multiple values"),
                   "right", options = list(container = "body")),
                 numericInput('pca_nrgenes', label = 'Nr of (most variant) genes:', value = 300,min = 50,max = 20000),
                 shinyBS::bsTooltip(
                   "pca_nrgenes", paste0("Number of genes to select for computing the principal components. The top n genes are",
                                         " selected ranked by their variance inter-samples"),
                   "right", options = list(container = "body")),
                 numericInput('pca_point_alpha', label = 'Alpha: ', value = 1,min = 0,max = 1,step = 0.01),
                 shinyBS::bsTooltip(
                   "pca_point_alpha", paste0("Color transparency for the plots. Can assume values from 0 (transparent) ",
                                             "to 1 (opaque)"),
                   "right", options = list(container = "body")),
                 numericInput('pca_label_size', label = 'Labels size: ', value = 2,min = 1,max = 8),
                 shinyBS::bsTooltip(
                   "pca_label_size", paste0("Size of the labels for the samples in the principal components plots"),
                   "right", options = list(container = "body")),
                 numericInput('pca_point_size', label = 'Points size: ', value = 2,min = 1,max = 8),
                 shinyBS::bsTooltip(
                   "pca_point_size", paste0("Size of the points to be plotted in the principal components plots"),
                   "right", options = list(container = "body")),
                 numericInput('pca_varname_size', label = 'Variable name size: ', value = 4,min = 1,max = 8),
                 shinyBS::bsTooltip(
                   "pca_varname_size", paste0("Size of the labels for the genes PCA - correspond to the samples names"),
                   "right", options = list(container = "body")),
                 numericInput('pca_scale_arrow', label = 'Scaling factor : ', value = 1,min = 0.01,max = 10),
                 shinyBS::bsTooltip(
                   "pca_scale_arrow", paste0("Scale value for resizing the arrow corresponding to the variables in the ",
                                             "PCA for the genes. It should be used for mere visualization purposes"),
                   "right", options = list(container = "body")),
                 selectInput("col_palette","Color palette",choices = list("hue","set1","rainbow")),
                 shinyBS::bsTooltip(
                   "col_palette", paste0("Select the color palette to be used in the principal components plots. The number of ",
                                         "colors is selected automatically according to the number of samples and to the levels ",
                                         "of the factors of interest and their interactions"),
                   "right", options = list(container = "body"))
        ),
        menuItem("Plot export settings", icon = icon("paint-brush"),

                 numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                 shinyBS::bsTooltip(
                   "export_width", paste0("Width of the figures to export, expressed in cm"),
                   "right", options = list(container = "body")),
                 numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2),
                 shinyBS::bsTooltip(
                   "export_height", paste0("Height of the figures to export, expressed in cm"),
                   "right", options = list(container = "body"))

        )
      ),

      dashboardBody(

        ## Define output size and style of error messages
        tags$head(
          tags$style(HTML("
                          .shiny-output-error-validation {
                          font-size: 15px;
                          color: forestgreen;
                          text-align: center;
                          }
                          "))
          ),

        tabBox(
          width=12,

          tabPanel(
            "About",
            includeMarkdown(system.file("extdata", "about.md",package = "pcaExplorer")),
            hr(),
            #             shiny::verbatimTextOutput("showuploaded1"),
            #             shiny::verbatimTextOutput("showuploaded2"),
            #             shiny::verbatimTextOutput("showuploaded3"),
            #             shiny::verbatimTextOutput("showuploaded4"),

            h4("Session Info"),
            verbatimTextOutput("sessioninfo"),
            footer()
          ),

          tabPanel(
            "Instructions",
            includeMarkdown(system.file("extdata", "instructions.md",package = "pcaExplorer"))
          ),

          tabPanel(
            "Data Preview",
            h1("Sneak peek in the data"),

            h3("General information on the provided SummarizedExperiment/DESeqDataSet"),
            shiny::verbatimTextOutput("showdata"),
            h3("Available metadata"),
            DT::dataTableOutput("showcoldata"),
            h3("Number of million of reads per sample"),
            plotOutput("reads_barplot"),
            h3("Basic summary for the counts"),
            verbatimTextOutput("reads_summary"),

            footer()
          ),

          tabPanel(
            "Samples View",
            p(h1('Principal Component Analysis on the samples'), "PCA projections of sample expression profiles onto any pair of components."),
            fluidRow(checkboxInput("sample_labels","Display sample labels",value = TRUE)),
            fluidRow(
              column(
                width = 6,
                plotOutput('samples_pca',brush = "pca_brush"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_samplesPca", "Download Plot"),
                    textInput("filename_samplesPca",label = "Save as...",value = "samplesPca.pdf"))
              ),
              column(
                width= 6,
                plotOutput("samples_scree"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_samplesScree", "Download Plot"),
                    textInput("filename_samplesScree",label = "Save as...",value = "samplesScree.pdf")),
                fluidRow(
                  column(
                    width = 6,
                    radioButtons("scree_type","Scree plot type:",choices=list("Proportion of explained variance"="pev","Cumulative proportion of explained variance"="cev"),"pev")
                  ),
                  column(
                    width = 6,
                    numericInput("scree_pcnr","Number of PCs to display",value=8,min=2)
                  )
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 6,
                plotOutput("samples_pca_zoom")
              ),
              column(
                width = 6,
                numericInput("ntophiload", "Nr of genes to display (top & bottom)",value = 10, min = 1, max=40),
                plotOutput("geneshiload")
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 6,
                p(h4('Outlier Identification'), "Toggle which samples to remove - suspected to be considered as outliers"),
                uiOutput("ui_outliersamples"),
                plotOutput("samples_outliersremoved")
              )

            )
          ),

          tabPanel(
            "Genes View",
            p(h1('Principal Component Analysis on the genes'), "PCA projections of genes abundances onto any pair of components."),

            # shinyURL.ui(),

            fluidRow(checkboxInput("variable_labels","Display variable labels",value = TRUE)),
            fluidRow(
              column(
                width = 4,
                h4("Main Plot - interact!"),
                plotOutput('genes_biplot',brush = 'pcagenes_brush',click="pcagenes_click"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesPca", "Download Plot"),
                    textInput("filename_genesPca",label = "Save as...",value = "genesPca.pdf"))),
              column(
                width = 4,
                h4("Zoomed window"),
                plotOutput("genes_biplot_zoom",click="pcagenes_zoom_click"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesZoom", "Download Plot"),
                    textInput("filename_genesZoom",label = "Save as...",value = "genesPca_zoomed.pdf"))),
              column(
                width = 4,
                h4("Boxplot of selected gene"),
                plotOutput("genes_biplot_boxplot"))),

            fluidRow(
              column(
                width = 6,
                h4("Zoomed heatmap"),
                plotOutput("heatzoom"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesHeatmap","Download Plot"),
                    textInput("filename_genesHeatmap",label = "Save as...",value = "genesHeatmap.pdf"))),
              column(
                width = 6,
                h4("Zoomed interactive heatmap"),
                fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
                fluidRow(d3heatmapOutput("heatzoomd3")))),
            hr(),
            fluidRow(
              column(
                width = 6,
                h4("Points selected by brushing - clicking and dragging:"),
                DT::dataTableOutput("pca_brush_out"),
                downloadButton('downloadData_brush', 'Download brushed points'),
                textInput("brushedPoints_filename","File name...")),
              column(
                width = 6,
                h4("Points selected by clicking:"),
                DT::dataTableOutput("pca_click_out"),
                textInput("clickedPoints_filename","File name..."),
                downloadButton('downloadData_click', 'Download clicked (or nearby) points'))
            )
          ),

          tabPanel(
            "Gene finder",
            fluidRow(
              h1("GeneFinder"),
              textInput("genefinder",label = "type in the name of the gene to search",value = NULL),
              shinyBS::bsTooltip(
                "genefinder", paste0("Type in the name of the gene to search. If no annotation is ",
                                     "provided, you need to use IDs that are the row names of the ",
                                     "objects you are using - count matrix, SummarizedExperiments ",
                                     "or similar. If an annotation is provided, that also contains ",
                                     "gene symbols or similar, the gene finder tries to find the ",
                                     "name and the ID, and it suggests if some characters are in a ",
                                     "different case"),
                "right", options = list(container = "body")),

              #               fluidRow(
              #                 column(
              #                   width = 6,
              #                   uiOutput("ui_selectID")
              #                 ),
              #                 column(
              #                   width = 6,
              #                   uiOutput("ui_selectName")
              #                 )
              #               ),
              # verbatimTextOutput("debugf"),


              verbatimTextOutput("searchresult"),
              verbatimTextOutput("debuggene"),
              checkboxInput("ylimZero","Set y axis limit to 0",value=TRUE),
              # plotOutput("newgenefinder_plot"),
              plotOutput("genefinder_plot"))
          ),

          tabPanel(
            "PCA2GO",
            h1("pca2go - Functional annotation of PC"),
            h4("Functions enriched in the genes with high loadings on the selected principal components"),
            # verbatimTextOutput("enrichinfo"),


            column(
              width = 6,
              uiOutput("ui_selectspecies")
            ),
            column(
              width = 6,
              uiOutput("ui_inputtype")
            ),

            shinyBS::bsTooltip(
              "ui_selectspecies", paste0("Select the species for the functional enrichment analysis, ",
                                         "choosing among the ones currently supported by limma::goana. ",
                                         "Alternatively, for other species, it can be possible to use one ",
                                         "of the available annotation packages in Bioconductor, and pre-",
                                         "computing the pca2go object in advance"),
              "bottom", options = list(container = "body")),
            verbatimTextOutput("speciespkg"),
            checkboxInput("compact_pca2go","Display compact tables",value=FALSE),
            shinyBS::bsTooltip(
              "compact_pca2go", paste0("Should I display all the columns? If the information content of the ",
                                       "tables is somehow too much for the screen width, as it can be for ",
                                       "objects generated by pca2go with the topGO routines, the app can ",
                                       "display just an essential subset of the columns"),
              "bottom", options = list(container = "body")),

            uiOutput("ui_computePCA2GO"),
            shinyBS::bsTooltip(
              "ui_computePCA2GO", paste0("Compute a pca2go object, using the limma::goana function, ",
                                         "after selecting the species of the experiment under investigation"),
              "bottom", options = list(container = "body")),

            fluidRow(
              column(width = 3),
              column(
                width = 6,
                DT::dataTableOutput("dt_pcver_pos")),
              column(width = 3)
            ),

            fluidRow(
              column(4,
                     DT::dataTableOutput("dt_pchor_neg")),
              column(4,
                     plotOutput("pca2go")),
              column(4,
                     DT::dataTableOutput("dt_pchor_pos"))
            ),
            fluidRow(
              column(width = 3),
              column(
                width = 6,
                DT::dataTableOutput("dt_pcver_neg")),
              column(width = 3)
            )
          ),



          tabPanel(
            "Multifactor Exploration",
            h1("Multifactor exploration of datasets with >= 2 experimental factors"),

            verbatimTextOutput("intro_multifac"),

            fluidRow(
              column(
                width = 6,
                uiOutput("covar1")
              ),
              column(
                width = 6,
                uiOutput("covar2")
              )
            ),
            fluidRow(
              column(
                width = 6,
                uiOutput("c1levels")
              ),
              column(
                width = 6,
                uiOutput("c2levels")
              )
            ),
            fluidRow(
              column(
                width = 6,
                uiOutput("colnames1")
              ),
              column(
                width = 6,
                uiOutput("colnames2")
              )
            ),


            shinyBS::bsTooltip(
              "covar1", paste0("Select the first experimental factor"),
              "bottom", options = list(container = "body")),
            shinyBS::bsTooltip(
              "covar2", paste0("Select the second experimental factor"),
              "bottom", options = list(container = "body")),
            shinyBS::bsTooltip(
              "c1levels", paste0("For factor 1, select two levels to contrast"),
              "bottom", options = list(container = "body")),
            shinyBS::bsTooltip(
              "c2levels", paste0("For factor 2, select two or more levels to contrast"),
              "bottom", options = list(container = "body")),
            shinyBS::bsTooltip(
              "colnames1", paste0("Combine samples belonging to Factor1-Level1 samples for each level in Factor 2"),
              "bottom", options = list(container = "body")),
            shinyBS::bsTooltip(
              "colnames2", paste0("Combine samples belonging to Factor1-Level2 samples for each level in Factor 2"),
              "bottom", options = list(container = "body")),



            actionButton("composemat","Compose the matrix",icon=icon("spinner")),
            shinyBS::bsTooltip(
              "composemat", paste0("Select first two different experimental factors, for example ",
                                   "condition and tissue. For each factor, select two or more ",
                                   "levels. The corresponding samples which can be used are then displayed ",
                                   "in the select boxes. Select an equal number of samples for each of ",
                                   "the levels in factor 1, and then click the button to compute the ",
                                   "new matrix which will be used for the visualizations below"),
              "bottom", options = list(container = "body")),


            fluidRow(
              column(4,
                     selectInput('pc_x_multifac', label = 'x-axis PC: ', choices = 1:8,
                                 selected = 1)
              ),
              column(4,
                     selectInput('pc_y_multifac', label = 'y-axis PC: ', choices = 1:8,
                                 selected = 2)
              )),

            # fluidRow(verbatimTextOutput("multifacdebug")),

            fluidRow(
              column(6,
                     plotOutput('pcamultifac',brush = 'pcamultifac_brush')),
              column(6,
                     plotOutput("multifaczoom"))
            ),
            fluidRow(downloadButton('downloadData_brush_multifac', 'Download brushed points'),
                     textInput("brushedPoints_filename_multifac","File name..."),
                     DT::dataTableOutput('pcamultifac_out'))


          )

        )

          ),

      skin="blue"

          )

  ## ------------------------------------------------------------------ ##
  ##                          Define server                             ##
  ## ------------------------------------------------------------------ ##

  myserver <- shinyServer(function(input, output, session) {

    exportPlots <- reactiveValues(
      samplesPca=NULL,
      samplesZoom=NULL,
      samplesScree=NULL,
      genesPca=NULL,
      genesZoom=NULL,
      genesBoxplot=NULL,
      genesHeatmap=NULL,
      genefinder=NULL
    )

    values <- reactiveValues()
    values$mydds <- dds
    values$myrlt <- rlt
    values$mycountmatrix <- countmatrix
    values$mymetadata <- coldata
    values$mypca2go <- pca2go
    values$myannotation <- annotation

    user_settings <- reactiveValues(save_width = 45, save_height = 11)

    # shinyURL.server()

    if(!is.null(dds)){
      if(!is(dds,"DESeqDataSet"))
        stop("dds must be a DESeqDataSet object. If it is a simple counts matrix, provide it to the countmatrix parameter!")
    }
    if(!is.null(rlt)){
      if(!is(rlt,"DESeqTransform"))
        stop("dds must be a DESeqTransform object")
    }

    # compute only rlt if dds is provided but not cm&coldata
    if(!is.null(dds) & (is.null(countmatrix) & is.null(coldata)) & is.null(rlt))
      withProgress(message = "computing rlog transformed values...",
                   value = 0,
                   {
                     values$myrlt <- rlogTransformation(dds)
                   })


    output$color_by <- renderUI({
      if(is.null(values$mydds))
        return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('color_by', label = 'Group/color by: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })



    ## Render the UI element to upload the count matrix
    output$upload_count_matrix <- renderUI({
      if (!is.null(dds) | !is.null(countmatrix)) {
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
      cm <- utils::read.delim(input$uploadcmfile$datapath, header = TRUE,
                              as.is = TRUE, sep = "\t", quote = "",
                              check.names = FALSE)

      return(cm)
    })




    output$upload_metadata <- renderUI({
      if (!is.null(dds) | !is.null(coldata)) {
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
      coldata <- utils::read.delim(input$uploadmetadatafile$datapath, header = TRUE,
                                   as.is = TRUE, sep = "\t", quote = "",
                                   check.names = FALSE)

      return(coldata)
    })


    output$upload_annotation <- renderUI({
      if (!is.null(annotation)) {
        NULL
      } else {
        return(fileInput(inputId = "uploadannotationfile",
                         label = "Upload an annotation file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readAnnotation <- reactive({
      if (is.null(input$uploadannotationfile))
        return(NULL)
      annodata <- utils::read.delim(input$uploadannotationfile$datapath, header = TRUE,
                                    as.is = TRUE, sep = "\t", quote = "",
                                    check.names = FALSE)

      return(annodata)
    })


    createDDS <- reactive({
      if(is.null(countmatrix) | is.null(coldata))
        return(NULL)

      dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                                    colData = coldata,
                                    design=~1)

      return(dds)


    })

    createRLT <- reactive({
      if(is.null(countmatrix) | is.null(coldata))
        return(NULL)

      rlt <- rlogTransformation(values$mydds)

      return(rlt)


    })

    observeEvent(createDDS,
                 {
                   if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                     values$mydds <- createDDS()
                 })

    observeEvent(createRLT,
                 {
                   if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                     values$myrlt <- createRLT()
                 })




    # as in http://stackoverflow.com/questions/29716868/r-shiny-how-to-get-an-reactive-data-frame-updated-each-time-pressing-an-actionb
    observeEvent(input$uploadcmfile,
                 {
                   values$mycountmatrix <- readCountmatrix()
                   if(!is.null(values$mymetadata)){
                     withProgress(message="Computing the objects...",value = 0,{

                       values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                              colData = values$mymetadata,
                                                              design=~1)
                       values$myrlt <- rlogTransformation(values$mydds)})
                   }
                 })

    observeEvent(input$uploadmetadatafile,
                 {
                   values$mymetadata <- readMetadata()
                   if(!is.null(values$mycountmatrix)){
                     withProgress(message="Computing the objects...",value = 0,{

                       values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                              colData = values$mymetadata,
                                                              design=~1)
                       values$myrlt <- rlogTransformation(values$mydds)})
                   }
                 })

    observeEvent(input$uploadannotationfile,
                 {
                   values$myannotation <- readAnnotation()
                 })



    output$showuploaded1 <- renderPrint({
      head(values$mycountmatrix)
    })
    output$showuploaded2 <- renderPrint({
      values$mymetadata
    })
    output$showuploaded3 <- renderPrint({
      values$mydds
    })
    output$showuploaded4 <- renderPrint({
      values$myrlt
    })


    colSel <- reactive({
      # find out how many colors to generate: if no factor is selected, either
      # return all say steelblue or all different

      if(!is.null(input$color_by)) {
        expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
        expgroups <- interaction(expgroups)
      } else {
        expgroups <- factor(colnames(values$myrlt))
        # return(rep("steelblue",ncol(values$myrlt))) # to return all same
      }

      nrgroups <- length(levels(expgroups))

      # hue_pal()(ncol(values$myrlt)/6) # or somewhat other way

      if(input$col_palette=="hue"){
        return(hue_pal()(nrgroups))
      }
      # hue_pal()(ncol(values$myrlt)/2) # or somewhat other way
      if(input$col_palette=="set1"){
        if(nrgroups <= 9) { # max color nr allowed for set1
          return(brewer_pal(palette = "Set1")(nrgroups))
        } else {
          return(hue_pal()(nrgroups)) # plus print message?
        }
      }
      # (ncol(values$myrlt)/2) # or somewhat other way
      if(input$col_palette=="rainbow"){
        return(rainbow(nrgroups))
      }
      # hue_pal()(nrgroups) # or somewhat other way
    })


    output$sessioninfo <- renderPrint({
      sessionInfo()
    })




    output$showdata <- renderPrint({
      values$mydds
    })

    output$showcoldata <- DT::renderDataTable({
      datatable(as.data.frame(colData(values$mydds)))
    })

    output$reads_barplot <- renderPlot({
      rr <- colSums(counts(values$mydds))/1e6
      if(is.null(names(rr)))
        names(rr) <- paste0("sample_",1:length(rr))
      rrdf <- data.frame(Reads=rr,Sample=names(rr),stringsAsFactors = FALSE)
      if (!is.null(input$color_by)) {
        rrdf$Group <- colData(values$mydds)[input$color_by][[1]]
        p <- ggplot(rrdf,aes_string("Sample",weight="Reads")) + geom_bar(aes_string(fill="Group"))
        p
      } else {
        p <- ggplot(rrdf,aes_string("Sample",weight="Reads")) + geom_bar()
        p
      }
    })

    output$reads_summary <- renderPrint({
      summary(colSums(counts(values$mydds))/1e6)
    })








    output$samples_pca <- renderPlot({
      res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA"
      )
      res <- res + theme_bw()
      exportPlots$samplesPca <- res
      res
    })

    output$samples_pca_zoom <- renderPlot({

      shiny::validate(
        need(!is.null(input$pca_brush),
             "Zoom in by brushing in the main plot panel above"
        )
      )
      # if(is.null(input$pca_brush))
      # return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

      res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA - zoom in"
      )
      res <- res + xlim(input$pca_brush$xmin,input$pca_brush$xmax) + ylim(input$pca_brush$ymin,input$pca_brush$ymax)
      res <- res + theme_bw()
      exportPlots$samplesZoom <- res
      res
    })

    output$samples_scree <- renderPlot({
      rv <- rowVars(assay(values$myrlt))
      select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
      pca <- prcomp(t(assay(values$myrlt)[select, ]))

      res <- pcascree(pca,type = input$scree_type, pc_nr = input$scree_pcnr, title="Scree plot for the samples PCA")
      res <- res + theme_bw()
      exportPlots$samplesScree <- res
      res
    })


    output$geneshiload <- renderPlot({
      rv <- rowVars(assay(values$myrlt))
      select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
      pca <- prcomp(t(assay(values$myrlt)[select, ]))

      par(mfrow=c(2,1))
      hi_loadings(pca,whichpc = as.integer(input$pc_x),topN = input$ntophiload,annotation = values$myannotation)
      hi_loadings(pca,whichpc = as.integer(input$pc_y),topN = input$ntophiload,annotation = values$myannotation)

    })


    output$ui_outliersamples <- renderUI({
      available_samples <- c("",colnames(values$myrlt))

      selectInput("outlierselection",label = "Select which sample(s) to remove - suspected outliers",choices = available_samples,multiple = TRUE)

    })

    output$samples_outliersremoved <- renderPlot({

      shiny::validate(
        need(input$outlierselection!="",
             message = "Select at least one sample to plot the new PCA where the selection is removed")
      )

      currentrlt <- values$myrlt
      allsamples <- colnames(currentrlt)

      outliersamples <- input$outlierselection
      currentrlt <- currentrlt[,setdiff(allsamples,outliersamples)]

      res <- pcaplot(currentrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA"
      )
      res <- res + theme_bw()
      # exportPlots$samplesPca <- res
      res

    })





    output$genes_biplot <- renderPlot({
      if(!is.null(input$color_by)) {
        expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
        expgroups <- interaction(expgroups)
      } else {
        expgroups <- colnames(values$myrlt)
      }
      colGroups <- colSel()[factor(expgroups)]

      res <- genespca(values$myrlt,
                      ntop = input$pca_nrgenes,
                      choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                      biplot = TRUE,
                      arrowColors = factor(colGroups),
                      groupNames = expgroups,
                      alpha=input$pca_point_alpha,coordEqual=FALSE,useRownamesAsLabels=FALSE,labels.size=input$pca_label_size,
                      point_size=input$pca_point_size,varname.size=input$pca_varname_size, scaleArrow = input$pca_scale_arrow,annotation=values$myannotation)
      exportPlots$genesPca <- res
      res
    })


    output$genes_biplot_zoom <- renderPlot({
      # if(is.null(input$pcagenes_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Zoom in by brushing in the main panel - this will also allow displaying the gene names"
        )
      )

      if(!is.null(input$color_by)) {
        expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
        expgroups <- interaction(expgroups)
      } else {
        expgroups <- colnames(values$myrlt)
      }
      colGroups <- colSel()[factor(expgroups)]

      res <- genespca(values$myrlt,
                      ntop = input$pca_nrgenes,
                      choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                      biplot = TRUE,
                      arrowColors = factor(colGroups),
                      groupNames = expgroups,
                      alpha=input$pca_point_alpha,coordEqual=FALSE,
                      var.axes=input$variable_labels, # workaround for a ggplot2 bug/missing thing: here details: https://github.com/hadley/ggplot2/issues/905
                      labels.size=input$pca_label_size,varname.size=input$pca_varname_size,
                      scaleArrow = input$pca_scale_arrow,point_size=input$pca_point_size,annotation=values$myannotation)

      res <- res +
        xlim(input$pcagenes_brush$xmin,input$pcagenes_brush$xmax) +
        ylim(input$pcagenes_brush$ymin,input$pcagenes_brush$ymax)
      exportPlots$genesZoom <- res
      res
    })


    output$genes_biplot_boxplot <- renderPlot({
      # if(length(input$color_by)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
      # if(is.null(input$pcagenes_zoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          !is.null(input$pcagenes_zoom_click),
          "Click to generate the boxplot for the selected gene"
        )
      )

      selectedGene <- curData_zoomClick()$ids
      selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
      # plotCounts(dds_cleaner,)

      shiny::validate(
        need(nrow(curData_zoomClick()) >0,message = "Click closer to a gene to get the boxplot")

      )

      genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

      onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
      genedata$plotby <- interaction(onlyfactors)

      res <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) +
        geom_boxplot(outlier.shape = NA) + scale_y_log10(name="Normalized counts") +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes_string(x="plotby",y="count"),position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions")
      exportPlots$genesBoxplot <- res
      res
    })


    # for reading in the brushed/clicked points
    curData_brush <- reactive({
      df2 <- genespca(values$myrlt,
                      ntop = input$pca_nrgenes,
                      choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                      biplot = TRUE,
                      # arrowColors = colGroups,
                      alpha=input$pca_point_alpha,
                      returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_click <- reactive({
      df2 <- genespca(values$myrlt,
                      ntop = input$pca_nrgenes,
                      choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                      biplot = TRUE,
                      # arrowColors = colGroups,
                      alpha=input$pca_point_alpha,
                      returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- nearPoints(df2, input$pcagenes_click,
                        threshold = 20, maxpoints = 3,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_zoomClick <- reactive({
      df2 <- genespca(values$myrlt,
                      ntop = input$pca_nrgenes,
                      choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                      biplot = TRUE,
                      # arrowColors = colGroups,
                      alpha=input$pca_point_alpha,
                      returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- nearPoints(df2, input$pcagenes_zoom_click,
                        threshold = 20, maxpoints = 1,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })



    output$pca_brush_out <- DT::renderDataTable({
      datatable(curData_brush(),options = list(pageLength = 50))
    })

    output$pca_click_out <- DT::renderDataTable({
      datatable(curData_click(),options = list(pageLength = 50))
    })




    output$heatzoomd3 <- renderD3heatmap({
      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Brush the main panel above to generate a heatmap"
        )
      )

      # if(is.null(input$pcagenes_brush)) return(NULL)

      brushedObject <- curData_brush()
      shiny::validate(
        need(
          nrow(brushedObject) > 1,
          "Brush to include at least two genes"
        )
      )

      selectedGenes <- brushedObject$ids
      toplot <- assay(values$myrlt)[selectedGenes,]
      rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]

      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding

      d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
    })


    output$heatzoom <- renderPlot({
      # if(is.null(input$pcagenes_brush)) return(NULL)
      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Brush the main panel above to generate a heatmap"
        )
      )

      brushedObject <- curData_brush()
      shiny::validate(
        need(
          nrow(brushedObject) > 1,
          "Brush to include at least two genes"
        )
      )
      selectedGenes <- brushedObject$ids
      toplot <- assay(values$myrlt)[selectedGenes,]
      rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
      # pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
      NMF::aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
      ## aheatmap is actually consistent in displaying the clusters with most of other heatmap packages
      ## keep in mind: pheatmap does somehow a better job if scaling/centering
    })





    #     displayed by default, with possibility to select from gene id provided as row names of the objects
    #     output$ui_selectID <- renderUI({
    #       allIDs <- withProgress(message = "loading the names in the UI",value = 0,
    #                              {
    #                                rownames(values$myrlt)
    #                              }
    #
    #
    #       )
    #       selectInput("selectID",label = "Select ID",choices = c("",allIDs),selected=NULL)
    #     })
    #     # additionally displayed if an annotation is provided
    #     output$ui_selectName <- renderUI({
    #       shiny::validate(
    #         need(
    #           !is.null(values$myannotation),
    #           "If you provide an annotation table, you could search by the corresponding name/ID"
    #         )
    #       )
    #
    #       selectInput("selectName",label = "Select gene name",choices = c("",values$myannotation$gene_name),selected=NULL)
    #       # selectInput("selectName")
    #     })
    #
    #     output$debugf <- renderPrint({
    #       input$selectID
    #     })


    #     output$newgenefinder_plot <- renderPlot({
    #       if (input$selectID == "")
    #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
    #
    #
    #       if(!is.null(values$myannotation)){
    #         if(input$selectName != "") {
    #           selectedGeneName <- input$selectName
    #           selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$selectName)]
    #         } else {
    #           if (input$selectID == "") {
    #             return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
    #           } else {
    #             selectedGene <- input$selectID
    #             selectedGeneName <- ifelse(!is.null(values$myannotation),
    #                                        values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))],
    #                                        "")
    #           }
    #         }
    #       } else if (input$selectID != ""){
    #         selectedGene <- input$selectID
    #         selectedGeneName <- ifelse(!is.null(values$myannotation),
    #                                    values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))],
    #                                    "")
    #       } else {
    #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
    #       }
    #
    #
    #       anno_id <- rownames(values$myannotation)
    #       anno_gene <- values$myannotation$gene_name
    #
    #       if(is.null(input$color_by))
    #         return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
    # #       if(is.null(input$color_by) & (input$selectName=="" | input$selectID ==""))
    # #         return(ggplot() + annotate("text",label="Select a gene and a factor to plot gene",0,0) + theme_bw())
    # #       if((input$selectName=="" | input$selectID ==""))
    # #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
    #       # if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
    #         # return(ggplot() + annotate("text",label="Gene not found...",0,0) + theme_bw())
    #
    # #       if (input$genefinder %in% anno_id) {
    # #         selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
    # #         selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
    # #       }
    # #       if (input$genefinder %in% anno_gene) {
    # #         selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
    # #         if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
    # #         selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
    # #       }
    #
    #
    #
    #
    #
    #       genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)
    #
    #       onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
    #       genedata$plotby <- interaction(onlyfactors)
    #
    #       p <- ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + labs(title=paste0("Normalized counts for ",selectedGeneName," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
    #
    #       if(input$ylimZero)
    #       {
    #         p <- p + scale_y_log10(name="Normalized counts - log10 scale",limits=c(1,max(genedata$count)))
    #       } else {
    #         p <- p + scale_y_log10(name="Normalized counts - log10 scale")
    #       }
    #       exportPlots$genefinder <- p
    #
    #       p
    #     })
    #

    output$searchresult <- renderPrint({

      if(is.null(input$color_by)) return("Select a factor to plot your gene")
      if(input$genefinder=="")
        return("Type in the gene name/id you want to plot")

      foundGeneID <- input$genefinder %in% rownames(values$myrlt)
      foundGeneName <- input$genefinder %in% values$myannotation$gene_name
      if(!foundGeneID){
        foundGeneID <- toupper(input$genefinder) %in% toupper(rownames(values$myrlt))
        if(foundGeneID){
          return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                        unique(rownames(values$myannotation)[which(toupper(input$genefinder)==toupper(rownames(values$myannotation)))]),"?"))
        } else {
          foundGeneNAME <- input$genefinder %in% values$myannotation$gene_name
          if(!foundGeneNAME){
            foundGeneNAME <- toupper(input$genefinder) %in% toupper(values$myannotation$gene_name)
            if(foundGeneNAME){
              return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                            unique(values$myannotation$gene_name[which(toupper(input$genefinder)==toupper(values$myannotation$gene_name))]),"?"))
            } else {return("Could not find the gene you typed!")}
          } else {
            fgn <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
            if (length(fgn) > 1) return(paste0("Found more than one gene with the selected gene name. Select one of the following: ",paste(selectedGene,collapse=", ")))
            selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]

            fg <- rownames(values$myannotation)[match(fgn,values$myannotation$gene_name)]
            return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

          }}
      } else {
        fg <- rownames(values$myannotation)[match(input$genefinder,rownames(values$myrlt))]
        return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

      }
    })



    output$genefinder_plot <- renderPlot({
      anno_id <- rownames(values$myrlt)
      anno_gene <- values$myannotation$gene_name

      if(is.null(input$color_by) & input$genefinder!="")
        return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
      if(is.null(input$color_by) & input$genefinder=="")
        return(ggplot() + annotate("text",label="Select a gene and a factor to plot gene",0,0) + theme_bw())
      if(input$genefinder=="")
        return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
      if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
        return(ggplot() + annotate("text",label="Gene not found...",0,0) + theme_bw())

      if (input$genefinder %in% anno_id) {
        selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
        selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
      }
      if (input$genefinder %in% anno_gene) {
        selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
        if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
        selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
      }
      genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

      onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
      genedata$plotby <- interaction(onlyfactors)

      p <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) + geom_boxplot() + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes_string(x="plotby",y="count"),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")

      if(input$ylimZero)
      {
        p <- p + scale_y_log10(name="Normalized counts - log10 scale",limits=c(1,max(genedata$count)))
      } else {
        p <- p + scale_y_log10(name="Normalized counts - log10 scale")
      }
      exportPlots$genefinder <- p

      p
    })


    output$ui_computePCA2GO <- renderUI({
      if(is.null(pca2go))
        actionButton("computepca2go","Compute the PCA2GO object",icon=icon("spinner"))
    })

    annoSpecies_df <- data.frame(species=c("","Anopheles","Arabidopsis","Bovine","Worm",
                                           "Canine","Fly","Zebrafish","E coli strain K12",
                                           "E coli strain Sakai","Chicken","Human","Mouse",
                                           "Rhesus","Malaria","Chimp","Rat",
                                           "Yeast","Streptomyces coelicolor", "Pig","Toxoplasma gondii",
                                           "Xenopus"),
                                 pkg=c("","org.Ag.eg.db",	"org.At.tair.db", "org.Bt.eg.db",	"org.Ce.eg.db",
                                       "org.Cf.eg.db",	"org.Dm.eg.db", "org.Dr.eg.db",	"org.EcK12.eg.db",
                                       "org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db",	"org.Mm.eg.db",
                                       "org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db",
                                       "org.Sc.sgd.db","org.Sco.eg.db",	"org.Ss.eg.db","org.Tgondii.eg.db",
                                       "org.Xl.eg.db"),
                                 stringsAsFactors = FALSE)
    annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species),]
    annoSpecies_df <- annoSpecies_df[annoSpecies_df$species %in% c("","Human", "Mouse", "Rat", "Fly", "Chimp"),]

    output$ui_selectspecies <- renderUI({
      if(is.null(values$mypca2go)) {
        selectInput("speciesSelect",label = "Select the species of your samples",choices = annoSpecies_df$species,selected="")
      }
    })
    output$ui_inputtype <- renderUI({
      if(is.null(values$mypca2go)) {
        selectInput("idtype",label = "Select the input type of your identifiers",
                    choices = c("ENSEMBL","SYMBOL","REFSEQ","ENTREZID"), selected = "ENSEMBL")
      }
    })


    output$speciespkg <- renderText({

      if(!is.null(values$mypca2go))
        return("pca2go object provided")

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




    computedPCA2GO <- eventReactive( input$computepca2go, {
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
      withProgress(message = "Computing the PCA2GO object...",
                   value = 0,
                   {
                     pcpc <- limmaquickpca2go(values$myrlt,background_genes = rownames(values$mydds),
                                              inputType = input$idtype,
                                              organism = gsub(".eg.db","",gsub("org.","",annopkg)))
                   })
      pcpc
    })


    observeEvent(input$computepca2go,
                 {
                   values$mypca2go <- computedPCA2GO()
                 })



    output$pca2go <- renderPlot({
      shiny::validate(
        need(
          !is.null(values$mypca2go),
          "Please provide a pca2go object to the app or alternatively click on the action button - could take some time to compute live!"
        )
      )
      # if(is.null(pca2go))
      # return(ggplot() + annotate("text",label="Provide a pca2go object to the app",0,0) + theme_bw())
      res <- pcaplot(values$myrlt,intgroup = input$color_by,
                     ntop = attr(values$mypca2go,"n_genesforpca"),
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title=paste0("PCA on the samples - ",attr(values$mypca2go,"n_genesforpca"), " genes used")

      )
      res
    })


    output$dt_pchor_pos <- DT::renderDataTable({
      if(is.null(values$mypca2go)) return(datatable(NULL))
      goe <- values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pchor_neg <- DT::renderDataTable({
      if(is.null(values$mypca2go)) return(datatable(NULL))
      goe <- values$mypca2go[[paste0("PC",input$pc_x)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pcver_pos <- DT::renderDataTable({
      if(is.null(values$mypca2go)) return(datatable(NULL))
      goe <- values$mypca2go[[paste0("PC",input$pc_y)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pcver_neg <- DT::renderDataTable({
      if(is.null(values$mypca2go)) return(datatable(NULL))
      goe <- values$mypca2go[[paste0("PC",input$pc_y)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$enrichinfo <- renderPrint({
      cat("enrich info:\n")
      # str(goEnrichs)
      class(input$pc_x)
      head(values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]])
      class(datatable(values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]]))
    })






    ## from here on, multifac APP
    output$intro_multifac <- renderText({
      if(!is.null(values$mydds))
        shiny::validate(
          need(ncol(colData(values$mydds)) > 1,
               message = "To use this section, you need a dataset where more than one experimental factor is available.")
        )
      return("Refer to the Instructions section if you need help on using this section")
    })


    output$covar1 <- renderUI({
      # if(is.null(values$myrlt))
      # return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('covar1', label = 'Select factor 1: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = FALSE)
    })

    output$covar2 <- renderUI({
      # if(is.null(values$myrlt))
      # return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('covar2', label = 'Select factor 2: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = FALSE)
    })


    output$c1levels <- renderUI({
      if(is.null(input$covar1))
        return(NULL)
      fac1lev <- levels(colData(values$myrlt)[[input$covar1]])
      selectInput('covar1levels', label = 'Factor 1 available levels: ',
                  choices = c(NULL, fac1lev), selected = NULL,multiple = TRUE) # actually 2
    })

    output$c2levels <- renderUI({
      if(is.null(input$covar2))
        return(NULL)
      fac2lev <- levels(colData(values$myrlt)[[input$covar2]])
      selectInput('covar2levels', label = 'Factor 2 available levels: ',
                  choices = c(NULL, fac2lev), selected = NULL,multiple = TRUE) # 2 or more are allowed!
    })

    output$colnames1 <- renderUI({
      if(is.null(values$myrlt))
        return(NULL)
      if(is.null(input$covar1))
        return(NULL)
      if(is.null(input$covar2))
        return(NULL)

      fac1 <- input$covar1
      fac2 <- input$covar2

      fac1_touse <- input$covar1levels
      fac2_touse <- input$covar2levels

      preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
      preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
      presel <- intersect(preselected_fac1,preselected_fac2)
      mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

      presel1 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[1]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]

      selectInput('picksamples1', label = 'Combine samples from Factor1-Level1 in the selected order: ',
                  choices = c(NULL, presel1), selected = NULL,multiple = TRUE)
    })


    output$colnames2 <- renderUI({
      if(is.null(values$myrlt))
        return(NULL)
      if(is.null(input$covar1))
        return(NULL)
      if(is.null(input$covar2))
        return(NULL)

      fac1 <- input$covar1
      fac2 <- input$covar2

      fac1_touse <- input$covar1levels
      fac2_touse <- input$covar2levels

      preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
      preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
      presel <- intersect(preselected_fac1,preselected_fac2)
      mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

      presel2 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[2]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]

      selectInput('picksamples2', label = 'Combine samples from Factor1-Level2 in the selected order: ',
                  choices = c(NULL, presel2), selected = NULL,multiple = TRUE)
    })



    composedMat <- eventReactive( input$composemat, {
      exprmat <- t(assay(values$myrlt))
      exprmat <- exprmat[,rowSums(counts(values$mydds) > 5)>2]

      withProgress(message = "Composing the matrix...",
                   value = 0,
                   {
                     pcmat <- cbind(exprmat[input$picksamples1,],
                                    exprmat[input$picksamples2,])
                   })
      pcmat
    })


    obj3 <- reactive({

      pcmat <- composedMat()
      aval <- 0.3
      fac2pal <- alpha(c("green","red","blue","orange","violet"),aval) # 5 are enough

      # colData(values$myrlt)[input$covar2][rownames(pcmat),]
      max.type <- apply(pcmat[,1:(ncol(pcmat)/2)],2,which.max)

      fac2_col <- factor(colData(values$myrlt)[input$covar2][rownames(pcmat),],
                         levels=unique(as.character(colData(values$myrlt)[input$covar2][rownames(pcmat),])))
      tcol.justMax <- fac2pal[fac2_col][max.type]
      # tcol.justMax <- ifelse(max.type <= 4,"green",ifelse(max.type <= 8,"red",ifelse(max.type <= 12,"blue","orange")))

      max.type2 <- apply(pcmat[,((ncol(pcmat)/2)+1):ncol(pcmat)],2,which.max)
      # tcol2.justMax <- ifelse(max.type2 <= 4,alpha("green",aval),ifelse(max.type2 <= 8,alpha("red",aval),ifelse(max.type2 <= 12,alpha("blue",aval),alpha("orange",aval))))

      tcol2.justMax <- fac2pal[fac2_col][max.type2]

      # using the median across replicates
      celltypes <- gsub("_R.","",rownames(pcmat))

      tcol <- tcol.justMax
      tcol2 <- tcol2.justMax
      # pcmat
      return(list(pcmat,tcol,tcol2))
    })


    output$pcamultifac <- renderPlot({
      pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]
      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x
      # set.seed(11)
      # for (i in 1:ncol(pcx)) {
      #   pcx[,i] <- pcx[,i] + rnorm(nrow(pcx),sd=diff(range(pcx[,i]))/100)
      # }
      plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],xlim=range(pcx[,plot.index[1]]),ylim=range(pcx[,plot.index[2]]),pch=20,col=tcol,cex=0.3)#,type="n")
      #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
      lcol <- ifelse(tcol != tcol2,"black","grey")
      for (i in 1:gene.no) {
        lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
      }
      points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
      points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)

    })


    output$multifaczoom <- renderPlot({
      if(is.null(input$pcamultifac_brush)) return(NULL)
      pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]
      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x

      plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],
           pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],
           xlim=c(input$pcamultifac_brush$xmin,input$pcamultifac_brush$xmax),
           ylim=c(input$pcamultifac_brush$ymin,input$pcamultifac_brush$ymax),
           pch=20,col=tcol,cex=0.3)#,type="n")
      #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
      lcol <- ifelse(tcol != tcol2,"black","grey")
      for (i in 1:gene.no) {
        lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
      }
      points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
      points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)

    })


    #       output$plot_brushinfo <- renderPrint({
    #         cat("input$pcagenes_brush:\n")
    #         str(input$pcagenes_brush)
    #       })

    curData_brush_multifac <- reactive({
      pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]

      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x


      firstPCselected <- c(
        pcx[1:offset,plot.index[1]][1:gene.no],
        pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no])

      secondPCselected <- c(
        pcx[1:offset,plot.index[2]][1:gene.no],
        pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no]
      )

      pcspcs <- data.frame(firstPC=firstPCselected,secondPC=secondPCselected,geneID=colnames(pcmat))
      rownames(pcspcs) <- c(paste0(colnames(pcmat)[1:gene.no],"_WT"),
                            paste0(colnames(pcmat)[(gene.no+1):(2*gene.no)],"_G37"))

      if(!is.null(values$myannotation))
        pcspcs$geneName <- values$myannotation$gene_name[match(pcspcs$geneID,rownames(values$myannotation))]


      res <- brushedPoints(pcspcs, input$pcamultifac_brush,xvar="firstPC",yvar="secondPC",)
      res
    })



    output$pcamultifac_out <- DT::renderDataTable({
      datatable(curData_brush_multifac())
    })


    output$downloadData_brush_multifac <- downloadHandler(
      filename = function() { paste(input$brushedPoints_filename_multifac, '.csv', sep='') },
      content = function(file) {
        if(length(input$pcamultifac_out_rows_selected)){
          data <- curData_brush_multifac()[input$pcamultifac_out_rows_selected,]
        } else {
          data <- curData_brush_multifac()
        }
        write.csv(data, file, quote=FALSE)
      }
    )






    ## all download handlers
    output$downloadData_brush <- downloadHandler(
      filename = function() { paste(input$brushedPoints_filename, '.csv', sep='') },
      content = function(file) {
        if(length(input$pca_brush_out_rows_selected)){
          data <- curData_brush()[input$pca_brush_out_rows_selected,]
        } else {
          data <- curData_brush()
        }
        write.csv(data, file, quote=FALSE)
      }
    )

    output$downloadData_click <- downloadHandler(
      filename = function() { paste(input$clickedPoints_filename, '.csv', sep='') },
      content = function(file) {
        write.csv(curData_click(), file, quote=FALSE)
      }
    )

    output$download_samplesPca <- downloadHandler(
      filename = function() { input$filename_samplesPca },
      content = function(file) {
        ggsave(file, exportPlots$samplesPca, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_samplesScree <- downloadHandler(
      filename = function() { input$filename_samplesScree },
      content = function(file) {
        ggsave(file, exportPlots$samplesScree, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesPca <- downloadHandler(
      filename = function() { input$filename_genesPca },
      content = function(file) {
        ggsave(file, exportPlots$genesPca, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesZoom <- downloadHandler(
      filename = function() { input$filename_genesZoom },
      content = function(file) {
        ggsave(file, exportPlots$genesZoom, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesHeatmap <- downloadHandler(
      filename=function(){
        input$filename_genesHeatmap
      },
      content = function(file){
        pdf(file)
        brushedObject <- curData_brush()

        selectedGenes <- brushedObject$ids
        toplot <- assay(values$myrlt)[selectedGenes,]
        rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
        aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
        dev.off()
      }
    )


  }) # end of pcaExplorer(dds,rlt,countmatrix,coldata,pca2go,annotation)
  shinyApp(ui = myui, server = myserver)

}


footer <- function(){
  tags$div(
    class = "footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        hr(),
        "ideal is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
        ". ",br(),
        "Development of the ideal package is on ",
        tags$a(href="https://github.com/federicomarini/ideal", "GitHub")
      )
    )
  )
}
