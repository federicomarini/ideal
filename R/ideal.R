### commented out on 24.11.2016


# # ideal.R
#
# #' Title
# #'
# #' @param dds
# #' @param rlt
# #' @param countmatrix
# #' @param coldata
# #' @param pca2go
# #' @param annotation
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #'
# #'
# #'
# #'
# #'
# #'
# #'
# #'
# ideal<- function(
#   res_obj = NULL,
#   dds_obj = NULL,
#   annotation_obj = NULL){
#
#   if ( !requireNamespace('shiny',quietly = TRUE) ) {
#     stop("ideal requires 'shiny'. Please install it using
#          install.packages('shiny')")
#   }
#
#   # get modes and themes for the ace editor
#   modes <- shinyAce::getAceModes()
#   themes <- shinyAce::getAceThemes()
#
#   ## upload max 300mb files - can be changed if necessary
#   options(shiny.maxRequestSize=300*1024^2)
#
#   ## ------------------------------------------------------------------ ##
#   ##                          Define UI                                 ##
#   ## ------------------------------------------------------------------ ##
#
#
#   ideal_ui <- shinydashboard::dashboardPage(
#     dashboardHeader(
#       title = paste0("ideal - Interactive Differential Expression AnaLysis ",
#                      packageVersion("ideal")),
#       titleWidth = 900,
#
#       # task menu for saving state to environment or binary data
#       shinydashboard::dropdownMenu(type = "tasks",icon = icon("cog"),badgeStatus = "success",
#                                    notificationItem(
#                                      text = actionButton("task_exit_and_save","Exit ideal & save",
#                                                          class = "btn_no_border",
#                                                          onclick = "setTimeout(function(){window.close();}, 100); "),
#                                      icon = icon("sign-out"),status = "primary"),
#                                    menuItem(
#                                      text = downloadButton("task_state_save","Save State as .RData"))
#       )
#     ),
#
#     dashboardSidebar(
#       width = 280,
#
#       menuItem("App settings",icon = icon("cogs"),
#                uiOutput("color_by"),
#                uiOutput("available_genes")
#                ),
#       menuItem("Plot export settings", icon = icon("paint-brush"))
#     ),
#
#     dashboardBody(
#       ## Define output size and style of error messages
#       tags$head(
#         tags$style(HTML("
#                         .shiny-output-error-validation {
#                         font-size: 15px;
#                         color: forestgreen;
#                         text-align: center;
#                         }
#                         "))
#         ),
#
#       ## main structure of the body for the dashboard
#       tabBox(
#         width=12,
#
#         tabPanel(
#           "Loaded Data",icon = icon("upload"),
#           p("Preview on the uploaded data"),
#           verbatimTextOutput("printdds"),
#           verbatimTextOutput("printres")
#           ),
#         tabPanel(
#           "Instructions",  icon = icon("info-circle"),
#           includeMarkdown(system.file("extdata", "instructions.md",package = "ideal")),
#           footer()
#         ),
#         tabPanel(
#           "Overview - Tabular", icon = icon("table"),
#           DT::dataTableOutput("table_res"),
#           plotOutput("pvals_hist"),
#           plotOutput("logfc_hist")
#           ),
#
#         tabPanel(
#           "Data Overview", icon = icon("eye"),
#           uiOutput("choose_fac"),
#           uiOutput("fac1"),
#           uiOutput("fac2"),
#           verbatimTextOutput("diyres_summary"),
#           verbatimTextOutput("diyres")
#           ),
#         tabPanel(
#           "MA Plot", icon = icon("photo"),
#           headerPanel("MA plot interactive exploration"),
#           fluidRow(verbatimTextOutput("deb")),
#           fluidRow(column(4,
#                           h4("MA plot - Interactive!"),
#                           plotOutput('plotma', brush = 'ma_brush')),
#                    column(4,
#                           h4("Zoomed section"),
#                           plotOutput("mazoom",click= 'mazoom_click')),
#                    column(4,
#                           h4("Boxplot for the selected gene"),
#                           plotOutput("geneplot")
#                    )
#           ),
#           plotOutput("genefinder_plot"),
#           fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
#           fluidRow(
#             column(4,
#                    checkboxInput("rowscale",label = "Scale by rows",value = TRUE)),
#             column(4,
#                    checkboxInput("pseudocounts","use log2(1+counts)",value = TRUE))
#           ),
#           fluidRow(
#             column(6,
#                    plotOutput("heatbrush")
#             ),
#             column(6,
#                    d3heatmapOutput("heatbrushD3"))
#           )
#           ,
#           # fluidRow(dataTableOutput('ma_brush_out')),
#           fluidRow(dataTableOutput("ma_brush_out"))
#         ),
#         tabPanel(
#           "Gene Finder", icon = icon("crosshairs"),
#
#           fluidRow(
#             column(6,
#                    plotOutput("bp1")
#             ),
#             column(6,
#                    plotOutput("bp2"))
#           ),
#           fluidRow(
#             column(6,
#                    plotOutput("bp3")
#             ),
#             column(6,
#                    plotOutput("bp4"))
#           ),
#
#           plotOutput("ma_highlight"),
#           # verbatimTextOutput("d1"),
#           DT::dataTableOutput("table_combi")),
#         tabPanel(
#           "Report Editor",
#           icon = icon("pencil"),
#
#
#           fluidRow(
#             column(
#               width = 6,
#               box(
#                 title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
#                 radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = T),
#                 textInput("report_title", "Title: "),
#                 textInput("report_author", "Author: "),
#                 radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false")),
#                 radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
#                 selectInput("report_theme", "Theme", choices = list("Default" = "default", "Cerulean" = "cerulean",
#                                                                     "Journal" = "journal", "Flatly" = "flatly",
#                                                                     "Readable" = "readable", "Spacelab" = "spacelab",
#                                                                     "United" = "united", "Cosmo" = "cosmo")),
#                 radioButtons("report_echo", "Echo the commands in the output", choices = list("Yes" = "TRUE", "No" = "FALSE")))),
#             column(
#               width = 6,
#               box(
#                 title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
#                 checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
#                 conditionalPanel(
#                   "input.enableAutocomplete",
#                   wellPanel(
#                     checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
#                     checkboxInput("enableRCompletion", "R code completion", TRUE)
#                   )
#                 ),
#
#                 selectInput("mode", "Mode: ", choices=modes, selected="markdown"),
#                 selectInput("theme", "Theme: ", choices=themes, selected="solarized_light"))
#             )
#             # ,
#             # column( # kept for debugging purposes!
#             #   width = 6,
#             #   verbatimTextOutput("loadedRmd")
#             # )
#           ),
#           fluidRow(
#             column(3,
#                    actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
#             ),
#             column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success"))
#           ),
#
#           tabBox(
#             width = NULL,
#             id="report_tabbox",
#             tabPanel("Report preview",
#                      icon = icon("file-text"),
#                      htmlOutput("knitDoc")
#             ),
#
#             tabPanel("Edit report",
#                      icon = icon("pencil-square-o"),
#                      aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
#                                value=readLines(system.file("extdata", "irt.Rmd",package = "ideal")),
#                                height="800px"))
#           )
#         ),
#         tabPanel(
#           "About", icon = icon("institution"),
#           includeMarkdown(system.file("extdata", "about.md",package = "ideal")),
#           hr(),
#           #             shiny::verbatimTextOutput("showuploaded1"),
#           #             shiny::verbatimTextOutput("showuploaded2"),
#           #             shiny::verbatimTextOutput("showuploaded3"),
#           #             shiny::verbatimTextOutput("showuploaded4"),
#
#           h4("Session Info"),
#           verbatimTextOutput("sessioninfo"),
#           footer()
#         )
#         ,tabPanel(
#           "devel", icon = icon("github"),
#           verbatimTextOutput("debugihw"),
#
#           plotOutput("ihwp1"),
#           plotOutput("ihwp2"),
#           plotOutput("ihwp3"),
#           plotOutput("ihwp4"),
#           plotOutput("ihwp5"),
#           plotOutput("ihwp6"),
#           plotOutput("ihwp7"),
#           plotOutput("ihwp8")
#
#           )
#
#       )
#
#     ),
#     skin="blue"
#
#   )
#
#
#   ideal_server <- shinyServer(function(input, output, session) {
#
#     ## placeholder for the figures to export
#     # exportPlots <- reactiveValues(
#       # expfig_fig1 <- NULL
#     # )
#
#     values <- reactiveValues()
#     values$res_obj <- res_obj
#     values$dds_obj <- dds_obj
#     values$annotation_obj <- annotation_obj
#
#
#     # if i want to focus a little more on the ihw object
#     values$ihwres <- NULL
#
#
#     ihwres <- reactive({
#       de_res <- as.data.frame(res_obj)
#       ihw_res <- ihw(pvalue ~ baseMean,  data = de_res, alpha = 0.1)
#       ihw_res
#     })
#
#     output$debugihw <- renderPrint({
#       ihw_res # replace then here with ihwres()
#     })
#
#     output$ihwp1 <- renderPlot({
#       plot(ihw_res)
#     })
#
#     output$ihwp2 <- renderPlot({
#       plot(ihw_res, what = "decisionboundary")
#     })
#
#     output$ihwp3 <- renderPlot({
#       gg <- ggplot(as.data.frame(ihw_res), aes(x = pvalue, y = adj_pvalue, col = group)) +
#         geom_point(size = 0.25) +
#         scale_colour_hue(l = 70, c = 150, drop = FALSE)
#
#       gg
#     })
#
#     output$ihwp4 <- renderPlot({
#       de_res <- res_airway ## CHANGE!!
#       de_res <- na.omit(de_res)
#       de_res$geneid <- as.numeric(gsub("ENSG[+]*", "", rownames(de_res)))
#
#       # set up data frame for plotting
#       df <- rbind(data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$baseMean)/nrow(de_res),
#                              covariate_type="base mean"),
#                   data.frame(pvalue = de_res$pvalue, covariate = rank(de_res$geneid)/nrow(de_res),
#                              covariate_type="gene id"))
#
#       ggplot(df, aes(x=covariate, y = -log10(pvalue))) +
#         geom_hex(bins = 100) +
#         facet_grid( . ~ covariate_type)
#
#     })
#
#
#
#     output$ihwp5 <- renderPlot({
#       de_res <- as.data.frame(res_airway) ## CHANGE!!
#       de_res <- na.omit(de_res)
#       ggplot(de_res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
#     })
#
#     output$ihwp6 <- renderPlot({
#       de_res <- as.data.frame(res_airway) ## CHANGE!!
#       de_res <- na.omit(de_res)
#       de_res$baseMean_group <- groups_by_filter(de_res$baseMean, 8)
#
#       ggplot(de_res, aes(x=pvalue)) +
#         geom_histogram(binwidth = 0.025, boundary = 0) +
#         facet_wrap( ~ baseMean_group, nrow = 2)
#     })
#
#     output$ihwp7 <- renderPlot({
#       de_res <- as.data.frame(res_airway) ## CHANGE!!
#       de_res <- na.omit(de_res)
#       de_res$baseMean_group <- groups_by_filter(de_res$baseMean, 8)
#       ggplot(de_res, aes(x = pvalue, col = baseMean_group)) + stat_ecdf(geom = "step")
#     })
#
#     output$ihwp8 <- renderPlot({
#       de_res <- as.data.frame(res_airway) ## CHANGE!!
#       de_res <- na.omit(de_res)
#
#       de_res$lfc_group <- groups_by_filter(abs(de_res$log2FoldChange),8)
#
#       ggplot(de_res, aes(x = pvalue)) +
#         geom_histogram(binwidth = 0.025, boundary = 0) +
#         facet_wrap( ~ lfc_group, nrow=2)
#     })
#
#
#
#
#     object <- res_obj
#     obj2 <- dds_obj
#     res_object <- res_obj
#
#
#     if(!is.null(res_obj) & !is.null(dds_obj)) {
#
#       normCounts <- as.data.frame(counts(estimateSizeFactors(dds_obj),normalized=T))
#       normCounts$id <- rownames(normCounts)
#       res_df <- FMmisc::deseqresult2tbl(res_obj)
#
#       combi_obj <- dplyr::inner_join(res_df,normCounts,by="id")
#       combi_obj$symbol <- annotation_obj$gene_name[match(combi_obj$id,annotation_obj$gene_id)]
#     }
#
#
#     output$color_by <- renderUI({
#       if(is.null(obj2))
#         return(NULL)
#       poss_covars <- names(colData(obj2))
#       selectInput('color_by', label = 'Group/color by: ',
#                   choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
#     })
#
#
#     output$available_genes <- renderUI({
#       if("symbol" %in% names(res_obj)) {
#         selectInput("avail_symbols", label = "Select the gene(s) of interest",
#                     choices = as.character(res_obj$symbol), selected = NULL, multiple = TRUE)
#       } else { # else use the rownames as identifiers
#         selectInput("avail_ids", label = "Select the gene(s) of interest - ids",
#                     choices = rownames(res_obj), selected = NULL, multiple = TRUE)
#       }
#     })
#
#
#
#     design_factors <- rev(attributes(terms.formula(design(dds_obj)))$term.labels)
#
#     output$choose_fac <- renderUI({
#       selectInput("choose_expfac",label = "choose the experimental factor to build the contrast upon",
#                   choices = c("",design_factors), selected = "")
#     })
#
#     output$fac1 <- renderUI({
#       shiny::validate(
#         need(input$choose_expfac!="",
#              "Please select the factor to build the contrast upon - contrast1"
#         )
#       )
#       fac1 <- input$choose_expfac
#       fac1_vals <- colData(dds_obj)[,fac1]
#       fac1_levels <- levels(fac1_vals)
#
#       selectInput("fac1_c1","c1",choices = c("",fac1_levels), selected = "")
#       # selectInput("fac1_c2","c2",choices = fac1_levels)
#     })
#
#     output$fac2 <- renderUI({
#       shiny::validate(
#         need(input$choose_expfac!="",
#              "Please select the factor to build the contrast upon - contrast2"
#         )
#       )
#       fac1 <- input$choose_expfac
#       fac1_vals <- colData(dds_obj)[,fac1]
#       fac1_levels <- levels(fac1_vals)
#
#       # selectInput("fac1_c1","c1",choices = fac1_levels)
#       selectInput("fac1_c2","c2",choices = c("",fac1_levels), selected = "")
#     })
#
#
#
#     output$diyres <- renderPrint({
#       shiny::validate(
#         need(input$choose_expfac!="" & input$fac1_c1 != "" & input$fac1_c2 != "" & input$fac1_c1 != input$fac1_c2 ,
#              "Please select the factor to build the contrast upon, and two different levels to build the contrast"
#         )
#       )
#
#       results(dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2))
#     })
#
#     output$diyres_summary <- renderPrint({
#       shiny::validate(
#         need(input$choose_expfac!="" & input$fac1_c1 != "" & input$fac1_c2 != "" & input$fac1_c1 != input$fac1_c2 ,
#              "Please select the factor to build the contrast upon, and two different levels to build the contrast"
#         )
#       )
#       summary(results(dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2)))
#     })
#
#
#
#     output$printdds <- renderPrint({
#
#       shiny::validate(
#         need(!is.null(dds_obj),
#              "Please provide a count matrix/dds object"
#         )
#       )
#
#       dds_obj
#
#     })
#
#     output$printres <- renderPrint({
#
#       shiny::validate(
#         need(!is.null(res_obj),
#              "Please provide a DESeqResults object"
#         )
#       )
#
#       print(sub(".*p-value: (.*)","\\1",mcols(res_obj, use.names=TRUE)["pvalue","description"]))
#       summary(res_obj,alpha = 0.05) # use fdr shiny widget
#
#     })
#
#
#     output$table_res <- DT::renderDataTable({
#       as.data.frame(res_obj[order(res_obj$padj),])[1:500,]
#
#     })
#
#
#     output$pvals_hist <- renderPlot({
#
#       res_df <- as.data.frame(res_obj)
#       p <- ggplot(res_df, aes(pvalue)) +
#         geom_histogram(binwidth = 0.01) + theme_bw()
#
#       p
#
#     })
#
#     output$logfc_hist <- renderPlot({
#
#       res_df <- as.data.frame(res_obj)
#       p <- ggplot(res_df, aes(log2FoldChange)) +
#         geom_histogram(binwidth = 0.1) + theme_bw()
#
#       p
#
#     })
#
#
#     output$dds_design <- renderPrint({
#       design(dds_obj)
#     })
#
#     output$res_names <- renderPrint({
#       resultsNames(dds_obj)
#     })
#
#
#
#     output$explore_res <- renderPrint({
#       expfac <- attributes(terms.formula(design(dds_obj)))$term.labels
#       expfac # plus, support up to four factors that are either there or not according to the length
#     })
#
#
#
#
#     output$plotma <- renderPlot({
#       plot_ma(object,annotation_obj = annotation_obj)
#     })
#
#     output$mazoom <- renderPlot({
#       if(is.null(input$ma_brush)) return(ggplot() + annotate("text",label="click and drag to zoom in",0,0) + theme_bw())
#
#       plot_ma(object,annotation_obj = annotation_obj) + xlim(input$ma_brush$xmin,input$ma_brush$xmax) + ylim(input$ma_brush$ymin,input$ma_brush$ymax) + geom_text(aes(label=genename),size=3,hjust=0.25, vjust=-0.75)
#     })
#
#
#     output$ma_highlight <- renderPlot({
#       if("symbol" %in% names(res_obj)) {
#         plot_ma_highlight(object,
#                         intgenes = input$avail_symbols,annotation_obj = annotation_obj)
#       } else {
#         plot_ma_highlight(object,
#                           intgenes = input$avail_ids,annotation_obj = annotation_obj)
#       }
#     })
#
#
#     curData <- reactive({
#       mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < 0.10),ID=rownames(object))
#       mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
#       # mama$yesorno <- ifelse(mama$isDE,"yes","no")
#       mama$yesorno <- ifelse(mama$isDE,"red","black")
#       mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
#
#
#       res <- brushedPoints(mama, input$ma_brush,xvar="logmean",yvar="lfc")
#       res
#     })
#
#
#     curDataClick <- reactive({
#       mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < 0.10),ID=rownames(object))
#       mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
#       # mama$yesorno <- ifelse(mama$isDE,"yes","no")
#       mama$yesorno <- ifelse(mama$isDE,"red","black")
#       mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
#
#
#       res <- nearPoints(mama, input$mazoom_click,threshold = 20, maxpoints = 1,
#                         addDist = TRUE)
#       res
#     })
#
#
#
#
#     output$ma_brush_out <- renderDataTable({
#       datatable(curData(),options=list(pageLength=100))
#
#       # alternative:
#
#       # curData()
#       # and...
#       # options=list(pageLength=100)) goes to renderDataTable
#       #         if (!is.null(input$ma_brush)) {
#       #           res <- enclosed_brush(mama, input$ma_brush)
#       #           # rv_ma$highlight_vars <- res
#       #         }  else {
#       #           res <- NULL
#       #         }
#       #
#       #         # TODO: total hack -- fix this correctly eventually
#       # #         if (is(res, 'data.frame')) {
#       # #           res <- dplyr::rename(res,
#       # #                                mean = mean_obs,
#       # #                                var = var_obs,
#       # #                                tech_var = sigma_q_sq,
#       # #                                final_sigma_sq = smooth_sigma_sq_pmax)
#       # #         }
#       #
#       #         res
#     })
#
#
#     output$heatbrush <- renderPlot({
#       if((is.null(input$ma_brush))|is.null(obj2)) return(NULL)
#
#       brushedObject <- curData()
#
#       selectedGenes <- brushedObject$ID
#       toplot <- assay(obj2)[selectedGenes,]
#       rownames(toplot) <- annotation_obj$gene_name[match(rownames(toplot),rownames(annotation_obj))]
#
#       if(input$pseudocounts) toplot <- log2(1+toplot)
#
#       if(input$rowscale) toplot <- pheatmap:::scale_mat(toplot,"row")
#
#       pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
#
#
#     })
#
#
#     output$heatbrushD3 <- renderD3heatmap({
#       if((is.null(input$ma_brush))|is.null(obj2)) return(NULL)
#
#       brushedObject <- curData()
#
#       selectedGenes <- brushedObject$ID
#       toplot <- assay(obj2)[selectedGenes,]
#       rownames(toplot) <- annotation_obj$gene_name[match(rownames(toplot),rownames(annotation_obj))]
#       mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding
#       if(input$pseudocounts) toplot <- log2(1+toplot)
#       if(input$rowscale) toplot <- pheatmap:::scale_mat(toplot,"row")
#
#       d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
#
#
#     })
#
#
#     output$deb <- renderPrint({
#       # curDataClick()
#       selectedGene <- curDataClick()$ID
#       #         selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
#       #         # plotCounts(dds_cleaner,)
#       #         genedata <- plotCounts(dds_cleaner,gene=selectedGene,intgroup = "condition",returnData = T)
#       #         genedata
#       # str(as.character(selectedGene))
#       selectedGene
#     })
#
#
#     output$geneplot <- renderPlot({
#
#       # if(length(input$color_by_G)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
#       if(is.null(input$ma_brush)) return(NULL)
#
#       if(is.null(input$mazoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())
#
#       selectedGene <- as.character(curDataClick()$ID)
#       selectedGeneSymbol <- annotation_obj$gene_name[match(selectedGene,rownames(annotation_obj))]
#       # plotCounts(dds_cleaner,)
#       genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by,returnData = T)
#
#       # onlyfactors <- genedata[match(input$color_by_G,colnames(genedata))]
#       # onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]
#       onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
#
#       ## intgroup can be a vector of factors. then i need interactions of the two factors
#       # plotCounts(ddsmf_global,gene="ENSMUSG00000026080",intgroup = c("tissue","condition"),returnData = T) -> dh
#       # dh$f1f2 <- interaction(dh$tissue,dh$condition)
#       # dh  %>% ggplot(aes(x=f1f2,y=count,fill=f1f2)) + geom_boxplot()
#
#
#       ## TODO: make the intgroup/colr by also somehow multiple-selectable?
#
#       # genedata$plotby <- lapply(1:ncol(onlyfactors),function(arg) onlyfactors[,arg]) %>% interaction()
#       genedata$plotby <- interaction(onlyfactors)
#
#
#
#       ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + scale_y_log10(name="Normalized counts") + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
#       # exportPlots$genesZoom <- res
#       # res
#     })
#
#
#     output$genefinder_plot <- renderPlot({
#
#       shiny::validate(
#         need(
#           length(input$color_by)>0,
#           "Select an experimental factor in the Group/color by element in the sidebar"
#         )
#       )
#
#       if(is.null(input$ma_brush)) return(NULL)
#
#       if(is.null(input$mazoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())
#
#
#       selectedGene <- as.character(curDataClick()$ID)
#       selectedGeneSymbol <- annotation_obj$gene_name[match(selectedGene,annotation_obj$gene_id)]
#
#       p <- ggplotCounts(dds_obj, selectedGene, intgroup = input$color_by,annotation_obj=annotation_obj)
#
#       p
#
#
#
#     })
#
#
#     cur_combires <- reactive({
#
#
#       if("symbol" %in% names(res_obj)) {
#         sel_genes <- input$avail_symbols
#         sel_genes_ids <- annotation_obj$gene_id[match(sel_genes,annotation_obj$gene_name)]
#       } else {
#         sel_genes_ids <- input$avail_ids
#       }
#
#       if(length(sel_genes) > 0) {
#         combi_obj[match(sel_genes_ids,combi_obj$id),]
#       } else {
#         combi_obj
#       }
#
#
#     })
#
#
#     # output$d1 <- renderPrint({
#     #   length(cur_combires())
#     # })
#
#     output$table_combi <- DT::renderDataTable({
#       datatable(cur_combires(),options = list(scrollX=TRUE))
#     })
#
#
#
#     output$bp1 <- renderPlot({
#       shiny::validate(
#         need(
#           length(input$color_by)>0,
#           "Select an experimental factor in the Group/color by element in the sidebar"
#         )
#       )
#       shiny::validate(
#         need(
#           (length(input$avail_symbols)>0 | length(input$avail_ids)>0),
#           "Select at least a gene to plot"
#         )
#       )
#       if(length(input$avail_symbols)>0) {
#         # got the symbol, look for the id
#         mysym <- input$avail_symbols[1]
#         myid <- annotation_obj$gene_id[match(mysym, annotation_obj$gene_name)]
#       } else {
#         myid <- input$avail_ids[1]
#         # make it optional if annot is available
#         if(!is.null(annotation_obj)) {
#           mysim <- annotation_obj$gene_name[match(myid, annotation_obj$gene_id)]
#         } else {
#           mysim <- ""
#         }
#       }
#       p <- ggplotCounts(dds_obj, myid, intgroup = input$color_by,annotation_obj=annotation_obj)
#       p
#     })
#
#     output$bp2 <- renderPlot({
#       shiny::validate(
#         need(
#           length(input$color_by)>0,
#           "Select an experimental factor in the Group/color by element in the sidebar"
#         )
#       )
#       shiny::validate(
#         need(
#           (length(input$avail_symbols)>1 | length(input$avail_ids)>1),
#           "Select at least a second gene to plot"
#         )
#       )
#       if(length(input$avail_symbols)>0) {
#         # got the symbol, look for the id
#         mysym <- input$avail_symbols[2]
#         myid <- annotation_obj$gene_id[match(mysym, annotation_obj$gene_name)]
#       } else {
#         myid <- input$avail_ids[2]
#         # make it optional if annot is available
#         if(!is.null(annotation_obj)) {
#           mysim <- annotation_obj$gene_name[match(myid, annotation_obj$gene_id)]
#         } else {
#           mysim <- ""
#         }
#       }
#       p <- ggplotCounts(dds_obj, myid, intgroup = input$color_by,annotation_obj=annotation_obj)
#       p
#     })
#
#     output$bp3 <- renderPlot({
#       shiny::validate(
#         need(
#           length(input$color_by)>0,
#           "Select an experimental factor in the Group/color by element in the sidebar"
#         )
#       )
#       shiny::validate(
#         need(
#           (length(input$avail_symbols)>2 | length(input$avail_ids)>2),
#           "Select at least a third gene to plot"
#         )
#       )
#       if(length(input$avail_symbols)>0) {
#         # got the symbol, look for the id
#         mysym <- input$avail_symbols[3]
#         myid <- annotation_obj$gene_id[match(mysym, annotation_obj$gene_name)]
#       } else {
#         myid <- input$avail_ids[3]
#         # make it optional if annot is available
#         if(!is.null(annotation_obj)) {
#           mysim <- annotation_obj$gene_name[match(myid, annotation_obj$gene_id)]
#         } else {
#           mysim <- ""
#         }
#       }
#       p <- ggplotCounts(dds_obj, myid, intgroup = input$color_by,annotation_obj=annotation_obj)
#       p
#     })
#
#     output$bp4 <- renderPlot({
#       shiny::validate(
#         need(
#           length(input$color_by)>0,
#           "Select an experimental factor in the Group/color by element in the sidebar"
#         )
#       )
#       shiny::validate(
#         need(
#           (length(input$avail_symbols)>3 | length(input$avail_ids)>3),
#           "Select at least a fourth gene to plot"
#         )
#       )
#       if(length(input$avail_symbols)>0) {
#         # got the symbol, look for the id
#         mysym <- input$avail_symbols[4]
#         myid <- annotation_obj$gene_id[match(mysym, annotation_obj$gene_name)]
#       } else {
#         myid <- input$avail_ids[4]
#         # make it optional if annot is available
#         if(!is.null(annotation_obj)) {
#           mysim <- annotation_obj$gene_name[match(myid, annotation_obj$gene_id)]
#         } else {
#           mysim <- ""
#         }
#       }
#       p <- ggplotCounts(dds_obj, myid, intgroup = input$color_by,annotation_obj=annotation_obj)
#       p
#     })
#
#
#
#
#
#
#     ## REPORT EDITOR
#     ### yaml generation
#     rmd_yaml <- reactive({
#       paste0("---",
#              "\ntitle: '", input$report_title,
#              "'\nauthor: '", input$report_author,
#              "'\ndate: '", Sys.Date(),
#              "'\noutput:\n  html_document:\n    toc: ", input$report_toc, "\n    number_sections: ", input$report_ns, "\n    theme: ", input$report_theme, "\n---\n\n",collapse = "\n")
#     })
#
#
#     # rmd_full <- reactive({
#     #   paste0(rmd_yaml(),"\n",
#     #          readLines("reportTemplate.Rmd"))
#     # })
#     # output$loadedRmd <- renderPrint({
#     #   # rmd_yaml() # or rmd_full()
#     #   paste0(
#     #     # rmd_yaml(),
#     #     paste0(readLines("reportTemplate.Rmd"),collapse = "\n"))
#     #   # head(paste0(rmd_yaml(),
#     #   # readLines("reportTemplate.Rmd")),collapse="\n")
#     # })
#
#     ### loading report template
#     # update aceEditor module
#     observe({
#       # loading rmd report from disk
#       inFile <- system.file("extdata", "irt.Rmd",package = "ideal")
#
#       isolate({
#         if(!is.null(inFile) && !is.na(inFile)) {
#
#           rmdfilecontent <- paste0(readLines(inFile),collapse="\n")
#
#           shinyAce::updateAceEditor(session, "acereport_rmd", value = rmdfilecontent)
#         }
#       })
#     })
#
#
#     ### ace editor options
#     observe({
#       autoComplete <- if(input$enableAutocomplete) {
#         if(input$enableLiveCompletion) "live" else "enabled"
#       } else {
#         "disabled"
#       }
#
#       updateAceEditor(session, "acereport_rmd", autoComplete = autoComplete,theme=input$theme, mode=input$mode)
#       # updateAceEditor(session, "plot", autoComplete = autoComplete)
#     })
#
#     #Enable/Disable R code completion
#     rmdOb <- aceAutocomplete("acereport_rmd")
#     observe({
#       if(input$enableRCompletion) {
#         rmdOb$resume()
#       } else {
#         rmdOb$suspend()
#       }
#     })
#
#     ## currently not working as I want with rmarkdown::render, but can leave it like this - the yaml will be taken in the final version only
#     output$knitDoc <- renderUI({
#       input$updatepreview_button
#       return(isolate(HTML(knit2html(text = input$acereport_rmd, fragment.only = TRUE, quiet = TRUE))))
#     })
#
#     # Generate and Download module
#     output$saveRmd <- downloadHandler(
#       filename = function() {
#         if(input$rmd_dl_format == "rmd") {
#           "report.Rmd"
#         } else {
#           "report.html"
#         }
#       },
#       content = function(file) {
#
#         # knit2html(text = input$rmd, fragment.only = TRUE, quiet = TRUE))
#
#         tmp_content <-
#           paste0(rmd_yaml(),
#                  input$acereport_rmd,collapse = "\n")
#         # input$acereport_rmd
#         if(input$rmd_dl_format == "rmd") {
#           cat(tmp_content,file=file,sep="\n")
#         } else {
#           # write it somewhere too keeping the source
#           # tmpfile <- tempfile()
#           # file.create(tmpfile)
#           # fileConn<- file(tempfile())
#           # writeLines(tmp_content, fileConn)
#           # close(fileConn)
#           if(input$rmd_dl_format == "html") {
#             cat(tmp_content,file="ideal_tempreport.Rmd",sep="\n")
#             rmarkdown::render(input = "ideal_tempreport.Rmd",
#                               output_file = file,
#                               # fragment.only = TRUE,
#                               quiet = TRUE)
#           }
#         }
#       })
#
#
#
#
#
#
#
#     ## STATE SAVING
#     ### to environment
#     observe({
#       if(is.null(input$task_exit_and_save) || input$task_exit_and_save ==0 ) return()
#
#       # quit R, unless you are running an interactive session
#       if(interactive()) {
#         # flush input and values to the environment in two distinct objects (to be reused later?)
#         isolate({
#           assign(paste0("ideal_inputs_",
#                         gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
#                  reactiveValuesToList(input), envir = .GlobalEnv)
#           assign(paste0("ideal_values_",
#                         gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
#                  reactiveValuesToList(values), envir = .GlobalEnv)
#           stopApp("ideal closed, state successfully saved to global R environment.")
#         })
#       } else {
#         stopApp("ideal closed")
#         q("no")
#       }
#     })
#
#     ### to binary data
#     saveState <- function(filename) {
#       isolate({
#         LiveInputs <- reactiveValuesToList(input)
#         # values[names(LiveInputs)] <- LiveInputs
#         r_data <- reactiveValuesToList(values)
#         save(LiveInputs, r_data , file = filename)
#       })
#     }
#
#     output$task_state_save <- downloadHandler(
#       filename = function() {
#         paste0("idealState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".RData")
#       },
#       content = function(file) {
#         saveState(file)
#       }
#     )
#
#
#
#
#   })
#
#   shinyApp(ui = ideal_ui, server = ideal_server)
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
# ############################# helper funcs #################################
#
# plot_ma <- function(object, which_beta, which_model = 'full',
#                     sig_level = 0.10,
#                     point_alpha = 0.2,
#                     sig_color = 'red',
#                     annotation_obj = NULL, # TODO: add a check, if not available skip this part
#                     hlines = NULL
#
# ) {
#   mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
#   mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
#   mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
#   # mama$yesorno <- ifelse(mama$isDE,"yes","no")
#   mama$yesorno <- ifelse(mama$isDE,"red","black")
#
#   p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))
#
#   if(!is.null(hlines)) {
#     p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) + geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
#   }
#
#
#   p <- p + geom_point(alpha = point_alpha)
#   p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-3,3)
#   # p <- p + xlab('mean( log( counts + 0.5 ) )')
#   # p <- p + ylab(paste0('beta: ', which_beta))
#
#   #   if (!is.null(highlight)) {
#   #     suppressWarnings({
#   #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
#   #     })
#   #     if (nrow(highlight) > 0) {
#   #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
#   #     } else {
#   #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
#   #     }
#   #   }
#
#
#   p <- p + theme_bw()
#
#
#   p
# }
#
#
#
# plot_ma_highlight <- function(object, which_beta, which_model = 'full',
#                     sig_level = 0.10,
#                     point_alpha = 0.2,
#                     sig_color = "red",
#                     intgenes = NULL,
#                     intgenes_color = "steelblue",
#                     labels_intgenes = TRUE,
#                     annotation_obj = annotation_obj # SAME AS IN FUNC ABOVE
#
# ) {
#   mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
#   mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
#   mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
#   # mama$yesorno <- ifelse(mama$isDE,"yes","no")
#   mama$yesorno <- ifelse(mama$isDE,"red","black")
#
#   p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))
#   p <- p + geom_point(alpha = point_alpha)
#   p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-4,4) # or even remove the y limit
#   # p <- p + xlab('mean( log( counts + 0.5 ) )')
#   # p <- p + ylab(paste0('beta: ', which_beta))
#
#   #   if (!is.null(highlight)) {
#   #     suppressWarnings({
#   #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
#   #     })
#   #     if (nrow(highlight) > 0) {
#   #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
#   #     } else {
#   #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
#   #     }
#   #   }
#
#   if(!is.null(intgenes)){
#     # check that genes are in the results object:
#     # will check the id as rownames or an explicit symbol column
#
#
#     # now here for the symbol
#     res_df <- as.data.frame(object)
#     res_df$logmean <- log10(res_df$baseMean)
#     df_intgenes <- res_df[res_df$symbol %in% intgenes,]
#
#     p <- p + geom_point(data = df_intgenes,aes(logmean, log2FoldChange), color = intgenes_color, size = 4)
#
#     if(labels_intgenes) {
#       p <- p + geom_text(data = df_intgenes,aes(logmean, log2FoldChange,label=symbol),
#                          color = intgenes_color, size=5,hjust = 0, nudge_x = 0.2)
#     }
#
#   }
#
#
#   p <- p + theme_bw()
#
#   p
# }
#
#
# ggplotCounts <- function(dds,gene,intgroup="condition",annotation_obj=annotation_obj,...){
#   df <- plotCounts(dds,gene,intgroup,returnData = T,...)
#   df$sampleID <- rownames(df)
#   genesymbol <- annotation_obj$gene_name[match(gene,annotation_obj$gene_id)]
#   jittered_df <- df
#   # jittered_df$conditionj <- jitter(as.numeric(factor(jittered_df$condition)))
#   jittered_df$countj <- jitter(jittered_df$count)
#
#   onlyfactors <- df[,match(intgroup,colnames(df))]
#   df$plotby <- interaction(onlyfactors)
#
#
#   p <- df %>%
#     ggplot(aes(x=plotby,y=count,col=plotby)) +
#     geom_boxplot(outlier.shape = NA) +
#     # geom_text(data = jittered_df,aes(x=conditionj,y=countj,label=sampleID)) +
#     geom_text(aes(label=sampleID),hjust=-.1,vjust=0) +
#     labs(title=paste0("Normalized counts for ",genesymbol," - ",gene)) +
#     scale_x_discrete(name="") +
#     geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) +
#     scale_color_discrete(name="Experimental\nconditions") +
#     scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.1,NA)) +
#     theme_bw()
#
#   p
#
# }
#
#
# # combineTogether <- function(normCounts,resuTable,anns) {
# #   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
# #   anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
# #   combinedCountsAndRes$Description <- anns2$description
# #   return(combinedCountsAndRes)
# # }
# #
# #
# # combine_resucounts <- function(normCounts,resuTable) {
# #   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
# #   # anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
# #   # combinedCountsAndRes$Description <- anns2$description
# #   return(combinedCountsAndRes)
# # }
# #
#
#
#
#
#
#
# footer <- function(){
#   tags$div(
#     class = "footer",
#     style = "text-align:center",
#     tags$div(
#       class = "foot-inner",
#       list(
#         hr(),
#         "ideal is a project developed by Federico Marini in the Bioinformatics division of the ",
#         tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
#         ". ",br(),
#         "Development of the ideal package is on ",
#         tags$a(href="https://github.com/federicomarini/ideal", "GitHub")
#       )
#     )
#   )
# }
