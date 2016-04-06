library("d3heatmap")
library("shiny")
library("ggvis")
library("DT")

load("../../linuxHome/032-ruf-macrophages/cm2.RData")

plot_ma <- function(object, which_beta, which_model = 'full',
                    sig_level = 0.10,
                    point_alpha = 0.2,
                    sig_color = 'red'

) {
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
  mama$genename <- cm2$fromgtf[match(mama$ID,rownames(cm2))]
  mama$logmean <- log(mama$mean) # TO ALLOW FOR BRUSHING!!
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")

  p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))
  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-3,3)
  # p <- p + xlab('mean( log( counts + 0.5 ) )')
  # p <- p + ylab(paste0('beta: ', which_beta))

  #   if (!is.null(highlight)) {
  #     suppressWarnings({
  #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
  #     })
  #     if (nrow(highlight) > 0) {
  #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
  #     } else {
  #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
  #     }
  #   }

  p
}

enclosed_brush <- function(df, brush) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  xvar <- brush$mapping$x
  yvar <- brush$mapping$y

  xbool <- brush$xmin <= df[[xvar]] & df[[xvar]] <= brush$xmax
  ybool <- brush$ymin <= df[[yvar]] & df[[yvar]] <= brush$ymax

  df[xbool & ybool,]
}

ggmama <- function(object){
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < 0.05),ID=rownames(object))
  mama$genename <- cm2$fromgtf[match(mama$ID,rownames(cm2))]
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")
  # mama$logmean <- log(mama$mean)

  library(ggvis)
  malabs <- function(data){
    if(is.null(data)) return(NULL)
    data$genename
  }
  mama = subset(mama, mean != 0)

  mama %>% ggvis( ~log(mean), ~lfc,key:=~genename,fill:=~yesorno) %>% layer_points(size:=20,opacity:=0.5,size.hover := 100,fill.hover:="steelblue") %>% add_tooltip(malabs, "hover") %>% scale_numeric("y",domain = c(-2,2),clamp=T)
}


## TODOS
# parameters instead of fixed values
# some other tweaks
# try one run from the scratch
# list whatever is required on top of the pca_SUPALIVE function
### EMBED THIS INTO INTERACTIVE HTML REPORTS? e.g. a la maplots!!!!



ma_SUPALIVE <- function(object,obj2){



  serser <-

    shinyServer(function(input, output, session) {


      reactive({
        ggmama(object)
      }) %>%
        bind_shiny("p")

      output$plotma <- renderPlot({
        plot_ma(object)
      })

      output$mazoom <- renderPlot({
        if(is.null(input$ma_brush)) return(ggplot() + annotate("text",label="click and drag to zoom in",0,0) + theme_bw())

        plot_ma(object) + xlim(input$ma_brush$xmin,input$ma_brush$xmax) + ylim(input$ma_brush$ymin,input$ma_brush$ymax) + geom_text(aes(label=genename),size=3,hjust=0.25, vjust=-0.75)
      })


      curData <- reactive({
        mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < 0.10),ID=rownames(object))
        mama$genename <- cm2$fromgtf[match(mama$ID,rownames(cm2))]
        # mama$yesorno <- ifelse(mama$isDE,"yes","no")
        mama$yesorno <- ifelse(mama$isDE,"red","black")
        mama$logmean <- log(mama$mean) # TO ALLOW FOR BRUSHING!!


        res <- brushedPoints(mama, input$ma_brush,xvar="logmean",yvar="lfc")
        res
      })


      curDataClick <- reactive({
        mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < 0.10),ID=rownames(object))
        mama$genename <- cm2$fromgtf[match(mama$ID,rownames(cm2))]
        # mama$yesorno <- ifelse(mama$isDE,"yes","no")
        mama$yesorno <- ifelse(mama$isDE,"red","black")
        mama$logmean <- log(mama$mean) # TO ALLOW FOR BRUSHING!!


        res <- nearPoints(mama, input$mazoom_click,threshold = 20, maxpoints = 1,
                          addDist = TRUE)
        res
      })




      output$ma_brush_out <- renderDataTable({
        datatable(curData(),options=list(pageLength=100))

        # alternative:

        # curData()
        # and...
        # options=list(pageLength=100)) goes to renderDataTable
        #         if (!is.null(input$ma_brush)) {
        #           res <- enclosed_brush(mama, input$ma_brush)
        #           # rv_ma$highlight_vars <- res
        #         }  else {
        #           res <- NULL
        #         }
        #
        #         # TODO: total hack -- fix this correctly eventually
        # #         if (is(res, 'data.frame')) {
        # #           res <- dplyr::rename(res,
        # #                                mean = mean_obs,
        # #                                var = var_obs,
        # #                                tech_var = sigma_q_sq,
        # #                                final_sigma_sq = smooth_sigma_sq_pmax)
        # #         }
        #
        #         res
      })


      output$heatbrush <- renderPlot({
        if((is.null(input$ma_brush))|is.null(obj2)) return(NULL)

        brushedObject <- curData()

        selectedGenes <- brushedObject$ID
        toplot <- assay(obj2)[selectedGenes,]
        rownames(toplot) <- cm2$fromgtf[match(rownames(toplot),rownames(cm2))]

        if(input$pseudocounts) toplot <- log2(1+toplot)

        if(input$rowscale) toplot <- pheatmap:::scale_mat(toplot,"row")

        pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))


      })


      output$heatbrushD3 <- renderD3heatmap({
        if((is.null(input$ma_brush))|is.null(obj2)) return(NULL)

        brushedObject <- curData()

        selectedGenes <- brushedObject$ID
        toplot <- assay(obj2)[selectedGenes,]
        rownames(toplot) <- cm2$fromgtf[match(rownames(toplot),rownames(cm2))]
        mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding
        if(input$pseudocounts) toplot <- log2(1+toplot)
        if(input$rowscale) toplot <- pheatmap:::scale_mat(toplot,"row")

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


      output$geneplot <- renderPlot({

        # if(length(input$color_by_G)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
        if(is.null(input$ma_brush)) return(NULL)

        if(is.null(input$mazoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

        selectedGene <- as.character(curDataClick()$ID)
        selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
        # plotCounts(dds_cleaner,)
        genedata <- plotCounts(obj2,gene=selectedGene,intgroup = "condition",returnData = T)

        # onlyfactors <- genedata[match(input$color_by_G,colnames(genedata))]
        # onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]
        onlyfactors <- genedata[,match("condition",colnames(genedata))]

        ## intgroup can be a vector of factors. then i need interactions of the two factors
        # plotCounts(ddsmf_global,gene="ENSMUSG00000026080",intgroup = c("tissue","condition"),returnData = T) -> dh
        # dh$f1f2 <- interaction(dh$tissue,dh$condition)
        # dh  %>% ggplot(aes(x=f1f2,y=count,fill=f1f2)) + geom_boxplot()


        ## TODO: make the intgroup/colr by also somehow multiple-selectable?

        # genedata$plotby <- lapply(1:ncol(onlyfactors),function(arg) onlyfactors[,arg]) %>% interaction()
        genedata$plotby <- interaction(onlyfactors)



        ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + scale_y_log10(name="Normalized counts") + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
        # exportPlots$genesZoom <- res
        # res
      })









    })

  uiui <- (
    # var_range <- function(id, label, variable) {
    #   rng <- range(variable, na.rm = TRUE)
    #   sliderInput(id, label, rng[1], rng[2], rng)
    # }

    navbarPage("supaLiveMA",

               tabPanel("Overview",
                        fluidRow(
                          div(h3('maLive'), align = 'center')
                        ),
                        fluidRow("some text")
               ),
               navbarMenu("MA plot with...",
                          tabPanel("gene names visible on hovering mouse",
                                   ggvisOutput("p")
                          )
                          ,
                          tabPanel("brushable surface",headerPanel("MA plot interactive exploration"),
                                   fluidRow(verbatimTextOutput("deb")),
                                   fluidRow(column(4,
                                                   h4("MA plot - Interactive!"),
                                                   plotOutput('plotma', brush = 'ma_brush')),
                                            column(4,
                                                   h4("Zoomed section"),
                                                   plotOutput("mazoom",click= 'mazoom_click')),
                                            column(4,
                                                   h4("Boxplot for the selected gene"),
                                                   plotOutput("geneplot")
                                            )
                                   ),
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
                                   fluidRow(dataTableOutput("ma_brush_out"))
                          )

               )



    )
  )

  shinyApp(ui = uiui, server = serser)
}




