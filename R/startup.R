# startup
source('/Volumes/users$/marinif/Development/ideal/R/ideal.R')
library(shinydashboard)
library(shiny)
library(d3heatmap)
ls()
library(pcaExplorer)
example(pcaExplorer)
library(DESeq2)
library(ggplot2)
library(shinyAce)
library(DT)
library(knitr)
library(rmarkdown)
library(pheatmap)
dds_airway <- DESeq2::DESeq(dds_airway)
res_airway <- results(dds_airway,contrast=c("dex","trt","untrt"))
anno_df <- get_annotation_orgdb(dds_airway,"org.Hs.eg.db","ENSEMBL")
load("/Volumes/marinif/032-ruf-macrophages/cm2.RData")
res_airway$symbol <- anno_df$gene_name[match(rownames(res_airway),anno_df$gene_id)]
minires <- res_airway[1:200,]
res_airway
# library(FMmisc)
library(IHW)

ccmm <- counts(dds_airway)
eedd <- read.csv("design_commas.txt")
source('/Volumes/users$/marinif/Development/ideal/R/ideal_snapshot_oct14.R')

ideal(res_obj = minires,dds_obj = dds_airway,annotation_obj = anno_df)
