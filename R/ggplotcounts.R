



ggplotCounts <- function(dds,gene,intgroup="condition",annotation_obj=annotation_obj,...){
  df <- plotCounts(dds,gene,intgroup,returnData = TRUE,...)
  df$sampleID <- rownames(df)
  genesymbol <- annotation_obj$gene_name[match(gene,annotation_obj$gene_id)]
  jittered_df <- df
  # jittered_df$conditionj <- jitter(as.numeric(factor(jittered_df$condition)))
  jittered_df$countj <- jitter(jittered_df$count)

  onlyfactors <- df[,match(intgroup,colnames(df))]
  df$plotby <- interaction(onlyfactors)


  p <-
    ggplot(df, aes(x=plotby,y=count,col=plotby)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_text(data = jittered_df,aes(x=conditionj,y=countj,label=sampleID)) +
    geom_text(aes(label=sampleID),hjust=-.1,vjust=0) +
    labs(title=paste0("Normalized counts for ",genesymbol," - ",gene)) +
    scale_x_discrete(name="") +
    geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) +
    scale_color_discrete(name="Experimental\nconditions") +
    scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.1,NA)) +
    theme_bw()

  p

}
