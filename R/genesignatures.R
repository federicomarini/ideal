
read_gmt <- function(gmtfile){
  # TODO: some checks on the gmt file format?
  
  input_lines <- strsplit(readLines(gmtfile), "\t")
  # the two first elements in each signature would be the id and e.g. its source
  pathways <- lapply(input_lines, tail, -2)
  # pick the name from the provided info in the gmt file
  names(pathways) <- sapply(input_lines, head, 1)
  return(pathways)
}


sig_heatmap <- function(vst_data, my_signature,
                        de_only = FALSE, annovec, title = "",
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        center_mean = TRUE, scale_row = FALSE) {
  
  mydata <- assay(vst_data)
  
  
  signature_original_ids <- names(annovec)[match(my_signature,annovec)]
  
  mydata_sig <- mydata[signature_original_ids,]
  
  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove,]
  
  if(center_mean)
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)
  
  pheatmap(mydata_sig,
           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
           scale = ifelse(scale_row,"row","none"),main = title,
           labels_row = my_signature[!to_remove])
}
