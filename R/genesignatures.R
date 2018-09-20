#' Read in a GMT file
#' 
#' Returns a list of pathways from a GMT file.
#' 
#' @param gmtfile A character value, containing the location of the GMT formatted
#' file. It can also be a file found online
#'
#' @return A list of vectors, one for each pathway in the GMT file.
#' @export
#'
#' @examples
#' # this example reads in the freely available pathways from wikipathways
#' mysigs <- read_gmt(
#'   "http://data.wikipathways.org/current/gmt/wikipathways-20180910-gmt-Homo_sapiens.gmt")
#' head(mysigs)
#' # see how the gene identifiers are encoded as ENTREZ id
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
                        res_data = NULL, FDR = 0.05,
                        de_only = FALSE, annovec, title = "",
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        center_mean = TRUE, scale_row = FALSE
                        # ,
                        #anno_colData
                        ) {
  
  mydata <- assay(vst_data)
  
  signature_original_ids <- names(annovec)[match(my_signature,annovec)]
  
  sig_to_keep <- (signature_original_ids %in% rownames(mydata))#
  my_signature <- my_signature[sig_to_keep]
  signature_original_ids <- signature_original_ids[sig_to_keep]
  
  mydata_sig <- mydata[signature_original_ids,]
  
  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove,]
  
  if(center_mean)
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)
  
  if(de_only) {
    de_res <- deseqresult2DEgenes(res_data,FDR)
    de_genes <- de_res$id
    de_to_keep <- rownames(mydata_sig) %in% de_genes
    mydata_sig <- mydata_sig[de_to_keep,]
  }
  
  pheatmap(mydata_sig,
           # annotation_col = anno_colData,
           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
           scale = ifelse(scale_row,"row","none"),main = title,
           labels_row = annovec[rownames(mydata_sig)])
}
