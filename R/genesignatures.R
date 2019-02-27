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
#'   "http://data.wikipathways.org/20180910/gmt/wikipathways-20180910-gmt-Homo_sapiens.gmt")
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


#' Plot a heatmap of the gene signature on the data
#' 
#' Plot a heatmap for the selected gene signature on the provided data, with the 
#' possibility to compactly display also DE only genes
#'
#' @param vst_data A \code{\link{DESeqTransform}} object - usually the variance
#' stabilized transformed data, which will be used to extract the expression values 
#' @param my_signature A character vector, usually named, containing the genes 
#' which compose the gene signature
#' @param res_data A \code{\link{DESeqResults}} object. If not provided, it can
#' be computed during the execution of the application
#' @param FDR Numeric value between 0 and 1, the False Discovery Rate
#' @param de_only Logical, whether to display only DE genes belonging to the pathway - 
#' defaults to FALSE
#' @param annovec A named character vector, with the corresponding annotation across IDs
#' @param title Character, title for the heatmap
#' @param cluster_rows Logical, whether to cluster rows - defaults to TRUE
#' @param cluster_cols Logical, whether to cluster column - defaults to FALSE. 
#' Recommended to be set to TRUE if de_only is also set to TRUE
#' @param center_mean Logical, whether to perform mean centering on the expression 
#' values. Defaults to TRUE, as it improves the general readability of the heatmap
#' @param scale_row Logical, whether to perform row-based standardization of the 
#' expression values
#'
#' @return A plot based on the \code{pheatmap} function
#' @export
#'
#' @examples
#' # with the well known airway package...
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' \dontrun{
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#' vst_airway <- DESeq2::vst(dds_airway)
#' library(org.Hs.eg.db)
#' annovec <- mapIds(org.Hs.eg.db, rownames(dds_airway),"ENTREZID","ENSEMBL")
#' mysignatures <- read_gmt(
#'   "http://data.wikipathways.org/20190210/gmt/wikipathways-20190210-gmt-Homo_sapiens.gmt")
#' mysignature_name <- "Lung fibrosis%WikiPathways_20190210%WP3624%Homo sapiens"
#' library(pheatmap)
#' sig_heatmap(vst_airway,
#'             mysignatures[[mysignature_name]],
#'             res_data = res_airway,
#'             de_only = TRUE,
#'             annovec = annovec,
#'             title = mysignature_name,
#'             cluster_cols = TRUE
#'             )
#' }
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
