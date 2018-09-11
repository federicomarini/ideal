
read_gmt <- function(gmtfile){
  # TODO: some checks on the gmt file format?
  
  input_lines <- strsplit(readLines(gmtfile), "\t")
  # the two first elements in each signature would be the id and e.g. its source
  pathways <- lapply(input_lines, tail, -2)
  # pick the name from the provided info in the gmt file
  names(pathways) <- sapply(input_lines, head, 1)
  return(pathways)
}


sig_heatmap <- function(vst_data, signature, de_only = TRUE) {
  
}
