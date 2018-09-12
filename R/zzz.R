#' @importFrom shiny addResourcePath
.onLoad <- function(...) {
  # Create link to logo
  shiny::addResourcePath("ideal", system.file("www", package="ideal"))
}

