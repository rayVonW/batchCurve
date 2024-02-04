#' Get path to readr example
#'
#' batchCurve comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' readr_example()
#' readr_example("TRno4850.csv")
batchCurve_example <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "batchCurve"))
  } else {
    system.file("extdata", file, package = "batchCurve", mustWork = TRUE)
  }
}
