#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom magrittr %>%
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".data"))
