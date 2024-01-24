
#' Checks if numeric
#'
#' @param x object
#' @param arg variable name.
#' @param call environment.
#' @noRd
check_if_not_numeric <- function(x,
                                 arg = rlang::caller_arg(x),
                                 call = rlang::caller_env()) {
  if (!is.numeric(x)) {
    cli::cli_warn(
      c(
        "Provided vector, {.arg {arg}}, must be a numeric vector",
        "x" = " You've supplied a {.cls {class(x)} vector.}"
      ), call = call
    )
  }

}

#' Check if file exists and is not empty
#'
#' @param x variable name.
#' @param call environment.
#' @noRd
abort_not_found <- function(x,
                            call = rlang::caller_env()) {

  if (!file.exists(x)) {
    cli::cli_abort(c('File not found.',
                     "x" = "`{x}` does not exist",
                     "i" = "Try checking file paths or name"
    ), call = call)
  }

  if (file.size(x) < 3) {
    cli::cli_abort(c('Can not read empty files.',
                     "x" = "`{x}` is empty"
    ), call = call)
  }

}

#' checks arguments types
#'
#' @param arg variable name.
#' @param must required type.
#' @param not object.
#' @param call environment.
#' @noRd
abort_bad_argument <- function(arg, must, not = NULL,
                               call = rlang::caller_env()) {

  msg <- glue::glue("`{arg}` must {must}")
  if (!is.null(not)) {
    not <- typeof(not)
    msg2 <- glue::glue("`{arg}` supplied is a {not}.")
  }

  cli::cli_abort(c('wrong argument type',
                   "i" = msg,
                   "x" = msg2
  ), call = call)
}
