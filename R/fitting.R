

#' drc fitting function
#'
#' @param data
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
fit.LL4 <- function(data, assay) {
  #fitting function
  control <- drc::drmc(errorm = TRUE,noMessage = TRUE, warnVal = -1)

  possibleError <- tryCatch(
    LL.4 <- drc::drm(data = data,value~dose,fct = drc::LL.4(fixed = c(NA,NA, NA, NA)),
                na.action = na.omit, control = control),
    return(LL.4),
    error = function(e) e,
    silent = TRUE
  )

  if (inherits(possibleError, "error")) {
    cat('    | Convergence failed. The model was not fit!\n')

    fail <- c(assay,'Failed',NA,NA,NA,NA,NA,NA,NA,NA)
    fail_p <- failed_plots(data , assay)
    return(list(fail,fail_p))
    c <- c + 1
    next
  }
}

#' adds model coefficients to data
#'
#' @param model
#' @param assays
#'
#' @return
#' @export
#'
#' @examples
retrieve_results <- function(model, assays){
  #get coefs and calculate GI50
  if (length(model$coefficients) > 3){
    coefs <- setNames(
      model$coefficients, c("hill", "min_value", "max_value", "IC50"))
  }else{
    coefs <- setNames(
      model$coefficients,
      c("hill", "min_value", "IC50")
    )
    coefs$max_value = 100}
  print(model$coefficients)


  GI50 <- with(
    as.list(coefs),
    exp(
      log(IC50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )
  #Collect IC50 and IC90 with estimated standard errors
  SE<- data.frame(ED(model, c(50, 90), display = FALSE))

  assays$IC50<- as.numeric(SE[1,1])

  assays$IC50_SE<- as.numeric(SE[1,2])
  assays$GI50<- GI50
  assays$IC90<- as.numeric(SE[2,1])
  assays$IC90_SE<- as.numeric(SE[2,2])

  assays$min <- as.numeric(coefs[2])
  assays$max <- as.numeric(coefs[3])
  assays$hill_slope <- as.numeric(coefs[1])

  return(assays)
}
