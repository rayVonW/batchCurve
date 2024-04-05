

#' drc fitting function
#'
#' @param data Normalised measurements from a drug assay.
#'
#' @return An LL4 fitted model or error to handle
#' @noRd
fit.LL4 <- function(data) {
  # control arguments
  control <- drc::drmc(errorm = T,noMessage = TRUE, warnVal = -1)
  #fitting function
  possibleError <- tryCatch(

    drc::drm(data = data,
                  formula =  value~dose,
                  fct = drc::LL.4(fixed = c(NA,NA, NA, NA),
                                  names = c("hill",
                                            "min",
                                            "max",
                                            "IC50")),
                  na.action = stats::na.omit,
                  control = control),

    error = function(e) e,
    silent = TRUE
  )
  return(possibleError)

}

#' adds model coefficients to data
#'
#' @param model An LL4 model fitted to the data.
#' @param assays A data frame of assay meta data.
#' @param locate A data frame of coordinates.
#'
#' @return Assay meta data frame with LL4 coefficients appended
#' @noRd
retrieve_results <- function(model, assays, locate){

  coefs <- model$coefficients
  names(coefs) <- gsub(":.*","",names(coefs))
  if ( length(coefs) < 3) {coefs['max'] <- 100}

  GI50 <- with(
    as.list(coefs),
    exp(
      log(IC50) + (1 / hill) * log(max / (max - 2 * min))
    )
  )
  #Collect IC50 and IC90 with estimated standard errors
  SE <- data.frame(drc::ED(model, c(50, 90), display = FALSE))

  assays$IC50 <- as.numeric(SE[1,1])

  assays$IC50_SE <- as.numeric(SE[1,2])
  assays$GI50 <- GI50
  assays$IC90 <- as.numeric(SE[2,1])
  assays$IC90_SE <- as.numeric(SE[2,2])

  assays$min <- as.numeric(coefs[2])
  assays$max <- as.numeric(coefs[3])
  assays$hill_slope <- as.numeric(coefs[1])

  assays$Zprime <- round(1 - ((3 * (
    locate$no_drug_sd + locate$blank_sd
  )) / abs(locate$no_drug - locate$blank)), 2)

  return(assays)
}

#' Generate curves from LL4 coefficients for drug assay plots.
#'
#' @param data A data frame of model coefficients and drug assay data.
#'
#' @return A data frame of predicted values using the model and dose range.
#' @noRd
generate_curve <- function(data){
  #filter out NA where n = 2 or failed fit
  data <- data[stats::complete.cases(data),]
  data$IC50 <- as.numeric(as.character(data$IC50))

  curves <- data.frame()
  keys <- unique(data$IC50_key)

  for (k in  keys) {
    data2 <-  subset(data, IC50_key == k)


    # create log drug concentration sequence
    x <- exp(seq(log(min(data2$dose)),
                 log(max(data2$dose)),
                 length.out = 50))

    data2 <-  unique(subset(data2, select = c(
                    'IC50',
                    'min',
                    'max',
                    'hill_slope',
                    'IC50_key',
                    'cell',
                    'date',
                    'index',
                    'compound')))

    # LL.4 function with models coefficients
    f <- function(x) data2[1,2] + (data2[1,3] - data2[1,2]) / (1 + exp(data2[1,4]*(log(x) - log(data2[1,1]))))
    #MIN + ((MAX-MIN) / (1 + exp(HILLSLOPE*(log(x)-log(IC50)))))

    # predict the values for log drug
     y0 <- f(x)

    curve <- data.frame(dose = x, predicted = y0)

    curve$IC50_key <- k
    curve$cell <- unique(data2$cell)
    curve$compound <- unique(data2$compound)
    curve$index <- unique(data2$index)
    curve$date <- unique(data2$date)

    curves <- rbind(curves, curve)

  }
  return(curves)
}
