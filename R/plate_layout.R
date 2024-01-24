

#' finds drug assay data on a 96 or 384 well plate
#'
#' @param assay
#' @param plate
#'
#' @return
#' @export
#'
#' @examples
get_locations <- function(assay, plate){
  #The four formats are stored in sysdata.rda and are lazy loaded:
  #triplicate_384, triplicate_96, duplicate_384, duplicate_96
  #see github for layouts

  locate <- NULL

  if (assay$format == '384w') {

    if (assay$replicates == 3) {
      locate <- filter(triplicate_384, position_name == assay$position_id)
    }
    if (assay$replicates == 2) {
      locate <- filter(duplicate_384, position_name == assay$position_id)
    }
    locate$blank <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 13])))
    locate$blank_sd <-
      sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 13])))
    locate$no_drug <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$no_drug_sd <-
      sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
  }

  if (assay$format == '96w') {
    if (assay$replicates == 3) {
      locate <- filter(triplicate_96, position_name == assay$position_id)
    }
    if (assay$replicates == 2) {
      locate <- filter(duplicate_96, position_name == assay$position_id)
    }
    locate$blank <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$blank_sd <-
      sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$no_drug <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 11])))
    locate$no_drug_sd <-
      sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 11])))
  }

  return(locate)
}


#' subsets plate data per assay, converts values to percentages then adds the dose curve
#'
#' @param data
#' @param location
#' @param assays_info
#'
#' @return
#' @export
#'
#' @examples
get_assay_data <- function(data, location, assays_info){
  #retrieve assay specific data
  data_set <- as.data.frame(t(data[location$set_rows_t:location$set_rows_b,
                                   location$set_cols_l:location$set_cols_r]))

  #add no drug and blank data
  data_set[nrow(data_set) + 1,] = c( rep(location$no_drug,assays_info$replicates))
  data_set[] <- lapply(data_set, function(x) as.numeric(as.character(x)))
  data_set_blank <- data_set %>% mutate_all(.funs = list(~. - location$blank))

  #Create percentage dataframe

  data_pct <- data.frame(sapply(names(data_set_blank),
                                function(x) { data_set_blank[paste0(x, "_pct")] <<- ((data_set_blank[x] / data_set_blank[nrow(data_set_blank),x])*100) }))


  #create dose series and add to data
  dilutionF <- as.numeric(as.character(assays_info$dilution_factor))
  startDose <- as.numeric(as.character(assays_info$starting_uM))
  dose <- c()

  for (q in 1:nrow(data_pct )){
    dose[q] <- startDose
    startDose <- startDose/dilutionF

  }

  data_pct  <- cbind(dose = dose, data_pct )

  #rename columns
  if (ncol(data_pct ) == 3){
    colnames(data_pct ) <- c("dose", "R1", "R2")
  }
  else {
    colnames(data_pct ) <- c("dose", "R1", "R2", "R3")
  }

  #convert to long format
  data_pct  <- reshape2::melt(data_pct , id.vars = "dose")

  data_pct$value <- ifelse(data_pct$value > 100, 100, data_pct$value)
  data_pct$value <- ifelse(data_pct$value < 0, 0, data_pct$value)

  return(data_pct)

}
