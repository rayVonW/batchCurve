

#' finds drug assay data on a 96 or 384 well plate
#'
#' @param assay A dataframe of a single meta entry
#' @param plate A dataframe of plate reader measurements
#'
#' @return A dataframe of plate position coordinates and control values
#' @noRd
get_locations <- function(assay, plate){
  #The four formats are stored in sysdata.rda and are lazy loaded:
  #triplicate_384, triplicate_96, duplicate_384, duplicate_96
  #see github for layouts

  locate <- NULL

  if (assay$format == '384w') {
    if (!assay$position_id %in% c('pos_1','pos_2','pos_3','pos_4',
                                       'pos_5','pos_6','pos_7','pos_8',
                                       'pos_9','pos_10')) {
      cli::cli_abort(c(
                     'x' = 'Position ids format fail. ',
                     'i' = 'Expected pattern: pos_[1:10].'))
    }
    if (assay$replicates == 3) {
      locate <- triplicate_384 %>%
        dplyr::filter(.data$position_name == assay$position_id)
    }
    if (assay$replicates == 2) {
      locate <- duplicate_384 %>%
        dplyr::filter(.data$position_name == assay$position_id)
    }
    locate$blank <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 13])))
    locate$blank_sd <-
      stats::sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 13])))
    locate$no_drug <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$no_drug_sd <-
      stats::sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
  }

  if (assay$format == '96w') {
    if (assay$replicates == 3) {
      if (!assay$position_id %in% c('pos_1','pos_2')) {
        cli::cli_abort(c('x' = 'Position ids format fail. ',
                       'i' = 'Expected pattern: pos_[1:2].'))
      }
      locate <- triplicate_96 %>%
        dplyr::filter(.data$position_name == assay$position_id)
    }
    if (assay$replicates == 2) {
      if (!assay$position_id %in% c('pos_1','pos_2', 'pos_3')) {
        cli::cli_abort(c('x' = 'Position ids format fail. ',
                       'i' = 'Expected pattern: pos_[1:3].'))
      }
      locate <- duplicate_96 %>%
        dplyr::filter(.data$position_name == assay$position_id)
    }
    locate$blank <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$blank_sd <-
      stats::sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 12])))
    locate$no_drug <-
      mean(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 11])))
    locate$no_drug_sd <-
      stats::sd(as.numeric(
        as.character(plate[locate$set_rows_t:locate$set_rows_b, 11])))
  }

  return(locate)
}


#' subsets plate data per assay, converts values to percentages then adds the dose curve
#'
#' @param data A dataframe of measurement values per
#' @param location A dataframe of position coordinates
#' @param assays_info A dataframe of meta data
#'
#' @return A dataframe of normalised values for a dose-response assay
#' @noRd
#'
get_assay_data <- function(data, location, assays_info){
  #retrieve assay specific data
  data_set <- as.data.frame(t(data[location$set_rows_t:location$set_rows_b,
                                   location$set_cols_l:location$set_cols_r]))

  #add no drug and blank data
  data_set[nrow(data_set) + 1,] = c( rep(location$no_drug,assays_info$replicates))
  data_set[] <- lapply(data_set, function(x) as.numeric(as.character(x)))
  data_set_blank <- data_set %>% dplyr::mutate_all(.funs = list(~. - location$blank))

  #Create percentage dataframe

  data_pct <- data.frame(sapply(names(data_set_blank),
                                function(x) { data_set_blank[paste0(x, "_pct")] <<- (
                                  (data_set_blank[x] / data_set_blank[nrow(data_set_blank),x])*100) }))


  #create dose series and add to data
  dilutionF <- as.numeric(as.character(assays_info$dilution_factor))
  startDose <- as.numeric(as.character(assays_info$starting_uM))
  dose <- c()

  for (q in 1:nrow(data_pct )) {
    dose[q] <- startDose
    startDose <- startDose/dilutionF

  }

  data_pct  <- cbind(dose = dose, data_pct )

  #rename columns
  if (ncol(data_pct ) == 3) {
    rep_names <- c("dose", "replicate_1", "replicate_2")
  }
  else {
    rep_names <- c("dose", "replicate_1", "replicate_2", "replicate_3")
  }
  colnames(data_pct ) <- rep_names

  #convert to long format
data_pct  <- data.frame(stats::reshape(data = data_pct,
                                         direction = "long",
                                         v.names = "value",
                                         varying = rep_names[-1],
                                         idvar = "dose",
                                         timevar = "replicate"))

  rownames(data_pct) <- NULL

  data_pct$value <- ifelse(data_pct$value > 100, 100, data_pct$value)
  data_pct$value <- ifelse(data_pct$value < 0, 0, data_pct$value)

  return(data_pct)

}

