


#' reads in csv of meta data for drug assays
#'
#' @param file_path A path and name of the meta file for analysis.
#'
#' @return A dataframe of user supplied assay meta data.
#' @export
read_csv_file <- function(file_path) {
  abort_not_found(file_path)

  df <-  try(utils::read.csv(file = file_path,
                    header = T,
                    check.names = F)) %>%
    dplyr::mutate_if(is.character, stringr::str_trim) %>%
    stats::na.omit()

  return(df)
}

#' validates the meta file for assays and summarise contents
#'
#' @param df from the user csv file
#'
#' @return A dataframe of user supplied assay meta data.
#' @export
#'
validate_meta <- function(df) {
  nms <- c('plate_id','position_id',
           'format',	'replicates',
           'index',	'compound','cell',
           'starting_uM','dilution_factor')

  if (!length(colnames(df)) == 9) {
    cli::cli_abort(c('Error: Expected 9 meta data columns:',
                     "x" = " `{file_path}` has {length(colnames(df))}",
                     "i" = "See test data meta file for example"
    ))
  }

  colnames(df) <- nms
  df$plate_id <- gsub("-","_", df$plate_id)
  df$IC50_key <- stringi::stri_rand_strings(nrow(df), 14)
  check_if_not_numeric(df$starting_uM)
  check_if_not_numeric(df$dilution_factor)
  n <- df %>% dplyr::group_by(.data$plate_id, .data$position_id) %>%
    dplyr::filter(n() > 1)
  if (dim(n)[1] > 0) {
    cli::cli_abort(c('duplicate positions used on same plate:',
                     "i" = "check meta file"
    ))
  }
  z <- df %>% dplyr::select(.data$plate_id, .data$compound, .data$cell) %>%
    dplyr::summarise_all(dplyr::n_distinct)
  cli::cli_inform(c(
    "v" = "Meta csv read successfully",
    "i" = "It contains {length(df$plate_id)} assay{?s}.",
    "i" = "These are in {z$plate_id} plate{?s}.",
    "i" = "Testing {z$compound} compound{?s} with {z$cell} cell line{?s}."
  ))

  return(df)
}

#' Checks for the existence of raw plate reader files, prefixed with either 'TRno' or 'automated'
#'
#' @param file_path A path to the directory where raw data is stored with  meta file
#' @export
#' @return A list of raw files for analysis.
check_raw_files <- function(file_path) {
  file.names <- dir(dirname(file_path),
                    pattern = "TRno|automated",
                    full.names = T)

  if (length(file.names) == 0) {
    cli::cli_abort(c('Raw data files not found.',
                     "x" = "No files starting with (TRno*) or (automated*) found.",
                     "i" = "check file names conform or path"))
  }

  cli::cli_inform(c(
    "v" = "{length(file.names)} raw data file{?s} found"))

  return(file.names)
}


#' checks raw data files
#'
#' @param file List of raw data files matching prefix
#' @param meta The meta dataframe.
#'
#' @return A list of dataframes of raw data for each file/plate
#' @export
import_plate <-  function(file, meta) {
  j <- 1
  plates <- list()
  for (i in file) {

    abort_not_found(i)
    data <- as.data.frame(utils::read.csv(i,
                                   sep = ",",
                                   header = FALSE,
                                   skip = 10,
                                   stringsAsFactors = FALSE,
                                   blank.lines.skip = FALSE))

    tag <- as.data.frame(utils::read.csv(i,
                                  sep = ",",
                                  header = FALSE,
                                  nrows = 3,
                                  stringsAsFactors = FALSE))[3,1]

    if (nchar(tag) < 5) {
      cli::cli_alert(c('Plate ID not found.',
                       "x" = "Raw data missing ID in A3",
                       "i" = "check file format"))
      next
    }
    id <- gsub(" ","", strsplit(tag, " ")[[1]][2], fixed = TRUE)
    id <- gsub("-","_", id)
    if (!any(id == meta$plate_id)  ) {
      cli::cli_alert(c('Plate ID does not match meta data.',
                       "x" = "ID: {id} not known",
                       "i" = "check file format"))
      next
    }

    #get date of plate read from file
    date <- as.data.frame(read.csv(i,
                                   sep = ",",
                                   header = FALSE,
                                   nrows = 3,
                                   stringsAsFactors = FALSE))[2,1]

    date <- strsplit(date, " ")[[1]][2]
    plates[[j]] <- list(data,id, date, basename(i))
    j <- j + 1
  }
  return(plates)
}







