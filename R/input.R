


#' reads in csv of meta data for drug assays
#'
#' @param file_path A path and name of the meta file for analysis.
#'
#' @return A dataframe of user supplied assay meta data.
#' @noRd
read_csv_file <- function(file_path) {

  if (!file.exists(file_path)) {
    cli::cli_abort(c('File not found.',
                     "x" = "`{file_path}` does not exist",
                     "i" = "Try checking file paths or name"))
  }

  if (file.size(file_path) < 200) {
    cli::cli_abort(c('Can not read empty files.',
                     "i" = "`{file_path}` is empty"))
  }

  df <-  try(utils::read.csv(file = file_path,
                    header = T,
                    check.names = F)) %>%
    stats::na.omit()

  return(df)
}

#' validates the meta file for assays and summarise contents
#'
#' @param df A data frame of user supplied meta csv data.
#'
#' @return A data frame of user supplied assay meta data.
#' @noRd
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
  #generate a unique ID tag for each assay
  df <- df %>% dplyr::mutate(IC50_key = stringi::stri_rand_strings(nrow(df),
                                                             14),
                       .before = .data$plate_id)

  df <- as.data.frame(apply(df,2, trimws))
  df$starting_uM <- as.numeric(as.character(df$starting_uM))
  df$dilution_factor <- as.numeric(as.character(df$dilution_factor))
  if ((!is.numeric(df$starting_uM) | !is.numeric(df$dilution_factor))) {
    cli::cli_abort(
      c(
        "Provided dose vector must be a numeric.",
        "x" = " check meta starting_uM and dilution factor data type"
      ))}

  n <- df %>% dplyr::group_by(.data$plate_id, .data$position_id) %>%
    dplyr::filter(dplyr::n() > 1)
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
#' @noRd
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
#' @return A list of plate data: 1. raw data, 2. plate ID, 3.date 4. filename
#' @noRd
import_plate <-  function(file, meta) {
  j <- 1
  plates <- list()
  for (i in file) {

    if (!file.exists(i)) {
      cli::cli_inform(c("x" = "{basename(i)} does not exist",
                        "i" = "Skipping File."
      ))
      next
    }

    if (file.size(i) < 200) {
      cli::cli_inform(c("x" = "{basename(i)} is empty",
                       "i" = "Skipping file."
      ))
      next
    }
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

      cli::cli_inform(c("x" = "Plate ID missing from A3 {basename(i)}",
                      "i" = "Skipping file."))
      next
    }
    id <- gsub(" ","", strsplit(tag, " ")[[1]][2], fixed = TRUE)
    id <- gsub("-","_", id)
    if (!any(id == meta$plate_id)  ) {
      cli::cli_inform(c("x" = "ID: {id} from {basename(i)} not in meta file",
                        "i" = "Skipping file."))
      next
    }

    #get date of plate read from file
    date <- as.data.frame(utils::read.csv(i,
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



#' Merge multiple csv results or data files.
#' @param prefix output file prefix.
#' @param file_list A list of files to merge.
#'
#' @return A data frame of merged csv files
#' @noRd
merge_csv_files <- function(file_list, prefix = 'group') {
  # Read each CSV file into a list of data frames
  data_frames_list <- lapply(file_list, read.csv)

  # Combine all data frames into one
  merged_data_frame <- do.call(rbind, data_frames_list)

  write.csv(merged_data_frame, paste(prefix,'_merge.csv'), row.names = F)
  return(merged_data_frame)
}


#' Filter results data frame based on vector of IDs.
#'
#' @param dataframe A data frame of results
#' @param column_name A column you wish to filter by
#' @param filter_vector List of IDs to filter for
#'
#' @return A filtered data frame.
#' @noRd
filter_dataframe_by_vector <- function(dataframe, column_name = 'compound', filter_vector) {
  # Filter the dataframe
  filtered_df <- subset(dataframe, dataframe[[column_name]] %in% filter_vector)

  return(filtered_df)
}
