#' Dose response meta data example
#'
#' An example meta data file describing drug assays for analysis ...
#' @docType data
#' @usage data(example_meta)
#' @keywords dataset
#' @format ## `example_meta`
#' A data frame with 23 rows and 9 columns:
#' \describe{
#'   \item{plate_id}{ID given to each plates reader file}
#'   \item{position_id}{Position of an assay on plate, 'pos_x'}
#'   \item{format}{Plate format, either '384w' or '96w'}
#'   \item{replicates}{How many replicates per assay used, either '2' or '3'}
#'   \item{index}{label for compound, used for sorting output}
#'   \item{compound}{Supplier ID of compound}
#'   \item{cell}{Cell ID}
#'   \item{starting_uM}{Drug concentration in column 1 of plate, i.e. the highest concentration in uM}
#'   \item{dilution_factor}{The serial dilution factor used in the assay}
#'   ...
#' }
#' @source author
"example_meta"
