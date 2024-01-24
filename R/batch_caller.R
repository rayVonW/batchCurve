

#' iterates over assays for IC50 estimation
#'
#' @param file_path A file containing the batches meta data
#' @param batch_id A unique character value to identify the batch
#'
#' @return A dataframe of LL4 model coefficiants for dose-response assays
#' @export

run <- function(file_path, batch_id){

   r <- read_csv_file(file_path) %>% validate_meta()
   f <- check_raw_files(file_path = file_path)
   d <- import_plate(f, r)

  # To dump data
   all_Plots <- list()
   results_all <- data.frame(matrix(ncol = 21, nrow = 0))
   curves_all <- data.frame(matrix(ncol = 8, nrow = 0))
   data_all <- data.frame(matrix(ncol = 8, nrow = 0))
   c <- 1
  #
   for (i in 1:length(d)) {
  #   #process per plate/file
     plate_assays <-  r %>% dplyr::filter(.data$plate_id  == d[[i]][[2]])
     plate_assays$file_name <- d[[i]][[4]]
     plate_assays$date <- format(as.Date(d[[1]][[3]],
                                         format = "%d/%m/%Y"))
     plate <- as.data.frame(d[[1]][[1]])

     for (j in 1:nrow(plate_assays)) {
        #get location coordinates per assay
        locate <- get_locations(plate_assays[j,], plate)
        #retrieve data from plate per assay
        data <- get_assay_data(plate, locate, plate_assays[j,])


       # fit model
       m <- fit.LL4(data, plate_assays[j,])

       if (inherits(m, "error")) {
         cat('    | Convergence failed. The model was not fit!\n')

         fail <- c(plate_assays[j,],
                   'Failed',
                   NA,NA,NA,NA,NA,NA,NA,NA)
         fail_p <- ''#failed_plots(data , assay)
         c <- c + 1
         next
       }

       #model_results <- retrieve_results(m, plate_assays[j,] )
       return(data)
  }
  }
}
