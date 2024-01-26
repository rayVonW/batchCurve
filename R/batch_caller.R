

#' iterates over assays for IC50 estimation
#'
#' @param file_path A file containing the batches meta data
#' @param batch_id A unique character value to identify the batch
#'
#' @return A data frame of LL4 model coefficients for dose-response assays
#' @export

fit_data <- function(file_path, batch_id) {
  cli::cli_rule(left = "analysing batch: {batch_id}")
   # read and check [input.R]
   r <- read_csv_file(file_path) %>% validate_meta()
   cli::cli_rule(left = "Assessing raw data files")
   f <- check_raw_files(file_path = file_path)
   d <- import_plate(f, r)

  # To collect data
   results_all <- data.frame(matrix(ncol = 21, nrow = 0))
   cn <- NULL
   data_all <- data.frame(matrix(ncol = 5, nrow = 0))
   c <- 1

   for (i in 1:length(d)) {# loop over each plate:
     cli::cli_h2("Analysing  file: {d[[i]][[4]]}")
     #process per plate/file
     plate_assays <-  r %>% dplyr::filter(.data$plate_id  == d[[i]][[2]])
     plate_assays$file_name <- d[[i]][[4]]
     plate_assays$date <- format(as.Date(d[[1]][[3]],
                                         format = "%d/%m/%Y"))
     plate <- as.data.frame(d[[i]][[1]])

     for (j in 1:nrow(plate_assays)) {# functions per assay on plate:

        #get location coordinates per assay [plate_layout.R]
        locate <- get_locations(plate_assays[j,], plate)
        #retrieve normalised data per assay [plate_layout.R]
        data <- get_assay_data(plate, locate, plate_assays[j,])

        #format and append normalised data

        data_f <- data %>% dplyr::mutate(value = round(.data$value,2),
                                       dose = round(.data$dose, 5)) %>%
          tidyr::pivot_wider(names_from = replicate,
                             values_from = .data$value) %>%
          dplyr::mutate(key = plate_assays[j,'IC50_key'],
                        .before = .data$dose)

        if (!"replicate_3" %in% colnames(data_f)) { data_f$replicate_3 <- NA}

        data_all <- rbind(data_all, data_f)

       # fit model or collect error [fitting.R]
       m <- fit.LL4(data)

       if (inherits(m, "error")) {
         cli::cli_inform(c(
           "x" = "At {plate_assays[j,'position_id']}:
         {plate_assays[j,'index']} treated {plate_assays[j,'cell']},
         fit failed."))

         results_all[nrow(results_all) + 1,] <- c(plate_assays[j,],
                   'Failed',
                   NA,NA,NA,NA,NA,NA,NA,NA)

         c <- c + 1
         next
       }
       #  extract coefficients from model and add to results table
       cli::cli_inform(c(
         "v" = "At {plate_assays[j,'position_id']}:
         {plate_assays[j,'index']} treated {plate_assays[j,'cell']},
         fit successful."))

       # extract coefficient from model [fitting.R]
       rs <- retrieve_results(m,
                        plate_assays[j,],
                        locate)
       if (is.null(cn)) { colnames(results_all) <- colnames(rs)}
       results_all[nrow(results_all) + 1,] <- rs

       c <- c + 1
      }
   }

   #export all normalised data to csv file
   data_all %>% utils::write.csv(file = paste0(dirname(file_path),
                                               '/',batch_id,
                                               "_data.csv"),
                                 row.names = F)
   results_all %>% utils::write.csv(file = paste0(dirname(file_path),
                                                  '/',batch_id,
                                                  "_results.csv"),
                                 row.names = F)


   cli::cli_h2('Analysis complete.')
   cli::cli_inform(c(
     "i" = "Result and data csv files saved to {dirname(file_path)}"))
   return(list(results_all,data_all))
}



#' Produces group plots of dose response models from fit_data()
#'
#' @param results A data frame, from fit_data.
#' @param data A data frame, from fit_data.
#' @param batch_id An character ID for output file prefix.
#' @param plot.var A column name you wish to subset data per plot, default is "compound", you may wish specify "cell".
#' @param colour.var A column name you wish to differentiate data per plot by colour, default = "cell", you may wish specify "compound".
#' @param facet.var  A column name you wish to differentiate data by facet, default = NULL, you may wish to have "cell" on separate plots.
#' @param grid.var a numeric value specifying the grid size in the exported pdf, default = 3 i.e 3x3.
#' @return A list of plot objects per plot_var
#' @export

plot_fit <- function(results, data, batch_id,
                     plot.var = 'compound',
                     colour.var = 'cell',
                     facet.var = NULL,
                     grid.var = 3) {
  cli::cli_h2(left = "Generating plots")
  n <- data %>% tidyr::pivot_longer(cols = -c(.data$dose, .data$key),
                                 names_to = "replicate",
                                 values_to = "value") %>%
    dplyr::left_join(results, .data,by = c('IC50_key' = 'key')) %>%
    dplyr::arrange(.data$index)

  c <- generate_curve(n)

  g <- unique(n[[{{plot.var}}]])
  all_grouped_plots <- list()
  j <- 1
  for (i in g) {
    sub.c <- subset(c, c[[plot.var]] == i)
    sub.n <- subset(n, n[[plot.var]] == i)
    index <- unique(sub.n$index)
    suppressWarnings(MaxD <- max(sub.c$dose))
    suppressWarnings(MinD <- min(sub.c$dose))

    conc <- 'uM'
    if (MaxD <= 1) {
      sub.c$dose <- sub.c$dose*1000
      sub.n$dose <- sub.n$dose*1000
      conc <- 'nM'

      suppressWarnings(MinD <- min(sub.c$dose))
      suppressWarnings(MaxD <- max(sub.c$dose))
    }
    colour.var <- rlang::sym(colour.var)
    suppressWarnings(
      all_grouped_plots[[j]] <-
        ggplot2::ggplot(data = sub.n,
                        ggplot2::aes(x = .data$dose,
                                         y = .data$value,
                                         group = .data$IC50_key)) +
        ggplot2::geom_point(size = 0.4,
                        ggplot2::aes(colour = !!colour.var)) +
        ggplot2::geom_line(data = sub.c,
                        ggplot2::aes(x = .data$dose,
                                     y = .data$predicted,
                                     colour = !!colour.var,
                                     group = .data$IC50_key)) +
        ggplot2::scale_x_log10(labels = scales::comma,
                              limits = c(MinD,MaxD)) +
        ggplot2::annotation_logticks(outside = TRUE,
                                   side = 'b',
                                   color = 'black',
                       long = ggplot2::unit(0.35,"lines"),
                       mid = ggplot2::unit(0.2,"lines"),
                       short = ggplot2::unit(0.2,"lines")) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::ylim(-10,120) +
        ggplot2::labs(title = i, subtitle = index,
                    x = paste0('Drug[',conc,']'),
                    y = "Growth (%)") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = -1)))
     j <- j + 1
  }
  #generate and save pdf of all grouped plots
  group_name <- paste0(batch_id,'_grouped_plots.pdf')

  suppressWarnings(
    plots <- gridExtra::marrangeGrob(all_grouped_plots,
                        nrow = grid.var,
                        ncol = grid.var))
  suppressWarnings(
  ggplot2::ggsave(group_name, plots,
                  width = 11, height = 11,
                  units = "in",
                  dpi = 300))
  cli::cli_inform("v" = "Export complete")
  return(all_grouped_plots)
}

