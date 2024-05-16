
#' Wilcox rank sum of unpaired IC50 data
#'
#' @param data A data frame - exported object from fit_data().
#' @param prefix A string for output file names.
#'
#' @return A data frame object of mean IC50 data and p values
#' @export
compare <- function(data, prefix = 'group') {
  l <- unique(data$compound)

  if (colnames(data[13]) == 'IC50')  {
    names(data)[names(data) == 'IC50'] <- 'IC50_uM'
  }

  ev <- data.frame()

  for (i in unique(data$compound)) {
    dsub <- subset(data, compound == i)
    #pull each mutant
    c <- unique(subset( dsub$cell, grepl("_", dsub$cell)))
    all <- data.frame()

    for (j in c) {
      #ensure background is matched
      wt <- sub("_.*", "", j)
      if (dim(subset(dsub, cell == wt))[[1]] < 1) {
        cli::cli_inform(c(
          "v" = "For compound {i}, cell {j} no wild-type assays for {wt} found. Skipping cell."))
        next}
      dsubc <- dsub[which(dsub$cell == j | dsub$cell == wt),]
      dsubc <- dsubc[, c("compound", "cell", "IC50_uM", "hill_slope", 'starting_uM')]
      dsubc$group <- factor(dsubc$cell)

      #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
      r <- stats::wilcox.test(dsubc$IC50_uM~dsubc$group)

      s <- dsubc %>%
        dplyr::group_by(.data$cell) %>%
        dplyr::reframe(compound = compound,
                       cell = cell,
                       starting_uM = max(starting_uM),
                       mean_IC50_uM = round(exp(mean(log(IC50_uM))),4),
                       stdev = round(stats::sd(IC50_uM),4),
                       n = dplyr::n()) %>% unique

      dsubm <- subset(s, cell == j)
      dsubm$background <- wt
      dsubwt.m <- subset(s, cell == wt) %>% unique
      colnames(dsubwt.m) <- c('compound',
                              'cell',
                              'max_starting_uM',
                              'wild_type_IC50',
                              'wt_stdev',
                              'wt_n')

      dsubm <- dplyr::bind_cols(dsubm, dsubwt.m[3:6] ) %>%
        dplyr::mutate(FC = round(mean_IC50_uM/wild_type_IC50,2))
      dsubm$w <- r$statistic
      dsubm$p.value <- r$p.value
      all <- rbind(all , dsubm)
    }

    ev <- rbind(ev,all) %>% dplyr::arrange(compound)

  }
  utils::write.csv(ev, paste0(prefix, '_compare_IC50.csv'), row.names = FALSE)
  return(ev)
}


#' Plots barchart of mean I50 values and displays p value
#'
#' @param data A data frame direct output from compare()
#' @param prefix A string values to append output files
#' @param grid.var An integer value describing the plot grid size in output file
#'
#' @return A list of plot objects
#' @export
plot_mean <- function(data, prefix = 'group', grid.var = 2) {

  #pull wild types out for plotting
  control <- data %>% dplyr::select(compound,
                                  background,
                                  wild_type_IC50,
                                  max_starting_uM,
                                  wt_stdev,
                                  wt_n) %>% unique %>%
    dplyr::mutate(FC = round(wild_type_IC50/wild_type_IC50,2)) %>%
    dplyr::rename(cell = background,
                  mean_IC50_uM = wild_type_IC50,
                  stdev = wt_stdev, n = wt_n) %>%
    dplyr::bind_rows(data)

  control <- control %>% dplyr::mutate(background = ifelse(is.na(background),
                                                    cell, background),
                                label = gsub("3D7_|Dd2_","",control$cell))

  control$label <- gsub("_","\n",control$label)
  control$label <- as.vector(sapply(control$label ,
                                    function(x) ifelse(grepl("_", x),
                                                       x,
                                                       paste0(x, "\n"))))

  #annotation labels

  control <- control %>% dplyr::group_by(label) %>%
    dplyr::mutate(FC_label = paste0(label,'\n',mean_IC50_uM,
                                   '\u00B5','M',
                                   '\nFC:',FC,
                                   '\nn=',n),
           symb = ifelse(p.value < 0.001, '***',
                         ifelse(p.value < 0.01, '**',
                                ifelse(p.value < 0.05, '*', " "))))
  #loop over compounds
  c <- unique(data$compound)
  j <- 1

  l <- list()
  for (i in c) {
    csub <- control %>% dplyr::filter(compound == i)
    b <- unique(csub$background)

    #loop over backgrounds
    for (k in b) {

      csub.b <- csub %>% dplyr::filter(background == k)
      if (dim(csub.b)[[1]] < 1) {next}

      # dummy data set used to set y to max starting concentration of all assays
      dummy <- csub.b %>%
        dplyr::mutate(mean_IC50_uM = max(mean_IC50_uM)) %>%
        dplyr::select(compound, label, mean_IC50_uM, FC_label) %>%
        dplyr::filter(label == k)

      # change units if max IC50 is low
      units <- bquote(IC[50] ~(mu*"M"))
      if (max(csub.b$mean_IC50_uM) < 0.1) {
        units <- bquote(IC[50]~ "(nM)")
        csub.b$mean_IC50_uM <- csub.b$mean_IC50_uM * 1000
        csub.b$stdev <- csub.b$stdev * 1000
      }
      csub.b <- csub.b %>%
        dplyr::mutate(label_p = cumsum(mean_IC50_uM) + stdev,
                      FC_label_f = factor(FC_label),
                      name_l = nchar(label))

      # scale y max using rsd
      mx <- max(csub.b$mean_IC50_uM) + max(csub.b$stdev)*1.1
      rsd <- stats::sd(csub.b$mean_IC50_uM)/mean(csub.b$mean_IC50_uM)*100
      if (rsd < 100) {mx <- (csub.b$mean_IC50_uM + csub.b$mean_IC50_uM) * 1.5 }
      mx <- max(mx)
      l[[j]] <- csub.b %>%
        ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(FC_label_f, name_l), y = mean_IC50_uM)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        #ggplot2::geom_blank(data = dummy) +
        ggplot2::geom_col(ggplot2::aes(fill = label), position = "dodge") +
        ggplot2::scale_y_continuous(limits = c(NA, mx))+
        ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_IC50_uM - stdev,
                                          ymax = mean_IC50_uM + stdev),
                                         width = 0.2,
                                      position = ggplot2::position_dodge(.9)) +
        ggplot2::geom_text(ggplot2::aes(y = label_p,
                                      label = symb,  vjust = 0.2, size = 2)) +
        ggplot2::theme_classic() +
        ggplot2::labs(title = i,
                     subtitle = paste0('Reference: ', k),
                     x = '', y = units) +
        ggplot2::theme(legend.position = 'none',
                          axis.text.x = ggplot2::element_text(angle = 0,
                                                 hjust = 0.55,
                                                 vjust = 0.5))

      j <- j + 1

    }
  }
  sig <- 'p < 0.001 =  ***, p < 0.01 = **, p < 0.05 = *'
  #generate and save pdf of all grouped plots
  group_name <- paste0(prefix,'_IC50_compare.pdf')

  suppressWarnings(
    plots <- gridExtra::marrangeGrob(l,
                                     nrow = grid.var,
                                     ncol = grid.var,
                                     top = prefix,
                                     bottom = sig ))
  suppressWarnings(
   f <- ggplot2::ggsave(group_name, plots,
                    width = 11, height = 11,
                    units = "in",
                    dpi = 600))
  cli::cli_inform(c("v" = "Export complete"))
  return(l)

}



