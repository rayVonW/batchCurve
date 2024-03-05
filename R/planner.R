

#' Designs dose-response assays from 96 well compounds plate for Tecan D300e
#' @param input_file A file name containing Compound assay parameters or guide IC50 values.
#' @return A data frame of completed dose ranges.
#' @export
plan_ranges <- function(input_file) {

 input <- read_csv_file(input_file)

 for (i in 1:nrow(input)) {

   if (is.na(input[i,'starting_uM'])) {

     # set default range/stock for unguided assays
     if (is.na(input[i,'guide_IC50'])) {
       cli::cli_inform(c("x" = "No reccomendations for {input[i,'compound']}",
                         "i" = "Setting 10uM with 1:2 dilution."))
       input[i,'starting_uM'] <- 10
       input[i,'dilution_factor'] <- 2
       input[i,'stock_mM'] <- 10
       next
     }

     # set ranges for assays with guides
     input[i,'dilution_factor'] <- 2

     input[i,'starting_uM'] <- ifelse(input[i,'guide_IC50'] >= 2.34, 50,
                                ifelse(input[i,'guide_IC50'] >= 0.937, 25,
                                  ifelse(input[i,'guide_IC50'] >= 0.1875, 10,
                                    ifelse(input[i,'guide_IC50'] >= 0.0375, 2,
                                      ifelse(input[i,'guide_IC50'] >= 0.0075, 0.4,
                                        ifelse(input[i,'guide_IC50'] >= 0.0015, 0.08,
                                          ifelse(input[i,'guide_IC50'] >= 0.0003,
                                                 0.016,0.0032)))))))
     input[i,'stock_mM'] <- ifelse(input[i,'guide_IC50'] >= 0.1875, 10,
                                  ifelse(input[i,'guide_IC50'] >= 0.0375, 1,
                                   ifelse(input[i,'guide_IC50'] >= 0.0075, 0.1,
                                    ifelse(input[i,'guide_IC50'] >= 0.0015, 0.02,
                                     ifelse(input[i,'guide_IC50'] >= 0.0003,
                                                   0.01,0.001)))))

   }
   #adjust stocks for set assays
   if (!is.na(input[i,'starting_uM'])) {
     input[i,'stock_mM'] <- ifelse(input[i,'starting_uM'] >= 10, 10,
                             ifelse(input[i,'starting_uM'] >= 2, 1,
                              ifelse(input[i,'starting_uM'] >= 0.4, 0.1,
                               ifelse(input[i,'starting_uM'] >= 0.08, 0.02,
                                ifelse(input[i,'starting_uM'] >= 0.016,
                                                               0.01,0.001)))))
   }
 }

 return(input)
}


#' Converts assay parameters to 384w plate layouts in excel for Tecan D300e printing
#'
#' @param data A data frame of Compound dilution curve information.
#' @param prefix A text string for output file names.
#' @return An export .xlsx file with each sheet representing a 384w plate for Tecan D300e.
#' @export
build_plates <- function(data, prefix = 'example') {
  meta <- data.frame(matrix(ncol = 9))
  colnames(meta) <- c('plate_id','position_id','format','replicates','index','compound','cell','starting_uM','dilution_factor')

  OUT <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(OUT, sheetName = 'ranges')
  for (i in 1:length(seq(1,nrow(data), 8))) {

    sheet_name <- paste("column", i)
    openxlsx::addWorksheet(OUT, sheetName = sheet_name)
  }
  openxlsx::addWorksheet(OUT, sheetName = 'meta')
  openxlsx::writeData(OUT, sheet = 'ranges', x = data)
  sh.n <- 1
  c <- 1
  for (i in seq(1,nrow(data), 8)) {
    s <- data[i:(i + 7),]
    s <- s[!is.na(s$starting_uM), ]

    start <- 4
    column <- 3
    counter <- 1
    infrow <- 1
    t <- 1

    p <- data.frame(matrix(ncol = 25))

    for (j in 1:nrow(s)) {

      if (counter == 5) {
        column <- 15
        start <- start - 12}
      p[start:(start + 2),column] <- s[j,'starting_uM']
      p[start:(start + 2),column+1] <- s[j,'starting_uM']/s[j,'dilution_factor']
      p[start:(start + 2),column+2] <- p[start:(start + 2),column+1]/s[j,'dilution_factor']
      p[start:(start + 2),column+3] <- p[start:(start + 2),column+2]/s[j,'dilution_factor']
      p[start:(start + 2),column+4] <- p[start:(start + 2),column+3]/s[j,'dilution_factor']
      p[start:(start + 2),column+5] <- p[start:(start + 2),column+4]/s[j,'dilution_factor']
      p[start:(start + 2),column+6] <- p[start:(start + 2),column+5]/s[j,'dilution_factor']
      p[start:(start + 2),column+7] <- p[start:(start + 2),column+6]/s[j,'dilution_factor']
      p[start:(start + 2),column+8] <- p[start:(start + 2),column+7]/s[j,'dilution_factor']
      p[start:(start + 2),column+9] <- p[start:(start + 2),column+8]/s[j,'dilution_factor']
      p[infrow,1] <- s[j,'well']
      p[infrow,2] <- s[j,'stock_mM']
      p[infrow,3] <- 'mM'
      p[infrow + 1,1] <- 'uM'
      p[(infrow + 2):(infrow + 17),1] <- LETTERS[1:16]
      start <- start + 22
      counter <- counter + 1
      infrow <- infrow + 19
      pos <- t
      if ( t >= 5) { pos <-  t + 1}

      meta[c, 'plate_id'] <- paste0('P-', sh.n)
      meta[c, 'position_id'] <- paste0('pos_', pos)
      meta[c, 'index'] <- s[j,'index']
      meta[c, 'compound'] <- s[j,'compound']
      meta[c, 'cell'] <- ''
      meta[c, 'starting_uM'] <- s[j,'starting_uM']
      meta[c, 'dilution_factor'] <- s[j,'dilution_factor']

      t <- t + 1
      c <- c + 1

    }
    meta['format'] <- '384w'
    meta['replicates'] <- 3

    openxlsx::writeData(OUT, sheet = 'meta', x = meta)

    sheet_name <- paste("column", sh.n)
    openxlsx::writeData(OUT, sheet = sheet_name, x = p, colNames = FALSE)
    sh.n <- sh.n + 1


  }
  openxlsx::saveWorkbook(OUT, paste0(prefix,"_print_plan.xlsx"), overwrite = TRUE)

}






