
# batchCurve

## Overview

batchCurve a is bespoke analytical workflow for batch analysis of
dose-response assays. Assays formats must have a 10 point serial
dilution as well as a positive and negative growth control. The accepted
layouts assume either duplicate or triplicate measurement in either 96
or 384 well culture plates.

A batch of drug assays can be a combination of 96 and/or 384 well
plates. It is expected that there will be one raw plate reader csv file
per assay plate in the analysis folder. Assays that are on the same
plate must have the same number of replicates. more details of setup can
be found in `vignette("dose response analysis")`.

## Installation

There is an install_github function in devtools package that allows
users to install R packages hosted on GitHub. It requests developer’s
name/package name.

``` r
#First install devtools from CRAN
install.packages("devtools")
# load devtools
library(devtools)
# use function to install a github package
install_github("rayVonW/batchCurve")
```

## Usage

You can load the package with:

``` r
library(batchCurve)
#> _______________________________________
#>   This is version 0.1.0 of batchCurve
#> ____________________________________________
#> Available functions:
#> batchCurve_example
#> fit_data
#> plot_fit
```

There are two main functions to use:

1.  `fit_data()`: Fits log logistic models to a batch of raw
    dose-response data

- This expects as an argument a **meta.csv** file name, this file is to
  be found in a directory with all raw data file.
- A **batch id** text name - to prefix results files with.

``` r
#all meta and raw data stored in ~/batch
fit_data(file_path = '~/batch/meta-data.csv', 
         batch_id = 'batch01')
```

`vignette("dose-response-analysis")` gives an expanded introduction to
fitting batches of data.

2.  `plot_fit()`: Accepts batch log-logistic coefficients and normalised
    raw data to produce visualisations of dose-response assay data and
    allow assessment of their fit.

- It accepts a data frame of model **results** and a data frame of
  normalised **data**, both are exported directly from `fit_data()`.
- User can include a batch_id for linking visualisation to results in
  file prefixes.

``` r
#all meta and raw data stored in ~/batch
plot_fit(results = df1, 
         data = df2,
         batch_id = 'batch01')
```

`vignette("visualise-models")` gives an expanded introduction to
plotting out results.

## Acknowledgements

- `fit_data()` is a wrapper function for the R package \[drc\]
  (<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819>).
