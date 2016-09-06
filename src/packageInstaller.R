#!/usr/bin/env Rscript
#---------------------------------------------------------------------------
# packages, dir and file locations, output preferences, table of inputs
#---------------------------------------------------------------------------

# Rtools will need to be installed for zipping of outputs.  Rtools can be
#  downloaded and installed manually from:
#  https://cran.r-project.org/bin/windows/Rtools/

# packages needed
my.packs <- c(
  'ggplot2', 'RANN', 'foreach', 'doParallel', 'scales', 'car', 'rgl', 
  'fields', 'data.table', 'dplyr', 'RODBC', 'XML', 'methods', 'readr'
)

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {
  install.packages(
    my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
    dependencies = TRUE,
    repos = "http://cran.us.r-project.org"
  )
}

print("All Done")