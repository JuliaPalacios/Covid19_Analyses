library(shiny)
library(phylodyn)
library(ape)
library(lubridate)
library(INLA)
library(shinycssloaders)
library(phangorn) # need for upgma(), called in function_serial

# All paths must be *relative* within the shinyapp so it will work in the 
# app when it's packaged & uploaded
source("util.R") # utils for UI and Server
data_dir <- get_latest_data_dir()

# Load data files and helper scripts and verify they all exist
function_serial_script <- file.path(data_dir, "function_serial.R")
#meta_fp <- file.path(data_dir, "all_meta.tsv")
cases_fp <- file.path(data_dir, "owid-covid-data.csv")
#trees_dir <- file.path(data_dir, "trees")
dists_dir <- file.path(data_dir, "distlists")
for (f in c(function_serial_script, cases_fp)) {
  stopifnot(file.exists(f))
}
stopifnot(dir.exists(dists_dir))
stopifnot(length(list.files(dists_dir)) > 0)
source(function_serial_script)

# ========== Load latest metadata ==========
#meta_fig <- read.delim(meta_fp, header = T, sep = "\t", as.is = T)
cases_all <- read.csv(file = cases_fp, header = T)


# ========== Load data (saved RData from initial_distances.R) ==========
countries <- list()
for (f in list.files(dists_dir)) {
  stopifnot(tools::file_ext(f) == 'RData')
  countries[parse_distlist_filename(f)] <- file.path(dists_dir, f)
}

# mutation_rate <- .001
