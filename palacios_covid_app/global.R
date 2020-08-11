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

start_date <- decimal_date(as.Date("2019-12-01"))
data_date <- decimal_date(parse_date_from_dir(data_dir))

# Load data files and helper scripts and verify they all exist
function_serial_script <- file.path(data_dir, "function_serial.R")
#meta_fp <- file.path(data_dir, "all_meta.tsv")
cases_fp <- file.path(data_dir, "owid-covid-data.csv")
cases_usa_fp <- file.path(data_dir, "cases.csv")
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
cases_all <- read.csv(cases_fp)
cases_usa <- read.csv(cases_usa_fp)

# Change column names to match with cases_all
names(cases_usa)[names(cases_usa)=="deaths"] <- "total_deaths"

# Delete entries in owid-covid-data where data were not recorded (NA values),
# because it confuses the creation of plots.
cases_all <- cases_all[-which(is.na(cases_all$new_cases)),]

# Calculate new cases per day in USA since it is not in the CSV
new_daily_cases <- integer(nrow(cases_usa))
for(i in 1:nrow(cases_usa)) {
  row <- cases_usa[i,]
  yest <- toString(as.Date(row$date) - 1)
  yest_row <- cases_usa[which((cases_usa$state == row$state) & (cases_usa$date == yest)),]
  yest_cases <- if (nrow(yest_row) > 0) as.numeric(yest_row$cases) else 0
  new_daily_cases[i] <- row$cases - yest_cases # today - yesterday cases
}
cases_usa$new_cases <- new_daily_cases

# ========== Load data (saved RData from initial_distances.R) ==========
countries <- list()
for (f in list.files(dists_dir)) {
  stopifnot(tools::file_ext(f) == 'RData')
  countries[parse_distlist_filename(f)] <- file.path(dists_dir, f)
}

# mutation_rate <- .001
