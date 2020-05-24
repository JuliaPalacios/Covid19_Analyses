library(shiny)
library(phylodyn)
library(ape)
library(lubridate)
library(INLA)
library(shinycssloaders)


# All paths must be *relative* within the shinyapp so it will work in the 
# app when it's packaged & uploaded
source("util.R") # utils for UI and Server
data_dir <- get_latest_data_dir()

# Load data files and helper scripts and verify they all exist
function_serial_script <- file.path(data_dir, "function_serial.R")
meta_fp <- file.path(data_dir, "all_meta.tsv")
cases_fp <- file.path(data_dir, "owid-covid-data.csv")
trees_dir <- file.path(data_dir, "trees")
for (f in c(function_serial_script, meta_fp, cases_fp)) {
  stopifnot(file.exists(f))
}
stopifnot(dir.exists(trees_dir))
stopifnot(length(list.files(trees_dir)) > 0)
source(function_serial_script)

# ========== Load latest metadata ==========
meta_fig <- read.delim(meta_fp, header = T, sep = "\t", as.is = T)
cases_all <- read.csv(file = cases_fp, header = T)

# ========== Load data (latest computed trees) ==========
trees <- list()
for (f in list.files(trees_dir)) {
  stopifnot(tools::file_ext(f) == "tre")
  tree_meta <- parse_tree_filename(f)
  tree <- read.tree(file.path(trees_dir, f))
  trees[[tree_meta$country]] <- list(tree = tree, lastdate = tree_meta$lastdate)
}

# Cache output of BNPR and BNPR_PS because they are slow. Per app session, 
# no country needs to be calculated twice. 
bnp_cache <- list()
bnp_ps_cache <- list()

# Optionally, fill in these caches when the app starts. This might be nice to
# do for published shinyapps.io, but will be annoying when playing with the 
# ShinyApp locally because it will cause it to take too long to start up.
# TODO: this fails because it takes too long (when first client connects). 
# I'm not sure how to fix, so must keep this FALSE for now. 
calc_all_bnpr <- F
if (calc_all_bnpr) {
  print("Populating BNPR and BNPR_PS caches so the plots load faster in-app.")
  remaining <- length(names(trees))
  for (country in names(trees)) {
    print(paste("Estimated time:", remaining*0.1, "minutes"))
    bnp <- BNPR(trees[[country]]$tree)
    bnp_cache[[country]] <- bnp
    bnp_ps <- BNPR_PS(trees[[country]]$tree)
    bnp_ps_cache[[country]] <- bnp_ps
    remaining <- remaining - 1
  }
  print("Finished populating caches.")
}
