library(shiny)
library(phylodyn)
library(ape)
library(lubridate)
library(INLA)

# All paths must be relative within the shinyapp so it will work when uploaded
source("util.R") # utils for UI and Server
source("function_serial.R") # TODO: Copied from JuliaPalacios GitHub -- may be stale
data_dir <- "./data"
meta_fp <- file.path(data_dir, "all_meta.tsv")
cases_fp <- file.path(data_dir, "owid-covid-data.csv")

# ========== Load metadata ==========
meta_fig <- read.delim(meta_fp, header = T, sep = "\t", as.is = T)
cases_all <- read.csv(file = cases_fp, header = T)

# ========== Load data ==========
# get most recent computed trees
trees_dir <- list.dirs(path = data_dir, recursive = F)
trees_dir <- sort(trees_dir, decreasing = T)[[1]]
trees <- list()
for (f in list.files(trees_dir)) {
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
