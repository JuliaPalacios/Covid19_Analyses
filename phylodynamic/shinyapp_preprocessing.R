library(phylodyn)
library(ape)
library(phangorn)
library(lubridate)
library(INLA)

# Env settings
base_dir <- "~/script_dev/juliapr/Covid19_Analyses"
write_tree_files <- TRUE
delete_fasta_subset_files <- TRUE

# These should not change between user envs unless git fs changes
data_dir <- file.path(base_dir, "alignment/data")
gisaid_aligned_fp <- file.path(data_dir, "all_seq.fasta")
meta_fp <- file.path(data_dir, "all_meta.tsv")
app_dir <- file.path(base_dir, "palacios_covid_app")

# Load other code from repo
source(file.path(base_dir, "phylodynamic/function_serial.R"))
source(file.path(base_dir, "alignment/code/subset_data.R"))


compute_tree <- function(country, division = NULL) {
  if (is.null(division)) {
    subs2 <- seq(1, nrow(meta_fig))[meta_fig$country == country]
  } else {
    subs2 <- seq(1, nrow(meta_fig))[meta_fig$country == country &&
      meta_fig$division == division]
  }
  subs <- c(subs2, idx_root)

  fasta_subset_fp <- file.path(data_dir, paste("fastasub_", format(Sys.time(),
    "%Y%m%d%H%M%S", tz = "UTC"), ".fasta", sep = ""))
  subset.fasta(paste(base_dir, "/", sep = ""), subs, fasta_subset_fp)
  stopifnot(file.exists(fasta_subset_fp))

  print("Importing fasta file into R")
  gisaid_aligned_country <- ape::read.FASTA(fasta_subset_fp)
  fastafile_1 <- phangorn::as.phyDat(gisaid_aligned_country) # includes root
  fastafile_2 <- fastafile_1[-length(fastafile_1)] # does not include root

  print("Parsing sampling times from fasta")
  seq_names <- names(fastafile_2)
  samp_times <- c()
  for (r in seq_names) {
    # the format of this line has changed in prev iterations of the fasta file
    samp_times <- c(samp_times, strsplit(r, "[|]")[[1]][3])
  }
  print(paste("e.g.", samp_times[[1]]))
  samp_times <- lubridate::decimal_date(lubridate::date(samp_times))
  lastdate <- max(samp_times)
  samp_times <- lastdate - samp_times # normalize by lastdate
  name_samp <- cbind(samp_times, seq_names) # 2-column matrix (date & seq name)

  mu <- mu_linear_reg(fastafile_1)
  print(paste("mu =", mu))
  # fastafile_2 <- fastafile_1[-length(fastafile_1)] # TODO is this any different?
  tree <- serial_upgma(fastafile_2, mu, samp_times, name_samp, model = "F81")
  stopifnot(length(tree) == 4) # verify tree has been constructed properly

  if (delete_fasta_subset_files) {
    file.remove(fasta_subset_fp)
  }

  return(list(tree = tree, lastdate = lastdate))
}

# ========== Load data ==========
meta_fig <- read.delim(meta_fp, header = T, sep = "\t", as.is = T)

idx_root <- which(meta_fig$strain == "Wuhan-Hu-1/2019")

# Currently there's no valid workflow where this flag could be false, but
# it could potentially be useful in the future for auto upload if this
# functionality is moved into the app.
if (write_tree_files) {
  trees_dir <- file.path(app_dir, "data", format(Sys.time(),
      "%Y%m%d%H%M%S", tz = "UTC"))
  dir.create(trees_dir)
}
for (country in sort(unique(meta_fig$country))) {
  if (length(which(meta_fig$country == country)) < 20) {
    print(paste("Skipping ", country, "(<20 seqs)"))
  } else {
    print(paste("Computing tree for", country))
    treedata <- compute_tree(country)
    if (write_tree_files) {
      filename <- paste(country, "_", treedata$lastdate, ".tre", sep = "")
      write.tree(treedata$tree, file = file.path(trees_dir, filename))
    }
  }
}

# ========== copy/replace source files & data in the app ==========
# Need to copy some files over to the ShinyApp.