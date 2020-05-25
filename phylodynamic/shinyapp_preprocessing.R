library(phylodyn)
library(ape)
library(phangorn)
library(lubridate)
library(INLA)

# Env settings
base_dir <- "~/script_dev/juliapr/Covid19_Analyses"
delete_fasta_subset_files <- TRUE

# Files and scripts from elsewhere in the repository. These should not change 
# between user envs unless git fs changes
copy_fasta_subset_script <- file.path(base_dir, "alignment/code/copy_lines.pl")
function_serial_script <- file.path(base_dir, "phylodynamic/function_serial.R")
data_dir <- file.path(base_dir, "alignment/data")
gisaid_aligned_fp <- file.path(data_dir, "all_seq.fasta")
meta_fp <- file.path(data_dir, "all_meta.tsv")
cases_fp <- file.path(data_dir, "owid-covid-data.csv")
app_dir <- file.path(base_dir, "palacios_covid_app")

# Make sure everything we need exists before proceeding
for (f in c(copy_fasta_subset_script, function_serial_script, gisaid_aligned_fp, 
            meta_fp, cases_fp)) {
  stopifnot(file.exists(f))
}
stopifnot(dir.exists(file.path(app_dir, "data")))
source(function_serial_script)

# ========== Local functions ==========
subset_fasta <- function(subset_fp, inds) {
  # Adapted from subset_data.R
  fasta_lines <- c(rbind(2*inds-1, 2*inds))

  fasta_lines_tmp_fp <- paste(subset_fp, "_lines.bak", sep="")
  write.table(fasta_lines, file=fasta_lines_tmp_fp, col.names=FALSE, 
      row.names=FALSE, quote=FALSE)

  script_call <- paste("perl", copy_fasta_subset_script, fasta_lines_tmp_fp, 
      gisaid_aligned_fp, ">", subset_fp)
  
  system(script_call)
  file.remove(fasta_lines_tmp_fp)
}

compute_tree <- function(country, division = NULL) {
  if (is.null(division)) {
    subs2 <- seq(1, nrow(meta_fig))[meta_fig$country == country]
  } else {
    subs2 <- seq(1, nrow(meta_fig))[meta_fig$country == country &&
      meta_fig$division == division]
  }
  subs <- c(subs2, idx_root)

  fasta_subset_fp <- file.path(data_dir, paste("fastasub_", 
      gsub(" ", "-", country), "_", format(Sys.time(), "%Y%m%d%H%M%S", 
      tz = "UTC"), ".fasta", sep = ""))
  print(paste("Creating subset fasta file", fasta_subset_fp))
  subset_fasta(fasta_subset_fp, subs)
  stopifnot(file.exists(fasta_subset_fp))

  print("Importing fasta file into R")
  gisaid_aligned_country <- ape::read.FASTA(fasta_subset_fp)
  fastafile <- phangorn::as.phyDat(gisaid_aligned_country) # includes root
  ref <- which(names(fastafile) == 
                 "hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31")
  
  # Note: important to include ref at the end, because the function to estimate
  # mutation rate expects this order. Also, note that you need to be careful
  # with how you copy/collate/produce these vars. For example, 
  # c(fastafile, fastafile[idx]) will yield a different dtype (even though it
  # looks the same when being viewed and has the same `typeof` value).
  fastafile <- fastafile[c(seq(1, length(fastafile))[-ref], ref)]
  fastafile_2 <- fastafile[-length(fastafile)] # exclude ref

  seq_names <- names(fastafile_2)
  samp_times <- c()
  for (r in seq_names) {
    # the format of this line has changed in prev iterations of the fasta file
    samp_times <- c(samp_times, strsplit(r, "[|]")[[1]][3])
  }
  print(paste("Parsed sampling times from fasta, e.g.", samp_times[[1]]))
  samp_times <- lubridate::decimal_date(lubridate::date(samp_times))
  lastdate <- max(samp_times)
  samp_times <- lastdate - samp_times # normalize by lastdate
  name_samp <- cbind(samp_times, seq_names) # 2-column matrix (date & seq name)
  
  print("Estimating mutation rate")
  mu <- mu_linear_reg(fastafile)
  print(paste("mu =", mu))
  tree <- serial_upgma(fastafile_2, mu, samp_times, name_samp, model = "F81")
  stopifnot(length(tree) == 4) # verify tree has been constructed properly

  if (delete_fasta_subset_files) {
    file.remove(fasta_subset_fp)
  }

  return(list(tree = tree, lastdate = lastdate))
}

summarize_metadata <- function() {
  print(paste("A total of ", nrow(meta_fig), " samples from ",
        length(unique(meta_fig$country)), "countries --", date(), collapse=""))
  par(mfrow=c(1,2))
  plot(sort(table(meta_fig$country)),cex.axis = 0.35,las=2,
       xlab="",ylab="",main="Samples per country")
  plot(sort(table(meta_fig$division[meta_fig$country=="USA"])), cex.axis=.35, 
       las=2, xlab="", ylab="", main="Samples per US state")
}

check_country_name_consistency <- function() {
  # Check for any instances where the country names do not match between the
  # GISAID metadata dataset and the "owid-covid-data" cases dataset.
  cases_all <- read.csv(file = cases_fp, header = T)
  unmatched <- c()
  meta_countries <- unique(meta_fig$country)
  cases_countries <- unique(cases_all$location)
  for (c in meta_countries) {
    if (!is.element(c, cases_countries)) {
      unmatched <- c(unmatched, c)
    }
  }
  print("Could not find exact match country names between datasets for:")
  print(unmatched)
  print(paste("Make sure all unmatched countries are properly handled in the",
              "app by util::reformat_country."))
}


# ========== Load data ==========
meta_fig <- read.delim(meta_fp, header = T, sep = "\t", as.is = T)
summarize_metadata()

idx_root <- which(meta_fig$strain == "Wuhan-Hu-1/2019")

app_data_dir <- file.path(app_dir, "data", format(Sys.time(),
    "%Y%m%d%H%M%S", tz = "UTC"))
dir.create(app_data_dir)
trees_dir <- file.path(app_data_dir, "trees")
dir.create(trees_dir)

for (c in countries) {
  nseq <- length(which(meta_fig$country == c))
  if (nseq < 20 || nseq > 10000) {
    print(paste0("----- Skipping ", c, " (", nseq, " seqs) -----"))
  } else {
    print(paste0("----- Computing tree for ", c, " (", nseq, " seqs) -----"))
    treedata <- compute_tree(c)
    filename <- paste(c, "_", treedata$lastdate, ".tre", sep = "")
    write.tree(treedata$tree, file = file.path(trees_dir, filename))
  }
}

check_country_name_consistency()

# ========== copy/replace source files & data in the app ==========
# The ShinyApp needs to be standalone so any local dependencies must be copied
# into the app's data dir. We won't track any of it though because it will all
# be redundant and will be overwritten every time this file is run. 
file.copy(from=c(function_serial_script, meta_fp, cases_fp),
          to=app_data_dir, overwrite=T, recursive=F, copy.mode=T)
