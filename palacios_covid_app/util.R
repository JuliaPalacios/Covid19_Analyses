
get_latest_data_dir <- function() {
  data_dirs <- list.dirs(path = "./data", recursive = F)
  latest_data_dir <- sort(data_dirs, decreasing = T)[[1]]
}

reformat_country <- function(country) {
  if (country == "USA") {
    return("United States")
  }
  if (country == "Democratic Republic of the Congo") {
    return("Democratic Republic of Congo")
  }
  if (startsWith(country, 'UK - ')) {
    return('United Kingdom')
  }
  if (startsWith(country, 'USA - ')) {
    state <- substr(country, 7, nchar(country))
    if (state == 'Washington DC'){
      return('District of Columbia')
    }
    return(state)
  }

  return(country)
}

is_state <- function(country) {
  # return(length(strsplit(country, "[ - ]")[[1]]) > 1)
  return((startsWith(country, 'USA - ')))
}


parse_tree_filename <- function(tree_filename) {
  split <- strsplit(tree_filename, "[_]")
  loc <- split[[1]][1]
  dateext <- split[[1]][2]
  datestr <- substr(dateext, 1, nchar(dateext) - 4) # '.tre' ext is 4 chars
  return(list(country = loc, lastdate = as.numeric(datestr)))
}

parse_distlist_filename <- function(filename) {
  # Naming convention should change so that this is more possible to do
  name <- substr(filename, 1, nchar(filename) - 10) # 'dist.RData' suffix 10 chrs
  split <- strsplit(name, "[+]")
  
  prefix <- ''
  loc <- ''
  if (length(split[[1]]) == 2) {
    loc <- split[[1]][2]
    # this is a state - for now, special case UK and US
    if (split[[1]][1] == 'UnitedKingdom') {
      prefix <- 'UK - '
    } else if (split[[1]][1] == 'USA') {
      prefix <- 'USA - '
    }
  } else {
    loc <- split[[1]][1]
  }
  return(paste0(prefix, gsub('_', ' ', loc)))
}

compute_tree <- function(country, mu) {
  # TODO no cache for now
  
  # TODO consider miceadds#load.rdata2
  load(countries[[country]])
  dists <- listout
  
  # dists$n is the sample size
  # dists$seq_names should be the names of the sequence
  # dists$hamming is a R dist object of hamming matrices
  # dists$distGen is a R dist object of a genetic distance of choice
  
  samp_times <- c()
  for (r in dists$seq_names) {
    # the format of this line has changed in prev iterations of the fasta file
    samp_times <- c(samp_times, strsplit(r, "[|]")[[1]][3])
  }
  print(paste("Parsed sampling times from fasta, e.g.", samp_times[[1]]))
  samp_times <- lubridate::decimal_date(lubridate::date(samp_times))
  lastdate <- max(samp_times)
  samp_times <- lastdate - samp_times # normalize by lastdate
  name_samp <- cbind(samp_times, dists$seq_names) # 2 cols (date & seq name)

  # Ensure that there are at least 2 seqs from the most recent sampling date.
  # While there are <2 seqs on last date, remove the seq.
  while (table(samp_times)[1] == 1) {
    idx <- which(samp_times == 0)
    # remove the sequence observed only once
    dists$seq_names <- dists$seq_names[-idx]
    samp_times <- samp_times[-idx]
    name_samp <- name_samp[-idx,]
    # make the last sequence the new 0
    samp_times<-samp_times-min(samp_times)
    dists$hamming <- dists$hamming[-idx]
    dists$distGen <- dists$distGen[-idx, -idx]
  }
  
  # TODO not computing mutation rate for now.
  # mu <- mu_linear_reg_inputDist(dists)
  
  tree <- serial_upgma_inputDist(dists, mu, samp_times, name_samp)
  
  return(list(tree = tree, lastdate = lastdate))
}
