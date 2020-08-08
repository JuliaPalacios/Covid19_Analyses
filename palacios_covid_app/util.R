
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

  return(country)
}


parse_tree_filename <- function(tree_filename) {
  split <- strsplit(tree_filename, "[_]")
  loc <- split[[1]][1]
  dateext <- split[[1]][2]
  datestr <- substr(dateext, 1, nchar(dateext) - 4) # '.tre' ext is 4 chars
  return(list(country = loc, lastdate = as.numeric(datestr)))
}
