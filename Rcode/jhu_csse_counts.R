library(tidyverse)

cases_timeseries_global <- 
  read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")


country_rename <- function(vector){
  # vector <- sapply(vector, switch,
  #                    "He is"=1,
  #                    "She is"=1,
  #                    "He has"=2,
  #                    "She has"=2)
  
  vector <- ifelse(vector == "US", "USA",
                   ifelse(vector == "Congo (Kinshasa)", "Democratic Republic of the Congo",
                          ifelse(vector == "Czechia", "Czech Republic",
                                 ifelse(vector == "Taiwan*", "Taiwan",
                                        ifelse(vector == "Korea, South", "South Korea", as.character(vector))))))
  
  vector
}

case_counts <- cases_timeseries_global %>%
  select(-c(Lat, Long, Province.State)) %>%
  group_by(Country.Region) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "date",
               values_to = "case_counts") %>%
  mutate(date = str_sub(date, 2)) %>%
  mutate(date = as.Date(date, format = "%m.%d.%y")) %>%
  rename(country := Country.Region) %>%
  mutate(country = country_rename(country)) %>%
  {.}




