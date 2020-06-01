library(tidyverse)

cases_timeseries_global <- 
  read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

cases_timeseries_us <- 
  read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")

recovered_timeseries_global <- 
  read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")


################## Global manipulation

country_rename <- function(vector){
  
  vector <- ifelse(vector == "US", "USA",
                   ifelse(vector == "Congo (Kinshasa)", "Democratic Republic of the Congo",
                          ifelse(vector == "Czechia", "Czech Republic",
                                 ifelse(vector == "Taiwan*", "Taiwan",
                                        ifelse(vector == "Korea, South", "South Korea", as.character(vector))))))
  
  vector
}

cases_global <- cases_timeseries_global %>%
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


case_counts_diff_global <- cases_timeseries_global %>%
  select(-c(Lat, Long, Province.State))

case_counts_diff_global[,-1] <- t(apply(cbind(0,case_counts_diff_global[,-1]), 1, diff))


cases_increase_global <- case_counts_diff_global %>%
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


cases_recovered_global <- recovered_timeseries_global %>%
  select(-c(Lat, Long, Province.State)) %>%
  group_by(Country.Region) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "date",
               values_to = "case_counts_recovered") %>%
  mutate(date = str_sub(date, 2)) %>%
  mutate(date = as.Date(date, format = "%m.%d.%y")) %>%
  rename(country := Country.Region) %>%
  mutate(country = country_rename(country)) %>%
  {.}

cases_active_global <- cases_global %>%
  left_join(cases_recovered_global, by = c("country", "date")) %>%
  mutate(case_counts_active = pmax(0, case_counts - case_counts_recovered)) %>%
  select(-c(case_counts, case_counts_recovered)) %>%
  rename(case_counts := case_counts_active) %>% 
  {.}

###################################### US manipulation

cases_us <- cases_timeseries_us %>%
  select(Province_State, starts_with("X")) %>%
  group_by(Province_State) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "date",
               values_to = "case_counts") %>%
  mutate(date = str_sub(date, 2)) %>%
  mutate(date = as.Date(date, format = "%m.%d.%y")) %>%
  rename(state := Province_State) %>%
  {.}


case_counts_diff_us <- cases_timeseries_us %>%
  select(Province_State, starts_with("X"))

case_counts_diff_us[,-1] <- t(apply(cbind(0,case_counts_diff_us[,-1]), 1, diff))


cases_increase_us <- case_counts_diff_us %>%
  group_by(Province_State) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "date",
               values_to = "case_counts") %>%
  mutate(date = str_sub(date, 2)) %>%
  mutate(date = as.Date(date, format = "%m.%d.%y")) %>%
  rename(state := Province_State) %>%
  {.}



###################################### for export

case_counts_each_country <- cases_global

case_counts_all_together <- cases_global %>%
  group_by(date) %>%
  summarise(case_counts = sum(case_counts))


case_counts_diff_each_country <- cases_increase_global

case_counts_diff_all_together <- cases_increase_global %>%
  group_by(date) %>%
  summarise(case_counts = sum(case_counts))


case_counts_active_each_country <- cases_active_global

case_counts_active_all_together <- cases_active_global %>%
  group_by(date) %>%
  summarise(case_counts = sum(case_counts))


###########

case_counts_us_each_state <- cases_us

case_counts_us_all_together <- cases_us %>%
  group_by(date) %>%
  summarise(case_counts = sum(case_counts))


case_counts_diff_us_each_state <- cases_increase_us

case_counts_diff_us_all_together <- cases_increase_us %>%
  group_by(date) %>%
  summarise(case_counts = sum(case_counts))

