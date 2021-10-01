library(tidyverse)
library(lubridate)

source("code/common/conditional_run.R")


get_mobility_data = function (average_mode=c("trailing", "centered"),
                              days_average=7, rerun=F,
                              top_categories = c(
                                  "Colleges, Universities, and Professional Schools"="colleges",
                                  "Drinking Places (Alcoholic Beverages)"="drinking",
                                  "Elementary and Secondary Schools"="schools",
                                  "General Medical and Surgical Hospitals"="medical",
                                  "Grocery Stores"="grocery",
                                  "Museums, Historical Sites, and Similar Institutions"="museums/parks",
                                  "Restaurants and Other Eating Places"="restaurants"),
                                  sd_cols = c("median_home_dwell_time_relative")
                              ) {
  #browser()
  if (length(average_mode) > 1)
    average_mode = average_mode[1]
  
  ## read latest mobility data
  if( grepl("tacc",system("hostname",intern=T)) ) {
    SG.dir = "/work2/projects/utprojections/safegraph_data/final-metrics/"
  } else {
    SG.dir = "data/safegraph_data/final_metrics/"
  }
  visits.f = paste0(SG.dir,"visits_by_category_bycounty_baselined.csv")
  sdmetric.f = paste0(SG.dir,"sd-metrics_bycounty_baselined.csv")
  
  res.f = "data/mobility_msa_pca.Rda"
  
  
  #### Check if data in files is newer than our latest data
  if( length(rerun)==0) {
      need.run =  check.files.for.rerun.needed( files=c( visits.f, sdmetric.f), data.file = res.f)
   }


  ## if not, just load it from before
  
  if( (length(rerun)==0 && !need.run) || (length(rerun)>0 && !rerun) ) {
      load( res.f)
      return( res)
  } 

  visits =    read_csv( visits.f   )[ -1]
  sdmetrics = read_csv( sdmetric.f )[,-1]
  
  ## Load in MSA definitions and population weights
  county_population = "https://raw.githubusercontent.com/JieYingWu/COVID-19_US_County-level_Summaries/master/data/counties.csv" %>% 
    data.table::fread() %>% 
    as_tibble() %>% 
    select(fips=FIPS, state_short=State, county=Area_Name, pop=POP_ESTIMATE_2018) %>% 
    mutate(fips=sprintf("%05d", fips))
  
  rural_urban = read_csv("./data/counties_basicdata.csv")[, -1] %>% 
    mutate(fips = sprintf("%05d", fips)) %>% 
    mutate(rural_urban = Rural.urban_Continuum.Code_2013 / 9) %>% 
    select(fips, rural_urban)
  
  msa = read_csv("./data/csa.csv") %>% 
    filter(`Metropolitan/Micropolitan Statistical Area` == "Metropolitan Statistical Area",
           `State Name` != "Puerto Rico") %>% 
    mutate(fips = paste0(`FIPS State Code`, `FIPS County Code`)) %>%
    select(msa = `CBSA Title`,
           county = `County/County Equivalent`,
           state = `State Name`,
           central.outlying = `Central/Outlying County`,
           fips) %>%
    mutate(county = county %>%
             str_remove(" County") %>%
             str_remove(" city") %>%
             str_remove(" Parish"))  %>% 
    left_join(select(county_population, -county), by="fips") %>% 
    left_join(rural_urban, by="fips") %>% 
    left_join(group_by(., msa) %>% 
              summarize(
                msa_pop=sum(pop),
                rural_urban_msa=weighted.mean(rural_urban, pop, na.rm=TRUE)) %>% 
              ungroup() %>% 
              mutate(pop_rank=rank(-msa_pop)),
      by="msa") %>%
    filter(pop_rank <= 300) %>%  # keep only largest 300 areas
    mutate(county_weight = pop / msa_pop)
  
  msa_totals = msa %>% 
    select(msa, state, state_short, msa_pop, rural_urban_msa) %>% 
    distinct()

  # write_csv(msa, "data/msainfo.csv")
  
  ## death data
  death_data = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv' %>% 
    read_csv() %>% 
    left_join(select(msa, fips, msa), by="fips") %>% 
    group_by(msa, date) %>% 
    summarize(cumulative_deaths = sum(deaths)) %>% 
    ungroup() %>% 
    group_by(msa) %>% 
    arrange(date) %>% 
    mutate(deaths = pmax(0, c(0, diff(cumulative_deaths)))) %>% 
    ungroup()
  
  cases_data = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" %>% 
    read_csv()  %>%
    pivot_longer(cols = ends_with("21") | ends_with("20"),   ### HARD coded for year!!
                 names_to = "date_mdy",
                 values_to = "cumulative_cases") %>%
    mutate(date_ymd = mdy(date_mdy) %>% ymd()) %>%
    select(date = date_ymd,
           county = Admin2,
           state = Province_State,
           fips = FIPS,
           cumulative_cases) %>% 
           filter(!is.na(fips)) %>%
           mutate( fips=as.numeric(fips)) %>% 
    mutate(fips=sprintf("%05d", fips)) %>%  
    left_join(select(msa, fips, msa), by="fips") %>% 
    group_by(msa, date) %>% 
    summarize(cumulative_cases = sum(cumulative_cases)) %>%
    ungroup() %>% 
    filter(msa %in% unique(.env$msa$msa)) %>% 
    group_by(msa) %>% 
    arrange(msa, date) %>% 
    mutate(cases = pmax(0, c(0, diff(cumulative_cases)))) %>% 
    ungroup()

  
  county_poi = visits %>% 
    mutate(date = ymd(date)) %>%
    rename(fips = fips5d,
           state_short = region,
           `county_full` = County) %>% 
    filter(fips %in% unique(msa$fips)) %>% 
    left_join(select(msa, fips, msa), by="fips") %>%
    mutate(top_category = top_categories[top_category]) %>% 
    group_by(date, top_category, msa) %>% 
    summarize(pop = sum(pop),
              count = sum(count),
              count_baseline = sum(count_baseline)) %>% 
    mutate(count_pc = count / pop,
           count_relative = count / count_baseline) %>% 
    ungroup()
  
  sd = sdmetrics %>% 
    select(date, fips=fips5d, all_of(sd_cols)) %>% 
    left_join(county_population, by="fips") %>% 
    filter(fips %in% unique(msa$fips)) %>% 
    left_join(select(msa, fips, msa, county_weight), by="fips") %>% 
    group_by(date, msa) %>% 
    summarize(median_home_dwell_time_relative = sum(median_home_dwell_time_relative * county_weight))
    # summarize(median_home_dwell_time = sum(median_home_dwell_time * county_weight),
    #           full_time_work_behavior_devices = sum(full_time_work_behavior_devices))
  
  cases_msa = cases_data %>% 
    filter(msa=="Austin-Round Rock-Georgetown, TX"  )
  
  msa_data = county_poi %>%
#    filter(msa=="Austin-Round Rock-Georgetown, TX"  ) %>%
#    filter( date > ymd("2020-09-27")) %>% filter( date < ymd("2020-10-04")) %>% 
    select(date, msa, top_category, count_pc, count_relative) %>% 
    pivot_wider(names_from = top_category, values_from = c(count_pc, count_relative)) %>%
    left_join(sd, by=c("date", "msa")) %>% 
    # na.omit() %>% 
    group_by(msa) %>% 
    group_modify(~ { # add 7-week average
      df = arrange(.x, date)
      # add T-day running averages
      cur_cols = c(paste0("count_relative_", unname(top_categories)), sd_cols)
      new_cols = paste0(cur_cols, "_av")
      for (c in new_cols)
        df[[c]] = NA_real_  # df[[c]] = NA gave errors - could not add doubles later.
      offset = ifelse(average_mode == "trailing", 0, days_average %/% 2 + 1)
      for (i in days_average:nrow(df)) {
        start_index = i - days_average
        end_index = i
        df[i - offset, new_cols] = t(apply(df[start_index:end_index, cur_cols, drop=FALSE], 2, mean)) # added t() to avoid errors.
      }
      # output of group modify
      df
    }) %>% 
    ungroup() %>% 
    left_join(msa_totals, by="msa") %>% 
    left_join(death_data, by=c("msa", "date")) %>% 
    right_join(cases_data, by=c("msa", "date")) %>% 
    mutate(deaths_per_cap = deaths / msa_pop,
           cumulative_deaths_per_cap = cumulative_deaths / msa_pop) %>% 
    left_join(
      group_by(., msa) %>%
      summarize(threshold_day = date[min(which(cumulative_deaths_per_cap >= 3 / 1e7))]),
      by = "msa") %>% 
    mutate(days_since_thresh = as.integer(date - threshold_day)) %>% 
    mutate(weekend = chron::is.weekend(date)) # %>%
#    left_join( cases_msa, by=date)
    
        
#    list(msa_data=msa_data,
#         pca=pca,cases = cases_data %>% filter(msa=="Austin-Round Rock-Georgetown, TX"  ))
    res = list(msa_data=msa_data
               )
    save( res, file=res.f)
    return(res)
}
