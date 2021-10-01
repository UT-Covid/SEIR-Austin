library(tidyverse)

scenario.dir = "data/scenarios/"

source("code/common/small_func.R")

# Three scenarios:
# 1. No delta
growth1=c( shift=1000, mult=10)
# 2. growth as alpha. so 18 and 280
growth2=c( shift=280, mult=18)
# 3. Best fit. 10.5 and 256.5
growth3=c( shift=256.6, mult=10.5)

# The days are for reaching 50%, counting from 10.1

logit <- qlogis
logistic <- plogis

delta=function( x, growth=c(shift=256.6, mult=10.5)) {
  if( is.Date(x)) {
    x = as.numeric( x-da(10.1) )
  }
  logistic( ( x - growth["shift"] ) / growth["mult"])
}


school.day = function( d ) {
  ifelse( !chron::is.weekend( d ) & ( d > da( 8.21, y=21)), 1, 0 )
}


expand_state = function(names, dims) {
  dims=lapply(dims,seq.int)
  res = names
  for(l in dims)
    res = outer(res,l,paste,sep="_")
  c(res)
}

get.scenario = function( d, 
    par = list(
      vacc = "trend",
      dose = "doses",
      school = T,
      delta = "fit",
      vmult=1,
      dZ=0
    )
) {
    if( file.exists( paste0(scenario.dir,par$vacc) )) {
      vac =read_csv( paste0(scenario.dir,par$vacc) )[,-1]
    } else {
      if( par$vacc == "optimistic" & par$dose=="one") 
        vac =read_csv("data/scenarios/scenario_optimistic_one_dose.csv")[,-1]
      if( par$vacc == "trend" & par$dose=="one") 
        vac =read_csv("data/scenarios/scenario_trend_one_dose.csv")[,-1]
      if( par$vacc == "optimistic" & par$dose=="doses") 
        vac =read_csv("data/scenarios/scenario_optimistic_total_doses.csv")[,-1]
      if( par$vacc == "trend" & par$dose=="doses") 
        vac =read_csv("data/scenarios/scenario_trend_total_doses_jan1_start.csv")[,-1]
      if( par$vacc == "optimistic" & par$dose == "full") 
        vac =read_csv("data/scenarios/scenario_optimistic_fully_vacc.csv")[,-1]
      if( par$vacc == "trend" & par$dose == "full") 
        vac =read_csv("data/scenarios/scenario_trend_fully_vacc.csv")[,-1]
    }

    vac = vac %>% mutate( across( !date, as.integer) ) 
print(par)


    # diff starting from 0 - keeps vector length
    diff0 = function(x) {
      x=c(0,x)
      diff(x) %>%  ifelse( .<0, 0, .)
    }

    
    # Add empty rows before first vaccination
    min.vac.date = as.Date(min( vac$date))
    min.d = min(d)
    if( min.d < min.vac.date) {
      empty_vac = data.frame( date = seq( from=min.d, to=min.vac.date-1, by=1))
      empty_vac[,colnames(vac)[-1]] = 0
      vac = as_tibble( empty_vac ) %>% add_row( vac )
    }

    max.vac.date = as.Date(max( vac$date))
    max.d = max(d)

    if( max.d >  max.vac.date) {
      future_vac_row = tail( vac, 1)
      new_d = d[ d>max.vac.date ]
      vac = vac %>%
        bind_rows( purrr::map( new_d, ~ {
          D = future_vac_row
          D$date = .x
          D
        } ))
    }


    res = 
      vac %>% #
      mutate( date = as.Date( date ), .keep="unused") %>%  #
      mutate( across( !date, as.integer) ) %>%
      mutate( across( !date, diff0) ) %>%
       filter( date %in% d) 


    # If we don't have enough points, continue vaccinmations at same rate
    max.vac.date = as.Date(max( vac$date))
    max.d = max(d)

    if( max.d >  max.vac.date) {
      future_res_row = tail( res, 1)
      new_d = d[ d>max.vac.date ]
      res = res %>%
        bind_rows( purrr::map( new_d, ~ {
          D = future_res_row
          D$date = .x
          D
        } ))
    }

  
    res[res$date > da(7.2,y=21),-1] =     res[res$date > da(7.2,y=21),-1] * as.numeric(par$vmult)



    expand_state("new_V",c(5,2)) %>% order %>% order ->o 

    colnames(res)[(1:10)+1] = expand_state("new_V",c(2,5))[o]


    if( par$delta == "fit") growth.par = growth3
    if( par$delta == "alpha") growth.par = growth2
    if( par$delta == "none") growth.par = growth1

    
    res = res %>%
      bind_cols( variant_p =delta( res$date, growth=growth.par ) ) %>%
        mutate( dZ = par$dZ ) 

    if( "delta_IHR" %in% names(par) ) {
      res = res %>%
        mutate( var_YHR_mult = par$delta_IHR )
    } else {
      mutate( var_YHR_mult = 1 )
    }

    res = res %>%
      bind_cols( school    =school.day( res$date) * par$school ) 

    res
}



if( !exists( "par.index")) {
    par.index = 1
}

parameters.for.run = function( par.list = list(
      vacc = c(
        "data/scenarios/scenario_increased_total_doses_jan1_start.csv",
        "data/scenarios/scenario_trend_total_doses_jan1_start.csv"
      ), # trend / optimistic
      dose = c("doses"),   # one / full
      school = c(T,F),      
      delta = c("none", "fit","alpha"),  # fit: US fit / alpha: US alpha growth for delta / none: no delta
      vmult=c(1,2,4),
      dZ=c(0,-0.6) # effect of masking -0.6 is 45% reduction (log(0.55))
    ) ) {
    s = ""
    for( n in names( par.list) ) {
        v = paste( n, par.list[[n]], sep="=")
        s = c(outer( s, v, paste, sep=","))
    }
    x = gsub("^[,]","",s)
    y=strsplit(x,"[,]")
    sl=lapply(y,function(s) {
        r=strsplit(s,"=");
        res=lapply(r,function(x) x[[2]] )
        names(res)=sapply(r,function(x) x[[1]] )
        for( s in names(res))  res[[s]] = as(res[[s]], class( par.list[[s]])) 
        res
    } )
    list(pars=sl,names=x)
}

