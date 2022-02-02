# Main file for running Austin projections

#################################
#read command line arguments
args = commandArgs(trailingOnly=TRUE)
arg.list=c()

if( length(args)==1 && args[1] == "help") { #@note help
      cat("
Run Austin SEIR projection
## Commannd line arguments

First argment can be a modifier, or missing.
Allowed modifiers are `help`, `history`, `continue`, `short`

### Regular run: 
`Rscript run-seir.R Admit_discharge_datafile.xlsx ICY_datafile.xlsx`

### Continuing previous run:
`Rscript run-seir.R continue Rdata_session_file.Rda`
stage of run to continue is determined by name of file before/after stage.

`Rscript run-seir.R  continue Rdata_session_file.Rda after:stage`
continue run after stage. Possible stages: prepare, mf2, smooth, sim, analysis.

### History runs for performance analysis
`Rscript run-seir.R history Admit_discharge_datafile.xlsx [ndays]`
Load Admit_discharge_datafile.xlsx data file, but cut to first [ndays].

### Short run
`Rscript run-seir.R short Admit_discharge_datafile.xlsx ICY_datafile.xlsx`
Don't start run from beginning of hospitalization data, but instead load initial state from a previous run for a later date.
That run has to be carefully prepared top have enough initial variation and started at a relatively quiet part of the epidemic.
The start session is currently hard coded.
Short runs were not used for the PNAS paper analysis.

### Additional arguments
Additional aguments come in the form of keyword:argument

`no:dash` don't produce dashboard csv files
`no:pdf` don't produce pdf report
`no:mob` run without mobility data

`name.add:string` add string to output file name
")
}


timing=list()
timing$start=Sys.time()


library(session)

library(foreach)
library(doParallel)


library(tidyverse)
library(readxl)
library(pomp, quietly = T)
library(lubridate, quietly = T)
library(abind)

# Set size of cluster depending on computer name: clust.n. This will affect the number of parallel mf2 runs done.

if( grepl("frontera",system("hostname",intern=T)) ) {
  clust.n = 50
  clust.n.fork = 10
dash_dir = "/work2/07327/lachmann/frontera/austin-dash/data/"

} else if( grepl("stampede",system("hostname",intern=T)) ) {
  clust.n = 90
  clust.n.fork = 20
dash_dir = "/work2/07327/lachmann/frontera/austin-dash/data/"
} else if( grepl("mini",system("hostname",intern=T)) ) {
  clust.n = 6
  clust.n.fork = 3
dash_dir = "/Users/michael/cv19/austin-dash/data/"
} else {
  clust.n = 14
  clust.n.fork = 20
dash_dir = "/home/michael/cv19/austin-dash/data/"
}


# @note mif parameters
mif2.N.steps=300
mif2.iter = 10 # steps per single loop. Does not affect how many steps are taken, just monitoring.
mif2.particles = 3500

calc.r0=F


N_smooth_draws = 1000 #@note N_smooth_draws
N_particles_per_draw = 2500

hosp4csv_from_dashboard = T # take hosp for csv files from dashboard, not from admit file



# Fill in the data --------------------------------------------------------
# Read in covariates

source("code/common/sdmetrics-msas.R")

source("code/common/nextgen.R")
source("code/common/NGM_R_with_S.R")
source("code/common/data-processing.R")
source("code/models/mobility+agegroups+ar1+mu+pres+risk-simple.ml.R.cpp")
source("code/common/small_func.R")



days_average = 1
#days_average = 7

num_forecast_days = 60
normalization_date = "2020-02-21"
num_days_before_first_death = 22

mitigation_trans_start = ymd("2020-04-20")
mitigation_trans_end = ymd("2020-04-27")




# @note start args
# Manually give commend line arguments for debugging
if( length(args)==0) {
#args=c("continue","RDA/before.plot.20201220.Rda","after:after.smooth","dZfile:dZ.test1.csv","name.add:dZ1")
#args=c("continue","RDA/before.plot.20201219.Rda","after:after.smooth","name.add:dZ.flat","dZfile:dZ.test3.csv")
#args=c("continue","RDA/before.plot.20201219.Rda","after:after.smooth","name.add:dZ.bump","dZfile:dZ.test2.csv")
#args=c("continue","RDA/after.analysis.20210305.Rda","after:sim","name.add:mask2.0-15","dZfile:dZfile:dZ.mask.1.vac.0.4136.b117.-1.6.spr.1.csv","no:dash","no:pdf")
    #args=c("data/CV19AdmitDisch_UT_20210125.xlsx","data/Daily Hospitalization Dashboard Update 20210125.xlsx","Reg")
    args=c("history", "data/CV19AdmitDisch_UT_20210210.xlsx","30", "no:mob" ) 
 #  args=c("continue","after.mf2.20210308.Rda")
}


analysis.steps = c( "prepare", "mf2", "smooth","sim","analysis","saving")
analysis.continue.after ="all" # In case we want to continue analysis from a later step

do.this.step = function( analysis.continue, current.step ) {
    if( is.null( analysis.continue ) || (class(analysis.continue) != "character") ) {
        return(TRUE)
    }
    cat( "analysis.continue=",analysis.continue," current.step=",current.step,"\n")
    if( grepl( pattern = "^after.", analysis.continue)) {
        analysis.continue = gsub( pat="^after.",rep="",analysis.continue)
        after = T
    } else {
       after=F
    }
    step.done.i = which( analysis.continue == analysis.steps )
    next.step.i = which( current.step == analysis.steps )
    
    cat( "i:",step.done.i,",",next.step.i,ifelse(after,",after\n","\n"))

#    browser()
    if( after) {
        ret = (length(next.step.i)==0) || (length(step.done.i)==0) || ( next.step.i > step.done.i )
    } else {
        ret =(length(next.step.i)==0) || (length(step.done.i)==0) || ( next.step.i >= step.done.i )
    }
    ret
}

n.max.days=0
short.run=F

# Somewhat convoluted command line argument parsing

if( length(args)>1) {
  if( args[1] == "help") { #@note help
      cat("
Run Austin SEIR projection
## Commannd line arguments

First argment can be a modifier, or missing.
Allowed modifiers are `help`, `history`, `continue`, `short`

### Regular run: 
`Rscript run-seir.R Admit_discharge_datafile.xlsx ICY_datafile.xlsx`

### Continuing previous run:
`Rscript run-seir.R continue Rdata_session_file.Rda`
stage of run to continue is determined by name of file before/after stage.

`Rscript run-seir.R  continue Rdata_session_file.Rda after:stage`
continue run after stage. Possible stages: prepare, mf2, smooth, sim, analysis.

### History runs for performance analysis
`Rscript run-seir.R history Admit_discharge_datafile.xlsx [ndays]`
Load Admit_discharge_datafile.xlsx data file, but cut to first [ndays].

### Short run
`Rscript run-seir.R short Admit_discharge_datafile.xlsx ICY_datafile.xlsx`
Don't start run from beginning of hospitalization data, but instead load initial state from a previous run for a later date.
That run has to be carefully prepared top have enough initial variation and started at a relatively quiet part of the epidemic.
The start session is currently hard coded.
Short runs were not used for the PNAS paper analysis.

### Additional arguments
Additional aguments come in the form of keyword:argument

`no:dash` don't produce dashboard csv files
`no:pdf` don't produce pdf report
`no:mob` run without mobility data

`name.add:string` add string to output file name
")
  } else if( args[1] == "history") { #@note history
    if( length(args) < 3) {
      stop("history run requires 2 args: table_file n.max.weeks/days")
    }
    table.file = args[2]
    icu.file = "data/Daily Hospitalization Dashboard Update 10.26.2020.xlsx"
    run.name = table.file %>% 
        gsub( pat= "^.*[^0-9][^0-9][^0-9]([0-9][0-9.]*)[.]xlsx",rep="\\1") 
    
    n.max.days = as.numeric(  args[3] ) 
    cat("history run of only",n.max.days,"days\n")

    # @todo added no:agrs to history. Test!
    if( length(args) >3) {
          argl = strsplit(args[-(1:3)],":")
          arg.list=lapply(argl,function(z) setNames( z[2],z[1] )) %>% unlist   # @note generate arg.list
          
      } else {
          arg.list=c()
      }
  } else if( args[1] == "continue") { #@note continue
      if( length(args) < 2) {
          stop("continue requires a file to load as a second argument")
      }
      session.file = args[2]
      for( s in analysis.steps) {
        if( grepl(pattern = s, session.file)) {
            
            if( grepl(pat=paste0("after.",s),session.file)) {
                analysis.continue.after = paste0("after.", s)
            } else {
                analysis.continue.after = s
            }
            run.name = gsub(session.file,pattern = paste0("^(RDA/)*",s,"[.](.*)[.]Rda$"), repl="\\2")
            break
        }          
      }

      if( length(args) >2) {
          argl = strsplit(args[-(1:2)],":")
          arg.list=lapply(argl,function(z) setNames( z[2],z[1] )) %>% unlist   # @note generate arg.list
          if( "after" %in% names(arg.list)) {
              analysis.continue.after = arg.list["after"]
          }
          
      } else {
          arg.list=c()
      }

      
      cat("Continuing session \'",session.file,"\' after step \'", analysis.continue.after,"\' run.name=\'",run.name,"\'\n",sep="")
  }  else {# not history or continue
    if( args[1]=="short") { #@note short run
        short.run=T
        args = args[-1]
    }
    if( length(args) == 2) {
        print(2)
      table.file = args[1]
      run.name = table.file %>% 
        gsub( pat= "^.*[^0-9][^0-9][^0-9]([0-9][0-9.]*)[.]xlsx",rep="\\1") 
      icu.file = args[2]
    } else if( length(args) == 3) {
        print(3)
      table.file = args[1]
      run.name = table.file %>% 
        gsub( pat= "^.*[^0-9][^0-9][^0-9]([0-9][0-9.]*)[.]xlsx",rep="\\1") 
      run.name = paste0( run.name,args[3])
      icu.file = args[2]
    }  

  }
} else{ # no arguments
      run.name="test.run"
      table.file = "data/UT Template for COVID Hospital Admissions + Discharges 12.06.2020.xlsx"
    #  icu.file = "data/ICU Admit Discharge Meyers 10.21.2020.xlsx"
      icu.file = "data/Daily Hospitalization Dashboard Update 12.06.2020.xlsx"
}


# @todo init.states settings
if( short.run ) {
    init.states=list( 
        file="before.plot.11.23.2020B.Rds", 
        date0=ymd("2020-10-01") 
    )
}


if( ( gsub(pat="after.",rep="",analysis.continue.after) %in% analysis.steps)  )  {
    cat("run name:",run.name,"\n")
    save.vars=list(arg.list=arg.list, analysis.continue.after=analysis.continue.after, analysis.steps=analysis.steps,
                    do.this.step=do.this.step )
    load( session.file )
    #browser() #@audit-ok browser()
    for( v in names( save.vars )) {
        assign(v,save.vars[[v]])
    }
    rm(save.vars)
    if( "name.add" %in% names(arg.list) ) {
        run.name = paste0(run.name,".",arg.list["name.add"])
    }
    cat("After load step is ",analysis.continue.after)
    cat("run name after:",run.name,"\n")
} else {

    if( grepl("ICU Admit Discharge",icu.file)) {
        icu.data.file.type = "line data"
    } else if( grepl("Dashboard",icu.file) | grepl("KeyIndicators",icu.file)) {
        icu.data.file.type = "dashboard"
    } else {
        cat("Could not figure out if ICU file is dashboard or line data. Assuming line data\n")
        icu.data.file.type = "line data"
    }


    cat( "run name=",run.name)
    cat( "file=",table.file)


    cat("**** step 1: mif2\n")
    ##################################################################################
    # Read in new hospital data

    hospitalization_data = get_hospitalization_data( table.file )

    first_death = min(hospitalization_data$date[hospitalization_data$adm_total > 0])



    hospitalization_data %<>% 
            filter(date >= first_death - num_days_before_first_death) %>% 
            fill_dates(seq.Date(first_death - num_days_before_first_death, first_death - 1, by=1))

    hospitalization_data[ hospitalization_data$date < first_death,c("adm_total","recov_total","deaths_total")]=0
    
    hospitalization_data %<>% 
            mutate(day=0:(n() - 1),
            cum_adm_total = cumsum(adm_total),
            cum_recov_total = cumsum(recov_total),
            cum_deaths_total = cumsum(deaths_total),
            hospitalized = cum_adm_total - cum_recov_total - cum_deaths_total)


    cur.H.date = max( hospitalization_data$date, na.rm=T)
    min.H.date = min( hospitalization_data$date, na.rm=T)
    #
    ########################################################################
    # Read Mobility data

    top_categories <- c(
        "Colleges, Universities, and Professional Schools"="colleges",
        "Drinking Places (Alcoholic Beverages)"="drinking",
        "Elementary and Secondary Schools"="schools",
#       "General Medical and Surgical Hospitals"="medical",
        "Grocery Stores"="grocery",
        "Museums, Historical Sites, and Similar Institutions"="museums/parks",
        "Restaurants and Other Eating Places"="restaurants")

    sd_cols = c("median_home_dwell_time_relative")
    
    cat("reading mobility\n")
#    days_average=7
# BOOKMARK mobility_data

    source("code/common/sdmetrics-msas.R")
    mobility_data = get_mobility_data(days_average = days_average, 
                                    top_categories=top_categories,
                                    sd_cols = sd_cols, rerun=c()) # includes deaths, cases, #rerun=T



    # If we do a cut run, cut mobility data
    if( n.max.days > 0) {
        last.day = min.H.date + n.max.days
        mobility_data$msa_data %<>% filter( date <= last.day)
    }

    ### 
    # Add pca to msa_data 
    mobility_data = within( mobility_data, {
        # add pca
        cur_cols = c(paste0("count_relative_", unname(top_categories)), sd_cols)
        new_cols = paste0(cur_cols, "_av")
        pca = prcomp(data.matrix(na.omit(msa_data[ ,new_cols])), center = TRUE, scale=TRUE, rank=3)
        # comps = data.matrix(msa_data[ ,row.names(pca$rotation)]) %*% pca$rotation
        comps = predict(pca, newdata = msa_data)
        msa_data[ , paste0("PC", 1:3, "_av")] = comps
    } )


    #
    ########################################################################
    # Prepare covars
    mitigation_trans_dur = as.integer(mitigation_trans_end - mitigation_trans_start)

    covars = mobility_data$msa_data %>% 
        filter(msa == "Austin-Round Rock-Georgetown, TX") %>% 
#    filter( date < ymd("2020-12-18")) %>%
#       select(date, PC1=count_relative_drinking_av, PC2=median_home_dwell_time_relative_av, PC3=count_relative_schools_av, weekend=weekend, cases, cumulative_cases) %>%
        select(date, PC1=PC1_av, PC2=PC2_av, PC3=PC3_av, weekend=weekend, cases, cumulative_cases) %>%
        fill_dates(dates_from=hospitalization_data$date) %>% 
        mutate(weekend=as.numeric(chron::is.weekend(date))) %>% 
        right_join(select(hospitalization_data, date, day), by="date") %>% 
        mutate(maxT = n()) %>% 
        mutate(day = 0:(n() - 1)) %>% 
        mutate( dZ=0) %>%
        mutate( debug = 0) %>%
        mutate(awareness = case_when(date <= mitigation_trans_start ~ 0,
                                    date >  mitigation_trans_end ~ 1,
                                    TRUE ~ (1 / mitigation_trans_dur) * as.integer(date - mitigation_trans_start))) %>%
        mutate( future = 0 )


    x=covars[ covars$date=="2020-02-29",]; x$date = x$date+1
    covars[covars$date=="2020-03-01",] = x
    covars=covars[order(covars$date),]








    # if we do a cut run, cut hosp and covars
    if( n.max.days>0 ) {
        cat("reducing days\n")
        covars %<>% filter( date <= last.day ) 
        covars_last_date = max( covars$date)
        hospitalization_data %<>% filter( date <= covars_last_date )
        cur.H.date = max( hospitalization_data$date, na.rm=T)
    } else {
        covars_last_date = max( covars$date) # was max(mobility_data$msa_data$date)
    }


    # if (days_average == 7) { # if we already averaged, we don't need to average, just take last row
    # future_covars_row = tail(covars, 1)  
    # } else {
    # future_covars_row = tail(covars %>% filter(date <= covars_last_date), 7)  %>% 
    #     summarise_all(mean)
    # }

    #mobility_data$msa_data %>% 
    #  filter(msa == "Austin-Round Rock-Georgetown, TX")  %>% 
    #  select(date, count_relative_medical)

    covars_baseline = covars %>% 
        filter(date == normalization_date)

    # set PC's to be 0 exactly on the baseline day
    for (var in c("PC1", "PC2", "PC3")) {
        covars[ ,var] = covars[ , var] - as.numeric(covars_baseline[ , var])
    }



    # @todo changed here to remove mobility
    if( ( "mob" %in%  arg.list[ names(arg.list)=="no"])) {
        cat( "setting mobility PC to 1\n" )
        covars %<>% mutate( .keep="all",
                            PC1=1, PC2=1, PC3 = 1)      
        run.name = paste0( run.name,"-no:mob")  
    }


    # plot(covars$awareness)
    maxT = nrow(covars)

    # # Make sure covariates are fine
    # covars %>% 
    #   select(date, PC1, PC2, PC3) %>% 
    #   gather(variable, value, -date) %>% 
    #   ggplot(aes(x=date, y=value, color=variable)) +
    #   geom_line() +
    #   geom_hline(yintercept=0.0, lty=2)

    # grid.arrange(g1, g2, nrow=1)

    #############################################################################################
    # Create the pomp model ---------------------------------------------------

    if( F ) {
        load("last.coef.for.init.Rdata")
    }



    # Simulate from the prior pomp model and plot ------------------------------------
    # covid_mod %>%
    #   pomp::simulate(nsim = 1, format="data.frame") %>%
    #   View

    # a= covid_mod %>% 
    #     pomp::simulate(nsim = 1, format="data.frame") 


    # covid_mod %>% 
    # pomp::simulate(nsim = 100, format="data.frame") %>% 
    # ggplot(aes(day, total_new_H, group = .id)) + 
    # geom_line(alpha = .4) +
    # labs(x = "Day")


    #save.session(file=paste0("before.mf2.",run.name,".Rda"))

    # Conduct mif2 fitting ----------------------------------------------------
    timing$before.mf2=Sys.time()

    guess.init = function( init_pars, n=10, log_parameters, logit_parameters, sd=1) {
        m = length( init_pars)
        res = matrix( init_pars, n,m, byrow = T)
        colnames(res) = names(init_pars)
        log.i = names( init_pars) %in% log_parameters
        logit.i = names( init_pars) %in% logit_parameters
        neither.i = !logit.i & !log.i
        res[,neither.i] = rnorm( n * sum( neither.i), rep( init_pars[ neither.i], each=n), sd=sd )
        res[,log.i    ] = rnorm( n * sum( log.i    ), rep( log( init_pars[ log.i    ]), each=n), sd=sd ) %>% exp
        res[,logit.i  ] = rnorm( n * sum( logit.i  ), rep( logit( init_pars[ logit.i  ]), each=n), sd=sd ) %>% expit
        res
    }

# BOOKMARK optim.max()
###
# Find mf2 with highest likelihood by sampling till difference is bigger than 3 times standard error or the mean
optim.max=function( obs, N=10) {
  ll=c()
  while( length(obs) > 1) {
    ii = rep(seq_along(obs), N )
    l = foreach( i=ii , .packages = c("magrittr", "pomp") ) %dopar% {
      c( i, obs[[i]] %>% pfilter() %>% logLik )
    } %>% Reduce( rbind,.) %>% {tapply(.[,2],.[,1],c)} %>% Reduce(rbind,.) 

    ll=cbind(ll,l)
#    browser()
    n = dim(ll)[2]
    m = apply( ll, 1, mean)
    s = apply( ll, 1, sd)/ n^0.5
    o = order(m, decreasing = T)
    obs=obs[o]
    ll  = ll[o,]
    m  =  m[o]
    s  =  s[o]
    i = m > m[1]- 3*s[1]
    obs=obs[i]
    ll=ll[i,]
    cat("Finding best mf2. left:",dim(ll)[1]," ",head(m,2),"\n")
  }
  obs[[1]]
}

# @note find.mif function
find.mif = function( covid_mod, data, covars, N.steps=200, N.particles=2000, clust.n=14, dt = 1) {
  my.cl = makeCluster( clust.n, outfile="" )
  registerDoParallel( my.cl )
  
  mif2.N.steps= N.steps
  mif2.iter = min( 10, N.steps)
  # browser() ;

  mf2s = foreach( cpu=1:clust.n, 
                        .packages = c("magrittr", "pomp"), 
                        .export=c("covid_init_pars","covid_paramnames", "smoothed_states_t0"),
                        .verbose = T,
                        .errorhandling = "remove"
  ) %dopar% {
    mf2 <- covid_mod %>%
        mif2(
            Np= N.particles,
            Nmif= mif2.iter,
            cooling.fraction.50=0.6,
            paramnames= covid_paramnames,
            params = covid_init_pars,
            partrans=parameter_trans(log=c("sig_Z", "sig_R", "r", "sig_mu", "sig_R", "sig_delta_b1","sig_delta_b2"), #DLM
                                            logit=c("psi_R", "psi_mu", "psi")),
            rw.sd= rw.sd( log_beta_0 = 0.03, 
                        # b1=0.03, b12=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.02,  #DLM
                        r= 0.02, # if optimize r, uncomment the prior in dmeas
                        sig_Z=0.02, psi=0.0,
                        sig_delta_b1 = 0.01, sig_delta_b2 = 0.01, b1_0 = 0.03, b2_0 = 0.03, # DLM
                        sig_mu=0.02, psi_mu=0.02, log_mu_0=0.0,
                        sig_R=0.02, psi_R=0.0, log_gamma_H_0=0.0) #MICHAEL TRYING recovery that is more variable psi_R=0.02
        )


    for( i in (attributes(mf2)$Nmif/mif2.iter):(mif2.N.steps/mif2.iter) ) {
        print( c(cpu=cpu, mif2.iter=i, time=system.time( { mf2 %>% continue(Nmif= mif2.iter) ->mf2 } )[1] ) )
    }
    mf2
  }
#   ll = foreach( i=1:length(mf2s) , 
#                 .packages = c("magrittr", "pomp")
#                 ) %dopar% { 
#                   replicate( 100, mf2s[[i]] %>% pfilter() %>% logLik ) %>% mean
#                 }
  #  browser()
  m = optim.max( mf2s, N=20)
  stopCluster(my.cl)
  m
}



    mf2.par = function( clust.n, pop.n) {
        res.run=list()
        
        repop = function( mf) {
            m = length(mf)
            i = c( 1:m, sample( 1:m, pop.n-m, rep=T))
            mf[ i]
        }
        init.m = guess.init( covid_init_pars, n=pop.n, log_parameters, logit_parameters, sd=0.2 )
        
        my.cl = makeCluster(clust.n)
        registerDoParallel( my.cl )
        
        foreach(i=1:clust.n,
                .inorder=FALSE, .packages=c("pomp","magrittr"), 
                .export = c("covid_mod","log_parameters","logit_parameters","covid_paramnames"),
                .errorhandling = "remove"
        ) %dopar% {
            covid_mod %>%
            mif2(
                Np= mif2.particles,
                Nmif= 5,
                paramnames= covid_paramnames,
                params = init.m[i,],
                partrans=parameter_trans(log=log_parameters, logit=logit_parameters),
                #rw.sd=rw.sd(a = 0.03, b1=0.03, b12=0.02, b2=0.03, b3=0.03, b4=0.02, b5=0.02,
                rw.sd=rw.sd( log_beta_0 = 0.03, 
                            # b1=0.03, b12=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.02,  #DLM
                            r= 0.02, # if optimize r, uncomment the prior in dmeas
                            sig_Z=0.02, psi=0.0,
                            sig_delta_b1 = 0.01, sig_delta_b2 = 0.01, b1_0 = 0.03, b2_0 = 0.03, # DLM
                            sig_mu=0.02, psi_mu=0.02, log_mu_0=0.0,
                            sig_R=0.02, psi_R=0.0, log_gamma_H_0=0.0), #MICHAEL TRYING recovery that is more variable psi_R=0.02
                cooling.fraction.50=0.95
            ) -> m1
            ll <- replicate(10, logLik(pfilter(m1,Np=10000))) %>% mean(na.rm=T)
            list(mif=m1,ll=logmeanexp(ll,se=TRUE))
        } -> mf
        
        i= sapply(mf,function(i) {i$ll[1]}) %>% is.na
        mf1=mf[!i]
        
        res.run[[1]] = lapply( mf, function(i) i$ll)
        
        for( cool in rep(0.5,30)) {
            foreach(i=1:length(mf1),
                    .inorder=FALSE, .packages=c("pomp","magrittr"), .errorhandling = "remove"
            ) %dopar% {
            mf1[[i]]$mif %>% 
                continue(Nmif=10,cooling.fraction=cool) -> m1
            ll <- logLik(pfilter(m1,Np=3500))
            list(mif=m1,ll=logmeanexp(ll,se=TRUE))
            } -> mf1
            
            l = sapply(mf1,function(i) {i$ll[1]})
            res.run[[length(res.run)+1]]=l
            mf1=mf1[!is.na(l)]
            #    l=l[!is.na(l)]
            #    mf1=mf1[l>median(l)]
            mf1=repop(mf1)
        }
        stopCluster(my.cl)
        list(mf=mf1,res.run=res.run)
    }

    # @todo init.states: pomp
    if( exists("init.states") ) {
        init_res = readRDS( init.states$file )
        if( is.null( init.states$t0)) {
            init.states$t0 = which( init_res$new_covars$date == init.states$date0 )
        }
        smoothed_states_t0 = init_res$smoothed_states[ covid_statenames, init.states$t0, ]
        rm( "init_res")
        #browser()
        after.t0 = -(1:(init.states$t0-1))

        hospitalization_data = hospitalization_data[ after.t0,]
        covars = covars[ after.t0,]

        rinit_smooth <- function( smoothed_states_t0, ...  ) {
            n=dim(smoothed_states_t0)[2]
            i =sample.int(n,size=1)
            smoothed_states_t0[,i]
        }


        covid_mod = pomp(
            data=  hospitalization_data,
            times="day", 
            t0 = init.states$t0-1,
            rinit=rinit_smooth, smoothed_states_t0 = smoothed_states_t0,
            rmeasure=covid_rmeas,
            dmeasure=covid_dmeas_admitdischarge_nb, 
            covar=covariate_table(select(covars, -date), times="day"),
            rprocess = pomp::euler(covid_rprocess,1),
            statenames= covid_statenames,
            paramnames= covid_paramnames,
            params = covid_init_pars,
            accumvars = covid_accum_vars
        )
    } else {
        smoothed_states_t0=c()
        covid_mod = pomp(
            data=hospitalization_data,
            times="day", 
            t0 = 0,
            rinit=rinit, 
            rmeasure=covid_rmeas,
            dmeasure=covid_dmeas_admitdischarge_nb, 
            covar=covariate_table(select(covars, -date), times="day"),
            rprocess = pomp::euler(covid_rprocess,1),
            statenames= covid_statenames,
            paramnames= covid_paramnames,
            params = covid_init_pars,
            accumvars = covid_accum_vars
       )
    }

    mif.type="find.mif2"

    if( mif.type=="here" ) {
        mf2 = covid_mod %>%
            mif2(
                Np=3500,
                Nmif= mif2.iter,
                cooling.fraction.50=0.5,
                paramnames= covid_paramnames,
                params = covid_init_pars,
                partrans=parameter_trans(log=c("sig_Z", "sig_R", "r", "sig_mu", "sig_R", "sig_delta_b1","sig_delta_b2"), #DLM
                                        logit=c("psi_R", "psi_mu", "psi")),
                rw.sd= rw.sd( log_beta_0 = 0.03, 
                            # b1=0.03, b12=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.02,  #DLM
                            r= 0.02, # if optimize r, uncomment the prior in dmeas
                            sig_Z=0.02, psi=0.0,
                            sig_delta_b1 = 0.01, sig_delta_b2 = 0.01, b1_0 = 0.03, b2_0 = 0.03, # DLM
                            sig_mu=0.02, psi_mu=0.02, log_mu_0=0.0,
                            sig_R=0.02, psi_R=0.0, log_gamma_H_0=0.0) #MICHAEL TRYING recovery that is more variable psi_R=0.02
            )
        
        
        for( i in (attributes(mf2)$Nmif/mif2.iter):(mif2.N.steps/mif2.iter) ) {
            print( c(mif2.iter=i, time=system.time( { mf2 %>% continue(Nmif= mif2.iter) ->mf2 } )[1] ) )
        }
        
        #  for( i in 1:50 ) {
        #    print( c(mif2.iter=i, time=system.time( { mf2 %>% continue(Nmif= mif2.iter, cooling.fraction.50=0.5) ->mf2 } )[1] ) )
        #  }
    
    } else if( mif.type=="mf2.par") {
        
        
        mf = mf2.par( clust.n, clust.n )
        
        l=sapply(mf$mf,function(i) {i$ll[1]})
        mf2 = mf$mf[[which.max(l)]]$mif
    } else if( mif.type=="find.mif2") {
      mf2 = find.mif( covid_mod=covid_mod, data=hospitalization_data, covars=covars,  #@note call to find.mif2
                      N.steps = mif2.N.steps, N.particles = mif2.particles, clust.n = clust.n, dt=global.dt)
    }
    timing$after.mf2=Sys.time()



    #for( i in 2:10 ) {
    #  print( c(mif2.iter=i, time=system.time( { mf2 %>% continue(Nmif= mif2.iter) ->mf2 } )[1] ) )
    #}


    # # simulate from fitted model
    # mf2 %>% 
    # pomp::simulate(nsim = 100, format="data.frame") %>% 
    # ggplot(aes(day, total_new_H, group = .id)) + 
    # geom_line(alpha = .4) +
    # labs(x = "Day")

    # diagnostics and fit ------------------------------------------------------

    # likelihood
    replicate(10, mf2 %>% pfilter() %>% logLik()) %>% 
    logmeanexp(se=TRUE)

    # initial vs fitted parms
    cbind(covid_init_pars, coef(mf2)) %>% round(4)

    #plot(mf2)

#    cat("Now saving in parallel...\n")
#    mcparallel( detached=T, silent=T, 
#                expr = {
                    save.session(file=paste0("after.mf2.",run.name,".Rda"))
#                } )
    cat("done saving\n")
} # if not analysis.continue.after %in% analysis.steps


cat("Check if \'",analysis.continue.after,"\' is in ", paste( tail(analysis.steps,-1), collapse=", " ),"\n",sep="")

if( do.this.step( analysis.continue.after, current.step="smooth") ) { # If it is in it means we need to jump over this step
    cat("Not in\n")
    cat("*** step 2: smooth\n")
    # browser() # @audit-ok browser() smooth

    # @todo This is ugly. Better would be if covars and hospitalization data didn't change
    if( exists("init.states") ) {
        init_res = readRDS( init.states$file )
        if( is.null( init.states$t0)) {
            init.states$t0 = which( init_res$new_covars$date == init.states$date0 )
        }
        rm( "init_res")
        #browser()

        hospitalization_data %<>% filter( date >= init.states$date0 )
        covars %<>% filter( date >= init.states$date0 )
    }

    ## Estimate a smoothed distribution over beta --------------------------------

#########################
    #browser() # @audit-ok browser()

    #@note ****************** smooth
#    my.cl = makeForkCluster( clust.n.fork, outfile="")
    my.cl = makeCluster( clust.n.fork, outfile="")
    registerDoParallel( my.cl )

    # Filtered states  p(Z_t | Y_1:t)
    # Smoothed states  p(Z_t | Y_1:T)  t <= T
    # Forecasted states p(Z_t | Y_1:T) = int p(Z_t | Z_T)p(Z_T | Y_{1:T}) dZ_T t > T

    # construct the raw filtered states for the smoothed distribution
    # From each pfilter, one traj is generated
    # BOOKMARK smooth_dist
    smooth_dist = foreach(i=1:N_smooth_draws, 
                        .packages = c( "pomp")
    ) %dopar% {
        mf_filter = pfilter( mf2, Np = N_particles_per_draw, save.states=TRUE, filter.mean=TRUE, filter.traj=TRUE)
        #st=drop(filter.traj(mf_filter))
        #cur.i = length(covars$date)
        
        #     list(ll = logLik(mf_filter)  + 0*dpois( st["total_H",cur.i], hospitalization_data$hospitalized[cur.i] ,log=T ), states=drop(filter.traj(mf_filter)))
    #    print( c(cpu=cpu, smooth.iter=i) )
        list(ll = logLik(mf_filter) , states=drop(filter.traj(mf_filter)))
    }

#########################



    stopCluster(my.cl)
    
    n_states = nrow(smooth_dist[[1]]$states)

    # extract the weights for each filtered trajectory
    smoothing_weights = smooth_dist %>%
        lapply(function(x) x$ll) %>%
        unlist %>%
        (function(x) {exp(x - max(x))})

    # bind up the filtered trajectories for each state in an array

    smooth_dim = dim( smooth_dist[[1]]$state)
    filtered_states_unsampled = smooth_dist %>%
        lapply(function(x) x$states) %>%
        unlist %>%
        array(dim=c( smooth_dim, length(smooth_dist)))

        # resample indices and get smoothed states
  #  sample_ind = pomp::systematic_resample(smoothing_weights)
  #  filtered_states = filtered_states_unsampled[,,sample_ind]

    sample_ind = pomp::systematic_resample(smoothing_weights)
    filtered_states = filtered_states_unsampled[,,sample_ind]


}

#browser() #@audit browser() sim
cat(analysis.continue.after,"\n")
if( do.this.step( analysis.continue.after, current.step="sim") ) { # If it is in it means we need to jump over this step
    cat("Not in\n")
    cat("*** step 3: future sim\n")

    num_forecast_days = 60  # @audit added num_forecast_days=120 because it isn't loaded


#    browser() #@audit-ok browser()


    # @todo This is ugly. Better would be if covars and hospitalization data didn't change
    if( exists("init.states") ) {
        init_res = readRDS( init.states$file )
        if( is.null( init.states$t0)) {
            init.states$t0 = which( init_res$new_covars$date == init.states$date0 )
        }
        rm( "init_res")
        #browser()
        hospitalization_data %<>% filter( date >= init.states$date0 )
        covars %<>% filter( date >= init.states$date0 )
    }



#    sample_ind = pomp::systematic_resample(smoothing_weights)
#    filtered_states = filtered_states_unsampled[,,sample_ind]

    # ------------------------------------------------------------------------
    # Augment with forecasts counterfactual using sim

    max_dat_date = max(hospitalization_data$date)
    min_dat_date = min(hospitalization_data$date)

    # Modify for forecasting
    new_dates = seq(max_dat_date + 1, max_dat_date + num_forecast_days, by=1)


    # makes a counterfactual covariate matrix for the future

    type = "stay"
    covars_last_date = max( covars$date)
    
    if (days_average == 7) { # if we already averaged, we don't need to average, just take last row
        future_covars_row = tail(covars, 1)  
    } else {
        # argument "date:" uses past covars
        if( "date" %in% names(arg.list) ) {
            future_covars_row = tail(covars %>% filter(date <= ymd(arg.list["date"])), 7)  %>% 
                summarise_all(mean)
            # future_beta_states_i is the time from which to take the states
            # that affect beta.
            # We are using attributes(mf2)["covar"][[1]]@times here
            future_beta_states_i = which( covars$day[ covars$date ==ymd(arg.list["date"]) ] == attributes(mf2)["covar"][[1]]@times  ) 
        } else {
            future_covars_row = tail(covars %>% filter(date <= covars_last_date), 7)  %>% 
                summarise_all(mean)
            future_beta_states_i = which( covars$day[ covars$date == covars_last_date] ==  attributes(mf2)["covar"][[1]]@times  )
        }
    }


    # @todo new_covars
    new_covars = covars %>% 
        bind_rows(purrr::map(new_dates, ~ {
            D = future_covars_row
            D$date = .x
            D}) %>% 
            bind_rows()
        ) %>%
        mutate( day = (0:(n() - 1))+min(covars$day)   ) %>%
        mutate( dZ = 0) %>%    # should be set in covars, but just in case it wasn't
        mutate( weekend = as.numeric(chron::is.weekend(date))) %>%
        mutate( future = ifelse( date > max_dat_date, 1, 0))

    

    # @note if argument Z is provided, state Z is increased by that amount.
    if( "dZfile" %in% names(arg.list) ) {
        dZ = read.csv(arg.list["dZfile"],as.is=T)
        if( "date" %in% colnames(dZ)) {
            dZ$date %<>% ymd
            d = new_covars$date[ new_covars$date %in% dZ$date ]
            new_covars[ new_covars$date %in% dZ$date ,"dZ"] = log(dZ[ dZ$date %in% d ,"dZ"] )
        } else {
            new_covars[ dim(covars)[1]+(1:dim(dZ)[1]) ,"dZ"] = log(dZ[,"dZ"] )
        }
    } 

    future_covars = tail(new_covars, num_forecast_days)

    #browser() # @audit-ok browser()
    #@todo future.pomp
    future.pomp = pomp(
            data=tail( hospitalization_data,1 ),
            times="day", 
            t0 = maxT -1,
            rinit=rinit,
            rmeasure=covid_rmeas,
            dmeasure=covid_dmeas_admitdischarge_nb, 
            covar=covariate_table(select(new_covars, -date), times="day"),
            rprocess = pomp::euler(covid_rprocess,1),
            statenames= covid_statenames,
            paramnames= covid_paramnames,
            params = coef(mf2),
            accumvars = covid_accum_vars
        )

    # build an rinit function given parameter i
    # This function returns a function
    future_rinit = function(i) {
        rinit_f = function(...) {
                        s1 = setNames( filtered_states[ , dim(smooth_dist[[1]]$states)[2], i],
                                       dimnames(smooth_dist[[1]]$states)$variable
                                )
                        s2 = setNames( filtered_states[ , future_beta_states_i,            i],
                                       dimnames(smooth_dist[[1]]$states)$variable
                                )
                        s1[ beta_state_vars] = s2[ beta_state_vars]
                        s1
                    }
        rinit_f
    }

    sim_future_states =  sapply( 1:dim(filtered_states)[3], simplify="array",function(i) {
        future.pomp %>% 
        pomp:::simulate(nsim=1, 
                    format="arrays",
                    rinit=future_rinit(i),
                    statenames=covid_statenames,
                    times=future_covars$day ) %>% 
        .[["states"]] %>% 
        .[,1,]
#        aperm(c(1, 3, 2)) 
    } )

    filtered_states_new = abind::abind(filtered_states, sim_future_states, along=2)
   # filtered_states_new = abind::abind(filtered_states, sim_future_states, along=2)
    my_dim_names = dimnames(smooth_dist[[1]]$states)
    my_dim_names$sample = 1:length(smooth_dist)

    my_dim_names$time = c(my_dim_names$time, maxT:(maxT + num_forecast_days - 1))
    dimnames(filtered_states_new) = my_dim_names

    ## The following was dropped, now we do resampling after smoothing
    # # resample indices and get smoothed states
    # sample_ind = pomp::systematic_resample(smoothing_weights)
    #browser()
    smoothed_states = filtered_states_new


#    browser()
    # myfun = function(beta) get_r0(beta, phi, pop_totals, pop_ratio_mat, eta, tau, gamma_A, gamma_Y, omega_A, omega_Y, omega_P_overall)

    # now we have the smoothed posterior over all states
    ###########

 #   cat("Now saving in parallel...\n")
#    mcparallel( detached=T, silent=T, 
#                expr = {
    if( !("save.sim" %in%  arg.list[ names(arg.list)=="no"] ) ) {
                    save.session(file=paste0("after.sim.",run.name,".Rda"))
#                } )
        cat("done saving\n")
    }
    
#restore.session(file=paste0("after.smooth.",run.name,".Rda"))
} # if( !(analysis.continue.after %in% tail(analysis.steps,-1) ) ) { 


if( do.this.step( analysis.continue.after, current.step="analysis") )  { 
    cat("*** step 4: analysis\n")
    if( exists("init.states") ) {
#    browser()
        old = readRDS( init.states$file )
        i = old$new_covars$date < min( covars$date)
        old$new_covars$dZ=0

        nonexistant.cols.old = colnames(covars)[!(colnames( covars) %in% colnames( old$new_covars))]
        if( length(nonexistant.cols.old) > 0 ) {
            old$new_covars[, nonexistant.cols.old] = covars[ 1, nonexistant.cols.old]
        }
        covars     = rbind( old$new_covars[i,],     covars) # @note covars extended here
        i = old$hospitalization_data$date < min( hospitalization_data$date)
        hospitalization_data = rbind( old$hospitalization_data[i,], hospitalization_data)
        i = old$new_covars$date < min( new_covars$date)
        new_covars = rbind( old$new_covars[i,], new_covars)
        j = min( dim(old$smoothed_states)[3], dim(smoothed_states)[3] )
        j = 1:j
        s = dimnames(smoothed_states)[[1]]
        smoothed_states = abind( old$smoothed_states[s,i,j], smoothed_states[s,,j], along=2)
        rm(list=c("old","s"))
    }


    # helper function for extracting a confidence interval
    # states is an n_states X n_samples matrix
    getCI =   function(states, alpha=alpha, name="states") {
        states_mean = apply(states, 1, mean,na.rm=T)
        states_sd = apply(states, 1, sd,na.rm=T)
        states_lo = apply(states, 1, function(x){quantile(x,alpha/2,na.rm=T)})
        states_hi = apply(states, 1, function(x){quantile(x,1-alpha/2,na.rm=T)})
        states_med = apply(states, 1, function(x){quantile(x,0.5,na.rm=T)})
        names=paste( name,c("mean","sd","lo","hi","med"),sep="_")
        res = data.frame(states_mean, states_sd, states_lo, states_hi, states_med)
        colnames(res) = names
        res
    }


    rolling = function( a, dims=c(2), days=7 ) {
        if( is.null(dim(a))) {
            x = cumsum(a)
            c( rep(NA,days), tail(x,-days)-head(x,-days)) / days
        } else {
            apply(a, dims, function(x) {
            x = cumsum(x)
            c( rep(NA,days), tail(x,-days)-head(x,-days)) / days
            })
        }
    }


    alpha=0.05

    # First Beta is NA, fix that.
    smoothed_states['Beta', 1,] = smoothed_states['Beta', 2,]

    # compute r0, rt
    pop_totals = colSums(pop)
    pop_ratio_mat = get_pop_ratio(pop_totals)

    myfun = function(beta) get_r0(beta, (phi), pop_totals, pop_ratio_mat, sigma, tau, 
                                        gamma_A, gamma_Y, omega_A, omega_Y, omega_P_overall)
    myfun2 = function(beta,S_prop) get_rt(beta, (phi), pop, pop_ratio_mat, sigma, tau, rho_A, rho_Y, 
                                        gamma_A, gamma_Y, omega_A, omega_Y, omega_P, S_prop) 
    
    my.cl = makeCluster(clust.n) #######################################################################
    registerDoParallel( my.cl )


    if( calc.r0) {
        r0 = foreach(i=1:N_smooth_draws, .combine = cbind) %dopar% {
            r0_i = purrr::map_dbl(smoothed_states["Beta", , i], myfun)
        }
     
        dimnames(r0) = dimnames( smoothed_states)[2:3]
    }

        state.names = gsub("_1_1","",dimnames(smoothed_states)[[1]][grep("^[A-Z]+_1_1",dimnames(smoothed_states)[[1]])])
        all.states = dimnames(smoothed_states)[[1]][grep("^[A-Z]+_[0-9]+_[0-9]+",dimnames(smoothed_states)[[1]])]
        S.states = matrix( paste( "S",matrix(1:2,2,5), matrix(1:5,2,5,byrow = T), sep="_"), 2, 5)

        rt = foreach(i=1:dim(smoothed_states)[3], .combine = cbind, .packages = c("magrittr")) %dopar% {
            Age.Pop.N = sapply(1:5,function(g) smoothed_states[ all.states[grep(paste0("_",g,"$"),all.states)]  ,,i] %>% colSums ) 
            Age.S.N = sapply(1:5,function(g) smoothed_states[ S.states[,g]  ,,i] %>% colSums ) 
            Age.S.prop = Age.S.N / Age.Pop.N
            #if(sum(is.na(Age.S.prop))>0 | sum(is.infinite( Age.S.prop))>0) { na.i <<- c(na.i,i) }
            rt_i = sapply( seq_len(nrow(Age.S.N)), function(j) myfun2(beta=smoothed_states["Beta",j , i], S_prop=Age.S.prop[j,])   )
        }
        dimnames(rt) = dimnames( smoothed_states)[2:3]

    if( !calc.r0) r0=rt

    stopCluster(my.cl)  #######################################################################

    r0 = array( r0, dim = c(1,dim(smoothed_states)[2:3]))
    rt = array( rt, dim = c(1,dim(smoothed_states)[2:3]))
    dimnames(r0)=c(list("R0"), dimnames(smoothed_states)[2:3])
    dimnames(rt)=c(list("Rt"), dimnames(smoothed_states)[2:3])
    smoothed_states =abind( smoothed_states, r0, rt, along=1)


    betaCI = getCI(smoothed_states['Beta', ,], alpha=alpha)
    r0CI = getCI( smoothed_states['R0', ,], alpha=alpha)
    r0CI = cbind( r0CI, fracIncr=rowMeans(r0>1)*100)

    rtCI = getCI( smoothed_states['Rt', ,], alpha=alpha)
    rtCI = cbind( rtCI, fracIncr=rowMeans(rt>1)*100)

    rt7 = smoothed_states['Rt', ,] %>% rolling
    rt7CI = getCI( rt7, alpha=alpha)
    rt7CI = cbind( rt7CI, fracIncr=rowMeans(rt7>1)*100)


    newHCIlist = list()
    dyingHCIlist = list()
    recovHCIlist = list()
    for (i in 1:5) {
        newHCIlist[[i]]= getCI(smoothed_states[paste0("new_H_1_", i), ,] +
                                smoothed_states[paste0("new_H_2_", i), ,], alpha=alpha, name=paste0("new_H",i))
        dyingHCIlist[[i]]= getCI(smoothed_states[paste0("dying_H_1_", i), ,] +
                                    smoothed_states[paste0("dying_H_2_", i), ,], alpha=alpha , name=paste0("dying_H",i) )
        recovHCIlist[[i]]= getCI(smoothed_states[paste0("recovering_H_1_", i), ,] +
                                    smoothed_states[paste0("recovering_H_2_", i), ,], alpha=alpha, name=paste0("recovering_H",i))
    }
    newHCI = getCI(smoothed_states['total_new_H',,], alpha=alpha, name="new_H")
    newH7CI = getCI(smoothed_states['total_new_H',,] %>% rolling , alpha=alpha, name="new_H7")
    recoveringHCI = getCI(smoothed_states['total_recovering_H',,], alpha=alpha, name="recovering_H")
    dyingHCI = getCI(smoothed_states['total_dying_H',,], alpha=alpha, name="dying_H")
    SCI = getCI(smoothed_states['S_1_1',,] + smoothed_states['S_2_1',,] +
                smoothed_states['S_1_2',,] + smoothed_states['S_2_2',,] +
                smoothed_states['S_1_3',,] + smoothed_states['S_2_3',,] +
                smoothed_states['S_1_4',,] + smoothed_states['S_2_4',,] +
                smoothed_states['S_1_5',,] + smoothed_states['S_2_5',,]  , alpha=alpha)
    ECI = getCI(  smoothed_states['E_1_1',,] + smoothed_states['E_2_1',,] +
                    smoothed_states['E_1_2',,] + smoothed_states['E_2_2',,] +
                    smoothed_states['E_1_3',,] + smoothed_states['E_2_3',,] +
                    smoothed_states['E_1_4',,] + smoothed_states['E_2_4',,] +
                    smoothed_states['E_1_5',,] + smoothed_states['E_2_5',,]  , alpha=alpha)
    RCI = getCI(  smoothed_states['R_1_1',,] + smoothed_states['R_2_1',,] +
                    smoothed_states['R_1_2',,] + smoothed_states['R_2_2',,] +
                    smoothed_states['R_1_3',,] + smoothed_states['R_2_3',,] +
                    smoothed_states['R_1_4',,] + smoothed_states['R_2_4',,] +
                    smoothed_states['R_1_5',,] + smoothed_states['R_2_5',,]   , alpha=alpha)
    IYCI = getCI( smoothed_states['IY_2_1',,] + smoothed_states['IY_1_1',,] + 
                    smoothed_states['IY_2_2',,] + smoothed_states['IY_1_2',,] + 
                    smoothed_states['IY_2_3',,] + smoothed_states['IY_1_3',,] + 
                    smoothed_states['IY_2_4',,] + smoothed_states['IY_1_4',,] +
                    smoothed_states['IY_2_5',,] + smoothed_states['IY_1_5',,]  ,
                alpha=alpha)
    IACI = getCI(  smoothed_states['IA_2_1',,] + smoothed_states['IA_1_1',,] +
                    smoothed_states['IA_2_2',,] + smoothed_states['IA_1_2',,] + 
                    smoothed_states['IA_2_3',,] + smoothed_states['IA_1_3',,] + 
                    smoothed_states['IA_2_4',,] + smoothed_states['IA_1_4',,] +
                    smoothed_states['IA_2_5',,] + smoothed_states['IA_1_5',,]  ,
                alpha=alpha)
    PACI = getCI( smoothed_states['PA_2_1',,] + smoothed_states['PA_1_1',,] +
                    smoothed_states['PA_2_2',,] + smoothed_states['PA_1_2',,] + 
                    smoothed_states['PA_2_3',,] + smoothed_states['PA_1_3',,] + 
                    smoothed_states['PA_2_4',,] + smoothed_states['PA_1_4',,] +
                    smoothed_states['PA_2_5',,] + smoothed_states['PA_1_5',,]  ,
                alpha=alpha)
    PYCI = getCI( smoothed_states['PY_2_1',,] + smoothed_states['PY_1_1',,] +
                    smoothed_states['PY_2_2',,] + smoothed_states['PY_1_2',,] + 
                    smoothed_states['PY_2_3',,] + smoothed_states['PY_1_3',,] + 
                    smoothed_states['PY_2_4',,] + smoothed_states['PY_1_4',,] +
                    smoothed_states['PY_2_5',,] + smoothed_states['PY_1_5',,]  ,
                alpha=alpha)

    infectious = c("PA","PY","IA","IY")
    infectious = expand_state( infectious, c(2,5))

    smoothed_states_inf_total = apply( smoothed_states[infectious,,],2:3,sum)

    infCI = getCI( smoothed_states_inf_total, alpha=alpha)
    names(infCI) = gsub("states", "inf", names(infCI))


    alpha=0.05
    HCI   = getCI(smoothed_states['total_H',,], alpha=0.05)
    ZCI   = getCI(smoothed_states['Z',,], alpha=alpha)
    muCI  = getCI(1.0 / exp(coef(mf2)["log_mu_0"] + smoothed_states['Z_mu',,]), alpha=alpha)
    muRCI = getCI(1.0 / exp(coef(mf2)["log_gamma_H_0"] + smoothed_states['Z_R',,]), alpha=alpha)
    NIYCI = getCI(smoothed_states['NIY',,], alpha=alpha)
    NICI  = getCI(smoothed_states['NI',,], alpha=alpha)

    b1CI = getCI(smoothed_states['b1',,], alpha=alpha)
    b2CI = getCI(smoothed_states['b2',,], alpha=alpha)

    ## Percentage increase
    days4increase = 14 
    PercIncCI = smoothed_states[infectious,,] %>% {apply(.,c(2:3),sum)} %>% 
    { (tail(., -days4increase ) / head(.,-days4increase ) -1 )*100 } %>%
    {getCI( ., alpha=alpha) }


    names(PercIncCI) = gsub("states", "PercInc", names(PercIncCI))

    # gather for plotting
    names(betaCI) = gsub("states", "Beta", names(betaCI))
    names(r0CI) = gsub("states", "r0", names(r0CI))
    names(rtCI) = gsub("states", "rt", names(rtCI))
    names(rt7CI) = gsub("states", "rt", names(rtCI))
    names(NICI) = gsub("states", "NI", names(NICI))
    names(b1CI) = gsub("states", "b1", names(b1CI))
    names(b2CI) = gsub("states", "b2", names(b2CI))


    for (i in 1:5) {
        # names(r0CIlist[[i]]) = gsub("states", paste0("r0", i), names(r0CIlist[[i]]))
        names(newHCIlist[[i]]) = gsub("states", paste0("new_H", i), names(newHCIlist[[i]]))
        names(dyingHCIlist[[i]]) = gsub("states", paste0("dying_H", i), names(dyingHCIlist[[i]]))
        names(recovHCIlist[[i]]) = gsub("states", paste0("recovering_H", i), names(recovHCIlist[[i]]))
    }
    names(newHCI) = gsub("states", "newH", names(newHCI))
    names(newH7CI) = gsub("states", "newH7", names(newH7CI))
    names(recoveringHCI) = gsub("states", "recoveringH", names(recoveringHCI))
    names(dyingHCI) = gsub("states", "dyingH", names(dyingHCI))
    names(ZCI) = gsub("states", "Z", names(ZCI))
    names(HCI) = gsub("states", "H", names(HCI))
    names(muCI) = gsub("states", "mu", names(muCI))
    names(muRCI) = gsub("states", "muR", names(muRCI))
    names(NIYCI) = gsub("states", "NIY", names(NIYCI))

    names(SCI) = gsub("states", "S", names(SCI))
    names(ECI) = gsub("states", "E", names(ECI))
    names(RCI) = gsub("states", "R", names(RCI))
    names(PACI) = gsub("states", "PA", names(PACI))
    names(PYCI) = gsub("states", "PY", names(PYCI))
    names(IACI) = gsub("states", "IA", names(IACI))
    names(IYCI) = gsub("states", "IY", names(IYCI))






    ########################
    # Calculate doubling times

    slopes = function(A, days=14) {
        n = dim(A)[1]
        res = sapply(1:(n-days), function(i) {
            x=i:(i+days)
            
            apply(A[x,],2,function(y) lsfit(x,y)$coef["X"])
        })
        t(res)
    }



    a = log(
    smoothed_states["NI",,] %>% {apply(.,2,function(x) {
        res=c(0,diff(x))
        res[res<0]=0;res
    })}+1e-6
    )
    x = (log(2) / slopes(a))

    # Put unique/first of duplicate elements first
    i = duplicated( smoothed_states["Beta",dim(covars)[1],])
    smoothed_states = smoothed_states[,,c( which(!i), which(i))]


    if( icu.data.file.type == "dashboard") {
        icu.sheet = excel_sheets(icu.file) %>% grep(pattern = "ICU")
        ICU=read_excel( icu.file,sheet=icu.sheet,skip=1) # This was the old dashboard
        ICU$date=daS( 3.12, len=dim(ICU)[1])
        icu.df = data.frame( date=ICU$date, ICU[,c(2:4)])

        if( hosp4csv_from_dashboard) {
            hosp.sheet = excel_sheets(icu.file) %>% grep(pattern = "Hosp")
            hosp4csv = read_excel( icu.file,sheet=icu.sheet,skip=1)
            hosp4csv$date = daS( 3.12, len=dim(ICU)[1])
            hosp.col = "CV19 Inpatients"
            hosp7.col = "CV19 Inpatients 7 Day MA"
            hosp4csv.df = data.frame( date =  hosp4csv$date, hospitalized=pull(hosp4csv[,hosp.col]),hospitalized7=pull(hosp4csv,hosp7.col))
        }
        
    } else {
        ICU2 = read.icu.admit.table( icu.file2) 
        ICU2[,"CV19 ICU Inpatients"] = cumsum( ICU2$adm - ICU2$dis)
        icu.df=data.frame( date=ICU$date, ICU[,c(1,2,4)])
    }

    icu.last.i = which(new_covars$date == tail( icu.df$date,1))[1]
    icu.ratio = tail(icu.df$CV19.ICU.Inpatients,1) / HCI$H_med[ icu.last.i ] 

    total_icu = smoothed_states["total_H",,,drop=F] * icu.ratio
    dimnames(total_icu)[[1]]="total_ICU"

    smoothed_states =abind( smoothed_states, total_icu, along=1)
    
    ICUCI   = getCI(smoothed_states['total_ICU',,], alpha=0.05, name = "ICU")


    timing$end1=Sys.time()

    cat("Now saving in parallel...\n")
#    job.save = mcparallel( detached=T, silent=T, 
#                expr = {
                    save.session(file=paste0("RDA/after.analysis.",run.name,".Rda"))
#                    "save done!"
#                } )
    cat("done saving\n")
    #library(session); restore.session(file=paste0("RDA/before.plot.",run.name,".Rda"))
    timing$end2=Sys.time()
    
} # if( !(analysis.continue.after %in% tail(analysis.steps,-2) ) ) 
#browser()

if( do.this.step( analysis.continue.after, current.step="saving")) {
if( n.max.days == 0 ) { # not a history run

        cat("*** step 5: saving\n")
    
        if( ! ("dash" %in%  arg.list[ names(arg.list)=="no"] ) ) {
            cat("dash dir:", dash_dir,"\n")

            NId.2rate.CI = data.frame( date=tail(new_covars$date, -14), getCI(x, alpha=alpha , name = "NIx2") ) %>%
                mutate(observed = as.integer(date <= cur.H.date) )

            n.spagettis = min( 100, dim(smoothed_states)[3])
            write.csv( data.frame( NId.2rate.CI,NIx2=x[,1:n.spagettis]), file=paste0(dash_dir,"NIx2.csv" ),row.names=F )


            a = log(smoothed_states["total_new_H",,]+1e-6)
            x = log(2) / slopes(a)
            H.2rate.CI = data.frame( date=tail(new_covars$date, -14), getCI(x, alpha=alpha, name = "Hx2" ) ) %>%
                mutate(observed = as.integer(date <= cur.H.date) )

            write.csv( data.frame( H.2rate.CI ,Hx2=x[,1:n.spagettis]), file=paste0(dash_dir,"Hx2.csv") )




            ######




            # Calculate ratio of states. Can then be saved, to use as initial state ratio
            state_ratio = 
                smoothed_states[expand_state(c("E","PY","PA","IY","IA","H"),c(2,5)),,] %>% 
                apply(.,1:2,mean) %>% .[,40] 


            



        smoothed_states[expand_state("IY",c(2,5)),,1:n.spagettis] %>% apply(.,2:3,sum) %>% {data.frame(date=new_covars$date,IY=.)} ->a
        covars[,c("date","cases")] %>%
            right_join(a,by="date") %>%
            write.csv(file=paste0(dash_dir,"IY.csv"),row.names=F)

        hospitalization_data$hospitalized7 = hospitalization_data$hospitalized %>% rolling
        hospitalization_data$adm_total7 = hospitalization_data$adm_total %>% rolling
        


        a = data.frame(date=new_covars$date, HCI,H=smoothed_states["total_H",,]  ) 
        if( hosp4csv_from_dashboard) {
            i = hospitalization_data$date < min( hosp4csv.df$date)
            hosp4csv.df = rbind( hospitalization_data[i,c("date","hospitalized","hospitalized7")],
                                hosp4csv.df )
            hosp4csv.df$hospitalized7 = hosp4csv.df$hospitalized %>% rolling
            hosp4csv.df[,c("date","hospitalized","hospitalized7")] %>%
                right_join(a,by="date") %>% 
                write.csv(file=paste0(dash_dir,"H.csv"),row.names=F)
        } else {
            hospitalization_data[,c("date","hospitalized","hospitalized7")] %>%
                right_join(a,by="date") %>% 
                write.csv(file=paste0(dash_dir,"H.csv"),row.names=F)
        }

        a = data.frame(date=new_covars$date, ICUCI,H=smoothed_states["total_ICU",,]  )  %>% 
                filter( date > min(icu.df$date))

        icu.df %>%
            right_join( a,by="date") %>% 
            write.csv(file=paste0(dash_dir,"ICU.csv"),row.names=F)

        a = data.frame(date=new_covars$date, newH7CI,new_H=smoothed_states["total_new_H",,] %>% rolling  ) 
        hospitalization_data[,c("date","adm_total","adm_total7")] %>%
            right_join(a,by="date") %>% 
            write.csv(file=paste0(dash_dir,"new_H.csv"),row.names=F)


        write.csv(NICI, file=paste0(dash_dir,"new_E.csv"),row.names=F)


        smoothed_states[expand_state("dying_H",c(2,5)),,1:n.spagettis] %>% apply(.,2:3,sum) %>% data.frame(date=new_covars$date,new_D=.) ->a
        hospitalization_data[,c("date","deaths_total")] %>%
            right_join(a,by="date") %>% 
            write.csv(file=paste0(dash_dir,"D.csv"),row.names=F)







        r0CI %>% { write.csv( data.frame( date=new_covars$date, .),  file=paste0( dash_dir ,"R0.csv"),row.names=F) }
        rt7CI %>% { write.csv( data.frame( date=new_covars$date, .),  file=paste0( dash_dir ,"Rt.csv"),row.names=F) }
        PercIncCI %>% { write.csv( data.frame( date=tail(new_covars$date,-days4increase), .),  file=paste0( dash_dir ,"PercInc.csv"),row.names=F) }

        # make all the csv files world readable
        #system(paste0("chmod 777 ",dash_dir,"*") )
    } # if( no dash)
 

    if( !( "pdf" %in%  arg.list[ names(arg.list)=="no"])) {
        cat("output: ",paste0("Austin-SEIR-",run.name,".pdf"),"\n")
        cat("run.name: ",run.name,"\n")
       library(rmarkdown)

        # Note: to load data from scratch, delete "smoothed_states"
        rmarkdown::render("Austin-SEIR.Rmd",output_file=paste0("Austin-SEIR-",run.name,".pdf"))
    }
#    message = mccollect( job.save)
#    cat("Save job says ","\n")
#    print(message)
#    Sys.sleep(5)
} else { # n.max.days >0 # This is a history run 
#    smoothed_indexes_save = c("total_H","total_new_H","total_recovering_H","total_dying_H","Rt","Beta")
     smoothed_indexes_save = T
    smoothed_states = smoothed_states[smoothed_indexes_save,,]
    save(smoothed_states, hospitalization_data, new_covars,mf2, file=paste0("HIST/before.plot.",run.name,"-",n.max.days,".Rda"))
}
}
