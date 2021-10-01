library(pomp)


#// @todo 
#// 1. change contact matrix
#// 2. Change contact matrix using schools
#// 3. Add delta frequency to covars
#// 4. Modify Beta by delta frequency
#// 5. Modify IHR by delta
#// 6. Add tauV


# these aux functions to gengamma_e useful cpp code form r objects -----

# takes an r matrix and gengamma_es the corresponding manual cpp code
mat2cpp = function(x, name, type="double") {
  stopifnot(is.matrix(x))
  out = sprintf("%s %s[%s][%s];\n", type, name, nrow(x), ncol(x))
  for (j in 1:ncol(x))
    for (i in 1:nrow(x))
      out[length(out) + 1] = sprintf(
        "%s[%s][%s]=%s;\n", name, i - 1, j - 1, x[i, j])
  paste(out, collapse="")
}

# takes an r vector and gengamma_es the corresponding manual cpp code
vec2cpp = function(x, name, type="double") {
  stopifnot(is.vector(x))
  out = sprintf("%s %s[%s];\n", type, name, length(x))
  for (i in 1:length(x))
    out[length(out) + 1] = sprintf(
      "%s[%s]=%s;\n", name, i - 1, x[i])
  paste(out, collapse="")
}

# this function creates cpp indexed e.g. double S[ndims=d]; from S1,...Sd
enable_indexes = function(vars, dims, type="double") {
  out = c()
  for (v in vars) {
    if (length(dims) == 1) {
      x = paste0(v, 1:dims)
      out[length(out) + 1] = vec2cpp(x, v, type)
    } else {
      x = matrix("", dims[1], dims[2])
      for (i in 1:dims[1])
        for (j in 1:dims[2])
          x[i, j] = paste0(v, paste0("_", i, "_", j))
        out[length(out) + 1] = mat2cpp(x, v, type)
    }
  }
  paste(out, collapse="")
}


# this functions upsates the true globals Sd = S[d - 1]
register_indexed_changes = function(vars, dims) {
  out = c()
  for (v in vars)
    if (length(dims) == 1) {
      for (d in 1:dims)
        out[length(out) + 1] = sprintf(
          "%s_%s=%s[%s];\n", v, d, v, d - 1)
    } else {
      for (i in 1:dims[1])
        for (j in 1:dims[2])
          out[length(out) + 1] = sprintf(
            "%s_%s_%s=%s[%s][%s];\n", v, i, j, v, i - 1, j - 1)
    }
  paste(out, collapse="")
}



# define the global data --------


# this is the overall
# symp_hosp_ratio = c(0.00048721, 0.00048721, 0.03287572, 0.11337395, 0.17733063)
# YHR_overall = c(0.0007018, 0.0007018, 0.04735, 0.1633, 0.2544)  # symp_hosp_ratio =
YHR = matrix(c(
    0.004021 ,  0.003091,   0.1903,   0.4114,   0.4879,
    0.0004021,  0.0003091,  0.01903,  0.04114,  0.04879),
  2, 5, byrow = TRUE)


# prob_group = matrix(c(high_risk_ratio,
#                       1 - high_risk_ratio), 2, 5, byrow = TRUE)
pop = matrix(c(9350,37451,156209,108196,103763,
               128527,327148,915894,249273,132505),
             2, 5, byrow = TRUE)
totalpop = colSums(pop)

YHR_overall = (YHR[1, ] * pop[1, ] + YHR[2, ] * pop[2, ]) / colSums(pop)

# one row per risk group
high_risk_ratio = c(0.082825, 0.141121, 0.165298, 0.329912, 0.470568)
high_risk = round(high_risk_ratio * totalpop)

totalpop = colSums(pop)
relpop = totalpop / sum(totalpop)

# c("<5", "5-17", "18-49", "50-64", "65+"
# phi = matrix(c(2.1609408, 2.1641174,  4.119176, 0.8085660, 0.2808456,
#                0.5973413, 8.1469701,  5.406032, 0.7366175, 0.2261466,
#                0.3822025, 2.4313917, 10.236284, 1.7006885, 0.2102197,
#                0.3523967, 1.8851003,  6.706929, 3.0596349, 0.5003303,
#                0.1897561, 0.8929089,  2.390370, 1.2034597, 1.2322626),
#               nrow=5, ncol=5, byrow=TRUE)

phi = matrix(c(
  1.8763, 2.023 , 4.0143, 0.7924, 0.2806,
  0.5499, 7.0608, 5.0218, 0.6974, 0.2249,
  0.37  , 2.1869, 8.7179, 1.4471, 0.2072,
  0.3283, 1.6234, 5.786 , 2.7861, 0.4968,
  0.1879, 0.8797, 2.3552, 1.1865, 1.2246),
  5, 5,  byrow=TRUE)


P = 0.44  # proportion of pre-symptomatic transmission (%)
tau= 0.57  # symp_prop
eta = 1 / 5.9  # t_symp_onset_rate
#// HFR = c(0.04, 0.1236, 0.03122, 0.10745, 0.23158)  # death probability
HFR = c(0.04, 0.1236, 0.022, 0.16, 0.27)  # death probability, James estimate
#// HFR = c(0.04, 0.1236, 0.025, 0.13, 0.25)  # inthe middle
sigma =  1 / 2.9
omega_Y = 1
omega_A = 0.66
rho_Y = 1 / 2.3 # pre-symptomatic to symptomatic rate
rho_A = rho_Y
gamma_Y = 0.25  # recovery rate
gamma_A = gamma_Y
omega_P = P / (1.0 - P) *
  (tau * omega_Y * (YHR / eta + (1 - YHR) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
  (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))
omega_P_overall = P / (1.0 - P) *
  (tau * omega_Y * (YHR_overall / eta + (1 - YHR_overall) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
  (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))


covid_constants = Csnippet(paste0(
  mat2cpp( phi, "av_contacts"),
  mat2cpp( YHR, "YHR"),
  vec2cpp( totalpop, "N"),
  vec2cpp( pop, "pop"),
  vec2cpp( omega_P_overall, "omega_P"),
  vec2cpp( HFR, "HFR")
))

# --------------------

expand_state = function(names, dims) {
  dims=lapply(dims,seq.int)
  res = names
  for(l in dims)
    res = outer(res,l,paste,sep="_")
  c(res)
}


covid_statenames = c(
  expand_state("H", c(2, 5)),
  expand_state("new_H", c(2, 5)),
  expand_state("recovering_H", c(2, 5)),
  expand_state("dying_H", c(2, 5)),
  expand_state("leaving_H", c(2, 5)),
  expand_state("leaving_IY", c(2, 5)),
  expand_state("S", c(2, 5)),
  expand_state("PA", c(2, 5)),
  expand_state("PY", c(2, 5)),
  expand_state("E", c(2, 5)),
  expand_state("EV", c(2, 5)),
  expand_state("leaving_S_rate", c(2, 5)),
  expand_state("new_E", c(2, 5)),
  expand_state("new_EV", c(2, 5)),
  expand_state("leaving_E", c(2, 5)),
  expand_state("leaving_EV", c(2, 5)),
  expand_state("IA", c(2, 5)),
  expand_state("IY", c(2, 5)),
  expand_state("R", c(2, 5)),
  expand_state("D", c(2, 5)),
  expand_state("new_V_S", c(2, 5)),
  expand_state("V", c(2, 5)),
  "NI",
  "NH",
  "LH",
  "NIY",
  "Beta",
  "Z",
  "total_new_H",
  "total_recovering_H",
  "total_dying_H",
  "total_H",
  "Z_mu",
  "Z_R",
  "b1",  ## DLM
  "b2"  ## DLM
)


#   Beta = exp(   log_beta_0  
#                   + (b1 +b1_0)       * PC1 
# //                + 0.5 * b12 * (PC1 * PC1 - 1)
#                   + (b2 +b2_0)       * PC2
# //                + b3        * PC3 
# //                + b4        * weekend 
# //                + b5        * awareness 
#                   + Z                           );

# These are the variables above that change over time
beta_state_vars = c("b1","b2","Z")


# covid_accum_vars = c(
#   "NI",
#   "NH",
#   "LH",
#   "NIY"
# )
covid_accum_vars = c()

covid_indexed_states = c("S", "E", "PY", "PA", "IY", "IA",
                         "H", "R", "D", "leaving_S_rate",
                         "dying_H", "recovering_H", "new_H", "leaving_H",
                         "leaving_IY", "new_E", "leaving_E",
                         "V")
covid_indexed_covars = c("new_V")

# time t=1,..,T
# data Y_t   
# latent Z_t
# reconstruct Z_t from the likelihood



# Z_t = (S_t, E_t, I_t, ...)
# filtered states:  p(Z_T | Y_{1:T}) latent states today given the history
# smoothed states:  p(Z_1:T | Y_{1:T}) latent state at any point t given all history
# forecasted states:

# ---------------------

covid_rprocess = Csnippet(paste0(
  enable_indexes(covid_indexed_states, dims=c(2, 5)),
  enable_indexes( covid_indexed_covars , dims=c(2, 5)),
  covid_constants, #'
  '
  // Exposed
  //#define BIN_dt( N, rate) ( ( (N) > 0.5) ? rbinom( (N) , 1.0 - exp(         -(rate)  * dt) ) : 0.0 )
// The code below gives 1-(1-rate)^dt instead of 1-exp(-rate*dt) above.
  #define sigmoid(x) (2.0*expit(2.0*(x))-1.0)
  #define BIN_dt( N, rate) ( ( (N) > 0.5) ? rbinom( (N) , 1.0 - exp(  log(1.0-(sigmoid(rate))) * dt) ) : 0.0 )
  #define BIN(    N, rate) ( ( (N) > 0.5) ? rbinom( (N) , (rate)                   ) : 0.0 )

//  b1 = b1 *0.99     + rnorm( 0.0, sig_delta_b1 * sqrt(dt)) ; //DLM
//  b2 = b2 *0.99     + rnorm( 0.0, sig_delta_b2 * sqrt(dt)) ; //DLM

  b1 = b1      + rnorm( 0.0, sig_delta_b1 * sqrt(dt)) * (1-future) ; //DLM
  b2 = b2      + rnorm( 0.0, sig_delta_b2 * sqrt(dt)) * (1-future) ; //DLM

  
  Z = future==1 ? Z: psi * Z + rnorm( 0.0, sig_Z        * sqrt(dt));
  
  
  // Z = Z + rnorm(0.0, sig_Z * sqrt(dt));
// Any variable here needs to be accounted for in beta_state_vars
  Beta = exp(   log_beta_0  
                  + (b1 +b1_0)       * PC1 
//                + 0.5 * b12 * (PC1 * PC1 - 1)
                  + (b2 +b2_0)       * PC2
//                + b3        * PC3 
//                + b4        * weekend 
//                + b5        * awareness 
                  + Z    
                  +dZ
              );
  // rate for being infected
//Beta = 0.0107244365  ;
  
  for (int g=0; g<2; g++) {
    for (int i=0; i<5; i++) {
      // rate of infectiousness
      leaving_S_rate[g][i] = 0.0;
      for (int j=0; j<5; j++) {
        if( debug>0) {
          Rprintf(\"N[j]:%g %g  \\n\",  //"
                Beta * av_contacts[i][j] *              
                  (   omega_A * omega_P[j] * (PA[0][j] + PA[1][j]) +
                      omega_Y * omega_P[j] * (PY[0][j] + PY[1][j]) +
                      omega_A * (IA[0][j] + IA[1][j]) +
                      omega_Y * (IY[0][j] + IY[1][j])) ,N[j]) ; //"
        }
        leaving_S_rate[g][i] += Beta * av_contacts[i][j] *              // j is infecting i
                  (   omega_A * omega_P[j] * (PA[0][j] + PA[1][j]) +
                      omega_Y * omega_P[j] * (PY[0][j] + PY[1][j]) +
                      omega_A * (IA[0][j] + IA[1][j]) +
                      omega_Y * (IY[0][j] + IY[1][j])
                  ) / N[j]; // should subtract dead?
        if( debug>0) {
          Rprintf(\"leaving[%d][%d] %g  \\n\",  //"
                  g,i,leaving_S_rate[g][i] ) ;
        }
      }
      // newly infected
//      new_E[g][i] = rbinom(S[g][i], 1.0 - exp(- leaving_S_rate[g][i] * dt));   
        if( debug>0) {
          Rprintf(\"S[%d][%d] %g  \\n\",  //"
                    g,i,S[g][i]) ;
        }
        new_E[g][i] = BIN_dt( S[g][i], leaving_S_rate[g][i] );    
        new_EV[g][i] = BIN_dt( V[g][i], leaving_S_rate[g][i] );    

      }
  }
  
  // Infections
  
  for (int g=0; g<2; g++) {
    for (int i=0; i<5; i++) {
//      leaving_E[g][i] = (E[g][i] > 0) ? rbinom(E[g][i], 1.0 - exp(- sigma * dt)) : 0.0; // end of encubation period
        leaving_E[g][i]  = BIN_dt( E[g][i],  sigma ) ;
        leaving_EV[g][i] = BIN_dt( EV[g][i], sigma ) ;
    }
  }
  
  
  double new_PA[2][5];
  double new_PY[2][5];
  double new_IA[2][5];
  double new_IY[2][5];
  
  for (int g=0; g<2; g++) {
    for (int i=0; i < 5; i++) {
      new_PA[g][i] = BIN(leaving_E[g][i], 1.0 - tau ) + BIN(leaving_EV[g][i], 1.0 - tauV );
      new_PY[g][i] = leaving_E[g][i] + leaving_EV[g][i] - new_PA[g][i];
      new_IA[g][i] = BIN_dt( PA[g][i]       , rho_A ) ;
      new_IY[g][i] = BIN_dt( PY[g][i]       , rho_Y ) ;
    }
  }
  
  
  // Recovering asymptomatics 
  
  double recovering_IA[2][5];
  
  for (int i=0; i<5; i++) {
    for (int g=0; g<2; g++) {
         recovering_IA[g][i] = BIN_dt(IA[g][i], gamma_A ) ;
    }
  }
  
  // Hospitalization and recovering symptomatics
  
  double recovering_IY[2][5];
  
  for (int g=0; g<2; g++) {
    for (int i=0; i<5; i++) {
      double pi =  gamma_Y * YHR[g][i] / (eta + (gamma_Y - eta) * YHR[g][i]);
      double leaving_IY_rate = ((1.0 - pi) * gamma_Y + pi * eta);
      leaving_IY[g][i]    = BIN_dt( IY[g][i], leaving_IY_rate ) ;
      new_H[g][i]         = BIN(leaving_IY[g][i], pi * eta / leaving_IY_rate );
      recovering_IY[g][i] = leaving_IY[g][i] - new_H[g][i];
    }
  }
  
  // Dying and Recovering from Hospitals
  Z_mu = psi_mu * Z_mu  + rnorm(0.0, sig_mu * sqrt(dt));
  Z_R  = psi_R  * Z_R   + rnorm(0.0, sig_R  * sqrt(dt));
  
  double mu      = exp(log_mu_0 + Z_mu);
  double gamma_H = exp(log_gamma_H_0 + Z_R);
  
  for (int g=0; g<2; g++) {
    for (int i=0; i<5; i++) {
      double nu =  HFR[i] * gamma_H / (mu +  (gamma_H - mu) * HFR[i]);
      double leaving_H_rate = nu * mu + (1.0 - nu) * gamma_H;
      leaving_H[g][i] = BIN_dt(H[g][i], leaving_H_rate ) ;
      recovering_H[g][i] = BIN(leaving_H[g][i], (1 - nu) * gamma_H / leaving_H_rate) ;
      dying_H[g][i] = leaving_H[g][i] - recovering_H[g][i];
    }
  }
  
  // Update states
  
  total_H = 0.0;
  total_new_H = 0.0;
  total_recovering_H = 0.0;
  total_dying_H = 0.0;
  
  for (int g=0; g<2; g++) {
    for (int i=0; i<5; i++) {
      new_V_S[g][i] = rhyper( S[g][i], pop[g][i]-S[g][i], new_V[g][i]) ;
      S[g][i] += - new_E[g][i] -new_V_S[g][i];
      E[g][i]  += new_E[g][i]  - leaving_E[g][i];
      EV[g][i] += new_EV[g][i] - leaving_EV[g][i];
      PA[g][i] += new_PA[g][i] - new_IA[g][i];
      PY[g][i] += new_PY[g][i] - new_IY[g][i];
      IA[g][i] += new_IA[g][i] - recovering_IA[g][i];
      IY[g][i] += new_IY[g][i] - leaving_IY[g][i];
      H[g][i] += new_H[g][i] - leaving_H[g][i];
      R[g][i] += recovering_IA[g][i] + recovering_IY[g][i] + recovering_H[g][i];
      R[g][i] -= ( new_V[g][i] - new_V_S[g][i] ) ; #// move recovered into vaccinated
      D[g][i] += dying_H[g][i];
      NI += leaving_E[g][i];
      NH += new_H[g][i];
      LH += leaving_H[g][i];
      NIY += new_IY[g][i];
      V[g][i] += new_V[g][i] ;
      
      // Totals for the likelihood
      total_H            += H[g][i];
      total_new_H        += new_H[g][i];
      total_recovering_H += recovering_H[g][i];
      total_dying_H      += dying_H[g][i];
    }
  }
  ',
  #'
  register_indexed_changes(covid_indexed_states, dims=c(2, 5))
))

covid_dmeas_admitdischarge_nb = Csnippet(paste0(
  #"
  "
  double loglike = 0.0;
  
  // clamp r for stability
  double r_ = r + 1e-6;

  // negative binomial
  if(  !R_IsNA( adm_total ) )  
    loglike += dnbinom_mu( adm_total,    r_, total_new_H, 1);
  if( !R_IsNA( recov_total) )
    loglike += dnbinom_mu( recov_total,       r_, total_recovering_H, 1);
  if( !R_IsNA( deaths_total) )
    loglike += dnbinom_mu( deaths_total,      r_, total_dying_H,      1);
    
  if( !R_IsNA( hospitalized) )
    loglike +=  dnbinom_mu( hospitalized, r_, total_H,            1); // annealed cause not necessary // Michael: stronger anneal, was 0.5*
  
  // prior
  loglike += dexp(r, 1.0, 1) / maxT; //
  // loglike += dgamma(sig_Z, 1.1, 1.1, 1) / maxT; //
  // loglike += dgamma(sig_mu, 1.1, 1.1, 1) / maxT; //
  // loglike += dgamma(sig_R, 1.1, 1.1, 1) / maxT; //
  loglike += dgamma(sig_delta_b1, 1.1, 1/1.1, 1) / maxT; //## DLM
  loglike += dgamma(sig_delta_b2, 1.1, 1/1.1, 1) / maxT; //## DLM
  
  loglike += dnorm(b1_0, 0.0, 1.0, 1) / maxT; //
  loglike += dnorm(b2_0, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b2, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b3, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b5, 0.0, 1.0, 1) / maxT; //
  

  lik = (give_log) ? loglike : exp(loglike);
  "
  #"
))

covid_rmeas = Csnippet(
  #'
  '
  adm_total = rnbinom_mu(r - 1e-6, total_new_H) ;            
  ' #'
)

init_vec_var = function(values, name) {
  out = c()
  for (i in 1:length(values))
    out[length(out) + 1] = sprintf("%s_%s=%s;\n", name, i, values[i])
  paste(out, collapse="")
}

init_mat_var = function(values, name) {
  out = c()
  for (i in 1:nrow(values))
    for (j in 1:ncol(values))
      out[length(out) + 1] = sprintf("%s_%s_%s=%s;\n", name, i, j, values[i, j])
    paste(out, collapse="")
}

# c("<5", "5-17", "18-49", "50-64", "65+"
init_infected = matrix(c(0, 0, 0, 0, 0,
                         0, 0, 1, 0, 0), 2, 5, byrow = TRUE)
zero_mat = matrix(0, 2, 5)

rinit = Csnippet(paste0(
  init_mat_var(pop, name="S"),
  init_mat_var(zero_mat, name="E"),
  init_mat_var(zero_mat, name="PA"),
  init_mat_var(zero_mat, name="PY"),
  init_mat_var(zero_mat, name="IA"),
  init_mat_var(zero_mat, name="leaving_IY"),
  init_mat_var(init_infected, name="IY"),
  init_mat_var(zero_mat, name="H"),
  init_mat_var(zero_mat, name="new_H"),
  init_mat_var(zero_mat, name="leaving_H"),
  init_mat_var(zero_mat, name="recovering_H"),
  init_mat_var(zero_mat, name="dying_H"),
  init_mat_var(zero_mat, name="R"),
  init_mat_var(zero_mat, name="D"),
  init_mat_var(zero_mat, name="V"),
  #"
  "
  NI=1;
  NIY=1;
  NH=0;
  LH=0;
  Z=Z_0;
  Z_mu=0;
  Z_R=0;
  total_new_H=0.0;
  total_recovering_H=0.0;
  total_dying_H=0.0;
  total_H=0;
  b1 = 0; // DLM
  b2 = 0; // DLM
  " #"
))






covid_init_pars =
  c(
    r = 15.19,
    psi = 0.97,  # AR(1) coefficient
    sig_Z = 0.04,  # for the transmission prob
    psi_mu = 0.97, # AR(1) coefficient
    sig_mu = 0.53, # for the transmission prob
    psi_R = 0.99, #0.77, # AR(1) coefficient #MICHAEL TRYING recovery that is more variable
    sig_R = 0.01, # for the transmission prob
    log_mu_0 = log(1/12),
    log_gamma_H_0 = log(1/8.7),
    log_beta_0 = log(0.06),
    Z_0=0.0,
#//    b1 = 0.25,  ## coefficient PC1 for transmission rate // DLM
#//   b12 = 0.007,
#//    b2 = 0.22,  ## coefficient PC2 for transmission rate
#//    b3 = -0.24, ## coefficient PC3 for transmission rate
#//    b4 = 0.0,
#//    b5 = -0.34,
    b1_0 = 0.25,         ## DLM
    sig_delta_b1 = 0.02, ## DLM
    b2_0 = 0.22,         ## DLM
    sig_delta_b2 = 0.02, ## DLM
    sigma = sigma, ## Incubation period (no transmission)
    tau= tau, # symp_prop
    eta = eta,  # t_symp_onset_rate
    omega_Y = omega_Y,
    omega_A = omega_A,
    rho_Y = rho_Y,  # pre-symptomatic rate
    rho_A = rho_A,
    gamma_Y = gamma_Y,  # recovery rate
    gamma_A = gamma_Y
  )

covid_paramnames = names(covid_init_pars)
