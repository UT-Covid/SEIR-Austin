get_rtV = function(beta,
                  phi, # av contact matrix
                  pop,
                  pop_ratio_mat,
                  sigma,
                  tau,
                  rho_A,
                  rho_Y,
                  gamma_A,
                  gamma_Y,
                  omega_A,
                  omega_Y,
                  omega_P,
                  S_prop,
                  epsilon_V,
                  tau_V,
                  V_prop) {
  
  n_age = nrow(phi)
  n_c = 6  # [E,EV,PA,PY,IA,IY]
  
  # zero matrix useful recycle below
  Zm = matrix(0, n_age, n_age)
  
  # Proportion of susceptibles and vaccinated
  S_prop_mat = matrix( S_prop,nrow = 5, ncol = 5, byrow=FALSE)
  V_prop_mat = matrix( V_prop,nrow = 5, ncol = 5, byrow=FALSE)
  
  
  # tiled av contact matrix
  big_C = matrix(0, n_c * n_age, n_c * n_age)
  for (i in 1:n_c)
    for (j in 1:n_c) {
      start_row = n_age * (i - 1) + 1
      end_row = n_age * i
      start_col = n_age * (j - 1) + 1
      end_col = n_age * j
      big_C[start_row:end_row, start_col:end_col] = phi 
    }
  
  ## Transmission and Transition Matrices
  # Transmission
  # E rows
  block11 = Zm * S_prop_mat
  block12 = Zm * S_prop_mat
  block13 = beta * pop_ratio_mat * matrix(omega_P * omega_A, n_age, n_age, byrow = TRUE) * S_prop_mat
  block14 = beta * pop_ratio_mat * matrix(omega_P * omega_Y, n_age, n_age, byrow = TRUE) * S_prop_mat
  block15 = beta * pop_ratio_mat * matrix(omega_A, n_age, n_age) * S_prop_mat
  block16 = beta * pop_ratio_mat * matrix(omega_Y, n_age, n_age) * S_prop_mat
  row1 = cbind(block11, block12, block13, block14, block15, block16)
  
  # E^V rows
  beta_V = beta * epsilon_V # relative transmission risk for vaccinated
  block21 = Zm * V_prop_mat
  block22 = Zm * V_prop_mat
  block23 = beta_V * pop_ratio_mat * matrix(omega_P * omega_A, n_age, n_age, byrow = TRUE) * V_prop_mat
  block24 = beta_V * pop_ratio_mat * matrix(omega_P * omega_Y, n_age, n_age, byrow = TRUE) * V_prop_mat
  block25 = beta_V * pop_ratio_mat * matrix(omega_A, n_age, n_age) * V_prop_mat
  block26 = beta_V * pop_ratio_mat * matrix(omega_Y, n_age, n_age) * V_prop_mat
  row2 = cbind(block21, block22, block23, block24, block25, block26)
  
  row3 = matrix(0, 4 * n_age, n_c * n_age)
  T_mat = big_C * rbind(row1, row2, row3)
  
  # Transition
  row1 = cbind(diag(-(sigma), n_age),
               Zm, Zm, Zm, Zm, Zm)
  row2 = cbind(Zm,
               diag(-(sigma), n_age),
               Zm, Zm, Zm, Zm)
  row3 = cbind(diag((1 - tau) * (sigma), n_age),
               diag((1 - tau_V) * (sigma), n_age),
               diag(-rho_A, n_age),
               Zm, Zm, Zm)
  row4 = cbind(diag(tau * (sigma), n_age),
               diag(tau_V * (sigma), n_age),
               Zm,
               diag(-rho_Y, n_age),
               Zm, Zm)
  row5 = cbind(Zm, Zm,
               diag(rho_A, n_age),
               Zm,
               diag(-gamma_A, n_age),
               Zm)
  row6 = cbind(Zm, Zm, Zm,
               diag(rho_Y, n_age),
               Zm,
               diag(-gamma_Y, n_age))
  Sigma_mat = rbind(row1, row2, row3, row4, row5, row6)
  
  
  # R0 computation
  K_L = - T_mat %*% solve(Sigma_mat)
  r0 = max( as.numeric(eigen(K_L)$values))
  
  r0
}


r0.test=F 

if( r0.test) {

beta <- 0.053
pop <- c( 137877.,  364599., 1072103.,  357469.,  236268.)
pop_ratio_mat <- matrix(pop, nrow=5, ncol=5, byrow=FALSE) / matrix(pop, nrow=5, ncol=5, byrow=TRUE)
sigma <- 1.0 / 2.9
tau <- 0.57

dP <- 2.3
rho_A <- 1.0 / dP
rho_Y <- 1.0 / dP

dA <- 4.0
dY <- 4.0
gamma_A <- 1.0 / dA
gamma_Y <- 1.0 / dY
omega_A <- 2.0 / 3.0
omega_Y <- 1.0

preS <- 0.44
omega_P <- preS / (1-preS) * (tau * omega_Y * dY + (1-tau) * omega_A * dA) / 
  (dP * (tau * omega_Y + (1-tau) * omega_A))

S_prop <- rep(1.0,length(pop))

phi <- rbind(c(1.87626037, 2.0230411 , 4.01430722, 0.79240015, 0.28061969),
             c(0.54989386, 7.06081823, 5.02181545, 0.69740545, 0.22488401),
             c(0.37001112, 2.18689999, 8.71793613, 1.44712065, 0.20719312),
             c(0.3283189 , 1.62344686, 5.78598967, 2.78607592, 0.49678338),
             c(0.18794811, 0.87971799, 2.35522847, 1.18650896, 1.22458335))

# New parameters to include vaccination
epsilon_V <- 0.055 # Relative susceptibility of vaccinated individuals
tau_V <- 0.023 # Symptomatic proportion for vaccinated individuals
V_prop <- rep(0.0,length(pop)) # proportion of population vaccinated in each age group

rt0 <- get_rt(beta,
              phi, # av contact matrix
              pop,
              pop_ratio_mat,
              sigma,
              tau,
              rho_A,
              rho_Y,
              gamma_A,
              gamma_Y,
              omega_A,
              omega_Y,
              omega_P,
              S_prop,
              epsilon_V,
              tau_V,
              V_prop
              )

print(paste('Rt calculated:',rt0))


### Checks
## Direct scaling using symptomatic proportion and relative susceptibility
S_prop <- rep(1.0,length(pop))
V_prop <- rep(0.0,length(pop)) # proportion of population vaccinated in each age group
rt_all_S <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                  omega_Y,omega_P,S_prop,epsilon_V,tau_V,V_prop)

S_prop <- rep(0.0,length(pop))
V_prop <- rep(1.0,length(pop)) # proportion of population vaccinated in each age group
rt_all_V <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                   omega_Y,omega_P,S_prop,epsilon_V,tau_V,V_prop)

# This should equal rt_all_V
check1 <- rt_all_S * epsilon_V * (tau_V * omega_Y + (1-tau_V) * omega_A) / (tau * omega_Y + (1-tau) * omega_A)
print(paste('Check 1, the 2 following values should be the same:',rt_all_V,'==',check1))


## If proportion of symptomatic was the same for vaccinated
S_prop <- rep(0.0,length(pop))
V_prop <- rep(1.0,length(pop)) # proportion of population vaccinated in each age group
rt_all_V_same_tau <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                   omega_Y,omega_P,S_prop,epsilon_V,tau,V_prop)

# This should equal rt_all_V_same_tau
check2 <- rt_all_S * epsilon_V
print(paste('Check 2, the 2 following values should be the same:',rt_all_V_same_tau,'==',check2))


## If susceptibility was the same for vaccinated
S_prop <- rep(0.0,length(pop))
V_prop <- rep(1.0,length(pop)) # proportion of population vaccinated in each age group
rt_all_V_same_suscep <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                            omega_Y,omega_P,S_prop,1.0,tau_V,V_prop)

# This should equal rt_all_V_same_suscep
check3 <- rt_all_S * (tau_V * omega_Y + (1-tau_V) * omega_A) / (tau * omega_Y + (1-tau) * omega_A)
print(paste('Check 3, the 2 following values should be the same:',rt_all_V_same_suscep,'==',check3))


## Simplest: if everyone was vaccinated with same parameters as susceptible non-vaccinated
S_prop <- rep(0.0,length(pop))
V_prop <- rep(1.0,length(pop)) # proportion of population vaccinated in each age group
rt_all_V_same_params <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                               omega_Y,omega_P,S_prop,1.0,tau,V_prop)

# This should equal rt_all_V_same_params
check4 <- rt_all_S
print(paste('Check 4, the 2 following values should be the same:',rt_all_V_same_params,'==',check4))


##  If half the population is vaccinated and half is susceptible
S_prop <- rep(0.5,length(pop))
V_prop <- rep(0.5,length(pop)) # proportion of population vaccinated in each age group
rt_all_split <- get_rt(beta,phi,pop,pop_ratio_mat,sigma,tau,rho_A,rho_Y,gamma_A,gamma_Y,omega_A,
                   omega_Y,omega_P,S_prop,epsilon_V,tau_V,V_prop)

check5 <- 0.5 * (rt_all_S + rt_all_V)
print(paste('Check 5, the 2 following values should be the same:',rt_all_split,'==',check5))


}