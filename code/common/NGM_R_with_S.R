get_rt = function(beta,
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
                  S_prop) {
  
  n_age = nrow(phi)
  n_c = 5  # [E,PA,PY,IA,IY]
  
  # zero matrix useful recycle below
  Zm = matrix(0, n_age, n_age)
  
  # Proportion of susceptibles
  S_prop_mat = matrix(S_prop,nrow = 5,ncol = 5,byrow=FALSE)
  
  
  # tiled av contact matrix
  big_C = matrix(0, n_c * n_age, n_c * n_age)
  for (i in 1:n_c)
    for (j in 1:n_c) {
      start_row = n_c * (i - 1) + 1
      end_row = n_c * i
      start_col = n_c * (j - 1) + 1
      end_col = n_c * j
      big_C[start_row:end_row, start_col:end_col] = phi 
    }
  
  # Transmission and Transition Matrices
  
  block11 = Zm * S_prop_mat
  block12 = beta * pop_ratio_mat * matrix(omega_P * omega_A, n_age, n_age, byrow = TRUE) * S_prop_mat
  block13 = beta * pop_ratio_mat * matrix(omega_P * omega_Y, n_age, n_age, byrow = TRUE) * S_prop_mat
  block14 = beta * pop_ratio_mat * matrix(omega_A, n_age, n_age) * S_prop_mat
  block15 = beta * pop_ratio_mat * matrix(omega_Y, n_age, n_age) * S_prop_mat
  row1 = cbind(block11, block12, block13, block14, block15)
  row2 = matrix(0, 4 * n_age, 5 * n_age)
  T_mat = big_C * rbind(row1, row2)
  
  row1 = cbind(diag(-(sigma), n_age),
               Zm, Zm, Zm, Zm)
  row2 = cbind(diag((1 - tau) * (sigma), n_age),
               diag(-rho_A, n_age),
               Zm, Zm, Zm)
  row3 = cbind(diag(tau * (sigma), n_age),
               Zm,
               diag(-rho_Y, n_age),
               Zm, Zm)
  row4 = cbind(Zm,
               diag(rho_A, n_age),
               Zm,
               diag(-gamma_A, n_age),
               Zm)
  row5 = cbind(Zm, Zm,
               diag(rho_Y, n_age),
               Zm,
               diag(-gamma_Y, n_age))
  Sigma_mat = rbind(row1, row2, row3, row4, row5)
  
  
  # R0 computation
  K_L = - T_mat %*% solve(Sigma_mat)
  r0 = max(as.numeric(eigen(K_L)$values))
  
  r0
}


