#// define the global data --------


#// this is the overall
#// symp_hosp_ratio = c(0.00048721, 0.00048721, 0.03287572, 0.11337395, 0.17733063)
#// YHR_overall = c(0.0007018, 0.0007018, 0.04735, 0.1633, 0.2544)  # symp_hosp_ratio =
YHR = matrix(dimnames = list( risk=1:2,
  age=c(
           "<5",    "5-17",  "18-49",  "50-64", "65+")),
  c(  0.004021 ,  0.003091 ,  0.1903 ,  0.4114 ,  0.4879 ,
      0.0004021,  0.0003091,  0.01903,  0.04114,  0.04879),
  2, 5, byrow = TRUE)


#// prob_group = matrix(c(high_risk_ratio,
#//                       1 - high_risk_ratio), 2, 5, byrow = TRUE)
# pop = matrix(dimnames = list( risk=1:2,
#   age=c(
#            "<5",    "5-17",  "18-49",  "50-64", "65+")),
#       c(   9350,     37451,   156209,   108196, 103763,
#          128527,    327148,   915894,   249273, 132505),
#              2, 5, byrow = TRUE)
# totalpop = colSums(pop)

# YHR_overall = (YHR[1, ] * pop[1, ] + YHR[2, ] * pop[2, ]) / colSums(pop)

#// one row per risk group
high_risk_ratio = c(0.082825, 0.141121, 0.165298, 0.329912, 0.470568)
# high_risk = round(high_risk_ratio * totalpop)

# relpop = totalpop / sum(totalpop)

# c("<5", "5-17", "18-49", "50-64", "65+"
# phi = matrix(c(2.1609408, 2.1641174,  4.119176, 0.8085660, 0.2808456,
#                0.5973413, 8.1469701,  5.406032, 0.7366175, 0.2261466,
#                0.3822025, 2.4313917, 10.236284, 1.7006885, 0.2102197,
#                0.3523967, 1.8851003,  6.706929, 3.0596349, 0.5003303,
#                0.1897561, 0.8929089,  2.390370, 1.2034597, 1.2322626),
#               nrow=5, ncol=5, byrow=TRUE)

phi = matrix(dimnames = list(   age1=c(
        "<5", 
        "5-17",
        "18-49",
        "50-64",
        "65+"),
  age=c(
        "<5", "5-17","18-49","50-64", "65+")),
  c(  1.8763, 2.023 , 4.0143, 0.7924, 0.2806,
      0.5499, 7.0608, 5.0218, 0.6974, 0.2249,
      0.37  , 2.1869, 8.7179, 1.4471, 0.2072,
      0.3283, 1.6234, 5.786 , 2.7861, 0.4968,
      0.1879, 0.8797, 2.3552, 1.1865, 1.2246),
  5, 5,  byrow=TRUE)


P = 0.44  # proportion of pre-symptomatic transmission (%)
tau= 0.57  # symp_prop
eta = 1 / 5.9  # t_symp_onset_rate
# HFR = c(0.04, 0.1236, 0.03122, 0.10745, 0.23158)  # death probability
HFR = c(0.04, 0.1236, 0.022, 0.16, 0.27)  # death probability, //James' estimate
# HFR = c(0.04, 0.1236, 0.025, 0.13, 0.25)  # inthe middle
sigma =  1 / 2.9
omega_Y = 1
omega_A = 0.66
rho_Y = 1 / 2.3 # pre-symptomatic to symptomatic rate
rho_A = rho_Y
gamma_Y = 0.25  # recovery rate
gamma_A = gamma_Y

# omega_P = P / (1.0 - P) *
#   (tau * omega_Y * (YHR / eta + (1 - YHR) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
#   (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))
# omega_P_overall = P / (1.0 - P) *
#   (tau * omega_Y * (YHR_overall / eta + (1 - YHR_overall) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
#   (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))
