rm(list=ls())
library(coda)
library(readxl)
library(tidyverse)
library(patchwork)
##########full code################
dat=data.frame(setting="Vancouver",
               Sex_groups=c("0-4","4-22","22-100"),
               mean_rate=c(1.31,9.13,39.84),
               mean_rate_dy=c(1.31,9.13,39.84)/180,
               Prop=c(0.60,0.35,0.05),
               Pop=c(0.60,0.35,0.05)*26100)

contact <- dat$mean_rate / 180 ##convert the 6months sexual contact rates to daily contact rates
contact_dy<-dat$mean_rate_dy  ###daily contact rate for each of the three groups
Norm_contact_dy<-contact_dy/sum(contact_dy) ##Normalize to ensure prop of E reflect the relative risk of exposure based on contact rates 
sum(Norm_contact_dy)
#contact <- dat$mean_rate / 26 ##convert the 6months sexual contact rates to weekly contact rates
N <- c(15660, 9135, 1305) ##population sizes for the different groups
burn_in_timesteps <- 30
T=169+burn_in_timesteps
ww_flow=read_excel("Data/case_data_V2.xlsx",sheet="flowrate")
flow_all<-ww_flow %>%select(sumflow)
flow_alldaily<-as.numeric(unlist(flow_all))
flow_mlalldaily<-c(rep(0, burn_in_timesteps), flow_alldaily*1e3)
flow_mlalldaily

##1. case data ################
case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat$total_cases
Cumulat_inc=case_dat$Cumulat_Inc

##2. WW data #########################
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
ww_std = ww_dat %>% select(log10_cp_per_person_per_day) #####standardised WW data
ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))
obs_days <- ww_dat$Day

###########Test the code###############
# Example: draw one row
extract_params <- function(draw_row) {
  list(
    beta = draw_row["beta"],
    kappa = draw_row["kappa"],
    phi = draw_row["phi"],
    Vea = draw_row["Vea"],
    rsa = 1 - draw_row["Vea"],
    Veb = draw_row["Veb"],
    rsb = 1 - draw_row["Veb"],
    report_frac = draw_row["report_frac"],
    report_lag = 7,
    m = draw_row["m"],
    mva = draw_row["m"],
    mvb = draw_row["m"],
    delta = 1 / draw_row["delta_inv"],
    epsilon = 1 / 2,
    theta = 1 / (draw_row["theta_invall"] / 3),
    omega = 1 / (draw_row["omega_invall"] / 3),
    deltava=1 / draw_row["delta_inv"],
    epsilonva=1/2,
    omegava=1 / (draw_row["omega_invall"] / 3),
    thetava=1 / (draw_row["theta_invall"] / 3),
    deltavb=1 / draw_row["delta_inv"],
    epsilonvb=1/2,
    omegavb=1 / (draw_row["omega_invall"] / 3),
    thetavb=1 / (draw_row["theta_invall"] / 3),
    mu = 0.18,
    tau = 0.776,
    shed_I1 = 10^12.3,
    shed_I2 = 10^11.5,
    shed_I3 = 10^10.2,
    shed_A1 = 10^12.3 * 0.776,
    shed_A2 = 10^11.5 * 0.776,
    shed_A3 = 10^10.2 * 0.776,
    shed_alpha = 10^12.3 * 0.776,
    shed_P1a=10^12.3 * 0.776,
    shed_P1b=10^12.3 * 0.776,
    shed_A1a= 10^12.3 * 0.776,
    shed_A2a= 10^11.5 * 0.776,
    shed_A3a= 10^10.2 * 0.776,
    shed_A1b= 10^12.3 * 0.776,
    shed_A2b= 10^11.5 * 0.776,
    shed_A3b= 10^10.2 * 0.776,
    shed_I1a = 10^12.3,
    shed_I2a = 10^11.5,
    shed_I3a = 10^10.2,
    shed_I1b = 10^12.3,
    shed_I2b = 10^11.5,
    shed_I3b = 10^10.2,
    transit_mean = draw_row["transit_time_mean"],
    transit_cv = draw_row["transit_time_cv"],
    tau_ww = draw_row["tau_ww"],
    mult = exp(draw_row["log_mult"]),
    va_delay = 14,
    vb_delay = 14,
    flow_mlalldaily=flow_mlalldaily,
    wwtp_population = 1278020
  )
}

set.seed(123)
load("D:/mpox25output/Comb_finaltestf.RData")
mcmc_matrixallcom <- as.matrix(do.call(rbind, Comb_finaltestf))
one_draw <- mcmc_matrixallcom[1000, ]
params <- extract_params(one_draw)

# Simulation function with vaccination scenario flexibility
mysimulate_model <- function(params, init, contact, N, T, 
                             vax_scenario = "baseline", 
                             v1_frac = 1, 
                             v2_frac = 1) {
  G <- length(N)
  
  alpha_v <- rep(0, T)
  mu_v <- rep(0, T)
  
  if (vax_scenario == "baseline") {
    alpha_v <- init$alpha_v
    mu_v <- init$mu_v
    
    # Force vaccinated compartments to shed like unvaccinated
    params$shed_P1a <- params$shed_alpha
    params$shed_P1b <- params$shed_alpha
    
    params$shed_I1a <- params$shed_I1
    params$shed_I2a <- params$shed_I2
    params$shed_I3a <- params$shed_I3
    params$shed_I1b <- params$shed_I1
    params$shed_I2b <- params$shed_I2
    params$shed_I3b <- params$shed_I3
    
    params$shed_A1a <- params$shed_A1
    params$shed_A2a <- params$shed_A2
    params$shed_A3a <- params$shed_A3
    params$shed_A1b <- params$shed_A1
    params$shed_A2b <- params$shed_A2
    params$shed_A3b <- params$shed_A3
    
  } else if (vax_scenario == "reduced_shedding") {
    alpha_v <- init$alpha_v
    mu_v    <- init$mu_v
    
    # Apply dose-specific shedding multipliers (remaining fraction)
    # Dose 1
    params$shed_P1a <- params$shed_alpha * v1_frac
    params$shed_I1a <- params$shed_I1    * v1_frac
    params$shed_I2a <- params$shed_I2    * v1_frac
    params$shed_I3a <- params$shed_I3    * v1_frac
    params$shed_A1a <- params$shed_A1    * v1_frac
    params$shed_A2a <- params$shed_A2    * v1_frac
    params$shed_A3a <- params$shed_A3    * v1_frac
    
    # Dose 2
    params$shed_P1b <- params$shed_alpha * v2_frac
    params$shed_I1b <- params$shed_I1    * v2_frac
    params$shed_I2b <- params$shed_I2    * v2_frac
    params$shed_I3b <- params$shed_I3    * v2_frac
    params$shed_A1b <- params$shed_A1    * v2_frac
    params$shed_A2b <- params$shed_A2    * v2_frac
    params$shed_A3b <- params$shed_A3    * v2_frac
    
    message(sprintf("[Scenario: reduced shedding] v1_frac=%.2f, v2_frac=%.2f", 
                    v1_frac, v2_frac))
    
  } else {
    stop("vax_scenario must be 'baseline' or 'reduced_shedding'")
  }
  

  # Extract parameters
  va_delay <- params$va_delay
  vb_delay <- params$vb_delay
  beta <- params$beta
  kappa <- params$kappa
  phi <- params$phi
  rsa <- params$rsa
  rsb <- params$rsb
  report_frac <- params$report_frac
  report_lag <- params$report_lag
  m <- params$m
  mva <- params$mva
  mvb <- params$mvb
  delta <- params$delta
  epsilon <- params$epsilon
  theta <- params$theta
  omega <- params$omega
  deltava <- params$deltava
  epsilonva <- params$epsilonva
  thetava <- params$thetava
  omegava <- params$omegava
  deltavb <- params$deltavb
  epsilonvb <- params$epsilonvb
  thetavb <- params$thetavb
  omegavb <- params$omegavb
  mu <- params$mu
  tau <- params$tau
  shed_I1 <- params$shed_I1
  shed_I2 <- params$shed_I2
  shed_I3 <- params$shed_I3
  shed_A1 <- params$shed_A1
  shed_A2 <- params$shed_A2
  shed_A3 <- params$shed_A3
  shed_alpha <- params$shed_alpha
  shed_P1a<-params$shed_P1a
  shed_P1b<-params$shed_P1b
  shed_A1a <- params$shed_A1a
  shed_A2a <- params$shed_A2a
  shed_A3a <- params$shed_A3a
  shed_A1b <- params$shed_A1b
  shed_A2b <- params$shed_A2b
  shed_A3b <- params$shed_A3b
  shed_I1a <- params$shed_I1a
  shed_I2a <- params$shed_I2a
  shed_I3a <- params$shed_I3a
  shed_I1b <- params$shed_I1b
  shed_I2b <- params$shed_I2b
  shed_I3b <- params$shed_I3b
  transit_mean <- params$transit_mean
  transit_cv <- params$transit_cv
  tau_ww <- params$tau_ww
  mult <- params$mult
  flow <- params$flow_mlalldaily
  wwtp_population <- params$wwtp_population
  
  
  va_queue <- array(0, dim = c(G, va_delay, T))
  vb_queue <- array(0, dim = c(G, vb_delay, T))
  
  # For vaccine effectiveness kick-in
  Va_eff <- matrix(0, G, T)
  Vb_eff <- matrix(0, G, T)
  shed_P<-matrix(0, G, T)
  shed_A<-matrix(0, G, T)
  shed_I<-matrix(0, G, T)
  P_total<-matrix(0, G, T)
  A_total<-matrix(0, G, T)
  I_total<-matrix(0, G, T)
  # Compartment initialization function
  comp_matrix <- function() matrix(0, G, T)
  S <- comp_matrix(); E <- comp_matrix(); P <- comp_matrix()
  A1 <- comp_matrix(); A2 <- comp_matrix(); A3 <- comp_matrix()
  I1 <- comp_matrix(); I2 <- comp_matrix(); I3 <- comp_matrix()
  R <- comp_matrix(); Cuminc <- comp_matrix()
  Va <- comp_matrix(); Eva <- comp_matrix(); Pva <- comp_matrix()
  A1va <- comp_matrix(); A2va <- comp_matrix(); A3va <- comp_matrix()
  I1va <- comp_matrix(); I2va <- comp_matrix(); I3va <- comp_matrix()
  Vb <- comp_matrix(); Evb <- comp_matrix(); Pvb <- comp_matrix()
  A1vb <- comp_matrix(); A2vb <- comp_matrix(); A3vb <- comp_matrix()
  I1vb <- comp_matrix(); I2vb <- comp_matrix(); I3vb <- comp_matrix()
  prevalence <- comp_matrix(); lambda_det <- comp_matrix()
  
  new_reported_cases <- matrix(0, G, T)
  total_lambda <- numeric(T)
  active_infected <- numeric(T)
  total_new_cases <- numeric(T)
  total_Cuminc <- numeric(T)
  daily_shedding<-numeric(T)
  
  # Mixing matrix setup
  mix_ass <- diag(1, G)
  mix_prp <- outer(contact, contact * N) / sum(contact * N)
  mix <- array(0, dim = c(G, G, T))
  
  for (g in 1:G) {
    for (g_ in 1:G) {
      mix[g, g_, 1] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }
  
  # Initial states
  S[, 1] <- init$S
  E[, 1] <- init$E
  P[, 1] <- init$P
  A1[, 1] <- init$A1
  I1[, 1] <- init$I1
  
  # Calculate initial prevalence and lambda
  for (g in 1:G) {
    prevalence[g, 1] <- sum(c(P[g,1], A1[g,1], I1[g,1])) / N[g]
    lambda_det[g, 1] <- sum(beta * contact[g] * mix[g, , 1] * prevalence[, 1])
  }
  total_lambda[1] <- sum(lambda_det[, 1])
  
  # ---- Initialize transition matrices for each group ----
  new_Es     <- matrix(0, G, T)
  new_Ps     <- matrix(0, G, T)
  new_A1s    <- matrix(0, G, T)
  new_A2s    <- matrix(0, G, T)
  new_A3s    <- matrix(0, G, T)
  new_RAs    <- matrix(0, G, T)
  new_I1s    <- matrix(0, G, T)
  new_I2s    <- matrix(0, G, T)
  new_I3s    <- matrix(0, G, T)
  new_RIs    <- matrix(0, G, T)
  
  new_Vas    <- matrix(0, G, T)
  new_Evas   <- matrix(0, G, T)
  new_Pvas   <- matrix(0, G, T)
  new_A1vas  <- matrix(0, G, T)
  new_A2vas  <- matrix(0, G, T)
  new_A3vas  <- matrix(0, G, T)
  new_RAvas  <- matrix(0, G, T)
  new_I1vas  <- matrix(0, G, T)
  new_I2vas  <- matrix(0, G, T)
  new_I3vas  <- matrix(0, G, T)
  new_RIvas  <- matrix(0, G, T)
  
  new_Vbs    <- matrix(0, G, T)
  new_Evbs   <- matrix(0, G, T)
  new_Pvbs   <- matrix(0, G, T)
  new_A1vbs  <- matrix(0, G, T)
  new_A2vbs  <- matrix(0, G, T)
  new_A3vbs  <- matrix(0, G, T)
  new_RAvbs  <- matrix(0, G, T)
  new_I1vbs  <- matrix(0, G, T)
  new_I2vbs  <- matrix(0, G, T)
  new_I3vbs  <- matrix(0, G, T)
  new_RIvbs  <- matrix(0, G, T)
  
  
  # Dynamics of the model (deterministic-like)
  for (t in 2:T) {
    
    # Final mixing matrix for three groups only
    for (g in 1:3) {
      for (g_ in 1:3) {
        mix[g, g_, t] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
      }
    }
    
    # Calculate prevalence once per time step (deterministic)
    # Calculate prevalence once per time step (deterministic)
    for (g in 1:3) {
      prevalence[g, t] <- ((tau*P[g, t-1]) + I1[g, t-1] + I2[g, t-1] + I3[g, t-1]+
                             (tau*A1[g, t-1]) + (tau*A2[g, t-1]) + (tau*A3[g, t-1])+
                             (tau*Pva[g, t-1])+I1va[g, t-1] + I2va[g, t-1] + I3va[g, t-1]+
                             (tau*A1va[g, t-1])+ (tau*A2va[g, t-1])+ (tau*A3va[g, t-1])+
                             (tau*Pvb[g, t-1])+I1vb[g, t-1] + I2vb[g, t-1] + I3vb[g, t-1]+
                             (tau*A1vb[g, t-1])+ (tau*A2vb[g, t-1])+ (tau*A3vb[g, t-1]))/ N[g]
      
    }
    
    # Calculate force of infection (lambda_det) for each group
    for (g in 1:3) {
      lambda_det[g, t] <- sum(beta * contact[g] * mix[g, 1:3, t] * prevalence[1:3, t])
    }
    total_lambda[t] <- sum(lambda_det[1:3, t])
    
    # Deterministic transitions instead of binomial draws
    for (g in 1:3) {
      new_Es[g, t-1] <- lambda_det[g, t] * S[g, t-1]     # Deterministic S to E transition
      new_Ps[g, t-1] <- delta * E[g, t-1]             # Deterministic E to P transition
      new_A1s[g, t-1] <- epsilon * m * P[g, t-1]      # Deterministic P to A1
      new_A2s[g, t-1] <- omega * A1[g, t-1]           # Deterministic A1 to A2
      new_A3s[g, t-1] <- omega * A2[g, t-1]           # Deterministic A2 to A3
      new_RAs[g, t-1] <- omega * A3[g, t-1]           # Deterministic A3 to R
      new_I1s[g, t-1] <- epsilon * (1 - m) * P[g, t-1] # Deterministic P to I1
      new_I2s[g, t-1] <- theta * I1[g, t-1]           # Deterministic I1 to I2
      new_I3s[g, t-1] <- theta * I2[g, t-1]           # Deterministic I2 to I3
      new_RIs[g, t-1] <- theta * I3[g, t-1]           # Deterministic I3 to R
      
      
      
      # Deterministic transitions for 1st dose vaccination
      
      new_Vas[g, t-1] <- alpha_v[t] * S[g, t-1]           # Deterministic S to Va transition.alpha is vaccination rate
      # Shift the Sv_queue (buffer of people waiting for protection)
      va_queue[g, 1, t] <- new_Vas[g, t-1]
      for (d in 2:va_delay) {
        va_queue[g, d, t] <- va_queue[g, d-1, t-1]
      }
      
      
      Va_eff[g, t] <- va_queue[g, va_delay, t-1]       # Individuals whose protection kicks in today
      new_Evas[g, t-1] <- lambda_det[g, t] *rsa* Va[g, t-1]     # Deterministic Va to Eva transition.Lambda*1-vaccine effectiveness(ve)
      new_Pvas[g, t-1] <- deltava * Eva[g, t-1]             # Deterministic Eva to Pva transition
      new_A1vas[g, t-1] <- epsilonva * mva * Pva[g, t-1]      # Deterministic Pva to A1va
      new_A2vas[g, t-1] <- omegava * A1va[g, t-1]           # Deterministic A1va to A2va
      new_A3vas[g, t-1] <- omegava * A2va[g, t-1]           # Deterministic A2va to A3va
      new_RAvas[g, t-1] <- omegava * A3va[g, t-1]           # Deterministic A3va to Rva
      new_I1vas[g, t-1] <- epsilonva * (1 - mva) * Pva[g, t-1] # Deterministic Pva to I1va
      new_I2vas[g, t-1] <- thetava * I1va[g, t-1]           # Deterministic I1va to I2va
      new_I3vas[g, t-1] <- thetava * I2va[g, t-1]           # Deterministic I2va to I3va
      new_RIvas[g, t-1] <- thetava * I3va[g, t-1]           # Deterministic I3va to Rv
      
      # Deterministic transitions for 2nd dose vaccination
      
      new_Vbs[g, t-1] <- mu_v[t] * Va[g, t-1]              # Deterministic Va to Vb transition.mu is vaccination rate
      # Shift the Sv_queue (buffer of people waiting for protection)
      vb_queue[g, 1, t] <- new_Vbs[g, t-1]
      for (d in 2:vb_delay) {
        vb_queue[g, d, t] <- vb_queue[g, d-1, t-1]
      }
      
      Vb_eff[g, t] <- vb_queue[g, vb_delay, t-1]
      new_Evbs[g, t-1] <- lambda_det[g, t] *rsb* Vb[g, t-1]     # Deterministic Vb to Evb transition.Lambda*1-vaccine effectiveness(ve)
      new_Pvbs[g, t-1] <- deltavb * Evb[g, t-1]             # Deterministic Evb to Pvb transition
      new_A1vbs[g, t-1] <- epsilonvb * mvb * Pvb[g, t-1]      # Deterministic Pvb to A1vb
      new_A2vbs[g, t-1] <- omegavb * A1vb[g, t-1]           # Deterministic A1vb to A2vb
      new_A3vbs[g, t-1] <- omegavb * A2vb[g, t-1]           # Deterministic A2vb to A3vb
      new_RAvbs[g, t-1] <- omegavb * A3vb[g, t-1]           # Deterministic A3vb to RAvb
      new_I1vbs[g, t-1] <- epsilonvb * (1 - mvb) * Pvb[g, t-1] # Deterministic Pvb to I1vb
      new_I2vbs[g, t-1] <- thetavb * I1vb[g, t-1]           # Deterministic I1vb to I2vb
      new_I3vbs[g, t-1] <- thetavb * I2vb[g, t-1]           # Deterministic I2vb to I3vb
      new_RIvbs[g, t-1] <- thetavb * I3vb[g, t-1]           # Deterministic I3vb to RIvb
      
      
      
      # Update compartments deterministically at each time step
      
      S[g, t] <- max(0, S[g, t-1] - new_Es[g, t-1]-new_Vas[g, t-1])
      E[g, t] <- max(0, E[g, t-1] + new_Es[g, t-1] - new_Ps[g, t-1])
      P[g, t] <- max(0, P[g, t-1] + new_Ps[g, t-1] - new_A1s[g, t-1] - new_I1s[g, t-1])
      A1[g, t] <- max(0, A1[g, t-1] + new_A1s[g, t-1] - new_A2s[g, t-1])
      A2[g, t] <- max(0, A2[g, t-1] + new_A2s[g, t-1] - new_A3s[g, t-1])
      A3[g, t] <- max(0, A3[g, t-1] + new_A3s[g, t-1] - new_RAs[g, t-1])
      I1[g, t] <- max(0, I1[g, t-1] + new_I1s[g, t-1] - new_I2s[g, t-1])
      I2[g, t] <- max(0, I2[g, t-1] + new_I2s[g, t-1] - new_I3s[g, t-1])
      I3[g, t] <- max(0, I3[g, t-1] + new_I3s[g, t-1] - new_RIs[g, t-1])
      R[g, t] <- max(0, R[g, t-1] + new_RAs[g, t-1] + new_RIs[g, t-1]+new_RAvas[g, t-1]+new_RIvas[g, t-1]+
                       new_RAvbs[g, t-1]+new_RIvbs[g, t-1])
      
      Cuminc[g, t] <- Cuminc[g, t-1] + new_Es[g, t-1]+ new_Evas[g, t-1] + new_Evbs[g, t-1]
      
      # Add comps with 1st dose vaccination Va
      
      Va[g, t] <- max(0, Va[g, t-1] + Va_eff[g, t] - new_Evas[g, t-1] - new_Vbs[g, t-1])
      Eva[g, t] <- max(0, Eva[g, t-1] + new_Evas[g, t-1] - new_Pvas[g, t-1])
      Pva[g, t] <- max(0, Pva[g, t-1] + new_Pvas[g, t-1] - new_A1vas[g, t-1] - new_I1vas[g, t-1])
      A1va[g, t] <- max(0, A1va[g, t-1] + new_A1vas[g, t-1] - new_A2vas[g, t-1])
      A2va[g, t] <- max(0, A2va[g, t-1] + new_A2vas[g, t-1] - new_A3vas[g, t-1])
      A3va[g, t] <- max(0, A3va[g, t-1] + new_A3vas[g, t-1] - new_RAvas[g, t-1])
      I1va[g, t] <- max(0, I1va[g, t-1] + new_I1vas[g, t-1] - new_I2vas[g, t-1])
      I2va[g, t] <- max(0, I2va[g, t-1] + new_I2vas[g, t-1] - new_I3vas[g, t-1])
      I3va[g, t] <- max(0, I3va[g, t-1] + new_I3vas[g, t-1] - new_RIvas[g, t-1])
      
      
      # Add comps with 2nd dose vaccination Vb
      
      Vb[g, t] <- max(0, Vb[g, t-1] + Vb_eff[g, t]-new_Evbs[g, t-1])
      Evb[g, t] <- max(0, Evb[g, t-1] + new_Evbs[g, t-1] - new_Pvbs[g, t-1])
      Pvb[g, t] <- max(0, Pvb[g, t-1] + new_Pvbs[g, t-1] - new_A1vbs[g, t-1] - new_I1vbs[g, t-1])
      A1vb[g, t] <- max(0, A1vb[g, t-1] + new_A1vbs[g, t-1] - new_A2vbs[g, t-1])
      A2vb[g, t] <- max(0, A2vb[g, t-1] + new_A2vbs[g, t-1] - new_A3vbs[g, t-1])
      A3vb[g, t] <- max(0, A3vb[g, t-1] + new_A3vbs[g, t-1] - new_RAvbs[g, t-1])
      I1vb[g, t] <- max(0, I1vb[g, t-1] + new_I1vbs[g, t-1] - new_I2vbs[g, t-1])
      I2vb[g, t] <- max(0, I2vb[g, t-1] + new_I2vbs[g, t-1] - new_I3vbs[g, t-1])
      I3vb[g, t] <- max(0, I3vb[g, t-1] + new_I3vbs[g, t-1] - new_RIvbs[g, t-1])
      
      
      # Calculate reported cases based on all infectious compartments with reporting lag
      # Calculate reported cases based on new_I1s and new_I1vas only
      step <- function(x) {
        ifelse(x >= 1, 1, 0)
      }
      new_reported_cases[g, t] <- step(t - report_lag) * report_frac *
        (new_I1s[g, max(1, t - report_lag)]+new_I1vas[g, max(1, t - report_lag)]+
           new_I1vbs[g, max(1, t - report_lag)])
      
    }
    
    # Total new reported cases at time t(incidence at t)
    total_new_cases[t] <- sum(new_reported_cases[1:3, t])
    
    total_Cuminc[t] <- sum(Cuminc[1:3, t])
    
    # Update active infections
    active_infected[t] <-
      sum(E[1:3, t]) +
      sum(P[1:3, t]) +
      sum(A1[1:3, t]) + sum(A2[1:3, t]) + sum(A3[1:3, t]) +
      sum(I1[1:3, t]) + sum(I2[1:3, t]) + sum(I3[1:3, t]) +
      sum(Pva[1:3, t]) + sum(Pvb[1:3, t]) +
      sum(A1va[1:3, t]) + sum(A2va[1:3, t]) + sum(A3va[1:3, t]) +
      sum(I1va[1:3, t]) + sum(I2va[1:3, t]) + sum(I3va[1:3, t]) +
      sum(A1vb[1:3, t]) + sum(A2vb[1:3, t]) + sum(A3vb[1:3, t]) +
      sum(I1vb[1:3, t]) + sum(I2vb[1:3, t]) + sum(I3vb[1:3, t])
    
    #Shedding Parameters (Fixed log10 converted outside)
    ##Viral Shedding Calculation (Daily)
    
    # for (g in 1:3) {
    #   
    #   # 1. Shedding from P compartments (pre-symptomatic)
    #   shed_P[g, t] <- shed_alpha * (
    #     P[g, t] +
    #       Pva[g, t] +
    #       Pvb[g, t]
    #   )
    #   
    #   # 2. Shedding from A compartments (asymptomatic)
    #   shed_A[g, t] <- 
    #     shed_A1 * A1[g, t] + shed_A2 * A2[g, t] + shed_A3 * A3[g, t] +
    #     shed_A1a * A1va[g, t] + shed_A2a * A2va[g, t] + shed_A3a * A3va[g, t] +
    #     shed_A1b * A1vb[g, t] + shed_A2b * A2vb[g, t] + shed_A3b * A3vb[g, t]
    #   
    #   # 3. Shedding from I compartments (symptomatic)
    #   shed_I[g, t] <- 
    #     shed_I1 * I1[g, t] + shed_I2 * I2[g, t] + shed_I3 * I3[g, t] +
    #     shed_I1a * I1va[g, t] + shed_I2a * I2va[g, t] + shed_I3a * I3va[g, t] +
    #     shed_I1b * I1vb[g, t] + shed_I2b * I2vb[g, t] + shed_I3b * I3vb[g, t]
    # }
    
    # for (g in 1:3) {
    #   
    #   # 1. Shedding from P compartments (pre-symptomatic)
    #   shed_P[g, t] <- params$shed_alpha * P[g, t] +
    #       params$shed_P1a*Pva[g, t] +
    #     params$shed_P1b*Pvb[g, t]
    #   
    #   # 2. Shedding from A compartments (asymptomatic)
    #   shed_A[g, t] <-
    #     params$shed_A1 * A1[g, t] + params$shed_A2 * A2[g, t] + params$shed_A3 * A3[g, t] +
    #     params$shed_A1a * A1va[g, t] + params$shed_A2a * A2va[g, t] + params$shed_A3a * A3va[g, t] +
    #     params$shed_A1b * A1vb[g, t] + params$shed_A2b * A2vb[g, t] + params$shed_A3b * A3vb[g, t]
    #   
    #   # 3. Shedding from I compartments (symptomatic)
    #   shed_I[g, t] <-
    #     params$shed_I1 * I1[g, t] + params$shed_I2 * I2[g, t] + params$shed_I3 * I3[g, t] +
    #     params$shed_I1a * I1va[g, t] + params$shed_I2a * I2va[g, t] + params$shed_I3a * I3va[g, t] +
    #     params$shed_I1b * I1vb[g, t] + params$shed_I2b * I2vb[g, t] + params$shed_I3b * I3vb[g, t]
    # }
    # 
    
    for (g in 1:3) {
      # Shedding from P compartments (pre-symptomatic)
      shed_P[g, t] <- 
        params$shed_alpha * P[g, t] +
        params$shed_P1a * Pva[g, t] +
        params$shed_P1b * Pvb[g, t]
      
      # Shedding from A compartments (asymptomatic)
      shed_A[g, t] <- 
        params$shed_A1 * A1[g, t] + params$shed_A2 * A2[g, t] + params$shed_A3 * A3[g, t] +
        params$shed_A1a * A1va[g, t] + params$shed_A2a * A2va[g, t] + params$shed_A3a * A3va[g, t] +
        params$shed_A1b * A1vb[g, t] + params$shed_A2b * A2vb[g, t] + params$shed_A3b * A3vb[g, t]
      
      # Shedding from I compartments (symptomatic)
      shed_I[g, t] <- 
        params$shed_I1 * I1[g, t] + params$shed_I2 * I2[g, t] + params$shed_I3 * I3[g, t] +
        params$shed_I1a * I1va[g, t] + params$shed_I2a * I2va[g, t] + params$shed_I3a * I3va[g, t] +
        params$shed_I1b * I1vb[g, t] + params$shed_I2b * I2vb[g, t] + params$shed_I3b * I3vb[g, t]
    }
    
    
    # Total shedding across all groups (optional, if still needed)
    daily_shedding[t] <- sum(shed_P[1:3, t]) + sum(shed_A[1:3, t]) + sum(shed_I[1:3, t])
    
    
    for (g in 1:3) {
      
      # Pre-symptomatic (P)
      P_total[g, t] <- P[g, t] + Pva[g, t] + Pvb[g, t]
      
      # Asymptomatic (A)
      A_total[g, t] <- A1[g, t] + A2[g, t] + A3[g, t] +
        A1va[g, t] + A2va[g, t] + A3va[g, t] +
        A1vb[g, t] + A2vb[g, t] + A3vb[g, t]
      
      # Symptomatic (I)
      I_total[g, t] <- I1[g, t] + I2[g, t] + I3[g, t] +
        I1va[g, t] + I2va[g, t] + I3va[g, t] +
        I1vb[g, t] + I2vb[g, t] + I3vb[g, t]
    }
    
  }
  
  ## Advection-dispersion-delay process
  # Define advection-dispersion-decay using intermediate variables
  
  ## Advection-dispersion-delay process
  sigma <- params$transit_mean * params$transit_cv
  tmax <- ceiling(params$transit_mean + 4 * sigma)
  g <- numeric(tmax)
  for (i in 1:tmax) {
    g[i] <- (1 / (sigma * sqrt(2 * pi))) * exp(-(i - params$transit_mean)^2 / (2 * sigma^2))
  }
  
  contrib <- matrix(0, nrow = tmax, ncol = T)
  delayed_conc <- rep(0, T)
  cp_per_person_mL_all <- rep(0, T)
  log10_conc_all <- rep(0, T)
  
  for (t in 1:T) {
    for (i in 1:tmax) {
      idx <- t - i + 1
      if (idx >= 1) {
        contrib[i, t] <- daily_shedding[idx] * g[i] * exp(-params$mu * i)
      }
    }
    delayed_conc[t] <- sum(contrib[, t])
    cp_per_person_all <- delayed_conc[t] * params$mult / wwtp_population
    cp_per_person_mL_all[t] <- cp_per_person_all * flow_mlalldaily[t]
    log10_conc_all[t] <- log10(cp_per_person_mL_all[t] + 1)  # Avoid log(0)
  }
  
  return(list(
    cases = total_new_cases,
    prevalence = active_infected,
    viral_load = log10_conc_all,
    cumulative_incidence = total_Cuminc,
    FOI=total_lambda,
    P_total=P_total,
    A_total=A_total, 
    I_total=I_total,
    shed_P=shed_P,
    shed_A=shed_A,
    shed_I=shed_I
  ))
}
###We used the median estimates from the above values to generate new initial conditions
E_start <- 14 * Norm_contact_dy;P_start <- 5* Norm_contact_dy;A1_start <- 5* Norm_contact_dy;I1_start <- 5* Norm_contact_dy
A2_start=c(0,0,0);A3_start=c(0,0,0);I2_start=c(0,0,0);I3_start=c(0,0,0);R_start=c(0,0,0);Cuminc_start=c(0,0,0)
S_start <- N - (E_start + P_start + A1_start + I1_start)
#####initials for vaccinated grps va
Va_start=c(0,0,0);Pva_start=c(0,0,0);Eva_start=c(0,0,0);A1va_start=c(0,0,0);A2va_start=c(0,0,0);
A3va_start=c(0,0,0);I1va_start=c(0,0,0);I2va_start=c(0,0,0);I3va_start=c(0,0,0)
#####initials for vaccinated grps vb
Vb_start=c(0,0,0);Pvb_start=c(0,0,0);Evb_start=c(0,0,0);A1vb_start=c(0,0,0);A2vb_start=c(0,0,0)
A3vb_start=c(0,0,0);I1vb_start=c(0,0,0);I2vb_start=c(0,0,0);I3vb_start=c(0,0,0)

###incorporate alpha(vaccination coverage 1st dose)
case_dat2=read_excel("Data/case_data_V2.xlsx",sheet="Mergeddat")
prop_cov = case_dat2$cov_dose1
prop_cov2 = case_dat2$cov_dose2

######now convert proportional vaccinated to a rate
p <- prop_cov
v <- rep(NA,length(p))
v[1] <- -log(1 - p[1])
for(i in 2:length(p)){
  x <- (1 - sum(p[1:i]))*exp(sum(v[1:(i-1)])) 
  v[i] <- -log(x)
}

# Add 30 zeros at the beginning to match the length of case data
alpha_v<- c(rep(0, burn_in_timesteps), v)   ###alpha is the per capita vaccination rate of 1st dose at time t
#alpha<- rep(0, 199)
#alpha_test<- c(rep(0, 9),0.1,rep(0, 189))
####convert 2nd dose coverage into a rate
p2 <- prop_cov2
v2 <- rep(NA,length(p2))
v2[1] <- -log(1 - p2[1])
for(i in 2:length(p2)){
  x2 <- (1 - sum(p2[1:i]))*exp(sum(v2[1:(i-1)])) 
  v2[i] <- -log(x2)
}

# Add 30 zeros at the beginning to match the length of case data
mu_v<- c(rep(0, burn_in_timesteps), v2)   ###mu is the per capita vaccination rate of 2nd dose at time t
#mu_test<- c(rep(0, 9),0.1,rep(0, 189))
#mutest set to zero
#mu<- rep(0, 199)


init_vals <- list(
  S = S_start,
  E = E_start,
  P = P_start,
  A1 = A1_start,
  A2 = A2_start,
  A3 = A3_start,
  I1 = I1_start,
  I2 = I2_start,
  I3 = I3_start,
  R = R_start,
  Cuminc = Cuminc_start,
  
  # First dose vaccinated
  Va = Va_start,
  Eva = Eva_start,
  Pva = Pva_start,
  A1va = A1va_start,
  A2va = A2va_start,
  A3va = A3va_start,
  I1va = I1va_start,
  I2va = I2va_start,
  I3va = I3va_start,
  
  # Second dose vaccinated
  Vb = Vb_start,
  Evb = Evb_start,
  Pvb = Pvb_start,
  A1vb = A1vb_start,
  A2vb = A2vb_start,
  A3vb = A3vb_start,
  I1vb = I1vb_start,
  I2vb = I2vb_start,
  I3vb = I3vb_start,
  
  ######vaccination coverage
  mu_v=mu_v,
  alpha_v=alpha_v
)

remaining_fracs <- c(0.8, 0.6, 0.4, 0.2)

runs <- map(remaining_fracs, function(f) {
  out <- mysimulate_model(params, init_vals, contact, N, T,
                          vax_scenario = "reduced_shedding",
                          v1_frac = f, v2_frac = f)
  out$label <- paste0("Reduced shedding: ", round((1 - f) * 100), "%")
  out
})




sim_baseline <- mysimulate_model(
  params = params,
  init = init_vals,     # your list of initial conditions
  contact = contact,    # contact rates from your `dat`
  N = N,                # population sizes
  T = 169 + burn_in_timesteps,  # full time horizon
  vax_scenario = "baseline"
)

sim_baseline

##############Now we sample multiple posterior draws
set.seed(123)
n_draws <- 200    # e.g., 5000 for full run
G <- 3            # number of groups
T <- 199          # timepoints
draw_indices <- sample(1:nrow(mcmc_matrixallcom), size = n_draws, replace = FALSE)

# Fractions for reduced shedding scenarios
fractions <- c(0.8, 0.6, 0.4, 0.2)

# Storage: baseline + one list per fraction
results <- list(
  baseline = list(
    cases = matrix(NA, n_draws, T),
    vload = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev = matrix(NA, n_draws, T),
    FOI = matrix(NA, n_draws, T),
    shed_P = array(NA, dim = c(n_draws, G, T)),
    shed_A = array(NA, dim = c(n_draws, G, T)),
    shed_I = array(NA, dim = c(n_draws, G, T)),
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  ),
  reduced_shedding = vector("list", length(fractions))
)

# Pre-allocate each reduced_shedding fraction slot
for (f_idx in seq_along(fractions)) {
  results$reduced_shedding[[f_idx]] <- list(
    frac = fractions[f_idx],
    cases = matrix(NA, n_draws, T),
    vload = matrix(NA, n_draws, T),
    cuminc = matrix(NA, n_draws, T),
    prev = matrix(NA, n_draws, T),
    FOI = matrix(NA, n_draws, T),
    shed_P = array(NA, dim = c(n_draws, G, T)),
    shed_A = array(NA, dim = c(n_draws, G, T)),
    shed_I = array(NA, dim = c(n_draws, G, T)),
    P_total = array(NA, dim = c(n_draws, G, T)),
    A_total = array(NA, dim = c(n_draws, G, T)),
    I_total = array(NA, dim = c(n_draws, G, T))
  )
}

# Loop over posterior draws
for (i in seq_along(draw_indices)) {
  draw <- draw_indices[i]
  params <- extract_params(mcmc_matrixallcom[draw, ])
  
  # ---- Baseline ----
  sim_base <- mysimulate_model(params, init_vals, contact, N, T,
                               vax_scenario = "baseline")
  
  results$baseline$cases[i, ]    <- sim_base$cases
  results$baseline$vload[i, ]    <- sim_base$viral_load
  results$baseline$cuminc[i, ]   <- sim_base$cumulative_incidence
  results$baseline$prev[i, ]     <- sim_base$prevalence
  results$baseline$FOI[i, ]      <- sim_base$FOI
  results$baseline$shed_P[i, , ] <- sim_base$shed_P
  results$baseline$shed_A[i, , ] <- sim_base$shed_A
  results$baseline$shed_I[i, , ] <- sim_base$shed_I
  results$baseline$P_total[i, , ] <- sim_base$P_total
  results$baseline$A_total[i, , ] <- sim_base$A_total
  results$baseline$I_total[i, , ] <- sim_base$I_total
  
  # ---- Reduced shedding scenarios ----
  for (f_idx in seq_along(fractions)) {
    f <- fractions[f_idx]
    sim_rs <- mysimulate_model(params, init_vals, contact, N, T,
                               vax_scenario = "reduced_shedding",
                               v1_frac = f, v2_frac = f)
    
    rs <- results$reduced_shedding[[f_idx]]
    rs$cases[i, ]     <- sim_rs$cases
    rs$vload[i, ]     <- sim_rs$viral_load
    rs$cuminc[i, ]    <- sim_rs$cumulative_incidence
    rs$prev[i, ]      <- sim_rs$prevalence
    rs$FOI[i, ]       <- sim_rs$FOI
    rs$shed_P[i, , ]  <- sim_rs$shed_P
    rs$shed_A[i, , ]  <- sim_rs$shed_A
    rs$shed_I[i, , ]  <- sim_rs$shed_I
    rs$P_total[i, , ] <- sim_rs$P_total
    rs$A_total[i, , ] <- sim_rs$A_total
    rs$I_total[i, , ] <- sim_rs$I_total
    
    # Save back into results list
    results$reduced_shedding[[f_idx]] <- rs
  }
}

# Compute median and CrI
summary_quantiles <- function(mat) {
  apply(mat, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}

# Make tidy df for plotting
make_plot_df <- function(summary_matrix, burn_in_timesteps, T, scenario_label, date_vec) {
  tibble(
    time = seq_len(T - burn_in_timesteps),
    median = summary_matrix[2, -(1:burn_in_timesteps)],
    lower  = summary_matrix[1, -(1:burn_in_timesteps)],
    upper  = summary_matrix[3, -(1:burn_in_timesteps)],
    #Date   = date_vec[-(1:burn_in_timesteps)],
    Date   = date_vec,
    scenario = scenario_label
  )
}

# Match predicted values to observed sampling days
match_to_observed_days <- function(pred_matrix, sample_days) {
  pred_matrix[, sample_days, drop = FALSE]
}

make_sampled_vload_df <- function(summary_matrix, obs_days, scenario_label, date_vec) {
  sampled_matrix <- summary_matrix[, obs_days, drop = FALSE]
  tibble(
    time = seq_along(obs_days),
    Date = date_vec[obs_days],
    median = sampled_matrix[2, ],
    lower  = sampled_matrix[1, ],
    upper  = sampled_matrix[3, ],
    scenario = scenario_label
  )
}

# Plotting helper
plot_sim_output <- function(data, ylab, title, log_y = FALSE) {
  ggplot(data, aes(x = Date, y = median, color = scenario, fill = scenario)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    labs(x = "Date", y = ylab, title = title, color = "Scenario", fill = "Scenario") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_y_continuous(trans = if (log_y) "log10" else "identity")
}

cases_summary_baseline  <- summary_quantiles(results$baseline$cases)
vload_summary_baseline  <- summary_quantiles(results$baseline$vload)
prev_summary_baseline   <- summary_quantiles(results$baseline$prev)
cuminc_summary_baseline <- summary_quantiles(results$baseline$cuminc)

cases_df <- make_plot_df(cases_summary_baseline, burn_in_timesteps, T, "Baseline", case_dat$Date)
vload_df <- make_plot_df(vload_summary_baseline, burn_in_timesteps, T, "Baseline", case_dat$Date)
prev_df  <- make_plot_df(prev_summary_baseline, burn_in_timesteps, T, "Baseline", case_dat$Date)
cuminc_df <- make_plot_df(cuminc_summary_baseline, burn_in_timesteps, T, "Baseline", case_dat$Date)

vload_sampled_df <- make_sampled_vload_df(vload_summary_baseline, obs_days, "Baseline", case_dat$Date)

# --- Reduced shedding fractions ---
for (rs in results$reduced_shedding) {
  frac_label <- paste0("Reduced shedding: ", round((1 - rs$frac) * 100), "%")
  
  cases_summary  <- summary_quantiles(rs$cases)
  vload_summary  <- summary_quantiles(rs$vload)
  prev_summary   <- summary_quantiles(rs$prev)
  cuminc_summary <- summary_quantiles(rs$cuminc)
  
  cases_df  <- bind_rows(cases_df,  make_plot_df(cases_summary, burn_in_timesteps, T, frac_label, case_dat$Date))
  vload_df  <- bind_rows(vload_df,  make_plot_df(vload_summary, burn_in_timesteps, T, frac_label, case_dat$Date))
  prev_df   <- bind_rows(prev_df,   make_plot_df(prev_summary, burn_in_timesteps, T, frac_label, case_dat$Date))
  cuminc_df <- bind_rows(cuminc_df, make_plot_df(cuminc_summary, burn_in_timesteps, T, frac_label, case_dat$Date))
  
  vload_sampled_df <- bind_rows(vload_sampled_df, make_sampled_vload_df(vload_summary, obs_days, frac_label, case_dat$Date))
}




















































































saveRDS(all_sim_resultse$baseline, "D:/Mpoxoutputb/baseline_resultse.rds")
saveRDS(all_sim_resultse$no_vax, "D:/Mpoxoutputb/novax_resultse.rds")
saveRDS(all_sim_resultse$reduced_shedding, "D:/Mpoxoutputb/reducedshed_resultse.rds")


baseline <- readRDS("D:/Mpoxoutputb/baseline_resultse.rds")
novax <- readRDS("D:/Mpoxoutputb/novax_resultse.rds")
reducedshed <- readRDS("D:/Mpoxoutputb/reducedshed_resultse.rds")

##########we now compute credible intervals############
summary_quantiles <- function(mat) {
  apply(mat, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}

# Summarize both scenarios
cases_summary_novax        <- summary_quantiles(novax$cases)
cases_summary_reduced      <- summary_quantiles(reducedshed$cases)
cases_summary_baseline     <- summary_quantiles(baseline$cases)

vload_summary_novax        <- summary_quantiles(novax$vload)
vload_summary_reduced      <- summary_quantiles(reducedshed$vload)
vload_summary_baseline     <- summary_quantiles(baseline$vload)

cuminc_summary_novax       <- summary_quantiles(novax$cuminc)
cuminc_summary_reduced     <- summary_quantiles(reducedshed$cuminc)
cuminc_summary_baseline    <- summary_quantiles(baseline$cuminc)

prev_summary_novax         <- summary_quantiles(novax$prev)
prev_summary_reduced       <- summary_quantiles(reducedshed$prev)
prev_summary_baseline      <- summary_quantiles(baseline$prev)

FOI_summary_novax          <- summary_quantiles(novax$FOI)
FOI_summary_reduced        <- summary_quantiles(reducedshed$FOI)
FOI_summary_baseline       <- summary_quantiles(baseline$FOI)

### Example: plot-ready df for cases
make_plot_df <- function(summary_matrix, burn_in_timesteps, T, scenario_label) {
  df <- tibble(
    time = 1:T,
    median = summary_matrix[2, ],
    lower = summary_matrix[1, ],
    upper = summary_matrix[3, ],
    scenario = scenario_label
  ) %>%
    slice(-(1:burn_in_timesteps)) %>%
    mutate(time = 1:n(),
           Date=case_dat$Date)  # Reset time to start from 1
  return(df)
}

match_to_observed_days <- function(pred_matrix, sample_days) {
  # pred_matrix: [iterations × total_days], e.g. output of your JAGS model
  # sample_days: vector of day indices with observed WW (e.g., c(3, 5, 9, 14, ...))
  
  # Check that all sample days are within bounds
  if (any(sample_days > ncol(pred_matrix))) {
    stop("Some sampling days exceed number of model timepoints.")
  }
  
  matched_matrix <- pred_matrix[, sample_days, drop = FALSE]
  return(matched_matrix)  # [iterations × number_of_observations]
}

make_sampled_vload_df <- function(summary_matrix, obs_days, scenario_label) {
  sampled_matrix <- match_to_observed_days(summary_matrix, obs_days)
  
  df <- tibble(
    time = 1:length(obs_days),
    Date=ww_dat$date,
    median = sampled_matrix[2, ],
    lower = sampled_matrix[1, ],
    upper = sampled_matrix[3, ],
    scenario = scenario_label
  )
  
  return(df)
}

cases_df_novax <- make_plot_df(cases_summary_novax, burn_in_timesteps, T, "No vaccination")
cases_df_reduced <- make_plot_df(cases_summary_reduced, burn_in_timesteps, T, "Reduced shedding")
cases_df_baseline <- make_plot_df(cases_summary_baseline, burn_in_timesteps, T, "Baseline")

vload_df_novax <- make_plot_df(vload_summary_novax, burn_in_timesteps, T, "No vaccination")
vload_df_reduced <- make_plot_df(vload_summary_reduced, burn_in_timesteps, T, "Reduced shedding")
vload_df_baseline <- make_plot_df(vload_summary_baseline, burn_in_timesteps, T, "Baseline")

obs_days <- ww_dat$Day + burn_in_timesteps  # double check this is correct if burn-in = 30
vload_sampled_novax <- make_sampled_vload_df(vload_summary_novax, obs_days, "No vaccination")
vload_sampled_reduced <- make_sampled_vload_df(vload_summary_reduced, obs_days, "Reduced shedding")
vload_sampled_baseline <- make_sampled_vload_df(vload_summary_baseline, obs_days, "Baseline")

prev_df_novax <- make_plot_df(prev_summary_novax, burn_in_timesteps, T, "No vaccination")
prev_df_reduced <- make_plot_df(prev_summary_reduced, burn_in_timesteps, T, "Reduced shedding")
prev_df_baseline <- make_plot_df(prev_summary_baseline, burn_in_timesteps, T, "Baseline")

cuminc_df_novax <- make_plot_df(cuminc_summary_novax, burn_in_timesteps, T, "No vaccination")
cuminc_df_reduced <- make_plot_df(cuminc_summary_reduced, burn_in_timesteps, T, "Reduced shedding")
cuminc_df_baseline <- make_plot_df(cuminc_summary_baseline, burn_in_timesteps, T, "Baseline")

FOI_df_novax <- make_plot_df(FOI_summary_novax, burn_in_timesteps, T, "No vaccination")
FOI_df_reduced <- make_plot_df(FOI_summary_reduced, burn_in_timesteps, T, "Reduced shedding")
FOI_df_baseline <- make_plot_df(FOI_summary_baseline, burn_in_timesteps, T, "Baseline")


#########Now include the same output for the baseline######################3
# Combine for faceted/overlaid plotting
#plot_cases_df <- bind_rows(cases_df_novax, cases_df_reduced,cases_df_baseline)
#plot_vload_df <- bind_rows(vload_df_novax, vload_df_reduced,vload_df_baseline)
#plot_prev_df <- bind_rows(prev_df_novax, prev_df_reduced,prev_df_baseline )
#plot_cuminc_df <- bind_rows(cuminc_df_novax, cuminc_df_reduced,cuminc_df_baseline)
#plot_vload_sampled <- bind_rows(vload_sampled_novax, vload_sampled_reduced,vload_sampled_baseline)


#########Now include the same output for the baseline######################3
# Combine for faceted/overlaid plotting
plot_cases_df <- bind_rows(cases_df_novax, cases_df_baseline)
plot_vload_df <- bind_rows(vload_df_novax, vload_df_baseline)
plot_prev_df <- bind_rows(prev_df_novax, prev_df_baseline )
plot_cuminc_df <- bind_rows(cuminc_df_novax, cuminc_df_baseline)
plot_vload_sampled <- bind_rows(vload_sampled_novax,vload_sampled_baseline)

###########Now we generate plots for these plots
plot_sim_output <- function(data, ylab, title, log_y = FALSE) {
  ggplot(data, aes(x = Date, y = median, color = scenario, fill = scenario)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    labs(
      x = "Date",
      y = ylab,
      title = title,
      color = "Scenario",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(trans = if (log_y) "log10" else "identity")
}
p_cases <- plot_sim_output(plot_cases_df, ylab = "Reported cases", title = "Simulated Cases")
p_vload <- plot_sim_output(plot_vload_df, ylab = "WW viral load", title = "Simulated WW Viral Load", log_y = TRUE)
p_vload_sampled <- plot_sim_output(plot_vload_sampled, ylab = "WW viral load", title = "WW Load on Sampling Days", log_y = TRUE)
p_prev <- plot_sim_output(plot_prev_df, ylab = "Prevalence", title = "Simulated Prevalence")
p_cuminc <- plot_sim_output(plot_cuminc_df, ylab = "Cumulative incidence", title = "Simulated Cumulative Incidence")


bb=p_cases+p_vload_sampled+plot_layout(guides = "collect") & theme(legend.position = "bottom")
bb

ggsave(
  filename = "D:/Mpox25output/Figures/simfit.tiff",
  plot = bb,  # use your actual plot variable here
  width = 10,
  height = 5,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


#########Contribution to viral load by sexual activity group#########
# Total shedding per group: sum P, A, I
summarise_relative_shedding_array <- function(shed_P, shed_A, shed_I, scenario_name = "Scenario") {
  # Total shedding per group [draw, group, time]
  shed_total <- shed_P + shed_A + shed_I
  
  # Total shedding across groups [draw, time]
  shed_all <- apply(shed_total, c(1, 3), sum, na.rm = TRUE)
  
  # Relative shedding: [draw, group, time]
  rel_shed <- sweep(shed_total, c(1, 3), shed_all, "/")
  rel_shed[is.na(rel_shed)] <- 0
  
  # Add time labels before melting
  dimnames(rel_shed)[[3]] <- as.character(1:dim(rel_shed)[3])
  
  # Reshape to long format
  rel_df <- as.data.frame.table(rel_shed, responseName = "prop") %>%
    dplyr::rename(draw = Var1, group = Var2, time = Var3) %>%
    dplyr::mutate(
      draw = as.integer(draw),
      group = factor(as.numeric(group), levels = 1:3,
                     labels = c("Group 1 (low activity)",
                                "Group 2 (medium activity)",
                                "Group 3 (high activity)")),
      time = as.integer(time),
      scenario = scenario_name
    )
  
  # Summarize across draws
  summary_df <- rel_df %>%
    dplyr::group_by(time, group, scenario) %>%
    dplyr::summarise(
      median = median(prop, na.rm = TRUE),
      lower = quantile(prop, 0.025, na.rm = TRUE),
      upper = quantile(prop, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_df)
}

# Summarize relative shedding contributions by compartment for each scenario
baseline_df <- summarise_relative_shedding_array(
  baseline$shed_P, baseline$shed_A, baseline$shed_I, "Baseline"
)

novax_df <- summarise_relative_shedding_array(
  novax$shed_P, novax$shed_A, novax$shed_I, "No vaccination"
)

shedred_df <- summarise_relative_shedding_array(
  reducedshed$shed_P, reducedshed$shed_A, reducedshed$shed_I, "Reduced shedding"
)

#relg_all <- dplyr::bind_rows(baseline_df, novax_df, shedred_df)
relg_all <- dplyr::bind_rows(baseline_df, novax_df)

relg_all_trimmed <- relg_all %>%
  group_by(group, scenario) %>%
  arrange(time) %>%
  slice(-(1:30)) %>%  # Remove first 30 rows per group-scenario combo
  mutate(time = row_number(),
         Date=case_dat$Date) %>% 
  ungroup()

my_colors <- c(
  "Group 1 (low activity)" = "#F8766D",
  "Group 2 (medium activity)" = "#00BA38",
  "Group 3 (high activity)" = "#619CFF"
)

plot_shedding_contributions <- function(shed_summary_df) {
  ggplot(shed_summary_df, aes(x = Date, y = median, fill = group, color = group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    facet_wrap(~ scenario, ncol = 1) +
    labs(
      #title = "Relative Contribution to Shedding by Sexual Activity Group",
      x = "Time (days)",
      y = "Proportion of Total Viral Shedding"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# Plot
shedplot=plot_shedding_contributions(relg_all_trimmed)
shedplot

###########Number of infectious people per group#############
# Function to summarize total infectious individuals per group
summarise_total_infectious <- function(P_total, A_total, I_total, burn_in = 30) {
  # Total infectious = P + A + I
  infected <- P_total + A_total + I_total  # shape: draws x groups x time
  
  # Subset time after burn-in
  infected <- infected[, , (burn_in + 1):dim(infected)[3]]
  
  # Create summary per group
  summaries <- lapply(1:dim(infected)[2], function(g) {
    mat <- infected[, g, ]  # matrix: draws x time
    data.frame(
      time = 1:dim(mat)[2],  # renumbered to start from 1 after burn-in
      median = apply(mat, 2, median, na.rm = TRUE),
      lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      group = paste0("Group ", g,
                     ifelse(g == 1, " (low activity)",
                            ifelse(g == 2, " (medium activity)", " (high activity)")))
    )
  })
  
  dplyr::bind_rows(summaries)
}


infect_summary_baseline <- summarise_total_infectious(
  baseline$P_total, baseline$A_total, baseline$I_total
) %>%
  dplyr::mutate(scenario = "Baseline")

infect_summary_novax <- summarise_total_infectious(
  novax$P_total, novax$A_total, novax$I_total
) %>%
  dplyr::mutate(scenario = "No vaccination")

# infect_summary_reduced <- summarise_total_infectious(
#   reducedshed$P_total, reducedshed$A_total, reducedshed$I_total
# ) %>%
#   dplyr::mutate(scenario = "Reduced shedding")


# Combine all
infect_all <- bind_rows(infect_summary_baseline, infect_summary_novax)

Infectionsplot <- ggplot(infect_all, aes(x = time, y = median, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  facet_wrap(~scenario, ncol = 1) +
  labs(
    x = "Time (days)",
    y = "Number of Infectious Individuals",
    #title = "Number of Infectious Individuals by Sexual Activity Group"
  ) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ShedInf_plot <- plot_shedding_contributions(relg_all_trimmed) | Infectionsplot
pp=ShedInf_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "D:/Mpox25output/Figures/simshed.tiff",
  plot = pp,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

 #################Contribution to viral load by compartment##########################
# Combine P, A, I compartments and remove burn-in
summarise_rel_shedding <- function(shed_P, shed_A, shed_I, scenario_name = "Scenario", burn_in = 30, date = NULL) {
  # Check input shape
  stopifnot(length(dim(shed_P)) == 3)
  
  # Collapse across groups: [draw, time]
  shedP_total <- apply(shed_P, c(1, 3), sum)
  shedA_total <- apply(shed_A, c(1, 3), sum)
  shedI_total <- apply(shed_I, c(1, 3), sum)
  
  # Remove burn-in
  shedP_total <- shedP_total[, (burn_in + 1):ncol(shedP_total)]
  shedA_total <- shedA_total[, (burn_in + 1):ncol(shedA_total)]
  shedI_total <- shedI_total[, (burn_in + 1):ncol(shedI_total)]
  
  # Total shedding across all compartments
  shed_total <- shedP_total + shedA_total + shedI_total
  
  # Relative contributions
  relcom_P <- shedP_total / shed_total
  relcom_A <- shedA_total / shed_total
  relcom_I <- shedI_total / shed_total
  
  # Summary helper
  summary_rel <- function(mat, source_label) {
    df <- data.frame(
      time = seq_len(ncol(mat)),  # Start from 1 after burn-in
      median = apply(mat, 2, median, na.rm = TRUE),
      lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      source = source_label
    )
    if (!is.null(date)) {
      if (length(date) != ncol(mat)) stop("Length of `date` must match number of timepoints after burn-in")
      df$date <- date
    }
    return(df)
  }
  
  # Combine summaries
  summary_df <- dplyr::bind_rows(
    summary_rel(relcom_P, "Pre-symptomatic (P)"),
    summary_rel(relcom_A, "Asymptomatic (A)"),
    summary_rel(relcom_I, "Symptomatic (I)")
  ) %>%
    dplyr::mutate(scenario = scenario_name)
  
  return(summary_df)
}

# Relative shedding by clinical compartment
rel_shedding_baseline <- summarise_rel_shedding(
  baseline$shed_P, baseline$shed_A, baseline$shed_I, "Baseline",date = case_dat$Date
)

rel_shedding_novax <- summarise_rel_shedding(
  novax$shed_P, novax$shed_A, novax$shed_I, "No vaccination",date = case_dat$Date
)

rel_shedding_reduced <- summarise_rel_shedding(
  reducedshed$shed_P, reducedshed$shed_A, reducedshed$shed_I, "Reduced shedding",date = case_dat$Date
)

summarise_abs_infected <- function(P_mat, A_mat, I_mat, scenario_name = "Scenario", burn_in = 30, date = NULL) {
  stopifnot(length(dim(P_mat)) == 3)
  
  # Collapse across groups: [draw, time]
  P_total <- apply(P_mat, c(1, 3), sum)
  A_total <- apply(A_mat, c(1, 3), sum)
  I_total <- apply(I_mat, c(1, 3), sum)
  
  # Remove burn-in
  P_total <- P_total[, (burn_in + 1):ncol(P_total)]
  A_total <- A_total[, (burn_in + 1):ncol(A_total)]
  I_total <- I_total[, (burn_in + 1):ncol(I_total)]
  
  # Summary helper
  summarise_draws <- function(mat, comp_label) {
    df <- data.frame(
      time = seq_len(ncol(mat)),
      median = apply(mat, 2, median, na.rm = TRUE),
      lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      compartment = comp_label
    )
    if (!is.null(date)) {
      if (length(date) != ncol(mat)) stop("Length of `date` must match number of timepoints after burn-in")
      df$date <- date
    }
    return(df)
  }
  
  # Combine summaries
  summary_df <- dplyr::bind_rows(
    summarise_draws(P_total, "Pre-symptomatic (P)"),
    summarise_draws(A_total, "Asymptomatic (A)"),
    summarise_draws(I_total, "Symptomatic (I)")
  ) %>%
    dplyr::mutate(scenario = scenario_name)
  
  return(summary_df)
}

# Total number of infected individuals by clinical compartment
infected_baseline <- summarise_abs_infected(
  baseline$P_total, baseline$A_total, baseline$I_total, "Baseline",date = case_dat$Date
)

infected_novax <- summarise_abs_infected(
  novax$P_total, novax$A_total, novax$I_total, "No vaccination",date = case_dat$Date
)

infected_reduced <- summarise_abs_infected(
  reducedshed$P_total, reducedshed$A_total, reducedshed$I_total, "Reduced shedding",date = case_dat$Date
)

# Combine all scenarios
rel_shedding_all <- dplyr::bind_rows(
  rel_shedding_baseline, rel_shedding_novax)

infected_all <- dplyr::bind_rows(
  infected_baseline, infected_novax)

plotrelcom=ggplot(rel_shedding_all, aes(x = date, y = median*100, fill = source)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), alpha = 0.2, color = NA) +
  geom_line(aes(color = source), size = 1) +
  facet_wrap(~scenario)+
  labs(
    x = "Date", y = "Proportion of total viral load (%)",
    title = "Relative contribution to viral shedding by infection stage"
  ) +
  scale_fill_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                               "Asymptomatic (A)" = "orange",
                               "Symptomatic (I)" = "red")) +
  scale_color_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                                "Asymptomatic (A)" = "orange",
                                "Symptomatic (I)" = "red")) +
  theme_minimal(base_size = 10) +
  theme(legend.title = element_blank())

plotinfected <- ggplot(infected_all, aes(x = date, y = median, fill = compartment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(color = compartment), size = 1) +
  facet_wrap(~scenario) +
  labs(
    x = "Date",
    y = "Total number of infected individuals",
    title = "Number of infected individuals by clinical stage with 95% CrI"
  ) +
  scale_fill_manual(values = c(
    "Pre-symptomatic (P)" = "skyblue",
    "Asymptomatic (A)" = "orange",
    "Symptomatic (I)" = "red"
  )) +
  scale_color_manual(values = c(
    "Pre-symptomatic (P)" = "skyblue",
    "Asymptomatic (A)" = "orange",
    "Symptomatic (I)" = "red"
  )) +
  theme_minimal(base_size = 10) +
  theme(legend.title = element_blank())

combined_plot <- plotrelcom / plotinfected + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(
  filename = "D:/Mpox25output/Figures/simshedcomp.tiff",
  plot = combined_plot,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

############sliding window correlation plot####################
# 1. Drop burn-in if needed
burn_in <- 30
cases_postburn <- cases_all[, (burn_in + 1):ncol(cases_all)]
vload_postburn <- vload_all[, (burn_in + 1):ncol(vload_all)]

# 2. Define sliding window sizes
window_sizes <- seq(7, ncol(cases_postburn), by = 7)
n_draws <- nrow(cases_postburn)

# 3. Compute correlation per draw and window
correlation_list_sim <- vector("list", length(window_sizes))
names(correlation_list_sim) <- window_sizes

for (w in seq_along(window_sizes)) {
  size <- window_sizes[w]
  cor_per_draw <- numeric(n_draws)
  
  for (d in 1:n_draws) {
    ww_sub <- vload_postburn[d, 1:size]
    case_sub <- cases_postburn[d, 1:size]
    
    if (sd(ww_sub) > 0 && sd(case_sub) > 0) {
      cor_per_draw[d] <- cor(ww_sub, case_sub)
    } else {
      cor_per_draw[d] <- NA
    }
  }
  correlation_list_sim[[w]] <- cor_per_draw
}

# 4. Summarize across draws
summary_df_sim <- do.call(rbind, lapply(seq_along(correlation_list_sim), function(i) {
  vals <- na.omit(correlation_list_sim[[i]])
  data.frame(
    window_size = window_sizes[i],
    frac_time = window_sizes[i] / max(window_sizes),
    median = median(vals),
    lower = quantile(vals, 0.025),
    upper = quantile(vals, 0.975)
  )
}))

# 6. Estimate peak transmission time
# Use the median (50%) row from your FOI_summary matrix
FOI_median_sim <- FOI_summary["50%", (burn_in + 1):ncol(FOI_summary)]
# Find time of peak FOI
transmission_peak_day_sim <- which.max(FOI_median_sim)
# Normalize by number of time steps to get fractional time
total_days_sim <- length(FOI_median_sim)
transmission_peak_frac_sim <- transmission_peak_day_sim / total_days_sim


##############Time to peak##############
cases_novac <- tibble(
  time = 1:T,
  median = cases_summary[2, ],
  lower = cases_summary[1, ],
  upper = cases_summary[3, ]
)


##########summarise viral load
viral_novac <- tibble(
  time = 1:T,
  median = vload_summary[2, ],
  lower = vload_summary[1, ],
  upper = vload_summary[3, ]
)

cases_novacb <- cases_novac %>% slice(-(1:burn_in_timesteps))
viral_novacb <- viral_novac %>% slice(-(1:burn_in_timesteps))


# WW dataframe
ww_dfvac <- data.frame(
  time = 1:T_postburn,
  value = viral_novacb$median,
  lower = viral_novacb$lower,
  upper = viral_novacb$upper,
  series = "WW Viral Load"
)

# Cases dataframe
cases_dfvac <- data.frame(
  time = 1:T_postburn,
  value = cases_novacb$median,
  lower = cases_novacb$lower,
  upper = cases_novacb$upper,
  series = "Cases"
)

peak_day_cases_novax<- which.max(cases_dfvac$value)
peak_day_ww_novax<- which.max(ww_dfvac$value)
lag_novax <- peak_day_ww_novax - peak_day_cases_novax
max_y_novax <- max(c(cases_novac$upper, viral_novac$upper), na.rm = TRUE)

df_novac <- rbind(cases_dfvac,ww_dfvac)

###############generate summaries for baseline scenario########################
##########correlation plot##########
ww_mat <- mcmc_matrixallcom[, grep("^log10_conc_all\\[", colnames(mcmc_matrixallcom))]
case_mat <- mcmc_matrixallcom[, grep("^mu_nb\\[", colnames(mcmc_matrixallcom))]

# 2. Define window sizes to compute sliding correlation
window_sizes <- seq(7, ncol(ww_mat), by = 7)
n_draws <- nrow(ww_mat)

# 3. Initialize a list to hold results from each draw
correlation_list <- vector("list", length(window_sizes))
names(correlation_list) <- window_sizes

# 4. Loop through each draw and each window to compute correlations
for (w in seq_along(window_sizes)) {
  size <- window_sizes[w]
  cor_per_draw <- numeric(n_draws)
  
  for (d in 1:n_draws) {
    ww_sub <- ww_mat[d, 1:size]
    case_sub <- case_mat[d, 1:size]
    
    # Avoid correlation with no variance
    if (sd(ww_sub) > 0 && sd(case_sub) > 0) {
      cor_per_draw[d] <- cor(ww_sub, case_sub)
    } else {
      cor_per_draw[d] <- NA
    }
  }
  correlation_list[[w]] <- cor_per_draw
}

# 5. Summarize across draws: median and 95% credible intervals
summary_df <- do.call(rbind, lapply(seq_along(correlation_list), function(i) {
  vals <- na.omit(correlation_list[[i]])
  data.frame(
    window_size = window_sizes[i],
    frac_time = window_sizes[i] / max(window_sizes),
    median = median(vals),
    lower = quantile(vals, 0.025),
    upper = quantile(vals, 0.975)
  )
}))

# 6. Optional: Mark the peak transmission
burn_in_timesteps <- 30
FOI <- as.data.frame(posterior_summarycom[grep("^total_lambda\\[", rownames(posterior_summarycom)), "median"])
FOI_final <- FOI[(burn_in_timesteps + 1):nrow(FOI), ]
transmission_peak_day <- which.max(FOI_final)
total_days <- ncol(ww_mat)
transmission_peak_frac <- transmission_peak_day / total_days

#############peak timing code########################
# Median and credible intervals (95%) per timepoint
ww_summaryb <- apply(ww_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
cases_summaryb <- apply(case_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))

# Time vector
time <- 1:ncol(ww_mat)

# WW dataframe
ww_df <- data.frame(
  time = time,
  value = ww_summaryb["50%", ],
  lower = ww_summaryb["2.5%", ],
  upper = ww_summaryb["97.5%", ],
  series = "WW Viral Load"
)

# Cases dataframe
cases_df <- data.frame(
  time = time,
  value = cases_summaryb["50%", ],
  lower = cases_summaryb["2.5%", ],
  upper = cases_summaryb["97.5%", ],
  series = "Cases"
)

# Combine both
df_long <- rbind(ww_df, cases_df)

peak_day_cases_baseline <- which.max(cases_df$value)
peak_day_ww_baseline <- which.max(ww_df$value)
lag_baseline <- peak_day_ww_baseline - peak_day_cases_baseline
max_y_baseline <- max(c(cases_df$upper, ww_df$upper), na.rm = TRUE)

###########function to plot all the figures side by side##########
plot_comparison <- function(df_baseline, df_counterfactual, y_label,title) {
   df_baseline$scenario <- "Baseline"
   df_counterfactual$scenario <- "No vaccination"
   df_combined <- bind_rows(df_baseline, df_counterfactual)
   
   ggplot(df_combined, aes(x = time, group = scenario)) +
     geom_line(aes(y = median_fit, color = scenario), size = 1) +
     geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = scenario), alpha = 0.2) +
     scale_color_manual(values = c("Baseline" = "blue4", "No vaccination" = "firebrick")) +
     scale_fill_manual(values = c("Baseline" = "skyblue", "No vaccination" = "salmon")) +
     labs(
       x = "Time (days)",
       y = y_label,
       title = title
     ) +
     theme_minimal(base_size = 14) +
     theme(
       legend.title = element_blank(),
       plot.title = element_text(face = "bold", hjust = 0.5)
     )
 }
 
#################correlation plots###############
# Add scenario labels
summary_df$scenario <- "Baseline"
summary_df_sim$scenario <- "No vaccination"

# Combine the data
summary_df_combined <- bind_rows(summary_df, summary_df_sim)

rownames(summary_df_combined)<-NULL

# Plot
finalplt_combined <- ggplot(summary_df_combined, aes(x = frac_time, y = median, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  #geom_vline(xintercept = transmission_peak_frac, linetype = "dashed", color = "blue") +
  #geom_text(aes(x = transmission_peak_frac, y = 0.945), label = "Baseline peak", 
            #color = "blue", size = 3.5, hjust = -0.1) +
  #geom_vline(xintercept = transmission_peak_frac_sim, linetype = "dashed", color = "red") +
  #geom_text(aes(x = transmission_peak_frac_sim, y = 0.88), label = "No vax", 
            #color = "red", size = 3.5, hjust = -0.1) +
  labs(
    title = "Sliding window correlation between predicted cases and viral load",
    x = "Fraction of time series used",
    y = expression("Pearson correlation (" * log[10]~"viral load vs. predicted cases)"),
    caption = "Shaded region indicates 95% credible interval"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Baseline" = "blue", "No vaccination" = "red")) +
  scale_fill_manual(values = c("Baseline" = "blue", "No vaccination" = "red")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

finalplt_combined 


#######difference in peak timing
df_long$scenario <- "Baseline"
# No Vax
df_novac$scenario <- "No vaccination"

# Combine into one
df_combined <- rbind(df_long, df_novac)

peak_df <- data.frame(
  scenario = c("Baseline", "No vaccination"),
  peak_cases = c(peak_day_cases_baseline, peak_day_cases_novax),
  peak_ww = c(peak_day_ww_baseline, peak_day_ww_novax),
  lag_days = c(lag_baseline, lag_novax),
  max_y_cases = c(3.1, 4.0),
  max_y_ww = c(5.6, 5.7)
)


p_cases <- ggplot(subset(df_combined, series == "Cases"), aes(x = time, y = value, fill = scenario, color = scenario)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  facet_wrap(~scenario, ncol = 2) +
  geom_vline(data = peak_df, aes(xintercept = peak_cases, color = scenario), linetype = "dashed") +
  geom_text(data = peak_df, aes(x = peak_cases + 2, y = max_y_cases * 0.85,
                                label = paste0("Lag = ", lag_days, " days")),
            inherit.aes = FALSE, hjust = 0, size = 4) +
  scale_color_manual(values = c("Baseline" = "red", "No vaccination" = "darkred")) +
  scale_fill_manual(values = c("Baseline" = "red", "No vaccination" = "darkred")) +
  labs(y = "Predicted cases", x = "Time (days)") +
  theme_minimal()

# Plot for WW Viral Load (log scale)
p_ww <- ggplot(subset(df_combined, series == "WW Viral Load"), aes(x = time, y = value, fill = scenario, color = scenario)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  facet_wrap(~scenario, ncol = 2) +
  geom_vline(data = peak_df, aes(xintercept = peak_ww, color = scenario), linetype = "dashed") +
  geom_text(data = peak_df, aes(x = peak_ww + 2, y = max_y_ww * 0.85,
                                label = paste0("Lag = ", lag_days, " days")),
            inherit.aes = FALSE, hjust = 0, size = 4) +
  scale_y_log10() +
  scale_color_manual(values = c("Baseline" = "blue", "No vaccination" = "darkblue")) +
  scale_fill_manual(values = c("Baseline" = "blue", "No vaccination" = "darkblue")) +
  labs(y = "Predicted viral load (log)", x = "Time (days)") +
  theme_minimal()

# Combine using patchwork
p_cases / p_ww + plot_layout(heights = c(1, 1))














































  
  
  
  
  
  