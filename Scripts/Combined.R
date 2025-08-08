#rm(list=ls())
library(R2jags)
library(runjags)
library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)
options(scipen=999)

###how to interpret beta in my model
####In this I am now fitting combined case and WW datasets

textstring <- "
model {
  
  mu= 0.18           # Viral decay, fix to 0.18
  
  ##priors that worked for final case model 
   # log_mult ~ dnorm(log(1.5e-9), 400) ##final case
   # mult <- exp(log_mult)
   # transit_time_mean ~ dnorm(2.5, 25) T(1.5, 4)
   # transit_time_cv ~ dnorm(0.3, 36) T(0.2, 0.8)
   # tau_ww ~ dgamma(50, 60)  # mean ≈ 0.83, SD ≈ 0.12#newcombmod8

   #######To reun using similar priors for WW
   ####this was for modelb
   #log_mult ~ dnorm(log(7e-9), 200) #reduce sd of scaling factor
   #mult <- exp(log_mult)
   #tau_ww ~ dgamma(10, 1)
   #transit_time_mean ~ dnorm(2.5, 4) T(1, 5)     # mean = 2.5 days, SD ≈ 0.5
   #transit_time_cv ~ dnorm(0.3, 36) T(0.15, 0.6) #wwmod12
   
   #########this is for model c
  # log_mult ~ dnorm(log(3e-9), 300)
  # mult <- exp(log_mult)
  # tau_ww ~ dgamma(40, 48)
  # transit_time_mean ~ dnorm(2.5, 9) T(1.3, 4.5)
  # transit_time_cv ~ dnorm(0.3, 36) T(0.2, 0.6)
   
   #########this is for model h
   
   # log_mult ~ dnorm(log(3e-9), 40)
   # mult <- exp(log_mult)
   # tau_ww ~ dgamma(40, 48)
   # transit_time_mean ~ dnorm(2.5, 1) T(1, 5)
   # transit_time_cv ~ dnorm(0.3, 3) T(0.1, 1)
   # 
    #########priors for original covid paper
    # log_mult ~ dunif(log(0.001), log(0.005))
    # mult <- exp(log_mult)
    # tau_ww ~ dgamma(40, 48)
    # transit_time_mean ~ dunif(1, 5)
   
    ##New wider priors
  log_mult ~ dnorm(log(3e-9), 2.5) T(log(1e-9), log(1e-8))
  mult <- exp(log_mult)
  tau_ww ~ dgamma(40, 48)
  transit_time_mean ~ dnorm(2.5, 0.25) T(0.1, 10)
  transit_time_cv ~ dnorm(0.3, 3) T(0.1, 1)
   

   # Estimate both mean and CV
   ###other parameters
  
   beta~dnorm(0.8, 100) T(0,1)
   kappa ~ dbeta(40, 2)  ### Mean ~0.95, 95% CI ≈ [0.85, 0.995]
   phi ~ dgamma(2, 0.5)  # Mean = 4, SD = ~2.8


  E0 ~ dpois(5)    #  A Poisson prior with mean 5
  P0 ~ dpois(1)  # Total pre-symptomatic individuals
  A10 ~ dpois(1) # Total asymptomatic individuals
  I10 ~ dpois(1) # Total symptomatic individuals

 
  
  #delta_inv<-6.26 ## duration in E comp 6.26(5.55-6.97)
  delta_inv~dgamma(300.65298,47.98256) ## duration in E comp 6.26(5.55-6.97)
  delta <-1/delta_inv   ##daily rate of transition from E to P
  deltava <-delta   ##daily rate of transition from Eva to Pva
  deltavb <-delta  ##daily rate of transition from Evb to Pvb
  
  epsilon_inv<-2     ##duration in P comp 2days
  epsilon <- 1/epsilon_inv ##daily rate of transition from P to A and I
  epsilonva <- epsilon ##daily rate of transition from Pva to Ava and Iva
  epsilonvb <- epsilon ##daily rate of transition from Pvb to Avb and Ivb
  
  m~dbeta(49.09,94.88)   ##proportion of Asymptomatic individuals 34%(26.65%-43.08%)
  mva<-m   ##proportion of vaccinated(1st dose) Asymptomatic individuals 34%(26.65%-43.08%)
  mvb<-m   ##proportion of vaccinated (2nd dose) Asymptomatic individuals 34%(26.65%-43.08%)
  
  
  theta_invall~dgamma(42.3,1.99) ##duration in all I comps 21(7-28)
  theta_inv<-theta_invall/3    ##duration in each I comp
  theta <- 1/theta_inv    ##daily rate of transition I1-I2-I3-R
  thetava <- theta   ##daily rate of transition I1va-I2va-I3va-R
  thetavb <- theta    ##daily rate of transition I1vb-I2vb-I3vb-R
  
  
  omega_invall~ dgamma(20.6,1.45)     ###duration in all A comps 14(7-21)days
  omega_inv<-omega_invall/3     ###duration in each A comp
  omega<-1/omega_inv     ##daily rate of transition A1-A2-A3-R
  omegava<-omega     ##daily rate of transition A1va-A2va-A3va-R
  omegavb<-omega     ##daily rate of transition A1vb-A2vb-A3vb-R
  
  
   # Reporting fraction and lag
  #report_frac~dbeta(20,6)
  #report_frac ~ dbeta(5, 15)    # Mean = 0.25, 95% CI ≈ [0.09, 0.47] worked for 666
  report_frac ~ dbeta(10, 10)  # Mean = 0.50

  report_lag <- 7     # 7-day reporting lag
  
  #####vaccine effectiveness of 1st dose 35.8% (95% CI, 22.1 to 47.1)
 
  #Vea ~ dbeta(25.97, 45.6)     # Mean ~0.36
  #Vea=0.30
  Vea ~ dbeta(49.3, 87.4)   ###tight priors around 36
  rsa=1-Vea   ##rsa is the residual susceptibility after dose1
  
  
  #####vaccine effectiveness of 2nd dose 66% (95% CI, 47 to 78)

  #Veb ~ dbeta(31, 14)
  #Veb=0.50
  
  Veb ~ dbeta(69.3, 35.6)  # tight priors around 66
  rsb=1-Veb   ##rsb is the residual susceptibility after dose2

  
  # Initial conditions for each group
 for (g in 1:3) {
    E_start[g] <- E0 * Norm_contact_dy[g]
    P_start[g] <- P0 * Norm_contact_dy[g]
    A1_start[g] <- A10 * Norm_contact_dy[g]
    I1_start[g] <- I10 * Norm_contact_dy[g]
    E[g,1]<-E_start[g]
    P[g, 1] <- P_start[g]
    A1[g, 1] <- A1_start[g]
    I1[g, 1] <- I1_start[g]
    S[g, 1] <- max(0,N[g] - E_start[g]-P_start[g]-A1_start[g]-I1_start[g])
    A2[g, 1] <- A2_start[g]
    A3[g, 1] <- A3_start[g]
    I2[g, 1] <- I2_start[g]
    I3[g, 1] <- I3_start[g]
    R[g, 1] <- R_start[g]
    Cuminc[g, 1]<-Cuminc_start[g]
    
    #####Inits for vaccinated comps va
    Va[g,1]<-Va_start[g]
    Eva[g,1]<-Eva_start[g]
    Pva[g, 1] <- Pva_start[g]
    A1va[g, 1] <- A1va_start[g]
    A2va[g, 1] <- A2va_start[g]
    A3va[g, 1] <- A3va_start[g]
    I1va[g, 1] <- I1va_start[g]
    I2va[g, 1] <- I2va_start[g]
    I3va[g, 1] <- I3va_start[g]
    
      #####Inits for vaccinated comps vb
    Vb[g,1]<-Vb_start[g]
    Evb[g,1]<-Evb_start[g]
    Pvb[g, 1] <- Pvb_start[g]
    A1vb[g, 1] <- A1vb_start[g]
    A2vb[g, 1] <- A2vb_start[g]
    A3vb[g, 1] <- A3vb_start[g]
    I1vb[g, 1] <- I1vb_start[g]
    I2vb[g, 1] <- I2vb_start[g]
    I3vb[g, 1] <- I3vb_start[g]
    
    for (d in 1:va_delay) {
    va_queue[g,d,1] <- 0
    vb_queue[g,d,1] <- 0
  }
 }

  # Define the assortative mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_ass[g, g_] <- ifelse(g == g_, 1, 0)
    }
  }
  
  # Define the proportional mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_prp[g, g_] <- (contact[g_] * N[g_]) / sum(contact * N)
    }
  }
  
  # Initialize prevalence for t = 1
  for (g in 1:3) {
    prevalence[g, 1] <- ((tau*P[g, 1]) + I1[g, 1] + I2[g, 1] + I3[g, 1]+
                          (tau*A1[g, 1]) + (tau*A2[g, 1]) + (tau*A3[g, 1])+
                          (tau*Pva[g, 1])+I1va[g, 1] + I2va[g, 1] + I3va[g, 1]+
                          (tau*A1va[g, 1])+ (tau*A2va[g, 1])+ (tau*A3va[g, 1])+
                          (tau*Pvb[g, 1])+I1vb[g, 1] + I2vb[g, 1] + I3vb[g, 1]+
                          (tau*A1vb[g, 1])+ (tau*A2vb[g, 1])+ (tau*A3vb[g, 1]))/ N[g]
  }
  

  # Initialize mixing matrix (mix) for t = 1
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix[g, g_, 1] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }

 
 # Calculate initial lambda_det for t = 1 based on the initial conditions
for (g in 1:3) {
  lambda_det[g, 1] <- sum(beta * contact[g] * mix[g, 1:3, 1] * prevalence[1:3, 1])
}

total_lambda[1] <- sum(lambda_det[1:3, 1])

 ####calculating new cases to include in likelihood

    total_new_cases[1] <- 0                    # Initial people in Icomp
    total_Cuminc[1] <- sum(Cuminc[1:3, 1])     # Initial cumulative Incidence
    
     active_infected[1] <- 0                  # Initial prevalence
     daily_shedding[1]<-0                     #Initial shed viral load
     
     # Initial group-specific daily shedding at t = 1
for (g in 1:3) {
  daily_shedding_g[g, 1] <- 0
}

     

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
      
      for (g in 1:3) {

  # 1. Shedding from P compartments (pre-symptomatic)
  shed_P[g, t] <- shed_alpha * (
    P[g, t] +
    Pva[g, t] +
    Pvb[g, t]
  )

  # 2. Shedding from A compartments (asymptomatic)
  shed_A[g, t] <- 
    shed_A1 * A1[g, t] + shed_A2 * A2[g, t] + shed_A3 * A3[g, t] +
    shed_A1a * A1va[g, t] + shed_A2a * A2va[g, t] + shed_A3a * A3va[g, t] +
    shed_A1b * A1vb[g, t] + shed_A2b * A2vb[g, t] + shed_A3b * A3vb[g, t]

  # 3. Shedding from I compartments (symptomatic)
  shed_I[g, t] <- 
    shed_I1 * I1[g, t] + shed_I2 * I2[g, t] + shed_I3 * I3[g, t] +
    shed_I1a * I1va[g, t] + shed_I2a * I2va[g, t] + shed_I3a * I3va[g, t] +
    shed_I1b * I1vb[g, t] + shed_I2b * I2vb[g, t] + shed_I3b * I3vb[g, t]
}

# Total shedding across all groups (optional, if still needed)
daily_shedding[t] <- sum(shed_P[1:3, t]) + sum(shed_A[1:3, t]) + sum(shed_I[1:3, t])

# Group-specific daily shedding
for (g in 1:3) {
  daily_shedding_g[g, t] <- shed_P[g, t] + shed_A[g, t] + shed_I[g, t]
}


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

sigma <- transit_time_mean * transit_time_cv  

for (i in 1:tmax) {
  g[i] <- (1 / (sigma * sqrt(2 * 3.14159))) * exp(-pow(i - transit_time_mean, 2) / (2 * pow(sigma, 2)))
}

for (t in 1:T) {
  for (i in 1:tmax) {
    lag[i, t] <- t - i + 1
    is_valid[i, t] <- step(lag[i, t])
    safe_index[i, t] <- max(lag[i, t], 1)
    safe_shedding[i, t] <- is_valid[i, t] * daily_shedding[safe_index[i, t]]
    contrib[i, t] <- safe_shedding[i, t] * g[i] * exp(-mu * i)
  }
  delayed_conc[t] <- sum(contrib[1:tmax, t])
}

for (i in 1:tmax) {
  g_kernel[i] <- (1 / (sigma * sqrt(2 * 3.14159))) * exp(-pow(i - transit_time_mean, 2) / (2 * pow(sigma, 2)))
}

# Group-specific delayed concentration
for (g in 1:3) {
  for (t in 1:T) {
    for (i in 1:tmax) {
      lag_g[i, t, g] <- t - i + 1
      is_valid_g[i, t, g] <- step(lag_g[i, t, g])
      safe_index_g[i, t, g] <- max(lag_g[i, t, g], 1)
      safe_shedding_g[i, t, g] <- is_valid_g[i, t, g] * daily_shedding_g[g, safe_index_g[i, t, g]]
      contrib_g[i, t, g] <- safe_shedding_g[i, t, g] * g_kernel[i] * exp(-mu * i)
    }
    delayed_conc_g[g, t] <- sum(contrib_g[1:tmax, t, g])
  }
}


  ##Case likelihood
for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    mu_nb[t] <- total_new_cases[t]+ 1e-6 # Mean for the Negative Binomial to ensure its not neg
    p_nb[t]<-phi/(phi+mu_nb[t]) ##convert dispersion param to success prob
    cases_obs[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)  # Negative Binomial likelihood
      # Posterior predictive distribution
    cases_pred[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)
}

 ##Viral load likelihood
  ##Scaling, normalization and log transformation
  
 for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
  cp_total_all[t] <- delayed_conc[t] * mult
  cp_per_person_all[t] <- cp_total_all[t] / wwtp_population
  cp_per_person_mL_all[t] <- cp_per_person_all[t] * flow_mlalldaily[t - burn_in_timesteps]
  log10_conc_all[t] <- log(cp_per_person_mL_all[t] + 1) / log(10)
}
for (w in 1:T_WWobs) {
   cp_total[w] <- delayed_conc[ww_sample_days[w]] * mult
  cp_per_person[w] <- cp_total[w] / wwtp_population
  #cp_per_person_mL[w] <- cp_per_person[w] * flow_mL_daily[w]
  #log10_conc[w] <- log(cp_per_person_mL[w] + 1) / log(10)
  cp_raw[w] <- cp_per_person[w] * flow_mL_daily[w]
  cp_safe[w] <- max(cp_raw[w], 1e-6)
  log10_conc[w] <- log(cp_safe[w]) / log(10)
  ww_obs[w] ~ dnorm(log10_conc[w], tau_ww)
  ww_pred[w] ~ dnorm(log10_conc[w], tau_ww)

}

# All days (for plotting/PPD over full timeline)
for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    cp_total_all_g[g, t] <- delayed_conc_g[g, t] * mult
    cp_per_person_all_g[g, t] <- cp_total_all_g[g, t] / wwtp_population
    cp_per_person_mL_all_g[g, t] <- cp_per_person_all_g[g, t] * flow_mlalldaily[t - burn_in_timesteps]
    log10_conc_all_g[g, t] <- log(cp_per_person_mL_all_g[g, t] + 1) / log(10)
  }
}

for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    # Posterior predictive draw for full daily series (all timesteps)
    log10_conc_all_g_pred[g, t] ~ dnorm(log10_conc_all_g[g, t], tau_ww)
  }
}

for (g in 1:3) {
  for (t in (burn_in_timesteps + 1):(T_caseobs + burn_in_timesteps)) {
    mu_nb_g[g, t] <- new_reported_cases[g, t] + 1e-6
    p_nb_g[g, t] <- phi / (phi + mu_nb_g[g, t])
    cases_pred_g[g, t - burn_in_timesteps] ~ dnegbin(p_nb_g[g, t], phi)
  }}
}"


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
##for each grp
##new inits for each grp
A2_start=c(0,0,0)
A3_start=c(0,0,0)
I2_start=c(0,0,0)
I3_start=c(0,0,0)
R_start=c(0,0,0)
Cuminc_start=c(0,0,0)

#####initials for vaccinated grps va
Va_start=c(0,0,0)
Pva_start=c(0,0,0)
Eva_start=c(0,0,0)
A1va_start=c(0,0,0)
A2va_start=c(0,0,0)
A3va_start=c(0,0,0)
I1va_start=c(0,0,0)
I2va_start=c(0,0,0)
I3va_start=c(0,0,0)


#####initials for vaccinated grps vb
Vb_start=c(0,0,0)
Pvb_start=c(0,0,0)
Evb_start=c(0,0,0)
A1vb_start=c(0,0,0)
A2vb_start=c(0,0,0)
A3vb_start=c(0,0,0)
I1vb_start=c(0,0,0)
I2vb_start=c(0,0,0)
I3vb_start=c(0,0,0)

case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat$total_cases 
T_caseobs=length(cases_obsb)


#I1_start+I2_start+I3_start
# case_dat=read_excel("data_25/case_data_V2.xlsx",sheet="cases")
# cases_obsb = case_dat$total_cases 
#T_obsb=length(cases_obsb)

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

########1st convert all parameters to log scale
tau=0.776
shed.I1 = 10^12.3
shed.I2 = 10^11.5  
shed.I3 = 10^10.2 
shed.I1a<-shed.I1 
shed.I2a<-shed.I2 
shed.I3a<-shed.I3
shed.I1b<-shed.I1 
shed.I2b<-shed.I2 
shed.I3b<-shed.I3
shed.A1 = shed.I1*tau 
shed.A2 = shed.I2*tau 
shed.A3 = shed.I3*tau 
shed.A1a<-shed.A1 
shed.A2a<-shed.A2 
shed.A3a<-shed.A3 
shed.A1b<-shed.A1 
shed.A2b<-shed.A2 
shed.A3b<-shed.A3 
alpha_log=shed.I1*tau  

# Convert all shedding rates from log10(cp/mL) to cp/mL (linear scale)
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
names(ww_dat)
ww_std = ww_dat %>% select(log10_cp_per_person_per_day) #####standardised WW data
#ww_stdlinear = ww_dat %>% select(cp_per_person_per_day) #####standardised WW data
#ww_raw = ww_dat %>% select(log10_daily_avg_cp_ml) #####unstandardised WW data
flow_L_daily=ww_dat %>% select(agg_flow_dat) ###aggregated flow data
flow_L_daily <- as.numeric(unlist(flow_L_daily))
is.numeric(flow_L_daily)
#ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))
ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))
T_WWobs <- length(ww_obs)

# Convert flow to mL
flow_mL_daily <- flow_L_daily * 1e3

##########Flow rate for all data points
ww_flow=read_excel("Data/case_data_V2.xlsx",sheet="flowrate")
flow_all<-ww_flow %>%select(sumflow)
flow_alldaily<-as.numeric(unlist(flow_all))
flow_mlalldaily<-flow_alldaily*1e3

###########preposses start date and end date into my data block
# Define week indices (assuming 7-day weeks)

########Extract sample days
# Suppose these are the raw observation days (relative to sampling calendar)
raw_ww_days <- ww_dat$Day # in days from beginning of *observation*, not simulation
# Shift all by burn-in
ww_sample_days <- raw_ww_days + burn_in_timesteps

# Example: Suppose sampling occurred on days 35, 38, 42, ..., relative to model timeline
# Check result

# Final JAGS DataList
dataListcomb <- list(
  N = N, T = T,
  contact = contact,
  Norm_contact_dy = Norm_contact_dy,
  A2_start = A2_start, A3_start = A3_start,
  I2_start = I2_start, I3_start = I3_start, R_start = R_start,
  Va_start = Va_start,
  Pva_start = Pva_start, Eva_start = Eva_start, A1va_start = A1va_start,
  A2va_start = A2va_start, A3va_start = A3va_start, I1va_start = I1va_start,
  I2va_start = I2va_start, I3va_start = I3va_start,
  Vb_start = Vb_start,
  Pvb_start = Pvb_start, Evb_start = Evb_start, A1vb_start = A1vb_start,
  A2vb_start = A2vb_start, A3vb_start = A3vb_start, I1vb_start = I1vb_start,
  I2vb_start = I2vb_start, I3vb_start = I3vb_start,Cuminc_start=Cuminc_start,
  T_WWobs = T_WWobs,
  T_caseobs=T_caseobs,
  burn_in_timesteps=burn_in_timesteps,
  alpha_v = alpha_v,
  mu_v = mu_v,
  tau=tau,
  va_delay = 14,
  vb_delay = 14,
  # Shedding inputs
  shed_I1 = shed.I1, shed_I2 = shed.I2, shed_I3 = shed.I3,
  shed_I1a = shed.I1a, shed_I2a = shed.I2a, shed_I3a = shed.I3a,
  shed_I1b = shed.I1b, shed_I2b = shed.I2b, shed_I3b = shed.I3b,
  shed_A1 = shed.A1, shed_A2 = shed.A2, shed_A3 = shed.A3,
  shed_A1a = shed.A1a, shed_A2a = shed.A2a, shed_A3a = shed.A3a,
  shed_A1b = shed.A1b, shed_A2b = shed.A2b, shed_A3b = shed.A3b,
  shed_alpha = alpha_log,
  # WW model inputs
  flow_mL_daily = flow_mL_daily,
  flow_mlalldaily=flow_mlalldaily,
  wwtp_population = 1278020,
  ww_sample_days=ww_sample_days,
  # Observed WW data
  ww_obs = ww_obs,  # or ww_raw if unstandardized
  cases_obs=cases_obsb,
  ####precomputed g to include in the advection dispersion decay model
  #transit_time_cv=0.3,     #std dev transit time between shedding and sampling sites (in days)
  tmax=15)  #Max mean = 5,Max SD = 5 × 0.5 = 2.5,Max plausible delay ≈ mean + 5×SD = 5 + (5×2.5) = 17.5
###precomputed start and end dates of defining the epi weeks
inits_list <- list(
  list(
    beta = 0.8, kappa = 0.95, report_frac = 0.50,
    log_mult = log(3.5e-9),  # Based on fixed value that worked
    tau_ww = 0.4,          # Around posterior median (0.43)
    transit_time_mean = 2.5,  # Close to prior mean
    transit_time_cv = 0.3,    # Close to prior mean
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 42
  ),
  list(
    beta = 0.9, kappa = 0.92, report_frac = 0.55,
    log_mult = log(2.5e-9),  # Slight variation for chain independence
    tau_ww = 0.5,
    transit_time_mean = 2.8,
    transit_time_cv = 0.35,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 99
  )
)


# #Run the model with different initial values for each chain
system.time({
  Combined_finaltestfb<- run.jags(textstring, data = dataListcomb,
                     monitor = c("ww_pred","cases_pred","log10_conc_all",
                                 "cases_pred_g","mu_nb_g",  # group-level predictive + params
                                 "log10_conc_all_g_pred", # all-days per group
                                 "log10_conc_all_g",
                                 "E0","P0","A10","I10",
                                 "log10_conc","mu_nb",
                                 "P_total", "A_total", "I_total",
                                 "shed_P", "shed_A", "shed_I","mult","log_mult",
                                 "tau_ww","transit_time_mean","transit_time_cv",
                                 "beta","kappa","phi",
                                 "total_Cuminc", "active_infected","total_lambda",
                                 "report_frac","Vea","Veb","m",
                                 "delta_inv","theta_invall","omega_invall"),
                     method="parallel",
                     #sample = 2000, adapt =500, burnin = 500, thin = 1,
                     sample = 40000, adapt =4000, burnin = 4000, thin = 2,
                     n.chains = 2, inits = inits_list,
                     summarise = FALSE)
})

Comb_finaltestfb<- as.mcmc.list(  Combined_finaltestfb)
save(Comb_finaltestfb,file="U:/mpox25output/Comb_finaltestfb.RData")

############generate output
load(file="U:/mpox25output/Comb_finaltestc.RData")
###generate traceplots
traceplot(Comb_finaltestc[, "transit_time_mean"],main="Mean transit time in sewer")
traceplot(Comb_finaltestc[, "transit_time_cv"],main="Standard deviation of transit mean time")
traceplot(Comb_finaltestc[, "mult"],main="Scaling factor of viral load")
traceplot(Comb_finaltestc[, "tau_ww"],main="Precision of the dnorm likelihood")
traceplot(Comb_finaltestc[, "beta"],main="Transmission parameter")
traceplot(Comblist_finald[, "kappa"],main="Mixing probability")
traceplot(Comblist_finald[, "phi"],main="Negative binomial dispersion parmeter")
traceplot(Comb_finaltestc[, "report_frac"],main="reporting fraction")
traceplot(Comblist_finald[, "m"],main="Proportion of Asymptomatic fraction")

#####extract samples for plotting
chain_1_samples <- Comblist_mod7[[1]]
mcmc_matrixall<-as.matrix(Combined_finaltestb)
###randomly sample the list to generate summaries of predicted data
mcmc_matrix<-as.matrix(chain_1_samples)
total_samples <- nrow(mcmc_matrix)
# Randomly sample 1000 indices from the total number of samples
sample_indices <- sample(1:total_samples, size = 9000, replace = FALSE)
# Extract the sampled rows from the mcmc_matrix
sampled_mcmc <- mcmc_matrix[sample_indices, ]

# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

# Summarize the posterior distributions for all parameters
#posterior_summary <- summary_median_CI(as.matrix(sampled_mcmc))
posterior_summaryb <- summary_median_CI(mcmc_matrixall)
posterior_summaryb[grep("mult|tau_ww|transit_time_mean|transit_time_cv", rownames(posterior_summaryb)), ]

# If you want to focus on specific parameters, e.g., total_new_cases and new_I3_cases:
total_new_WW_summary <- as.data.frame(posterior_summaryb[grep("ww_pred", rownames(posterior_summaryb)), ])
total_delayed_summary <- as.data.frame(posterior_summaryb[grep("delayed_conc", rownames(posterior_summaryb)), ])

#total_new_case_summary <- as.data.frame(posterior_summary[grep("cases_pred", rownames(posterior_summary)), ])
#prev_summary <- as.data.frame(posterior_summaryb[grep("active_infected", rownames(posterior_summaryb)), ])
#CumInc_summary <- as.data.frame(posterior_summaryb[grep("total_Cuminc", rownames(posterior_summaryb)), ])

###read in the observed WW data
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
#names(ww_dat)
ww_std = ww_dat %>% select(log10_cp_per_person_per_day) #####standardised WW data
ww_stdlinear = ww_dat %>% select(cp_per_person_per_day) #####standardised WW data
ww_raw = ww_dat %>% select(log10_daily_avg_cp_ml) #####unstandardised WW data
ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))
#ww_obs = as.numeric(unlist(ww_stdlinear$cp_per_person_per_day))
###read in the observed case data
#total_new_cases_summaryb<-total_new_cases_summary %>% mutate(Day=1:169)
# case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
# cases_obsb = case_dat %>% select(total_cases)
# 
# ##########Generate case data fir
# ##Create a data frame for plotting
# plot_casedat <- data.frame(
#   time = 1:nrow(cases_obsb),                    # time index
#   Date=case_dat$Date,
#   observed = cases_obsb$total_cases,                    # observed cases
#   median_fit = total_new_case_summary$median,          # model median fit
#   lower_ci = total_new_case_summary$lower_95_CI,          # lower 95% CI
#   upper_ci = total_new_case_summary$upper_95_CI           # upper 95% CI
# )
# 
# # Plot the observed cases and model fit with 95% CI
# plot_casefit=ggplot(plot_casedat, aes(x = time)) +
#   geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
#   geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
#   labs(x = "Time", y = "Reported Mpox cases", title = "fit vs. Observed Cases") +
#   theme_minimal() +
#   theme(legend.position = "top")
# 
# plot_casefit

##Create a data frame for plotting
plot_wwdata <- data.frame(
  time = 1:nrow(ww_std),                    # time index
  #Time=case_dat$Date,
  observed = ww_obs,                    # observed cases
  median_fit = total_new_WW_summary$median,          # model median fit
  lower_ci = total_new_WW_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_new_WW_summary$upper_95_CI           # upper 95% CI
)

# Plot the observed cases and model fit with 95% CI
plot_wwfit=ggplot(plot_wwdata, aes(x = time)) +
  geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
  geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
  labs(x = "Time", y = "Reported viral load", title = "fit vs. Observed viral load") +
  theme_minimal() +
  theme(legend.position = "top")

plot_wwfit

##Create a data frame for plotting
obs_days <- ww_sample_days  # length 48
obs_viral_load <- ww_raw$log10_daily_avg_cp_ml# log10(cp/mL/person)
total_delayed_summary$day <- 1:nrow(total_delayed_summary)

plot_df <- total_delayed_summary %>%
  filter(day %in% obs_days) %>%
  mutate(observedraw = obs_viral_load,
         observedstd = ww_obs)

#names(plot_df)
ggplot(plot_df, aes(x = day)) +
  geom_ribbon(aes(ymin = log10(lower_95_CI), ymax = log10(upper_95_CI)), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = log10(median)), color = "blue", size = 1) +
  geom_point(aes(y = observedraw), color = "black", shape = 16, size = 2) +
  geom_point(aes(y = observedstd), color = "red", shape = 16, size = 2) +
  labs(
    title = "Predicted Delayed Viral Load vs. Observed raw(black) vs Observed std(red)",
    x = "Day",
    y = "log10 copies/mL/person"
  ) +
  theme_minimal()

























###########plot prevalence
burn_in_timesteps <- 30
prev_summaryb <- prev_summary[(burn_in_timesteps + 1):nrow(prev_summary), ]

plot_prev <- data.frame(
  time = 1:nrow(cases_obsb),                    # time index
  #Time=case_dat$Date,
  #observed = cases_obsb,                    # observed cases
  median_fit = prev_summaryb$median,          # model median fit
  lower_ci = prev_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_summaryb$upper_95_CI          # upper 95% CI
)

plotprev=ggplot(plot_prev, aes(x = time)) +
  #geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
  geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
  labs(x = "Time", y = "Total prevalence(Actively infected)", title = "Prevalence-Combined model") +
  theme_minimal() +
  theme(legend.position = "top")

plotprev

CumInc_summary <- CumInc_summary[(burn_in_timesteps + 1):nrow(CumInc_summary), ]

plot_CumInc <- data.frame(
  time = 1:nrow(cases_obsb),                    # time index
  #Time=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = CumInc_summary$median,          # model median fit
  lower_ci = CumInc_summary$lower_95_CI,          # lower 95% CI
  upper_ci = CumInc_summary$upper_95_CI          # upper 95% CI
)

plotCumInc=ggplot(plot_CumInc, aes(x = time)) +
  #geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
  geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
  labs(x = "Time", y = "Cumulative incidence", title = "Cumulative incidence-Combined model") +
  theme_minimal() +
  theme(legend.position = "top")

plotCumInc














