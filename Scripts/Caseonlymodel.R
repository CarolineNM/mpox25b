#rm(list=ls())
library(R2jags)
library(runjags)
library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)
options(scipen=999)

###how to interpret beta in my model
##In this one I have constrained beta to range between 0 and 1

textstring <- "
model {
  
   beta~dnorm(0.8, 25) T(0,1) ##casemod 2.SD=0.2

  E0 ~ dpois(5)    #  A Poisson prior with mean 5
  P0 ~ dpois(1)  # Total pre-symptomatic individuals
  A10 ~ dpois(1) # Total asymptomatic individuals
  I10 ~ dpois(1) # Total symptomatic individuals
  
  #kappa ~ dbeta(20, 5)  ##wider prior ranging from 0 to 99%
  kappa ~ dbeta(40, 2)  ### Mean ~0.95, 95% CI ≈ [0.85, 0.995]
  
  
  #phi ~ dgamma(0.1, 0.1) ##negb11 prior
  phi ~ dgamma(2, 0.5)  # Mean = 4, SD = ~2.8
 
  
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


    total_new_cases[1] <- 0 # Initial people in Icomp
    total_Cuminc[1] <- sum(Cuminc[1:3, 1])     # Initial cumulative Incidence
    
     active_infected[1] <- 0                  # Initial prevalence
    

  # Dynamics of the model (deterministic-like)
  for (t in 2:T) {
  
   # Final mixing matrix for three groups only
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix[g, g_, t] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }
    

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

        new_Vas[g, t-1] <- alpha[t] * S[g, t-1]           # Deterministic S to Va transition.alpha is vaccination rate
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
        
          new_Vbs[g, t-1] <- mu[t] * Va[g, t-1]              # Deterministic Va to Vb transition.mu is vaccination rate
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
    
       # Total new reported cases at time t
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
  }
  

   for (t in (burn_in_timesteps + 1):(T_obs + burn_in_timesteps)) {
    mu_nb[t] <- total_new_cases[t]+ 1e-6 # Mean for the Negative Binomial to ensure its not neg
    p_nb[t]<-phi/(phi+mu_nb[t]) ##convert dispersion param to success prob
    cases_obs[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)  # Negative Binomial likelihood
     # Posterior predictive distribution
    cases_pred[t - burn_in_timesteps] ~ dnegbin(p_nb[t], phi)
   }
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



#I1_start+I2_start+I3_start
case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat$total_cases 
T_obsb=length(cases_obsb)


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
alpha<- c(rep(0, burn_in_timesteps), v)   ###alpha is the per capita vaccination rate of 1st dose at time t
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
mu<- c(rep(0, burn_in_timesteps), v2)   ###mu is the per capita vaccination rate of 2nd dose at time t
#mu_test<- c(rep(0, 9),0.1,rep(0, 189))
#mutest set to zero
#mu<- rep(0, 199)
# Data list for JAGS
dataList <- list(N = N, T = T, contact = contact,Norm_contact_dy=Norm_contact_dy,
                 A2_start = A2_start,A3_start=A3_start,
                 I2_start = I2_start,I3_start=I3_start, R_start = R_start,
                 Va_start=Va_start,
                 Pva_start=Pva_start,Eva_start=Eva_start,A1va_start=A1va_start,
                 A2va_start=A2va_start,A3va_start=A3va_start,I1va_start=I1va_start,
                 I2va_start=I2va_start,I3va_start=I3va_start,
                 Vb_start=Vb_start,
                 Pvb_start=Pvb_start,Evb_start=Evb_start,A1vb_start=A1vb_start,
                 A2vb_start=A2vb_start,A3vb_start=A3vb_start,I1vb_start=I1vb_start,
                 I2vb_start=I2vb_start,I3vb_start=I3vb_start,Cuminc_start=Cuminc_start,
                 cases_obs=cases_obsb,
                 T_obs=T_obsb,
                 alpha=alpha,
                 mu=mu,
                 tau=0.8,
                 va_delay= 14,
                 vb_delay=14,
                 burn_in_timesteps=burn_in_timesteps)


inits_list <- list(
  list(beta = 0.8, kappa = 0.95, report_frac = 0.50, phi = 4,Vea=0.3,Veb=0.5,.RNG.name = "base::Wichmann-Hill", .RNG.seed = 69),
  list(beta = 0.9, kappa = 0.92, report_frac = 0.55, phi = 5,Vea=0.3,Veb=0.5,.RNG.name = "base::Wichmann-Hill", .RNG.seed = 99)
)


#Run the model with different initial values for each chain
system.time({
  Case_modfinal<- run.jags(textstring, data = dataList,
                            monitor = c("beta", "kappa","phi","cases_pred",
                                        "total_lambda","report_frac","Vea","Veb","m",
                                        "delta_inv","theta_invall","omega_invall",
                                        "total_Cuminc", "active_infected"),
                            method="parallel",
                            #sample = 50000, adapt =10000, burnin = 10000, thin = 2,
                            sample = 30000, adapt =4000, burnin = 4000, thin = 2,
                            n.chains = 2, inits = inits_list,
                            summarise = FALSE)
})

Case_modlstfinal<- as.mcmc.list(Case_modfinal)
save(Case_modlstfinal,file="Output/Case_modlstfinal.RData")

###############################################################
load("Output/Case_modlist3.RData")
###generate traceplots
traceplot(Case_modlist3[, "beta"])
traceplot(Case_modlist3[, "kappa"])
traceplot(Case_modlist3[, "phi"])
traceplot(Case_modlist3[, "report_frac"])
traceplot(Case_modlist3[, "m"])
#####extract samples for plotting
chain_1_samples <- Case_modlist3[[1]]  
mcmc_matrixallb<-as.matrix(Case_modlist3)
###randomly sample the list to generate summaries of predicted data
mcmc_matrix<-as.matrix(chain_1_samples)
total_samples <- nrow(mcmc_matrix)
# Randomly sample 1000 indices from the total number of samples
sample_indices <- sample(1:total_samples, size = 45000, replace = FALSE)
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
posterior_summary <- summary_median_CI(as.matrix(sampled_mcmc))
posterior_summaryc <- summary_median_CI(mcmc_matrixallb)
posterior_summaryc[grep("beta|kappa|phi|report_frac|Vea|Veb|delta_inv|theta_invall|omega_invall", rownames(posterior_summaryc)), ]

# If you want to focus on specific parameters, e.g., total_new_cases and new_I3_cases:
total_new_cases_summary <- as.data.frame(posterior_summaryc[grep("cases_pred", rownames(posterior_summaryc)), ])
prev_summary <- as.data.frame(posterior_summaryc[grep("active_infected", rownames(posterior_summaryc)), ])
CumInc_summary <- as.data.frame(posterior_summaryc[grep("total_Cuminc", rownames(posterior_summaryc)), ])

#nrow(total_new_cases_summaryb)
total_new_cases_summaryb<-total_new_cases_summary %>% mutate(Day=1:169)
case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat %>% select(total_cases)

##Create a data frame for plotting
plot_data <- data.frame(
  time = 1:nrow(cases_obsb),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb$total_cases,                    # observed cases
  median_fit = total_new_cases_summary$median,          # model median fit
  lower_ci = total_new_cases_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_new_cases_summary$upper_95_CI           # upper 95% CI
)

# Plot the observed cases and model fit with 95% CI
plotfit=ggplot(plot_data, aes(x = time)) +
  geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
  geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
  labs(x = "Time", y = "Reported Mpox cases", title = "fit vs. Observed Cases") +
  theme_minimal() +
  theme(legend.position = "top")

plotfit

##############plot prevalence##################
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



##########plotting a smoothed version
# Smooth each line separately
loess_median <- loess(median ~ Day, data = total_new_cases_summaryb, span = 0.40)
loess_upper <- loess(upper_95_CI ~ Day, data = total_new_cases_summaryb, span = 0.40)
loess_lower <- loess(lower_95_CI ~ Day, data = total_new_cases_summaryb, span = 0.40)
# Predict the smoothed values
total_new_cases_summaryb$median_smoothed <- pmax(0,predict(loess_median))
total_new_cases_summaryb$upper_smoothed <- pmax(0,predict(loess_upper))
total_new_cases_summaryb$lower_smoothed <- pmax(0,predict(loess_lower))

##Create a data frame for plotting
plot_datasmoothed <- data.frame(
  time = 1:nrow(cases_obsb),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb$total_cases,                    # observed cases
  median_fit = total_new_cases_summaryb$median_smoothed,          # model median fit
  lower_ci = total_new_cases_summaryb$lower_smoothed,          # lower 95% CI
  upper_ci = total_new_cases_summaryb$upper_smoothed          # upper 95% CI
)

# Plot the observed cases and model fit with 95% CI
plotfitsm=ggplot(plot_datasmoothed, aes(x = time)) +
  geom_point(aes(y = observed), color = "black", size = 1) +  # observed cases
  geom_line(aes(y = median_fit), color = "blue", size = 1) +                      # model median fit
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.5) +  # 95% CI
  labs(x = "Time", y = "Reported Mpox cases", title = "fit vs. Observed Cases") +
  theme_minimal() +
  theme(legend.position = "top")

plotfitsm








