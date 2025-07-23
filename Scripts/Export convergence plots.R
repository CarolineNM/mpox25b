#########Call the libraries
####This has model fits,prevalence and cumulative incidence,relative contribution
rm(list=ls())
library(R2jags)
library(runjags)
library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)
library(bayesplot)
library(truncnorm)
library(ggridges)
library(patchwork)
library(scales)
options(scipen=999)

#########load the raw datasets####
##1. case data ################
case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat$total_cases
Cumulat_inc=case_dat$Cumulat_Inc
##2. WW data #########################
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
ww_std = ww_dat %>% select(log10_cp_per_person_per_day) #####standardised WW data
ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))

#######load the three outputs
load("D:/mpox25output/Comblist_finalg.RData")
load("D:/mpox25output/Combined_WWd.RData")
load("D:/mpox25output/Combined_casd.RData")

######generate the model fit
mcmc_matrixallWW<-as.matrix(Combined_WWd)
mcmc_matrixallcas<-as.matrix(Combined_casd) ###most recent version
mcmc_matrixallcom<-as.matrix(Comblist_finalg)

# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summarycas <- summary_median_CI(mcmc_matrixallcas)
posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)
posterior_summaryww <- summary_median_CI(mcmc_matrixallWW)

total_cas_summary <- as.data.frame(posterior_summarycas[grep("cases_pred", rownames(posterior_summarycas)), ])
total_casb_summary <- as.data.frame(posterior_summarycas[grep("mu_nb", rownames(posterior_summarycas)), ])
total_cas_ww_summary <- as.data.frame(posterior_summarycas[grep("^log10_conc\\[", rownames(posterior_summarycas)), ])
total_cas_wwb_summary <- as.data.frame(posterior_summarycas[grep("ww_pred", rownames(posterior_summarycas)), ])

total_ww_summary <- as.data.frame(posterior_summaryww[grep("ww_pred", rownames(posterior_summaryww)), ])
total_ww_summaryb <- as.data.frame(posterior_summaryww[grep("^log10_conc\\[", rownames(posterior_summaryww)), ])
total_ww_cas_summary <- as.data.frame(posterior_summaryww[grep("mu_nb", rownames(posterior_summaryww)), ])
total_ww_casb_summary <- as.data.frame(posterior_summaryww[grep("cases_pred", rownames(posterior_summaryww)), ])


total_com_ww_summary <- as.data.frame(posterior_summarycom[grep("ww_pred", rownames(posterior_summarycom)), ])
total_com_ww_summaryb <- as.data.frame(posterior_summarycom[grep("^log10_conc\\[", rownames(posterior_summarycom)), ])
total_com_cas_summary <- as.data.frame(posterior_summarycom[grep("cases_pred", rownames(posterior_summarycom)), ])
total_com_cas_summaryb <- as.data.frame(posterior_summarycom[grep("mu_nb", rownames(posterior_summarycom)), ])


###########generate data from the combined ww model##############3
plot_dataww <- data.frame(
  time = 1:nrow(ww_std),
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_ww_summary$median,          # model median fit
  lower_ci = total_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_summary$upper_95_CI          # upper 95% CI
)

plot_datawwb <- data.frame(
  time = 1:nrow(ww_std),
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_ww_summaryb$median,          # model median fit
  lower_ci = total_ww_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_summaryb$upper_95_CI          # upper 95% CI
)


plot_dataww_cases <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_ww_cas_summary$median,         # model median fit
  lower_ci = total_ww_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_cas_summary$upper_95_CI         # upper 95% CI
)

plot_dataww_casesb <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_ww_casb_summary$median,         # model median fit
  lower_ci = total_ww_casb_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_casb_summary$upper_95_CI         # upper 95% CI
)

###############generate data from the combined case model##########################
plot_datcas <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_cas_summary$median,         # model median fit
  lower_ci = total_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_summary$upper_95_CI          # upper 95% CI
)

plot_datcasb <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_casb_summary$median,         # model median fit
  lower_ci = total_casb_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_casb_summary$upper_95_CI          # upper 95% CI
)


plot_datacas_ww <- data.frame(
  time = 1:nrow(ww_std), 
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_cas_ww_summary$median,          # model median fit
  lower_ci = total_cas_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_ww_summary$upper_95_CI          # upper 95% CI
)


plot_datacas_wwb <- data.frame(
  time = 1:nrow(ww_std), 
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_cas_wwb_summary$median,          # model median fit
  lower_ci = total_cas_wwb_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_wwb_summary$upper_95_CI          # upper 95% CI
)

####################generate data from the combined model############################

plot_datcomcas <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_com_cas_summary$median,         # model median fit
  lower_ci = total_com_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_cas_summary$upper_95_CI          # upper 95% CI
)

plot_datcomcasb <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_com_cas_summaryb$median,         # model median fit
  lower_ci = total_com_cas_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_cas_summaryb$upper_95_CI          # upper 95% CI
)

plot_datcomww <- data.frame(
  time = 1:nrow(ww_std), 
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_com_ww_summary$median,          # model median fit
  lower_ci = total_com_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_ww_summary$upper_95_CI          # upper 95% CI
)

plot_datcomwwb <- data.frame(
  time = 1:nrow(ww_std),  
  Date=ww_dat$date,
  observed = ww_obs,                    # observed cases
  median_fit = total_com_ww_summaryb$median,          # model median fit
  lower_ci = total_com_ww_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_ww_summaryb$upper_95_CI          # upper 95% CI
)


#################Generate posterior plots##########
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Shared plot elements
plot_geom <- list(
  geom_point(aes(y = observed), color = "black", size = 1.2, alpha = 0.7),
  geom_line(aes(y = median_fit), color = "blue4", size = 1),
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.3)
)

# Cases – Case-only model
plot_casefit <- ggplot(plot_datcas, aes(x = Date)) +
  plot_geom +
  ylim(0, 20)+
  #ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    title = "(Case-only model)"
  ) +
  custom_theme

# Cases – Case-only model
plot_casefitb <- ggplot(plot_datcasb, aes(x = Date)) +
  plot_geom +
  ylim(0, 11)+
  #ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    #title = "(Case-only model)"
  ) +
  custom_theme

# Cases – ww-only model
plot_casewwfit <- ggplot(plot_dataww_cases, aes(x = Date)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    #title = "(WW-only model)"
  ) +
  custom_theme

# Cases – ww-only model
plot_casewwfitb <- ggplot(plot_dataww_casesb, aes(x = Date)) +
  plot_geom +
  ylim(0, 20)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    title = "(WW-only model)"
  ) +
  custom_theme


# Viral load – WW-only model
plot_wwfit <- ggplot(plot_dataww, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(WW-only model)"
  ) +
  custom_theme

plot_wwfitb <- ggplot(plot_datawwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(WW-only model)"
  ) +
  custom_theme


# Viral load – Case only model
plot_wwcasefit <- ggplot(plot_datacas_ww, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Case-only model)"
  ) +
  custom_theme

plot_wwcasefitb <- ggplot(plot_datacas_wwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Case-only model)"
  ) +
  custom_theme


# Viral load – Combined model
plot_comwwfit <- ggplot(plot_datcomww, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Combined model)"
  ) +
  custom_theme

plot_comwwfitb <- ggplot(plot_datcomwwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 9)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Combined model)"
  ) +
  custom_theme


# Cases – Combined model
plot_comcasefit <- ggplot(plot_datcomcas, aes(x = Date)) +
  plot_geom +
  ylim(0,20)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    title = "(Combined model)"
  ) +
  custom_theme


plot_comcasefitb <- ggplot(plot_datcomcasb, aes(x = Date)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    #title = "(Combined model)"
  ) +
  custom_theme



##########Predictive fits

plot_comcasefitb
plot_comcasefit
plot_comwwfitb
plot_comwwfit
plot_wwcasefit
plot_wwcasefitb
plot_wwfitb
plot_wwfit
plot_casewwfit
plot_casewwfitb
plot_casefitb
plot_casefit


ZZ=(plot_casefit+plot_comcasefit+plot_casewwfitb)/(plot_casefitb+plot_comcasefitb+plot_casewwfit)
ZZ+plot_layout(guides = "collect") & theme(legend.position = "bottom")

mm=(plot_casefit+plot_comcasefit)/(plot_casefitb+plot_comcasefitb)
mm

vv=(plot_wwfit+plot_comwwfit)/(plot_wwfitb+plot_comwwfitb)
vv


yy=(plot_wwfit+plot_comwwfit+plot_wwcasefitb)/(plot_wwfitb+plot_comwwfitb+plot_wwcasefit)
yy

ggsave(
  filename = "D:/Mpox25output/Figures/modfitb.tiff",
  plot = yy,  # use your actual plot variable here
  width = 13,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

##########generate data for model predicted prevalence
prev_wwsummary <- as.data.frame(posterior_summaryww[grep("active_infected", rownames(posterior_summaryww)), ])
prev_cassummary <- as.data.frame(posterior_summarycas[grep("active_infected", rownames(posterior_summarycas)), ])
prev_comsummary <- as.data.frame(posterior_summarycom[grep("active_infected", rownames(posterior_summarycom)), ])

burn_in_timesteps <- 30
prev_wwsummaryb <- prev_wwsummary[(burn_in_timesteps + 1):nrow(prev_wwsummary), ]
prev_cassummaryb <- prev_cassummary[(burn_in_timesteps + 1):nrow(prev_cassummary), ]
prev_comsummaryb <- prev_comsummary[(burn_in_timesteps + 1):nrow(prev_comsummary), ]


plot_wwprev <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index
  median_fit = prev_wwsummaryb$median,          # model median fit
  lower_ci = prev_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_cassummaryb$median,          # model median fit
  lower_ci = prev_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_cassummaryb$upper_95_CI         # upper 95% CI
)


plot_comprev <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index
  median_fit = prev_comsummaryb$median,          # model median fit
  lower_ci = prev_comsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_comsummaryb$upper_95_CI         # upper 95% CI
)


#############Generate prevalence plots
# Add model labels
plot_wwprev$model <- "WW-only"
plot_casprev$model <- "Case-only"
plot_comprev$model <- "Combined"

# Combine into one dataframe
plot_allprev <- bind_rows(plot_wwprev, plot_casprev, plot_comprev)

plot_prev_all <- ggplot(plot_allprev, aes(x = Date, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Date",
    y = "Active prevalence",
    #title = "Model-predicted prevalence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5))



plot_prev_all


##########Generate incidence plots##########
Cuminc_wwsummary <- as.data.frame(posterior_summaryww[grep("total_Cuminc", rownames(posterior_summaryww)), ])
Cuminc_cassummary <- as.data.frame(posterior_summarycas[grep("total_Cuminc", rownames(posterior_summarycas)), ])
Cuminc_comsummary <- as.data.frame(posterior_summarycom[grep("total_Cuminc", rownames(posterior_summarycom)), ])

burn_in_timesteps <- 30
Cuminc_wwsummaryb <- Cuminc_wwsummary[(burn_in_timesteps + 1):nrow(Cuminc_wwsummary), ]
Cuminc_cassummaryb <- Cuminc_cassummary[(burn_in_timesteps + 1):nrow(Cuminc_cassummary), ]
Cuminc_comsummaryb <- Cuminc_comsummary[(burn_in_timesteps + 1):nrow(Cuminc_comsummary), ]


plot_wwCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_wwsummaryb$median,          # model median fit
  lower_ci = Cuminc_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casCuminc <- data.frame(
  time = 1:nrow(case_dat), 
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_cassummaryb$median,          # model median fit
  lower_ci = Cuminc_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_cassummaryb$upper_95_CI         # upper 95% CI
)

plot_comCuminc <- data.frame(
  time = 1:nrow(case_dat),
  Date=case_dat$Date,# time index# time index
  median_fit = Cuminc_comsummaryb$median,          # model median fit
  lower_ci = Cuminc_comsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_comsummaryb$upper_95_CI         # upper 95% CI
)


#########generate plots for cumulative incidence
plot_wwCuminc$model <- "WW-only"
plot_casCuminc$model <- "Case-only"
plot_comCuminc$model <- "Combined"

# Combine into one dataframe
plot_allCuminc <- bind_rows(plot_wwCuminc, plot_casCuminc,plot_comCuminc)


plot_Cuminc_all <- ggplot(plot_allCuminc, aes(x = Date, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Date",
    y = "Cumulative Incidence",
    #title = "Model-predicted cumulative incidence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5))



y=plot_prev_all+plot_Cuminc_all+plot_layout(guides = "collect") & theme(legend.position = "bottom")


ggsave(
  filename = "D:/Mpox25output/Figures/prevcum.tiff",
  plot = y,  # use your actual plot variable here
  width = 13,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

case_mod  <- Combined_casd
ww_mod   <- Combined_WWd
com_mod <- Comblist_finalg


param_model_map <- list(
  beta = c("Case-only", "Viral load-only", "Combined"),
  kappa = c("Case-only", "Viral load-only", "Combined"),
  delta_inv=c("Case-only", "Viral load-only", "Combined"),
  theta_invall=c("Case-only", "Viral load-only", "Combined"),
  omega_invall=c("Case-only", "Viral load-only", "Combined"),
  Vea=c("Case-only", "Viral load-only", "Combined"),
  Veb= c("Case-only", "Viral load-only", "Combined"),
  phi = c("Case-only","Viral load-only", "Combined"),
  report_frac = c("Case-only","Viral load-only","Combined"),
  m = c("Case-only","Viral load-only","Combined"),
  transit_time_mean = c("Case-only","Viral load-only", "Combined"),
  transit_time_cv = c("Case-only","Viral load-only", "Combined"),
  mult = c("Case-only","Viral load-only", "Combined"),
  tau_ww = c("Viral load-only","Case-only", "Combined")
)

###########Add parameter labels
param_labels <- c(
  beta = "Transmission probability, β",
  kappa = "Mixing probability, κ",
  phi = "Negative binomial dispersion, φ",
  report_frac = "Reporting fraction",
  m = "Proportion asymptomatic, m",
  delta_inv = "Latent period (1/δ), days",
  theta_invall = "Infectious period (1/θ), days",
  omega_invall = "Asymptomatic period (1/ω), days",
  Vea = "Vaccine effectiveness (1st dose)",
  Veb = "Vaccine effectiveness (2nd dose)",
  transit_time_mean = "Mean transit time (days)",
  transit_time_cv = "Sd of transit time",
  mult = "Scaling factor for viral load",
  tau_ww = "Precision of WW likelihood (τ)"
)

# Set prior definitions
prior_definitions <- list(
  beta = function(n) rtruncnorm(n, 0, 1, mean = 0.8, sd = 0.1),
  kappa = function(n) rbeta(n, 40, 2),
  phi = function(n) rgamma(n, shape = 2, rate = 0.5),
  report_frac = function(n) rbeta(n, 10, 10),
  m = function(n) rbeta(n, 49.09, 94.88),
  delta_inv=function(n) rgamma(n, shape = 300.6, rate = 47.98),
  theta_invall=function(n) rgamma(n, shape = 42.3, rate = 1.99),
  omega_invall=function(n) rgamma(n, shape = 20.6, rate = 1.45),
  Vea=function(n) rbeta(n, shape1 = 49.3, shape2 = 87.4),
  Veb=function(n) rbeta(n, shape1 = 69.3, shape2 = 35.6),
  transit_time_mean = function(n) rtruncnorm(n, 1.3, 4.5, mean = 2.5, sd = sqrt(1/9)),
  transit_time_cv = function(n) rtruncnorm(n, 0.2, 0.6, mean = 0.3, sd = sqrt(1/36)),
  mult = function(n) exp(rnorm(n, log(3e-9), sqrt(1/300))),
  tau_ww = function(n) rgamma(n, shape = 40, rate = 48)
)

extract_param <- function(mcmc_obj, param_name, model_label) {
  df <- as.data.frame(as.matrix(mcmc_obj)[, param_name])
  colnames(df) <- "value"
  df$source <- model_label
  df$parameter <- param_name
  return(df)
}

##########compile full dataset
n_samples <- 60000
models <- list("Case-only" = case_mod, "Viral load-only" = ww_mod, "Combined" = com_mod)

all_param_dfs <- list()

for (param in names(param_model_map)) {
  # Simulate prior
  prior_df <- data.frame(
    value = prior_definitions[[param]](n_samples),
    source = "Prior",
    parameter = param
  )
  
  # Posteriors
  posteriors <- list()
  for (model_name in param_model_map[[param]]) {
    posteriors[[model_name]] <- extract_param(models[[model_name]], param, model_name)
  }
  
  all_param_dfs[[param]] <- bind_rows(prior_df, do.call(bind_rows, posteriors))
}

beta_df_all <- bind_rows(all_param_dfs)
beta_df_all$source <- factor(beta_df_all$source, levels = c("Combined","Viral load-only","Case-only","Prior"))
beta_df_all$parameter_label <- param_labels[beta_df_all$parameter]

WW_params <- c("mult", "transit_time_mean","tau_ww","beta","transit_time_cv")

plota=ggplot(beta_df_all %>% filter(parameter %in% WW_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plota

















###############plot the relative contribution to the total viral load
###To repeat this with new runs#####
###extract first

# Example: Summarizing P shedding by sexual activity group
# You can repeat for shed_A and shed_I as needed

# Extract matrices for each group
shedP_1 <- as.matrix(Comblist_finalg[, grep("shed_P\\[1,", varnames(Comblist_finalg))])
shedP_2 <- as.matrix(Comblist_finalg[, grep("shed_P\\[2,", varnames(Comblist_finalg))])
shedP_3 <- as.matrix(Comblist_finalg[, grep("shed_P\\[3,", varnames(Comblist_finalg))])

shedA_1 <- as.matrix(Comblist_finalg[, grep("shed_A\\[1,", varnames(Comblist_finalg))])
shedA_2 <- as.matrix(Comblist_finalg[, grep("shed_A\\[2,", varnames(Comblist_finalg))])
shedA_3 <- as.matrix(Comblist_finalg[, grep("shed_A\\[3,", varnames(Comblist_finalg))])

shedI_1 <- as.matrix(Comblist_finalg[, grep("shed_I\\[1,", varnames(Comblist_finalg))])
shedI_2 <- as.matrix(Comblist_finalg[, grep("shed_I\\[2,", varnames(Comblist_finalg))])
shedI_3 <- as.matrix(Comblist_finalg[, grep("shed_I\\[3,", varnames(Comblist_finalg))])


# Total shedding for each group (denominator)

shed1_total <- shedP_1 + shedA_1 + shedI_1
shed2_total <- shedP_2 + shedA_2 + shedI_2
shed3_total <- shedP_3 + shedA_3 + shedI_3
shedall_total <- shed1_total + shed2_total + shed3_total


# Relative contribution of each group
relg_1 <- shed1_total / shedall_total
relg_2 <- shed2_total/ shedall_total
relg_3 <- shed3_total / shedall_total

# Function to summarize posterior samples over time
summary_rel_group <- function(mat, group_label) {
  data.frame(
    time = 1:ncol(mat),
    Date=case_dat$Date,# time index
    median = apply(mat, 2, median, na.rm = TRUE),
    lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    group = group_label
  )
}

burn_in <- 30

# Remove first 30 timesteps from each matrix
relg_1_burned <- relg_1[, (burn_in):ncol(relg_1)]
relg_2_burned <- relg_2[, (burn_in):ncol(relg_2)]
relg_3_burned <- relg_3[, (burn_in):ncol(relg_3)]

# Combine summaries for plotting
relg_all <- dplyr::bind_rows(
  summary_rel_group(relg_1_burned, "Group 1 (low activity)"),
  summary_rel_group(relg_2_burned, "Group 2 (medium activity)"),
  summary_rel_group(relg_3_burned, "Group 3 (high activity)")
)

#relg_all=relg_all %>% mutate(Date=case_dat$Date)

plotshed=ggplot(relg_all, aes(x = Date, y = median, fill = group, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, linetype = 0) +
  labs(
    #title = "Relative contribution to shedding by sexual activity group",
    x = "Time",
    y = "Proportion of viral load shed"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position = "bottom")

plotshed

###########We now do the same for the total number of infected individuals
# Example extraction for each group across time (adapt to your variable names)

P1_total <- as.matrix(Comblist_finalg[, grep("P_total\\[1,", varnames(Comblist_finalg))])
P2_total <- as.matrix(Comblist_finalg[, grep("P_total\\[2,", varnames(Comblist_finalg))])
P3_total <- as.matrix(Comblist_finalg[, grep("P_total\\[3,", varnames(Comblist_finalg))])

A1_total <- as.matrix(Comblist_finalg[, grep("A_total\\[1,", varnames(Comblist_finalg))])
A2_total <- as.matrix(Comblist_finalg[, grep("A_total\\[2,", varnames(Comblist_finalg))])
A3_total <- as.matrix(Comblist_finalg[, grep("A_total\\[3,", varnames(Comblist_finalg))])

I1_total <- as.matrix(Comblist_finalg[, grep("I_total\\[1,", varnames(Comblist_finalg))])
I2_total <- as.matrix(Comblist_finalg[, grep("I_total\\[2,", varnames(Comblist_finalg))])
I3_total <- as.matrix(Comblist_finalg[, grep("I_total\\[3,", varnames(Comblist_finalg))])

# Sum across compartments per group
infected_1 <- P1_total + A1_total + I1_total
infected_2 <- P2_total + A2_total + I2_total
infected_3 <- P3_total + A3_total + I3_total


burn_in <- 30  # Number of initial time steps to exclude

# Remove burn-in columns
infected_1_burned <- infected_1[, (burn_in):ncol(infected_1)]
infected_2_burned <- infected_2[, (burn_in):ncol(infected_2)]
infected_3_burned <- infected_3[, (burn_in):ncol(infected_3)]

# Summary function for CrI and median over time
summarise_draws <- function(mat, group_label) {
  data.frame(
    time = 1:ncol(mat),
    Date=case_dat$Date,# time index
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    group = group_label
  )
}


# Apply summarization AFTER trimming
infected_df <- bind_rows(
  summarise_draws(infected_1_burned, "Group 1 (low activity)"),
  summarise_draws(infected_2_burned, "Group 2 (medium activity)"),
  summarise_draws(infected_3_burned, "Group 3 (high activity)")
)


# Plotting
# Plot 1: Relative contribution to viral shedding
plot_shed <- ggplot(relg_all, aes(x = Date, y = median, fill = group, color = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    #title = "Relative contribution to viral shedding by sexual activity group",
    x = "Date",
    y = "Proportion of total viral load"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Plot 2: Number of infectious individuals by group
plot_epi <- ggplot(infected_df, aes(x = Date, y = median, color = group, fill = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    #title = "Number of infectious individuals by sexual activity group",
    x = "Date",
    y = "Number of infectious individuals"
  ) +
  scale_color_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Combine plots with shared legend
mm=plot_shed + plot_epi + plot_layout(guides = "collect") & theme(legend.position = "bottom")



ggsave(
  filename = "D:/Mpox25output/Figures/relativeC.tiff",
  plot = mm,  # use your actual plot variable here
  width = 13,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)



# Extract matrices for each group
shedP_1 <- as.matrix(Comblist_finalg[, grep("shed_P\\[1,", varnames(Comblist_finalg))])
shedP_2 <- as.matrix(Comblist_finalg[, grep("shed_P\\[2,", varnames(Comblist_finalg))])
shedP_3 <- as.matrix(Comblist_finalg[, grep("shed_P\\[3,", varnames(Comblist_finalg))])

shedA_1 <- as.matrix(Comblist_finalg[, grep("shed_A\\[1,", varnames(Comblist_finalg))])
shedA_2 <- as.matrix(Comblist_finalg[, grep("shed_A\\[2,", varnames(Comblist_finalg))])
shedA_3 <- as.matrix(Comblist_finalg[, grep("shed_A\\[3,", varnames(Comblist_finalg))])

shedI_1 <- as.matrix(Comblist_finalg[, grep("shed_I\\[1,", varnames(Comblist_finalg))])
shedI_2 <- as.matrix(Comblist_finalg[, grep("shed_I\\[2,", varnames(Comblist_finalg))])
shedI_3 <- as.matrix(Comblist_finalg[, grep("shed_I\\[3,", varnames(Comblist_finalg))])



# Total shedding across all groups (denominator)
shedP_total <- shedP_1 + shedP_2 + shedP_3
shedA_total <- shedA_1 + shedA_2 + shedA_3
shedI_total <- shedI_1 + shedI_2 + shedI_3
shed_total <- shedP_total + shedA_total + shedI_total

# Total shedding for each group (denominator)
summary_rel <- function(mat, source_label) {
  data.frame(
    time = 1:ncol(mat),
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    source = source_label
  )
}

relcom_P=shedP_total/shed_total
relcom_A=shedA_total/shed_total
relcom_I=shedI_total/shed_total

rel_com_all <- bind_rows(
  summary_rel(relcom_P, "Pre-symptomatic (P)"),
  summary_rel(relcom_A, "Asymptomatic (A)"),
  summary_rel(relcom_I, "Symptomatic (I)")
)


plotrelcom=ggplot(rel_com_all, aes(x = time, y = median*100, fill = source)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), alpha = 0.2, color = NA) +
  geom_line(aes(color = source), size = 1) +
  labs(
    x = "Time", y = "Proportion of total viral load",
    title = "Relative contribution to viral shedding with 95% CrI(Combined)"
  ) +
  scale_fill_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                               "Asymptomatic (A)" = "orange",
                               "Symptomatic (I)" = "red")) +
  scale_color_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                                "Asymptomatic (A)" = "orange",
                                "Symptomatic (I)" = "red")) +
  theme_minimal(base_size = 10) +
  theme(legend.title = element_blank())

plotrelcom

# Sum across compartments per group
infected_P <- P1_total + P2_total + P3_total
infected_A <- A1_total + A2_total + A3_total
infected_I <- I1_total + I2_total + I3_total

# Summary function for CrI and median over time
summarise_draws <- function(mat, group_label) {
  data.frame(
    time = 1:ncol(mat),
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    group = group_label
  )
}

# Apply summarization
infected_All <- bind_rows(
  summarise_draws(infected_P, "Presymptomatic"),
  summarise_draws(infected_A, "Asymptomatic"),
  summarise_draws(infected_I, "Symptomatic")
)

plot_epi2 <- ggplot(infected_All , aes(x = time, y = median, color = group, fill = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    title = "Number of infectious individuals by infection stage",
    x = "Time (days)",
    y = "Number of infectious individuals"
  ) +
  scale_color_manual(values = c("orange","skyblue", "red")) +
  scale_fill_manual(values = c("orange","skyblue", "red")) +
  theme_minimal(base_size = 11)+ 
  theme(legend.position = "none")

vv=plot_epi2+plotrelcom
vv


ggsave(
  filename = "D:/Mpox25output/Figures/relv2.tiff",
  plot = vv,  # use your actual plot variable here
  width = 13,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


#####################Now geenerate a panel with each##############
# 1. Function to summarize posterior samples over time
summary_rel_group_comp <- function(mat, group_label, comp_label) {
  data.frame(
    time = 1:ncol(mat),
    median = apply(mat, 2, median, na.rm = TRUE),
    lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    group = group_label,
    compartment = comp_label
  )
}

# 2. Create summaries (replace these with your real matrices)
# For example: rel_P1 <- shedP_1 / (shedP_1 + shedP_2 + shedP_3), etc.

rel_P1 <- shedP_1 / (shedP_1 + shedP_2 + shedP_3)
rel_P2 <- shedP_2 / (shedP_1 + shedP_2 + shedP_3)
rel_P3 <- shedP_3 / (shedP_1 + shedP_2 + shedP_3)

rel_A1 <- shedA_1 / (shedA_1 + shedA_2 + shedA_3)
rel_A2 <- shedA_2 / (shedA_1 + shedA_2 + shedA_3)
rel_A3 <- shedA_3 / (shedA_1 + shedA_2 + shedA_3)

rel_I1 <- shedI_1 / (shedI_1 + shedI_2 + shedI_3)
rel_I2 <- shedI_2 / (shedI_1 + shedI_2 + shedI_3)
rel_I3 <- shedI_3 / (shedI_1 + shedI_2 + shedI_3)


rel_P1_sum <- summary_rel_group_comp(rel_P1, "Group 1", "Pre-symptomatic (P)")
rel_P2_sum <- summary_rel_group_comp(rel_P2, "Group 2", "Pre-symptomatic (P)")
rel_P3_sum <- summary_rel_group_comp(rel_P3, "Group 3", "Pre-symptomatic (P)")

rel_A1_sum <- summary_rel_group_comp(rel_A1, "Group 1", "Asymptomatic (A)")
rel_A2_sum <- summary_rel_group_comp(rel_A2, "Group 2", "Asymptomatic (A)")
rel_A3_sum <- summary_rel_group_comp(rel_A3, "Group 3", "Asymptomatic (A)")

rel_I1_sum <- summary_rel_group_comp(rel_I1, "Group 1", "Symptomatic (I)")
rel_I2_sum <- summary_rel_group_comp(rel_I2, "Group 2", "Symptomatic (I)")
rel_I3_sum <- summary_rel_group_comp(rel_I3, "Group 3", "Symptomatic (I)")

# 3. Combine all summaries
rel_all_comp_df <- bind_rows(
  rel_P1_sum, rel_P2_sum, rel_P3_sum,
  rel_A1_sum, rel_A2_sum, rel_A3_sum,
  rel_I1_sum, rel_I2_sum, rel_I3_sum
)


rel_all_comp_df$compartment=factor(rel_all_comp_df$compartment,levels=c("Pre-symptomatic (P)","Asymptomatic (A)",
                                                                        "Symptomatic (I)"))


# 4. Plotting
plotrel<-ggplot(rel_all_comp_df, aes(x = time, y = median, fill = group, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  facet_wrap(~ compartment, ncol = 1, scales = "free_y") +
  labs(
    title = "Relative contribution to viral load by group and compartment",
    x = "Time",
    y = "Proportion of viral load shed"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

plotrel

ggsave(
  filename = "D:/Mpox25output/Figures/plotrel.tiff",
  plot = plotrel,  # use your actual plot variable here
  width = 9,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


##############model input data###############
##1. case data ################
case_dat=read_excel("Data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat %>% select(Date,total_cases)
#Cumulat_inc=case_dat$Cumulat_Inc

##2. WW data #########################
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
ww_std = ww_dat %>% select(date,log10_cp_per_person_per_day) %>% rename(Date=date)




# Merge the two with full_join
merged_df <- full_join(cases_obsb, ww_std, by = c("Date"))

# Reshape to long format
df_long <- merged_df %>%
  pivot_longer(cols = c(total_cases, log10_cp_per_person_per_day), names_to = "type", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(type) %>%
  mutate(value_norm = (value - min(value, na.rm = TRUE)) / (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ungroup()

# Label for legend
df_long$type <- factor(df_long$type, 
                       levels = c("total_cases", "log10_cp_per_person_per_day"),
                       labels = c("Case report", "Wastewater concentration"))



custom_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

Plot=ggplot(df_long, aes(x = Date, y = value_norm, color = type)) +
  geom_line(size = 0.2) +
  scale_color_manual(values = c("brown", "steelblue")) +
  labs(x = "Date", y = "Normalized value", color = "Data type") +
  custom_theme

Plot























