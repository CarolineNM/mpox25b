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
#Cumulat_inc=case_dat$Cumulat_Inc
##2. WW data #########################
ww_dat=read_excel("Data/case_data_V2.xlsx",sheet="dailyWW")
ww_std = ww_dat %>% select(log10_cp_per_person_per_day) #####standardised WW data
ww_obs = as.numeric(unlist(ww_std$log10_cp_per_person_per_day))

#######load the three outputs
load("D:/mpox25output/Comb_castestf.RData")
load("D:/mpox25output/Comb_finaltestf.RData")
load("D:/mpox25output/Comb_WWtestf.RData")

######generate the model fit
mcmc_matrixallWW<-as.matrix(Comb_WWtestf)
mcmc_matrixallcas<-as.matrix(Comb_castestf) ###most recent version
mcmc_matrixallcom<-as.matrix(Comb_finaltestf)

#E0=8(3-14);P0<-2(0-5);A10<-1(0-3);I10<-2(0-6)
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
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 14, face = "bold", colour = "black"),
    axis.text = element_text(size = 11, colour = "black")
  )

# Updated shared plot elements for better visibility
plot_geom <- list(
  geom_point(aes(y = observed), color = "black", size = 1.5, alpha = 0.8),
  geom_line(aes(y = median_fit), color = "blue4", size = 1.2), # thicker line
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.4) # slightly darker ribbon
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
    title = "(Case-only model)"
  ) +
  custom_theme

# Cases – ww-only model
plot_casewwfit <- ggplot(plot_dataww_cases, aes(x = Date)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Date", y = "Fit vs. observed cases",
    title = "(WW-only model)"
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
  ylim(0, 8)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(WW-only model)"
  ) +
  custom_theme

plot_wwfitb <- ggplot(plot_datawwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 8)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(WW-only model)"
  ) +
  custom_theme


# Viral load – Case only model
plot_wwcasefit <- ggplot(plot_datacas_ww, aes(x = Date)) +
  plot_geom +
  coord_cartesian(ylim = c(NA, 8))+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8))+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Case-only model)"
  ) +
  custom_theme

plot_wwcasefitb <- ggplot(plot_datacas_wwb, aes(x = Date)) +
  plot_geom +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8))+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Case-only model)"
  ) +
  custom_theme


# Viral load – Combined model
plot_comwwfit <- ggplot(plot_datcomww, aes(x = Date)) +
  plot_geom +
  ylim(0, 8)+
  labs(
    x = "Date", y = "Fit vs. observed viral load",
    title = "(Combined model)"
  ) +
  custom_theme

plot_comwwfitb <- ggplot(plot_datcomwwb, aes(x = Date)) +
  plot_geom +
  ylim(0, 8)+
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
    title = "(Combined model)"
  ) +
  custom_theme

ZZ=(plot_casefit+plot_comcasefit+plot_casewwfitb)/(plot_casefitb+plot_comcasefitb+plot_casewwfit)
mm=ZZ+plot_layout(guides = "collect") & theme(legend.position = "bottom")


yy=(plot_wwfit+plot_comwwfit+plot_wwcasefitb)/(plot_wwfitb+plot_comwwfitb+plot_wwcasefit)
yy


# Make top row: Case-only, Combined, WW-only
top_row <- plot_casefitb + plot_comcasefitb + plot_casewwfit

# Make bottom row: Case-only, Combined, WW-only (matching order)
bottom_row <- plot_wwcasefit + plot_comwwfitb + plot_wwfitb

# Combine with patchwork, keeping consistent layout
mainfigure1 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "D:/Mpox25output/Figures/mainfigure1.pdf",
  plot = mainfigure1,
  width = 12,       # adjust width for your manuscript column/page size
  height = 7,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

# ggsave(
#   filename = "D:/Mpox25output/Figures/predfit2.tiff",
#   plot = yy,  # use your actual plot variable here
#   width = 13,
#   height = 6,
#   dpi = 300,
#   units = "in",
#   device = "tiff",
#   compression = "lzw"
# )

########generate prevalence and cumulative incidence
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
  filename = "D:/Mpox25output/Figures/prevcum2.tiff",
  plot = y,  # use your actual plot variable here
  width = 13,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


###################Mraginal posterior plots##########
####################refine this######################
case_mod  <- Comb_castestf
ww_mod   <- Comb_WWtestf
com_mod <- Comb_finaltestf

###Define parameter to model mapping
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
  transit_time_cv = c("Viral load-only", "Combined"),
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
  transit_time_mean = function(n) rtruncnorm(n, a = 0.1, b = 10,mean = 2.5, sd = sqrt(1 / 0.25)),
  transit_time_cv = function(n) rtruncnorm(n, a = 0.1, b = 1, mean = 0.3, sd = sqrt(1 / 3)),
  mult = function(n) exp(rtruncnorm(n, a = log(1e-9), b = log(1e-8), mean = log(3e-9), sd = sqrt(1 / 2.5))),
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
n_samples <- 40000
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

WW_params <- c("mult", "transit_time_mean","transit_time_cv","tau_ww")

plota=ggplot(beta_df_all %>% filter(parameter %in% WW_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plota

ggsave(
  filename = "U:/Mpox25output/marginalA.tiff",
  plot = plota,  # use your actual plot variable here
  width = 12,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)


case_params <- c("beta", "kappa","report_frac","m")

plotb=ggplot(beta_df_all %>% filter(parameter %in% case_params), 
             aes(x = value, y = source, fill = source)) +
  geom_density_ridges(alpha = 0.7, 
                      quantile_lines = TRUE, 
                      quantiles = 0.5,
                      scale = 0.9,
                      rel_min_height = 0.0005) +
  facet_wrap(~ parameter_label, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )
plotb

ggsave(
  filename = "U:/Mpox25output/marginalB.tiff",
  plot = plotb,  # use your actual plot variable here
  width = 12,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

################Trace plots##################
case_mod  <- Combined_castestc
ww_mod   <- Comb_WWtestd
com_mod <- Comb_finaltestd

#WW_params <- c("mult", "transit_time_mean","transit_time_cv","tau_ww")
params_to_plot <- c("mult", "transit_time_mean","transit_time_cv","tau_ww")

# Convert mcmc.list to long dataframe
mcmc_df <- as.data.frame(do.call(rbind, lapply(1:length(com_mod), function(i) {
  df <- as.data.frame(com_mod[[i]][, params_to_plot])
  df$Chain <- paste0("Chain_", i)
  df$Iteration <- 1:nrow(df)
  df
})))

# Reshape to long format
mcmc_long <- mcmc_df %>%
  pivot_longer(cols = all_of(params_to_plot), names_to = "Parameter", values_to = "Value")


param_labels <- c(
  transit_time_mean = "Mean transit time (days)",
  transit_time_cv = "Sd of transit time",
  mult = "Scaling factor for viral load",
  tau_ww = "Precision of WW likelihood (τ)"
)


mcmc_long$ParameterLabel <- param_labels[mcmc_long$Parameter]

# Plot
ww_traceplot<-ggplot(mcmc_long, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.6) +
  #facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
  facet_wrap(~ ParameterLabel, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(title = "Posterior Distributions (Combined)")


case_densplot<-ggplot(mcmc_long, aes(x = Value, fill = Chain)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ ParameterLabel, scales = "free", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(title = "Posterior Distributions (Combined)", x = "Value", y = "Density") +
  theme(legend.position = "right")

case_densplot

mm=ww_traceplot+case_densplot
comb=mm+plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "U:/Mpox25output/traceplotsB.tiff",
  plot = comb,  # use your actual plot variable here
  width = 14,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

###############The Gelman rubins diagnostics and ESS###############
params_of_interest <- c("transit_time_cv","transit_time_mean", "log_mult", "tau_ww")

ww_subset <- com_mod[, params_of_interest]

gelman <- gelman.diag(ww_subset, multivariate = FALSE)
ess <- effectiveSize(ww_subset)


# Create summary table
WW_summary_table <- data.frame(
  Parameter = rownames(gelman$psrf),
  Rhat = round(gelman$psrf[, 1], 3),
  ESS = round(ess[rownames(gelman$psrf)], 0)
)

print(WW_summary_table)


>>>>>>> e42687c2552bdaa868433054e0b73647b88823a3






