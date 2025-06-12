#########Call the libraries
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
load("D:/mpox25output/Combined_cas.RData")
load("D:/mpox25output/Combined_WWb.RData")
load("D:/mpox25output/Comblist_finale.RData")
load("D:/mpox25output/WW_modlstfinale.RData")


######What did we save in each model?

#options(max.print = 10000)  # or higher depending on your needs
#print(colnames(as.matrix(Combined_cas)))
####From the combined cas model we are interest in both the case fit and also the viral load fit(log10_con)
####We try both (log10_con) and log_10conb
###Can we use cases to predict viral load?
###We also plot the case_pred from here
###From model2 we have;
##Can we use viral load to predict cases?
####we plot total_new_cases
##we also plot the ww_pred from here


######generate the model fit
mcmc_matrixallcas<-as.matrix(Combined_cas)
mcmc_matrixallWW<-as.matrix(Combined_WWb)
mcmc_matrixallcom<-as.matrix(Comblist_finale)

# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summaryww <- summary_median_CI(mcmc_matrixallWW)
posterior_summarywwb <- summary_median_CI(mcmc_matrixallWWb)###
posterior_summarycas <- summary_median_CI(mcmc_matrixallcas)
posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)

#posterior_summarycom[grep("mult|tau_ww|transit_time_mean|transit_time_cv", rownames(posterior_summarycom)), ]
####generate the posterior fits
total_ww_summary <- as.data.frame(posterior_summaryww[grep("ww_pred", rownames(posterior_summaryww)), ])
total_ww_summaryb <- as.data.frame(posterior_summarywwb[grep("ww_pred", rownames(posterior_summarywwb)), ])
total_ww_cas_summary <- as.data.frame(posterior_summaryww[grep("total_new_cases", rownames(posterior_summaryww)), ])
total_cas_summary <- as.data.frame(posterior_summarycas[grep("cases_pred", rownames(posterior_summarycas)), ])
total_cas_ww_summary <- as.data.frame(posterior_summarycas[grep("log10_concb", rownames(posterior_summarycas)), ])
total_com_ww_summary <- as.data.frame(posterior_summarycom[grep("ww_pred", rownames(posterior_summarycom)), ])
total_com_cas_summary <- as.data.frame(posterior_summarycom[grep("cases_pred", rownames(posterior_summarycom)), ])

###########generate data for plotting
plot_dataww <- data.frame(
  time = 1:nrow(ww_std),                    
  observed = ww_obs,                    # observed cases
  median_fit = total_ww_summary$median,          # model median fit
  lower_ci = total_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_summary$upper_95_CI          # upper 95% CI
)

plot_datawwb <- data.frame(
  time = 1:nrow(ww_std),                    
  observed = ww_obs,                    # observed cases
  median_fit = total_ww_summaryb$median,          # model median fit
  lower_ci = total_ww_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_summaryb$upper_95_CI          # upper 95% CI
)



burn_in_timesteps <- 30
total_ww_cas_summaryb <- total_ww_cas_summary[(burn_in_timesteps + 1):nrow(total_ww_cas_summary), ]

plot_dataww_cases <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_ww_cas_summaryb$median,         # model median fit
  lower_ci = total_ww_cas_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_cas_summaryb$upper_95_CI         # upper 95% CI
)


plot_datcas <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_cas_summary$median,         # model median fit
  lower_ci = total_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_summary$upper_95_CI          # upper 95% CI
)


plot_datacas_ww <- data.frame(
  time = 1:nrow(ww_std),                    
  observed = ww_obs,                    # observed cases
  median_fit = total_cas_ww_summary$median,          # model median fit
  lower_ci = total_cas_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_ww_summary$upper_95_CI          # upper 95% CI
)


plot_datcomcas <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  Date=case_dat$Date,
  observed = cases_obsb,                    # observed cases
  median_fit = total_com_cas_summary$median,         # model median fit
  lower_ci = total_com_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_cas_summary$upper_95_CI          # upper 95% CI
)


plot_datcomww <- data.frame(
  time = 1:nrow(ww_std),                    
  observed = ww_obs,                    # observed cases
  median_fit = total_com_ww_summary$median,          # model median fit
  lower_ci = total_com_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_ww_summary$upper_95_CI          # upper 95% CI
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

# Viral load – WW-only model
plot_wwfit <- ggplot(plot_dataww, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title = "Fit vs. observed viral load (WW-only model)"
  ) +
  custom_theme

plot_wwfitb <- ggplot(plot_datawwb, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title = "Fit vs. observed viral load (WW-only model)"
  ) +
  custom_theme


# Viral load – Case only model
plot_wwcasefit <- ggplot(plot_datacas_ww, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title = "Fit vs. observed viral load (Case-only model)"
  ) +
  custom_theme

# Viral load – Combined model
plot_comwwfit <- ggplot(plot_datcomww, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title = "Fit vs. observed viral load (Combined model)"
  ) +
  custom_theme


# Cases – Case-only model
plot_casefit <- ggplot(plot_datcas, aes(x = time)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Time", y = "Reported mpox cases",
    title = "Fit vs. observed cases (Case-only model)"
  ) +
  custom_theme


plot_casewwfit <- ggplot(plot_dataww_cases, aes(x = time)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Time", y = "Reported mpox cases",
    title = "Fit vs. observed cases (WW-only model)"
  ) +
  custom_theme


# Cases – Combined model
plot_comcasefit <- ggplot(plot_datcomcas, aes(x = time)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Time", y = "Reported mpox cases",
    title = "Fit vs. observed cases (Combined model)"
  ) +
  custom_theme

# Combine plots using patchwork
x=(plot_casefit|plot_comcasefit|plot_casewwfit)/(plot_wwfitb|plot_comwwfit|plot_wwcasefit)
x

ggsave(
  filename = "D:/Mpox25output/Figures/modfit.tiff",
  plot = x,  # use your actual plot variable here
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
  time = 1:nrow(case_dat),                    # time index
  median_fit = prev_wwsummaryb$median,          # model median fit
  lower_ci = prev_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casprev <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  median_fit = prev_cassummaryb$median,          # model median fit
  lower_ci = prev_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = prev_cassummaryb$upper_95_CI         # upper 95% CI
)


plot_comprev <- data.frame(
  time = 1:nrow(case_dat),                    # time index
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

plot_prev_all <- ggplot(plot_allprev, aes(x = time, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Time",
    y = "Active prevalence",
    title = "Model-predicted prevalence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
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
  time = 1:nrow(case_dat),                    # time index
  median_fit = Cuminc_wwsummaryb$median,          # model median fit
  lower_ci = Cuminc_wwsummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_wwsummaryb$upper_95_CI         # upper 95% CI
)


plot_casCuminc <- data.frame(
  time = 1:nrow(case_dat),                    # time index
  median_fit = Cuminc_cassummaryb$median,          # model median fit
  lower_ci = Cuminc_cassummaryb$lower_95_CI,          # lower 95% CI
  upper_ci = Cuminc_cassummaryb$upper_95_CI         # upper 95% CI
)

plot_comCuminc <- data.frame(
  time = 1:nrow(case_dat),                    # time index
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


plot_Cuminc_all <- ggplot(plot_allCuminc, aes(x = time, y = median_fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Time",
    y = "Cumulative Incidence",
    title = "Model-predicted cumulative incidence by model type"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  scale_fill_manual(values = c("WW-only" = "blue", "Case-only" = "darkgreen", "Combined" = "Orange")) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5))



y=plot_prev_all+plot_Cuminc_all



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



###############plot the relative contribution to the total viral load
###To repeat this with new runs#####
###extract first
shedP_com <- as.matrix(Comblist_finale[, grep("shed_P\\[", varnames(Comblist_finale))])
shedA_com <- as.matrix(Comblist_finale[, grep("shed_A\\[", varnames(Comblist_finale))])
shedI_com <- as.matrix(Comblist_finale[, grep("shed_I\\[", varnames(Comblist_finale))])
shedP_WW <- as.matrix(WW_modlstfinale[, grep("shed_P\\[", varnames(WW_modlstfinale))])
shedA_WW <- as.matrix(WW_modlstfinale[, grep("shed_A\\[", varnames(WW_modlstfinale))])
shedI_WW <- as.matrix(WW_modlstfinale[, grep("shed_I\\[", varnames(WW_modlstfinale))])

shedcom_total <- shedP_com + shedA_com + shedI_com

relcom_P <- shedP_com / shedcom_total
relcom_A <- shedA_com / shedcom_total
relcom_I <- shedI_com / shedcom_total


summary_rel <- function(mat, source_label) {
  data.frame(
    time = 1:ncol(mat),
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    source = source_label
  )
}

rel_com_all <- bind_rows(
  summary_rel(relcom_P, "Pre-symptomatic (P)"),
  summary_rel(relcom_A, "Asymptomatic (A)"),
  summary_rel(relcom_I, "Symptomatic (I)")
)

shedww_total <- shedP_WW + shedA_WW + shedI_WW

relww_P <- shedP_WW /shedww_total
relww_A <- shedA_WW / shedww_total
relww_I <- shedI_WW / shedww_total

rel_ww_all <- bind_rows(
  summary_rel(relww_P, "Pre-symptomatic (P)"),
  summary_rel(relww_A, "Asymptomatic (A)"),
  summary_rel(relww_I, "Symptomatic (I)")
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


plotrelww=ggplot(rel_ww_all, aes(x = time, y = median*100, fill = source)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), alpha = 0.2, color = NA) +
  geom_line(aes(color = source), size = 1) +
  labs(
    x = "Time", y = "Proportion of total viral load",
    title = "Relative contribution to viral shedding with 95% CrI(WW_only)"
  ) +
  scale_fill_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                               "Asymptomatic (A)" = "orange",
                               "Symptomatic (I)" = "red")) +
  scale_color_manual(values = c("Pre-symptomatic (P)" = "skyblue",
                                "Asymptomatic (A)" = "orange",
                                "Symptomatic (I)" = "red")) +
  theme_minimal(base_size = 10) +
  theme(legend.title = element_blank())

plotrelww


xx=plotrelcom+plotrelww


ggsave(
  filename = "D:/Mpox25output/Figures/relv.tiff",
  plot = xx,  # use your actual plot variable here
  width = 13,
  height = 6,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)



##################Summarising correlation across posteriors#############
mcmc_matrixallcom <- as.matrix(Comblist_finale)

# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)

# Function to compute sliding window correlation
sliding_correlation <- function(pred_ww, pred_cases, window_sizes = seq(7, length(pred_ww), by = 7)) {
  result <- data.frame(
    window_size = integer(),
    correlation = numeric()
  )
  
  for (w in window_sizes) {
    # Subset both series up to the current window size
    ww_subset <- pred_ww[1:w]
    cases_subset <- pred_cases[1:w]
    
    # Compute correlation if enough variance
    if (sd(ww_subset) > 0 && sd(cases_subset) > 0) {
      corr <- cor(ww_subset, cases_subset, method = "pearson")
    } else {
      corr <- NA
    }
    
    result <- rbind(result, data.frame(window_size = w, correlation = corr))
  }
  
  return(result)
}



# Example usage:
# Replace these with your actual model output
ww_pred_median <- posterior_summarycom[grep("ww_pred", rownames(posterior_summarycom)), "median"]
cases_pred_median <- posterior_summarycom[grep("cases_pred", rownames(posterior_summarycom)), "median"]

# Calculate correlation trajectory
correlation_df <- sliding_correlation(ww_pred_median, cases_pred_median)
correlation_df$frac_time <- correlation_df$window_size  / max(correlation_df$window_size)
# Plot
XX=ggplot(correlation_df, aes(x = frac_time, y = correlation)) +
  geom_line(color = "blue", size = 1) +
  geom_point() +
  #geom_vline(xintercept = transmission_peak_day / total_days, linetype = "dashed") +
  labs(
    title = "Correlation of predicted cases and viral load",
    x = "Fraction of time series used",
    y = "Pearson correlation "
  ) +
  theme_minimal()


ggsave(
  filename = "D:/Mpox25output/Figures/Corplot.tiff",
  plot = XX,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

summary(rbeta(1000,7,3))











