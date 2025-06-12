########Troubleshooting current outputs
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
options(max.print = 10000)

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
print(colnames(as.matrix(Combined_cas)))

######generate the model fit
mcmc_matrixallcas<-as.matrix(Combined_cas)

# This gets all columns for total_new_cases
cols <- grep("log10_concb", colnames(mcmc_matrixallcas))
# Extract full MCMC samples for those timepoints
v_test <- as.data.frame(mcmc_matrixallcas[, cols])


hist(v_test[[which(colnames(v_test) == "log10_concb[20]")]],
     breaks = 50,
     main = "Posterior of log10_concb",
     xlab = "Predicted vload",
     col = "skyblue")


# Function to compute the median and 95% credible interval
summary_median_CI <- function(samples) {
  med <- apply(samples, 2, median)
  lower <- apply(samples, 2, quantile, probs = 0.025)
  upper <- apply(samples, 2, quantile, probs = 0.975)
  summary_table <- cbind(median = med, lower_95_CI = lower, upper_95_CI = upper)
  return(summary_table)
}

posterior_summarycas <- summary_median_CI(mcmc_matrixallcas)
total_cas_summary <- as.data.frame(posterior_summarycas[grep("cases_pred", rownames(posterior_summarycas)), ])
total_cas_ww_summary <- as.data.frame(posterior_summarycas[grep("log10_concb", rownames(posterior_summarycas)), ])
total_cas_ww_summaryb <- as.data.frame(posterior_summarycas[grep("log10_conc", rownames(posterior_summarycas)), ])


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

plot_datacas_wwb <- data.frame(
  time = 1:nrow(ww_std),                    
  observed = ww_obs,                    # observed cases
  median_fit = total_cas_ww_summaryb$median,          # model median fit
  lower_ci = total_cas_ww_summaryb$lower_95_CI,          # lower 95% CI
  upper_ci = total_cas_ww_summaryb$upper_95_CI          # upper 95% CI
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


# Viral load – Case only model
plot_wwcasefit <- ggplot(plot_dataww_cases, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title =  "Fit vs. observed cases (WW-only model)"
  ) +
  custom_theme


plot_wwcasefitb <- ggplot(plot_dataww_casesb, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title =  "Fit vs. observed cases (WW-only model)"
  ) +
  custom_theme





ggsave(
  filename = "D:/Mpox25output/Figures/plot_wwcasefit.tiff",
  plot =plot_wwcasefit,  # use your actual plot variable here
  width = 13,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)






