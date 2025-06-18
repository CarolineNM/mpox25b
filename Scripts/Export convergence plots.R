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
load("D:/mpox25output/Combined_WWc.RData")
load("D:/mpox25output/Combined_casb.RData")
load("D:/mpox25output/Combined_cas.RData")
#load("D:/mpox25output/Case_modlstfinalc.RData")

######What did we save in each model?
options(max.print = 10000)  # or higher depending on your needs
print(colnames(as.matrix(Combined_casb)))
cat(colnames(as.matrix(Combined_casb)), sep = "\n")
View(as.data.frame(colnames(as.matrix(Comblist_finalg))))

######generate the model fit
mcmc_matrixallWW<-as.matrix(Combined_WWc)
mcmc_matrixallcas<-as.matrix(Combined_casb)
mcmc_matrixallcasb<-as.matrix(Combined_cas)
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
posterior_summarycasb <- summary_median_CI(mcmc_matrixallcasb)
posterior_summarycom <- summary_median_CI(mcmc_matrixallcom)
posterior_summaryww <- summary_median_CI(mcmc_matrixallWW)

total_cas_summary <- as.data.frame(posterior_summarycas[grep("cases_pred", rownames(posterior_summarycas)), ])
total_casb_summary <- as.data.frame(posterior_summarycas[grep("mu_nb", rownames(posterior_summarycas)), ])
total_cas_ww_summary <- as.data.frame(posterior_summarycasb[grep("log10_conc", rownames(posterior_summarycasb)), ])

total_ww_summary <- as.data.frame(posterior_summaryww[grep("ww_pred", rownames(posterior_summaryww)), ])
total_ww_summaryb <- as.data.frame(posterior_summaryww[grep("^log10_conc\\[", rownames(posterior_summaryww)), ])
total_ww_cas_summary <- as.data.frame(posterior_summaryww[grep("mu_nb", rownames(posterior_summaryww)), ])

total_com_ww_summary <- as.data.frame(posterior_summarycom[grep("ww_pred", rownames(posterior_summarycom)), ])
total_com_ww_summaryb <- as.data.frame(posterior_summaryww[grep("^log10_conc\\[", rownames(posterior_summaryww)), ])
total_com_cas_summary <- as.data.frame(posterior_summarycom[grep("cases_pred", rownames(posterior_summarycom)), ])
total_com_cas_summaryb <- as.data.frame(posterior_summarycom[grep("mu_nb", rownames(posterior_summarycom)), ])


###########generate data for plotting posterior predictive plots
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
  median_fit = total_ww_cas_summary$median,         # model median fit
  lower_ci = total_ww_cas_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_ww_cas_summary$upper_95_CI         # upper 95% CI
)


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
  observed = ww_obs,                    # observed cases
  median_fit = total_com_ww_summary$median,          # model median fit
  lower_ci = total_com_ww_summary$lower_95_CI,          # lower 95% CI
  upper_ci = total_com_ww_summary$upper_95_CI          # upper 95% CI
)

plot_datcomwwb <- data.frame(
  time = 1:nrow(ww_std),                    
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
plot_casefit <- ggplot(plot_datcas, aes(x = time)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Time", y = "Reported mpox cases",
    title = "Fit vs. observed cases (Case-only model)"
  ) +
  custom_theme

# Cases – Case-only model
plot_casefitb <- ggplot(plot_datcasb, aes(x = time)) +
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

plot_comwwfitb <- ggplot(plot_datcomwwb, aes(x = time)) +
  plot_geom +
  ylim(0, 15)+
  labs(
    x = "Time", y = "Reported viral load",
    title = "Fit vs. observed viral load (Combined model)"
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


plot_comcasefitb <- ggplot(plot_datcomcasb, aes(x = time)) +
  plot_geom +
  ylim(0, 11)+
  labs(
    x = "Time", y = "Reported mpox cases",
    title = "Fit vs. observed cases (Combined model)"
  ) +
  custom_theme



##########Predictive fits

plot_comcasefitb
plot_comcasefit
plot_comwwfitb
plot_comwwfit
plot_wwcasefit
plot_wwfitb
plot_wwfit
plot_casewwfit
plot_casefitb
plot_casefit


ZZ=(plot_casefit+plot_comcasefit+plot_casewwfit)/(plot_casefitb+plot_comcasefitb+plot_casewwfit)
ZZ

yy=(plot_wwfit+plot_comwwfit+plot_wwcasefit)/(plot_wwfitb+plot_comwwfitb+plot_wwcasefit)
yy


ggsave(
  filename = "D:/Mpox25output/Figures/modfit.tiff",
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


# Total shedding across all groups (denominator)
shedP_total <- shedP_1 + shedP_2 + shedP_3
shedA_total <- shedA_1 + shedA_2 + shedA_3
shedI_total <- shedI_1 + shedI_2 + shedI_3
shed_total <- shedP_total + shedA_total + shedI_total

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
    median = apply(mat, 2, median, na.rm = TRUE),
    lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    group = group_label
  )
}

# Combine summaries for plotting
relg_all <- dplyr::bind_rows(
  summary_rel_group(relg_1, "Group 1 (low activity)"),
  summary_rel_group(relg_2, "Group 2 (medium activity)"),
  summary_rel_group(relg_3, "Group 3 (high activity)")
)



plotshed=ggplot(relg_all, aes(x = time, y = median, fill = group, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, linetype = 0) +
  labs(
    title = "Relative contribution to shedding by sexual activity group",
    x = "Time",
    y = "Proportion of viral load shed"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position = "bottom")



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
infected_df <- bind_rows(
  summarise_draws(infected_1, "Group 1 (low activity)"),
  summarise_draws(infected_2, "Group 2 (medium activity)"),
  summarise_draws(infected_3, "Group 3 (high activity)")
)

# Plotting
# Plot 1: Relative contribution to viral shedding
plot_shed <- ggplot(relg_all, aes(x = time, y = median, fill = group, color = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    title = "Relative contribution to viral shedding by sexual activity group",
    x = "Time (days)",
    y = "Proportion of total viral load"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# Plot 2: Number of infectious individuals by group
plot_epi <- ggplot(infected_df, aes(x = time, y = median, color = group, fill = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(
    title = "Number of infectious individuals by sexual activity group",
    x = "Time (days)",
    y = "Number of infectious individuals"
  ) +
  scale_color_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a", "#377eb8")) +
  theme_minimal(base_size = 11) +
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


##########Now plot the contribution from each group
##################Summarising correlation across posteriors#############
mcmc_matrixallcom <- as.matrix(Comblist_finalg)

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
ww_pred_median <- posterior_summarycom[grep("^log10_conc_all\\[", rownames(posterior_summarycom)), "median"]
cases_pred_median <- posterior_summarycom[grep("mu_nb", rownames(posterior_summarycom)), "median"]

# Calculate correlation trajectory
correlation_df <- sliding_correlation(ww_pred_median, cases_pred_median)
correlation_df$frac_time <- correlation_df$window_size  / max(correlation_df$window_size)


# Estimate transmission peak day from predicted cases
transmission_peak_day <- which.max(cases_pred_median)
total_days <- length(cases_pred_median)
transmission_peak_frac <- transmission_peak_day / total_days

# Enhanced plot
finalplt=ggplot(correlation_df, aes(x = frac_time, y = correlation)) +
  # Optional shaded pre- and post-peak phases
  annotate("rect", xmin = 0, xmax = transmission_peak_frac, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "gray70") +
  annotate("rect", xmin = transmission_peak_frac, xmax = 1, ymin = -Inf, ymax = Inf,
           alpha = 0.05, fill = "gray90") +
  
  # Correlation line
  geom_line(color = "grey", size = 1.2) +
  geom_point(color = "grey", size = 2) +
  
  # Vertical line at transmission peak
  geom_vline(xintercept = transmission_peak_frac, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = transmission_peak_frac, y = max(correlation_df$correlation, na.rm = TRUE),
           label = "Peak transmission", vjust = -1, hjust = -0.1, size = 4) +
  
  # Axis and labels
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0.9, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Temporal correlation between viral load and predicted cases",
    subtitle = "Combined model output, with dashed line at peak infectiousness",
    x = "Fraction of time series used",
    y = "Pearson correlation (log10 viral load vs. predicted cases)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 10)),
    panel.grid.minor = element_blank()
  )

finalplt



####counterfactual overlays
# Example structure if you compute multiple correlations:
correlation_all <- bind_rows(
  sliding_correlation(ww_pred_baseline, case_pred_baseline) %>% mutate(scenario = "Baseline"),
  sliding_correlation(ww_pred_alt, case_pred_alt) %>% mutate(scenario = "Faster transit"),
  ...
)

ggplot(correlation_all, aes(x = frac_time, y = correlation, color = scenario)) +
  geom_line(size = 1.1) +
  labs(x = "Fraction of time series", y = "Pearson correlation", color = "Scenario")







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











