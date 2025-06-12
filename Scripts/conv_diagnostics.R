#Parameters
Parameters[,"transit_time_mean"] main="Mean transit time in sewer"
Parameters[, "transit_time_cv"] main="Standard deviation of transit mean time"
Parameters[, "mult"] main="Scaling factor of viral load"
Parameters[, "tau_ww"] main="Precision of the dnorm likelihood"
Parameters[, "beta"] main="Transmission parameter"
Parameters[, "kappa"] main="Mixing probability"
Parameters[, "phi"] main="Negative binomial dispersion parmeter"
Parameters[, "report_frac"] main="reporting fraction"
Parameters[, "m"] main="Proportion of Asymptomatic fraction"


######start with beta and kappa which is shared across all the three models
set.seed(123)
vloadparams <- c("beta","kappa")  # Add/remove as needed
# Sample size
n_samples <- 30000

#beta~dnorm(0.8, 100) T(0,1) ##to change this to sd=0.1
#kappa ~ dbeta(40, 2)  ### Mean ~0.95, 95% CI ≈ [0.85, 0.995]

# 2. Now create the list with all priors
prior_samples <- list(
  beta = rtruncnorm(n_samples, a = 0, b = 1, mean = 0.8, sd = sqrt(1 / 100)),
  kappa   = rbeta(n_samples, shape1 = 40, shape2 = 2)
)


prior_df <- bind_rows(lapply(names(prior_samples), function(p) {
  data.frame(value = prior_samples[[p]], source = "Prior", parameter = p)
}))


#Replace with actual mcmc.list objects
case_mod<- Case_modlstfinalb
ww_mod  <- WW_modlstfinalb
com_mod <- Comblist_final


extract_param <- function(mcmc_obj, param_name, model_label) {
  df <- as.data.frame(as.matrix(mcmc_obj)[, param_name])
  colnames(df) <- "value"
  df$source <- model_label
  df$parameter <- param_name
  return(df)
}

# Loop over each parameter and model
get_posteriors <- function(param_name) {
  bind_rows(
    extract_param(case_mod, param_name, "Case-only"),
    extract_param(ww_mod, param_name, "Viral load-only"),
    extract_param(com_mod, param_name, "Combined")
  )
}

posterior_df <- bind_rows(lapply(vloadparams, get_posteriors))
beta_df_all <- bind_rows(prior_df, posterior_df)

# Make sure 'source' and 'parameter' are factors with correct ordering
beta_df_all$source <- factor(beta_df_all$source, 
                             levels = c("Prior","Viral load-only","Case-only","Combined")
)

beta_df_all$parameter <- factor(beta_df_all$parameter, 
                                levels = c("beta","kappa"))  # Adjust as needed

beta_param <- "beta"
kappa_param <- "kappa"
# Set grayscale colors (light to dark)
gray_palette <- c("#cccccc", "#999999", "#666666", "#333333")  

plot1=ggplot(beta_df_all%>% filter(parameter %in% beta_param), aes(x = value, y = source, fill = source)) +
  geom_density_ridges(
    alpha = 0.7,
    quantile_lines = TRUE,
    quantiles = 0.5,
    scale = 0.9,
    rel_min_height = 0.0005
  ) +
  scale_fill_manual(values = gray_palette) +
  scale_y_discrete(limits = c("Combined","Viral load-only","Case-only","Prior")) +
  #scale_y_discrete(limits = c("Combined","Prior")) +
  labs(
    title = "Prior vs Posterior Distributions of Transmission Probability",
    #title = "Prior vs Posterior Distributions of mean transit time",
    x = "Transmission probability, β",
    y = NULL
  ) +
  #coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

plot1


plot2=ggplot(beta_df_all%>% filter(parameter %in% kappa_param), aes(x = value, y = source, fill = source)) +
  geom_density_ridges(
    alpha = 0.7,
    quantile_lines = TRUE,
    quantiles = 0.5,
    scale = 0.9,
    rel_min_height = 0.0005
  ) +
  scale_fill_manual(values = gray_palette) +
  scale_y_discrete(limits = c("Combined","Viral load-only","Case-only","Prior")) +
  #scale_y_discrete(limits = c("Combined","Prior")) +
  labs(
    title = "Prior vs Posterior Distributions of mixing probability",
    #title = "Prior vs Posterior Distributions of mean transit time",
    x = "Mixing probability, K",
    y = NULL
  ) +
  #coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

plot2


plot1+plot2

















































############Marginal posterior plots
###We start with those that are common in the three models
set.seed(123)

vloadparams <- c("transit_time_mean", "transit_time_cv", "mult","tau_ww")  # Add/remove as needed

# Sample size
n_samples <- 50000

log_mult ~ dnorm(log(1.5e-9), 400) ##final case
mult <- exp(log_mult)
transit_time_mean ~ dnorm(2.5, 25) T(1.5, 4)
transit_time_cv ~ dnorm(0.3, 36) T(0.2, 0.8)
tau_ww ~ dgamma(50, 60)  # mean ≈ 0.83, SD ≈ 0.12#newcombmod8


# 1. Simulate log_mult first
log_mult_samples <- rnorm(n_samples, mean = log(1.5e-9), sd = sqrt(1 / 400))

# 2. Now create the list with all priors
prior_samples <- list(
  transit_time_mean = rtruncnorm(n_samples, a = 1.5, b = 4, mean = 2.5, sd = sqrt(1 / 25)),
  transit_time_cv   = rtruncnorm(n_samples, a = 0.2, b = 0.8, mean = 0.3, sd = sqrt(1 / 36)),
  mult              = exp(log_mult_samples),
  tau_ww            = rgamma(n_samples, shape = 50, rate = 60)
)


prior_df <- bind_rows(lapply(names(prior_samples), function(p) {
  data.frame(value = prior_samples[[p]], source = "Prior", parameter = p)
}))

# Replace with actual mcmc.list objects
#case_mod     <- Case_modlstfinal
#ww_mod       <- WW_modlstfinalb
combined_mod <- Comblist_final

extract_param <- function(mcmc_obj, param_name, model_label) {
  df <- as.data.frame(as.matrix(mcmc_obj)[, param_name])
  colnames(df) <- "value"
  df$source <- model_label
  df$parameter <- param_name
  return(df)
}

# Loop over each parameter and model
get_posteriors <- function(param_name) {
  bind_rows(
    #extract_param(ww_mod, param_name, "WW-only"),
    extract_param(combined_mod, param_name, "Combined")
  )
}

posterior_df <- bind_rows(lapply(vloadparams, get_posteriors))
beta_df_all <- bind_rows(prior_df, posterior_df)

# Make sure 'source' and 'parameter' are factors with correct ordering
beta_df_all$source <- factor(beta_df_all$source, 
                             levels = c("Prior","Combined"))

beta_df_all$parameter <- factor(beta_df_all$parameter, 
                                levels = c("transit_time_mean", "transit_time_cv", "mult","tau_ww"))  # Adjust as needed


WW_params <- c("transit_time_mean", "transit_time_cv")

plotcom=ggplot(beta_df_all%>% filter(parameter %in% WW_params), aes(x = value, y = source, fill = source)) +
  geom_density_ridges(
    scale = 1.5,  # Slightly less than before
    jittered_points = TRUE,
    position = "identity",  # Prevent stacking
    quantile_lines = TRUE,
    quantiles = 0.5
  ) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("gray30", "gray60", "gray70", "gray90")) +
  #coord_cartesian(xlim = c(1.5, 3.5), clip = "off") +
  labs(
    title = "Prior vs Posterior Distributions by Parameter",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines")
  )

pretty_labels <- c(
  transit_time_mean = "Mean transit time (days)",
  transit_time_cv = "CV of transit time"
)

plotcom+facet_wrap(~ parameter, scales = "free_x", ncol = 2,
             labeller = labeller(parameter = pretty_labels))


























