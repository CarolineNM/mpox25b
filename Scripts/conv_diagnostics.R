###########This has the convergence diagnostics,traceplots and density plots and marginal distribution plots
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


####################refine this######################
load("D:/mpox25output/Comblist_finalg.RData")
load("D:/mpox25output/Combined_WWd.RData")
load("D:/mpox25output/Combined_casc.RData")

#mcmc_matrixallcas<-as.matrix(Combined_cas)
#mcmc_matrixallWW<-as.matrix(Combined_WW)
#mcmc_matrixallcom<-as.matrix(Comblist_finale)

case_mod  <- Combined_casc
ww_mod   <- Combined_WWd
com_mod <- Comblist_finalg

#case_mod  <- Case_modlstfinalb
#ww_mod   <- WW_modlstfinale
#com_mod <- Comblist_finale
#com_modb <- Comblist_final

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



















