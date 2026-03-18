# ============================================================
# BayGMST_R: Bayesian global temperature reconstruction pipeline 
#
# This script:
#   1. Reads user configuration from config.yml
#   2. Loads reduced proxies, forcings, and instrumental temperatures
#   3. Aligns all inputs onto a common annual time grid
#   4. Applies forcing transformations / normalization
#   5. Prepares the data list required by the Stan model
#   6. Fits the Bayesian hierarchical model with CmdStan
#   7. Saves posterior summaries
#   8. Produces a reconstruction figure and posterior histograms
# ============================================================

library(config)
library(yaml)
library(cmdstanr)
library(car)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

cfg <- yaml::read_yaml("config.yml")
#str(cfg)
#cfg$rp_method

### SET cmdstan PATH
# CmdStan is the command line interface to Stan
set_cmdstan_path(path = cfg$folder_paths$cmdstan_path)

### FUNCTIONS:
# ------------------------------------------------------------
# Helper function: load the reduced-proxy dataset specified
# in the config file. The proxy method determines which CSV
# is loaded from the reduced-proxy directory.
# ------------------------------------------------------------
load_proxies <- function(rp_method = c("LASSO", "PCR", "SIR", "SPLS"),
                         base_dir = cfg$folder_paths$barboza_rps_path) {
  rp_method <- match.arg(toupper(rp_method), c("LASSO", "PCR", "SIR", "SPLS"))
  f <- file.path(base_dir, sprintf("RP_new_All_%s.csv", rp_method))
  stopifnot(file.exists(f))
  read.csv(f)
}

### INPUTS
# ------------------------------------------------------------
# Load all model inputs:
#   - reduced proxies
#   - external forcings
#   - instrumental temperature observations
#
# Also pull the reconstruction / calibration window from the
# config file.
# ------------------------------------------------------------
rp_method        <- cfg$rp_method
Proxies.in       <- load_proxies(rp_method)
Forcings.in      <- read.csv(cfg$folder_paths$forcings_path)
Forcings.in$year <- as.integer(rownames(Forcings.in))
Temperatures.in  <- read.csv(cfg$folder_paths$instr_temp_path)
colnames(Temperatures.in) <- c("year","T","l95","u95")

t1 <- cfg$partition_years$t1
t2 <- cfg$partition_years$t2
t3 <- cfg$partition_years$t3

# validate ordering: t1 <= t2 <= t3
# ------------------------------------------------------------
# Validate that the reconstruction window is well defined.
# Required ordering is:
#   t1 = start of full reconstruction window
#   t2 = start of instrumental period
#   t3 = end of analysis window
# with t1 <= t2 <= t3.
# ------------------------------------------------------------
if (anyNA(c(t1, t2, t3))) {
  stop("t1/t2/t3 contains NA.")
}
if (!is.numeric(t1) || !is.numeric(t2) || !is.numeric(t3)) {
  stop("t1/t2/t3 must be numeric.")
}
if (!(t1 <= t2 && t2 <= t3)) {
  stop(sprintf("Invalid partition years: require t1 <= t2 <= t3, got t1=%s, t2=%s, t3=%s", t1, t2, t3))
}

vol_coef <- cfg$vol_params$vol_coef
co2_coef <- cfg$co2_params$co2_coef
co2_c0   <- cfg$co2_params$c0

### BUILD DATAFRAME
# ------------------------------------------------------------
# Build the main analysis dataframe by aligning all inputs
# onto a common annual sequence from t1 to t3.
#
# The resulting dataframe contains:
#   year = annual time index
#   S    = solar forcing
#   V    = volcanic forcing
#   G    = greenhouse gas forcing
#   T    = instrumental temperature anomaly
#   R    = reduced proxy series
# ------------------------------------------------------------
years <- t1:t3

stopifnot(
  !anyDuplicated(Forcings.in$year),
  !anyDuplicated(Temperatures.in$year),
  !anyDuplicated(Proxies.in$Year)
)

Temps_inst <- subset(Temperatures.in, year >= t2)
iF <- match(years, Forcings.in$year)
iT <- match(years, Temps_inst$year)
iP <- match(years, Proxies.in$Year)

df <- data.frame(
  year = years,
  S = Forcings.in$solar[iF],      # solar forcing
  V = Forcings.in$volcanic[iF],   # volcanism forcing
  G = Forcings.in$CO2[iF],        # greenhouse gas (CO2) forcing
  T = Temps_inst$T[iT],           # will be NA for years not in Temperatures.in (e.g., < t2)
  R = Proxies.in$RP1[iP]          # reduced proxy
)

tail(df)

## quick missingness check by column
colSums(is.na(df))

### NORMALIZE
# ------------------------------------------------------------
# Apply the forcing transformations used by the model:
#   - volcanic forcing is transformed to a negative saturating form
#   - CO2 forcing is log-transformed relative to a baseline
#   - solar forcing is centered
# ------------------------------------------------------------
df$V <- -abs(vol_coef)*(1-exp(-df$V))
df$G <- co2_coef*log(df$G/co2_c0)
df$S <- df$S - mean(df$S)


### PREPARE DATA FOR STAN (VARIBALES BELOW ARE CONSISTENT WITH STAN CODE)
# ------------------------------------------------------------
# Prepare the observed and missing temperature indices for Stan.
#
# In this setup:
#   y = instrumental temperature series (partially observed)
#   z = reduced proxy series
#
# Stan receives both the observed temperature values and the
# index locations of observed vs. missing entries, so that
# missing historical temperatures can be estimated.
# ------------------------------------------------------------
y <- df$T
z <- df$R

NT      <- length(z)
idx_obs <- which(!is.na(y))
idx_mis <- which(is.na(y))
y_obs   <- as.vector(y[idx_obs])

data_list <- list(
  NT = NT, 
  NT_obs = length(idx_obs), 
  NT_mis = length(idx_mis),
  idx_obs = as.integer(idx_obs),
  idx_mis = as.integer(idx_mis),
  G = as.vector(df$G),
  S = as.vector(df$S),
  V = as.vector(df$V),
  y_obs = y_obs,
  z = z
)

### FIT BHM with STAN
# ------------------------------------------------------------
# Fit the Bayesian hierarchical model in Stan.
# The model file is read from the config, compiled, and then
# sampled using HMC through cmdstanr.
# ------------------------------------------------------------
message("Running STAN model now...")
mod <- cmdstan_model(cfg$folder_paths$stan_code_path)
t <- system.time({
  fit <- mod$sample(data = data_list, chains = 4, parallel_chains = 2,
                    iter_warmup = 500, iter_sampling = 2500)
})
elapsed_sec <- unname(t["elapsed"])
elapsed_sec
message("Done.")

# posterior summaries for parameters
# ------------------------------------------------------------
# Save posterior summaries for key model parameters so they
# can be inspected outside R or reused in later analysis.
# ------------------------------------------------------------
summ <- fit$summary(variables = c(
  "alpha0","alpha1","phi_R","phi_T",
  "beta0","betaG","betaS","betaV",
  "sigma_y","sigma_z"
))
out <- summ[]
out_path <- file.path(cfg$folder_paths$reconstruction_dir, "fit_post_summaries.csv") # NEED TO FIX THIS!!
write.csv(out, file = out_path, row.names = FALSE) # NEED TO FIX THIS!!


# PLOTTING
# ------------------------------------------------------------
# Extract posterior draws of the missing temperature states,
# which correspond to the reconstructed temperature series
# outside the observed instrumental period.
# ------------------------------------------------------------
draws_mean <- fit$draws("y_mis")
idx_names  <- paste0("y_mis[", seq_along(idx_mis), "]")
mat        <- posterior::as_draws_matrix(draws_mean)[, idx_names, drop = FALSE]
y1_post    <- cbind(
  t = idx_mis,
  mean = apply(mat, 2, mean),
  lo = apply(mat, 2, quantile, 0.025),
  hi = apply(mat, 2, quantile, 0.975)
)

df_pred <- data.frame(
  year = as.numeric(idx_mis + t1),
  mean = as.numeric(y1_post[, "mean"]),
  lo   = as.numeric(y1_post[, "lo"]),
  hi   = as.numeric(y1_post[, "hi"])
)

df_obs <- data.frame(
  year = as.numeric(idx_obs + t1),
  T    = as.numeric(y_obs)
)

# time series of reconstructions
# ------------------------------------------------------------
# Create the main time-series reconstruction plot showing:
#   - posterior mean reconstruction
#   - 95% credible interval ribbon
#   - instrumental temperature observations
# ------------------------------------------------------------
x_lim <- range(c(df_pred$year, df_obs$year), na.rm = TRUE)
y_lim <- range(c(df_obs$T, df_pred$lo, df_pred$hi), na.rm = TRUE)

p_ts <- ggplot() +
  geom_ribbon(
    data = df_pred, fill="#005AB5",
    aes(x = year, ymin = lo, ymax = hi),
    alpha = 0.4
  ) +
  geom_line(
    data = df_pred, color="#005AB5",
    aes(x = year, y = lo),
    alpha = 0.2
  ) +
  geom_line(
    data = df_pred, color="#005AB5",
    aes(x = year, y = hi),
    alpha = 0.2
  ) +
  geom_line(
    data = df_pred,
    aes(x = year, y = mean, color = "Posterior mean (w/ 95% CrI)"),
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_line(
    data = df_obs,
    aes(x = year, y = T, color = "T Anomaly, HadCRUT5"),
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "",
    values = c("T Anomaly, HadCRUT5" = "black", "Posterior mean (w/ 95% CrI)" = "red")
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(x = "year", y = "T (deg C)") +
  theme_light(base_size = 12) +
  theme(legend.position = c(0.22,0.85),
        legend.background = element_rect(fill = NA, color = NA))

# histograms of posterior dist. of parameters
# ------------------------------------------------------------
# Prepare posterior draws for selected structural parameters
# and visualize their posterior distributions.
#
# betaG, betaV, betaS = forcing effects
# phi_R, phi_T        = AR(1) persistence parameters
# ------------------------------------------------------------
betas <- c("betaG","betaV","betaS","phi_R","phi_T")

draws_df <- fit$draws(variables = betas, format = "df")

df_hist <- draws_df %>%
  select(all_of(betas)) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value")

beta_params <- c("betaG","betaV","betaS")
phi_params  <- c("phi_R","phi_T")

# max over the first three parameters
beta_max <- df_hist %>%
  filter(parameter %in% beta_params) %>%
  summarise(mx = max(value, na.rm = TRUE)) %>%
  pull(mx)

beta_min <- df_hist %>%
  filter(parameter %in% beta_params) %>%
  summarise(mn = min(value, na.rm = TRUE)) %>%
  pull(mn)

base_hist <- function(dat) {
  ggplot(dat, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, linewidth = 0.2,
                   color = "white", fill = "#1A85FF") +
    geom_vline(xintercept = 0.0) +
    geom_hline(yintercept = 0.0) +
    labs(x = "Posterior dist.", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}

p_beta <- df_hist %>%
  filter(parameter %in% beta_params) %>%
  base_hist() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(beta_min, beta_max))

p_phi <- df_hist %>%
  filter(parameter %in% phi_params) %>%
  base_hist() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-0.1, 1))

p_hist <- p_beta / p_phi

# plots combine side by side
# ------------------------------------------------------------
# Combine the reconstruction panel and posterior histogram
# panels into one final figure with a descriptive subtitle.
# ------------------------------------------------------------
sub_txt <- sprintf(
  "Reconstruction window: (%s, %s);  RP computed via %s;  Estimation via Stan (HMC)",
  t1, t2, rp_method
)
p <- p_ts + p_hist + plot_layout(widths = c(4, 1)) + plot_annotation(
  title = "BHM with AR(1) structure for both the R and T equations",
  subtitle = sub_txt
) &
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p

ggsave(paste0(cfg$folder_paths$figures_dir,"/reconstruction_ts.png"), plot = p, width = 10, height = 5, units = "in", dpi = 300, bg = "white")
