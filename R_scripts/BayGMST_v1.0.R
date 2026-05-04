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
t4 <- cfg$partition_years$t4

# validate ordering: t1 <= t2 <= t3
# ------------------------------------------------------------
# Validate that the reconstruction window is well defined.
# Required ordering is:
#   t1 = start of full reconstruction window
#   t2 = start of instrumental period
#   t3 = end of analysis window
#   t4 = end of future projections window
# with t1 <= t2 <= t3 <= 4.
# ------------------------------------------------------------
check_year <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("%s is missing.", name), call. = FALSE)
  }
  if (
    !is.numeric(x) ||
    length(x) != 1 ||
    is.na(x) ||
    !is.finite(x) ||
    x <= 0 ||
    x != as.integer(x)
  ) {
    stop(sprintf("%s must be a single positive integer.", name), call. = FALSE)
  }
}
check_year(t1, "t1")
check_year(t2, "t2")
check_year(t3, "t3")

if (!(t1 <= t2 && t2 <= t3)) {
  stop(sprintf("Invalid partition years: require t1 <= t2 <= t3, got t1=%s, t2=%s, t3=%s", t1, t2, t3))
}
if (!is.null(t4)){
  check_year(t4, "t4")
  if (anyNA(c(t4))) {
    stop("t4 contains NA.")
  }
  if (!is.numeric(t4)) {
    stop("t4 must be numeric.")
  }
  if (!(t3 < t4 && t4 <= 2100)) {
    stop(sprintf("Invalid partition years: require t3 < t4 <= 2100, got t3=%s, t4=%s", t3, t4))
  }
}

instru_year_min <- min(Temperatures.in$year, na.rm = TRUE)
instru_year_max <- max(Temperatures.in$year, na.rm = TRUE)
if (!all(c(t2, t3) >= instru_year_min & c(t2, t3) <= instru_year_max)) {
  stop(
    sprintf(
      "t2 and/or t3 are outside the observed instrumental temperature year range [%s, %s]. t2 = %s, t3 = %s",
      instru_year_min, instru_year_max, t2, t3
    )
  )
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

#tail(df)
#plot(df$year, df$V, type='l')

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
iter_warmup   <- cfg$stan_params$iter_warmup
iter_sampling <- cfg$stan_params$iter_sampling

if (iter_sampling < 1000) {
  stop("cfg$stan_params$iter_sampling must be at least 1000.")
}

message("Running STAN model now...")
mod <- cmdstan_model(cfg$folder_paths$stan_code_path)
t <- system.time({
  fit <- mod$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling
  )
})
elapsed_sec <- unname(t["elapsed"])
elapsed_sec
message("Done.")


tail(fit$summary(variables = c("y_ins_fitted")))



obs_idx <- !is.na(df$T)
y_ins_summary <- fit$summary(variables = "y_ins_fitted")
y_ins_summary <- fit$summary(
  variables = "y_ins_fitted",
  ~posterior::quantile2(.x, probs = c(0.025, 0.160, 0.840, 0.975)),
  "mean"
)
df_ins <- data.frame(
  year = df$year[obs_idx],
  T.obs = df$T[obs_idx],
  T.mean = y_ins_summary$mean,
  T.lo = y_ins_summary$q16,
  T.hi = y_ins_summary$q84,
  T.lolo = y_ins_summary$q2.5,
  T.hihi = y_ins_summary$q97.5
)
err <- df_ins$T.mean - df_ins$T.obs
perf_stats <- data.frame(
  RMSE = sqrt(mean(err^2, na.rm = TRUE)),
  MAE = mean(abs(err), na.rm = TRUE),
  Bias = mean(err, na.rm = TRUE),
  Correlation = cor(df_ins$T.obs, df_ins$T.mean, use = "complete.obs"),
  R2 = cor(df_ins$T.obs, df_ins$T.mean, use = "complete.obs")^2
)
perf_stats


plot(df_ins$year, df_ins$T.mean, type='l', col='red', lwd=2.0, lty=1, ylim=c(-0.7, +0.5))
lines(df_ins$year, df_ins$T.lolo,  type='l', col='red', lwd=1.5, lty=5)
lines(df_ins$year, df_ins$T.hihi, type='l', col='red', lwd=1.5, lty=5)
lines(df_ins$year, df_ins$T.obs, Temps_inst$T, type='l', col='black', lwd=1.0, lty=1)

# residual analysis
df_ins$resid <- df_ins$T.obs - df_ins$T.mean
acf(df_ins$resid, na.action = na.pass, main = "ACF of fitted residuals")
pacf(df_ins$resid, na.action = na.pass, main = "PACF of fitted residuals")



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
  lo = apply(mat, 2, quantile, 0.160),
  hi = apply(mat, 2, quantile, 0.840),
  lolo = apply(mat, 2, quantile, 0.025),
  hihi = apply(mat, 2, quantile, 0.975)
)

df_pred <- data.frame(
  year = as.numeric(idx_mis + t1),
  mean = as.numeric(y1_post[, "mean"]),
  lo   = as.numeric(y1_post[, "lo"]),
  hi   = as.numeric(y1_post[, "hi"]),
  lolo   = as.numeric(y1_post[, "lolo"]),
  hihi   = as.numeric(y1_post[, "hihi"])
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
x_lim <- range(c(df_pred$year, df_obs$year, df_ins$year), na.rm = TRUE)
y_lim <- range(
  c(df_obs$T, df_pred$lolo, df_pred$hihi, df_ins$T.lolo, df_ins$T.hihi),
  na.rm = TRUE
)
r2_label <- sprintf("Instrumental Period R² = %.2f", perf_stats$R2)
annot_x <- x_lim[2] - 0.01 * diff(x_lim)
annot_y <- y_lim[1] + 0.05 * diff(y_lim)

alpha_dark = 0.90
alpha_lite = 0.25

p_ts <- ggplot() +
  geom_ribbon(
    data = df_pred,
    aes(
      x = year,
      ymin = lolo,
      ymax = hihi,
      fill = "95% Credible Band",
      alpha = "95% Credible Band"
    )
  ) +
  geom_ribbon(
    data = df_pred,
    aes(
      x = year,
      ymin = lo,
      ymax = hi,
      fill = "68% Credible Band",
      alpha = "68% Credible Band"
    )
  ) +
  geom_line(
    data = df_pred, color = "cyan3",
    aes(x = year, y = lolo),
    alpha = 0.2
  ) +
  geom_line(
    data = df_pred, color = "cyan3",
    aes(x = year, y = hihi),
    alpha = 0.2
  ) +
  geom_line(
    data = df_pred,
    aes(x = year, y = mean, color = "Pre-Instrumental Reconstruction (Post. Mean)"),
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_ribbon(
    data = df_ins,
    aes(x = year, ymin = T.lo, ymax = T.hi),
    fill = "darkorange3",
    alpha = alpha_dark,
    na.rm = TRUE
  ) +
  geom_ribbon(
    data = df_ins,
    aes(x = year, ymin = T.lolo, ymax = T.hihi),
    fill = "darkorange2",
    alpha = alpha_lite,
    na.rm = TRUE
  ) +
  geom_line(
    data = df_obs,
    aes(x = year, y = T, color = "HadCRUT5 (Instrumental Observations)"),
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_line(
    data = df_ins,
    aes(x = year, y = T.mean, color = "Instrumental Reconstruction (Post. Mean)"),
    linewidth = 0.55,
    alpha = 0.90,
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "GMST Anomaly",
    values = c(
      "HadCRUT5 (Instrumental Observations)" = "grey20",
      "Pre-Instrumental Reconstruction (Post. Mean)" = "darkorchid4",
      "Instrumental Reconstruction (Post. Mean)" = "limegreen"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "95% Credible Band" = "cyan3",
      "68% Credible Band" = "cyan3"
    )
  ) +
  scale_alpha_manual(
    name = "",
    values = c(
      "95% Credible Band" = alpha_lite,
      "68% Credible Band" = alpha_dark
    )
  )  +
  guides(
    alpha = "none",
    color = guide_legend(
      byrow = TRUE,
      keyheight = unit(0.55, "lines")
    ), 
    fill = guide_legend(
      override.aes = list(
        alpha = c(alpha_dark, alpha_lite)
      ), keyheight = unit(0.55, "lines")
    )
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  annotate(
    "text",
    x = annot_x,
    y = annot_y,
    label = r2_label,
    hjust = 1,
    vjust = 0,
    size = 3.5
  ) +
  labs(x = "year", y = "GMST Anomaly (°C)") +
  theme_light(base_size = 11) +
  theme(
    legend.position = c(0.01, 0.60),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.01, "lines")
  )

p_ts
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

# --- trace plots
trace_df <- draws_df %>%
  pivot_longer(
    cols = all_of(betas),
    names_to = "parameter",
    values_to = "value"
  )

p_trace <- ggplot(trace_df, aes(x = .iteration, y = value, group = .chain, color = factor(.chain))) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  labs(
    x = "Iteration",
    y = "Draw value",
    color = "Chain",
    title = "Trace plots"
  ) +
  theme_bw()

ggsave(paste0(cfg$folder_paths$figures_dir,"/trace_plots.png"), 
       plot = p_trace, width = 10, height = 5, units = "in", dpi = 300, bg = "white")
# ------ # 

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
                   bins = 80, linewidth = 0.15,
                   color = "white", fill = "gray40") +
    geom_vline(xintercept = 0.0) +
    geom_hline(yintercept = 0.0) +
    labs(x = "Posterior dist.", y = "Density") +
    scale_x_continuous(
      labels = function(x) signif(x, 1)
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 3)
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          panel.grid.minor = element_blank())
}

p_beta <- df_hist %>%
  filter(parameter %in% beta_params) %>%
  base_hist() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(beta_min, beta_max))
p_beta
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
  "Reconstruction window: (%s, %s);  RP computed via %s;  AR(1) structure in T and R equations",
  t1, t2, rp_method
)
p <- p_ts + p_hist + plot_layout(widths = c(4, 1)) + plot_annotation(
  title = "GMST Reconstruction using a Reduced Proxy",
  subtitle = sub_txt
) &
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p

ggsave(paste0(cfg$folder_paths$figures_dir,"/reconstruction_ts.png"), plot = p, width = 10, height = 5, units = "in", dpi = 300, bg = "white")
