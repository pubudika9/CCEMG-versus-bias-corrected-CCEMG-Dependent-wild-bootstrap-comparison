### Monte Carlo comparison: standard CCEMG versus ABC whole-sum DWB inference.
### The bootstrap bandwidth is fixed at l_N = 1, so independent cross-sectional
### multipliers are applied directly to the complete corrected unit-level slope
### contributions. The reported ABC standard error uses the exact conditional
### covariance, while multiplier draws are used for the two tail quantiles.

rm(list = ls())

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required. Install it with install.packages('MASS').")
}

read_int_env <- function(name, default_value, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default_value))))
  if (is.na(value) || value < minimum) as.integer(default_value) else value
}

read_num_env <- function(name, default_value) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default_value))))
  if (is.na(value) || !is.finite(value)) default_value else value
}

read_num_vector_env <- function(name, default_values, minimum = -Inf) {
  raw_value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw_value)) return(default_values)
  values <- suppressWarnings(as.numeric(trimws(strsplit(raw_value, ",", fixed = TRUE)[[1]])))
  if (anyNA(values) || any(!is.finite(values)) || any(values < minimum)) {
    stop(sprintf("%s must contain comma-separated numbers greater than or equal to %s.", name, minimum))
  }
  unique(values)
}

safe_inv <- function(x, eps = 1e-10) {
  if (!is.finite(x)) return(0)
  if (abs(x) < eps) return(sign(x + eps) / eps)
  1 / x
}

make_Mw <- function(y_bar_t, x_bar_t) {
  H <- cbind("(Intercept)" = 1, y_bar_t, x_bar_t)
  Tn <- nrow(H)
  diag(Tn) - H %*% MASS::ginv(crossprod(H)) %*% t(H)
}

# Pesaran's original mean-group covariance, preserved from the supplied code.
var_ccemg_pesaran <- function(betas_full, beta_hat) {
  N <- length(betas_full)
  if (N <= 1L) return(0)
  sum((betas_full - beta_hat)^2) / (N * (N - 1L))
}

build_factor_loadings <- function(N, rank_satisfied) {
  if (rank_satisfied) {
    list(
      gamma_i1 = rnorm(N, mean = 1.0, sd = 0.2),
      gamma_i2 = rnorm(N, mean = 0.0, sd = 0.2),
      Gamma_x_i1 = rnorm(N, mean = 0.5, sd = 0.5),
      Gamma_x_i2 = rnorm(N, mean = 0.5, sd = 0.5)
    )
  } else {
    list(
      gamma_i1 = rnorm(N, mean = 1.0, sd = 0.2),
      gamma_i2 = rnorm(N, mean = 1.0, sd = 0.2),
      Gamma_x_i1 = rnorm(N, mean = 0.5, sd = 0.5),
      Gamma_x_i2 = rnorm(N, mean = 0.5, sd = 0.5)
    )
  }
}

# Frisch-Waugh-Lovell computes exactly the slopes from
# lm(y ~ x + y_bar + x_bar), while fitting the full panel much faster.
fit_cce_panel <- function(Y, X, rows = seq_len(nrow(Y))) {
  Y_sub <- Y[rows, , drop = FALSE]
  X_sub <- X[rows, , drop = FALSE]
  Tn <- nrow(Y_sub)
  N <- ncol(Y_sub)
  if (Tn < 4L) stop("At least four periods are required for individual CCE regressions.")

  y_bar <- rowMeans(Y_sub)
  x_bar <- rowMeans(X_sub)
  H <- cbind("(Intercept)" = 1, y_bar = y_bar, x_bar = x_bar)
  qr_H <- qr(H, tol = 1e-10)

  X_tilde <- as.matrix(qr.resid(qr_H, X_sub))
  Y_tilde <- as.matrix(qr.resid(qr_H, Y_sub))
  denominators <- colSums(X_tilde^2)
  numerators <- colSums(X_tilde * Y_tilde)
  betas <- ifelse(is.finite(denominators) & denominators > 1e-10,
                  numerators / denominators, 0)

  gamma_controls <- as.matrix(qr.coef(qr_H, Y_sub - sweep(X_sub, 2L, betas, "*")))
  Gamma_controls <- as.matrix(qr.coef(qr_H, X_sub))
  gamma_controls[!is.finite(gamma_controls)] <- 0
  Gamma_controls[!is.finite(Gamma_controls)] <- 0

  residuals_y <- Y_tilde - sweep(X_tilde, 2L, betas, "*")
  q_values <- denominators / Tn

  list(
    betas = as.numeric(betas),
    gammas = t(gamma_controls[2:3, , drop = FALSE]),
    Gammas = t(Gamma_controls[2:3, , drop = FALSE]),
    residuals_y = residuals_y,
    residuals_x = X_tilde,
    q_values = as.numeric(q_values),
    y_bar = y_bar,
    x_bar = x_bar,
    T = Tn,
    N = N
  )
}

# This uses the subtraction sign in Equation (44): the estimated leading bias
# is subtracted from the mean-group estimator. All nuisance variances are
# estimated from CCE residuals, so the correction is feasible rather than
# using the latent simulated v_it and epsilon_it.
analytical_bias_correction <- function(fit) {
  N <- fit$N
  Tn <- fit$T
  betas_full <- fit$betas
  gamma_matrix <- fit$gammas
  Gamma_matrix <- fit$Gammas

  Cbar <- rbind(
    colMeans(gamma_matrix + sweep(Gamma_matrix, 1L, betas_full, "*")),
    colMeans(Gamma_matrix)
  )
  CplusT <- t(MASS::ginv(Cbar))
  denom <- fit$q_values
  inverse_denom <- vapply(denom, safe_inv, numeric(1))
  df_x <- max(1L, Tn - 3L)
  df_eps <- max(1L, Tn - 4L)
  sig_v_hat <- colSums(fit$residuals_x^2) / df_x
  sig_eps_hat <- colSums(fit$residuals_y^2) / df_eps

  b0_i <- numeric(N)
  b1_i <- numeric(N)
  for (i in seq_len(N)) {
    gamma_i_vec <- matrix(gamma_matrix[i, ], ncol = 1L)
    Gamma_i_vec <- matrix(Gamma_matrix[i, ], ncol = 1L)
    svec_i <- matrix(c(sig_eps_hat[i], 1), ncol = 1L)

    b0_i[i] <- (inverse_denom[i] * sig_v_hat[i]) *
      as.numeric(t(c(betas_full[i], 1)) %*% CplusT %*% gamma_i_vec)
    b1_i[i] <- inverse_denom[i] *
      as.numeric(t(Gamma_i_vec) %*% CplusT %*% svec_i)
  }

  b01_hat <- sqrt(Tn / N) * (mean(b0_i) + mean(b1_i))
  beta_mg <- mean(betas_full)

  beta_abc <- beta_mg - b01_hat / sqrt(N * Tn)
  unit_bias_adjustment <- -(b0_i + b1_i) / N
  unit_corrected_slopes <- betas_full + unit_bias_adjustment

  # The ABC estimator is exactly the average of these corrected unit slopes.
  if (!isTRUE(all.equal(mean(unit_corrected_slopes), beta_abc,
                        tolerance = 1e-9))) {
    stop("Unit-level ABC contributions do not reproduce the ABC estimate.")
  }

  list(beta_abc = beta_abc, bias_term = b01_hat,
       unit_corrected_slopes = unit_corrected_slopes,
       unit_bias_adjustment = unit_bias_adjustment,
       b0_i = b0_i, b1_i = b1_i,
       sig_v_hat = sig_v_hat, sig_eps_hat = sig_eps_hat)
}

# Bartlett multipliers can be generated as moving averages of independent N(0,1)
# innovations across units: xi_i = ell^(-1/2) sum_{j=0}^{ell-1} u_{i+j}.
# The input is the complete centered cluster influence for each individual.
generate_bartlett_dwb_draws <- function(unit_influence, bootstrap_reps, bandwidth) {
  unit_influence <- as.numeric(unit_influence)
  N <- length(unit_influence)
  ell <- as.integer(bandwidth)
  extended_length <- N + ell - 1L

  innovation_coefficients <- numeric(extended_length)
  for (offset in seq_len(ell)) {
    indices <- offset:(offset + N - 1L)
    innovation_coefficients[indices] <- innovation_coefficients[indices] + unit_influence
  }

  innovations <- matrix(rnorm(extended_length * bootstrap_reps),
                        nrow = extended_length, ncol = bootstrap_reps)
  as.numeric(crossprod(innovation_coefficients, innovations)) / sqrt(N * ell)
}

whole_sum_dwb_inference <- function(fit, abc, bootstrap_reps,
                                    confidence_level, interval_type,
                                    finite_sample_correction = TRUE) {
  N <- fit$N
  Tn <- fit$T
  beta_abc <- abc$beta_abc

  # sqrt(T)*(beta_ABC,i-beta_ABC) estimates the entire unit-level summand:
  # sqrt(T/N)*v_i + T^(-1/2)*sigma_v,i^(-1)*V_i'*epsilon_i.
  centered_slopes <- abc$unit_corrected_slopes - beta_abc
  unit_influence <- sqrt(Tn) * centered_slopes
  if (finite_sample_correction && N > 1L) {
    unit_influence <- sqrt(N / (N - 1L)) * unit_influence
  }

  bandwidth <- 1L
  statistic_draws <- generate_bartlett_dwb_draws(
    unit_influence, bootstrap_reps, bandwidth
  )

  # Exact conditional covariance implied by the Bartlett multipliers.
  unit_distance <- abs(outer(seq_len(N), seq_len(N), "-"))
  bartlett_weights <- pmax(1 - unit_distance / bandwidth, 0)
  conditional_sigma <- as.numeric(
    t(unit_influence) %*% bartlett_weights %*% unit_influence / N
  )

  sigma_boot_mc <- stats::var(statistic_draws)
  if (!is.finite(sigma_boot_mc)) sigma_boot_mc <- 0
  standard_error <- sqrt(max(conditional_sigma, 0) / (N * Tn))

  variance_denominator <- if (finite_sample_correction) {
    N * (N - 1L)
  } else {
    N^2
  }
  corrected_slope_variance <- sum(centered_slopes^2) / variance_denominator
  if (!isTRUE(all.equal(conditional_sigma / (N * Tn),
                        corrected_slope_variance, tolerance = 1e-12))) {
    stop("The l_N=1 DWB covariance identity failed.")
  }
  tail_probability <- (1 - confidence_level) / 2

  # Asymmetric basic-bootstrap interval.  Because statistic_draws approximate
  # sqrt(NT)*(beta_ABC-beta_0), the upper bootstrap quantile determines the
  # lower endpoint and the lower quantile determines the upper endpoint:
  # [beta_ABC-q*_(1-alpha/2)/sqrt(NT),
  #  beta_ABC-q*_(alpha/2)/sqrt(NT)].
  if (interval_type %in% c("asymmetric", "basic", "percentile")) {
    q_low <- as.numeric(stats::quantile(statistic_draws, probs = tail_probability,
                                        names = FALSE, type = 8L))
    q_high <- as.numeric(stats::quantile(statistic_draws, probs = 1 - tail_probability,
                                         names = FALSE, type = 8L))
    lower <- beta_abc - q_high / sqrt(N * Tn)
    upper <- beta_abc - q_low / sqrt(N * Tn)
  } else {
    critical_value <- as.numeric(stats::quantile(abs(statistic_draws),
                                                 probs = confidence_level,
                                                 names = FALSE, type = 8L))
    half_width <- critical_value / sqrt(N * Tn)
    lower <- beta_abc - half_width
    upper <- beta_abc + half_width
  }

  if (!is.finite(lower) || !is.finite(upper) || lower > upper) {
    stop("ABC bootstrap confidence interval endpoints are invalid.")
  }

  list(lower = lower, upper = upper, se = standard_error,
       sigma_boot_mc = sigma_boot_mc, conditional_sigma = conditional_sigma,
       estimator_variance = conditional_sigma / (N * Tn),
       unit_score_variance = stats::var(unit_influence),
       bandwidth = bandwidth, bootstrap_dimension = "cross_section")
}

generate_panel <- function(N, Tn, alpha, slope_mode, rank_satisfied,
                           beta_true, heterogeneity_sd, eps_sd) {
  f1 <- as.numeric(stats::arima.sim(model = list(ar = 0.5), n = Tn))
  f2 <- as.numeric(stats::arima.sim(model = list(ar = 0.5), n = Tn))

  if (slope_mode == "heterogeneous") {
    eta_i <- heterogeneity_sd * rnorm(N)
    beta_i <- beta_true + eta_i / (N^alpha)
  } else {
    beta_i <- rep(beta_true, N)
  }

  loads <- build_factor_loadings(N, rank_satisfied)
  a_i <- rnorm(N)
  A_i <- rnorm(N)
  X <- matrix(0, nrow = Tn, ncol = N)
  Y <- matrix(0, nrow = Tn, ncol = N)

  for (i in seq_len(N)) {
    v_it <- as.numeric(stats::arima.sim(model = list(ar = runif(1, 0.05, 0.95)), n = Tn))
    eps_it <- eps_sd * as.numeric(stats::arima.sim(
      model = list(ar = runif(1, 0.05, 0.95)), n = Tn
    ))

    X[, i] <- A_i[i] + loads$Gamma_x_i1[i] * f1 +
      loads$Gamma_x_i2[i] * f2 + v_it
    Y[, i] <- a_i[i] + beta_i[i] * X[, i] + loads$gamma_i1[i] * f1 +
      loads$gamma_i2[i] * f2 + eps_it
  }

  list(Y = Y, X = X, beta_i = beta_i)
}

run_single_cell <- function(N, Tn, alpha, slope_mode, rank_satisfied,
                            replications, beta_true, heterogeneity_sd, eps_sd,
                            power_null, power_null_grid, bootstrap_reps,
                            confidence_level, interval_type,
                            finite_sample_correction) {
  est_ccemg <- numeric(replications)
  est_abc <- numeric(replications)
  coverage_ccemg <- logical(replications)
  coverage_boot <- logical(replications)
  size_ccemg <- logical(replications)
  size_boot <- logical(replications)
  interval_length_ccemg <- numeric(replications)
  interval_length_boot <- numeric(replications)
  se_ccemg <- numeric(replications)
  se_boot <- numeric(replications)
  conditional_sigma_values <- numeric(replications)
  estimator_variance_values <- numeric(replications)
  bootstrap_sigma_mc_values <- numeric(replications)
  unit_score_variance_values <- numeric(replications)
  bandwidth_values <- integer(replications)
  power_ccemg <- matrix(FALSE, nrow = replications, ncol = length(power_null_grid))
  power_boot <- matrix(FALSE, nrow = replications, ncol = length(power_null_grid))
  z_critical <- stats::qnorm((1 + confidence_level) / 2)

  for (replication in seq_len(replications)) {
    panel <- generate_panel(N, Tn, alpha, slope_mode, rank_satisfied,
                            beta_true, heterogeneity_sd, eps_sd)
    fit <- fit_cce_panel(panel$Y, panel$X)
    beta_mg <- mean(fit$betas)
    abc <- analytical_bias_correction(fit)

    se_standard <- sqrt(max(var_ccemg_pesaran(fit$betas, beta_mg), 0))
    ci_standard <- c(beta_mg - z_critical * se_standard,
                     beta_mg + z_critical * se_standard)
    if (!is.finite(ci_standard[1L]) || !is.finite(ci_standard[2L]) ||
        ci_standard[1L] > ci_standard[2L]) {
      stop("CCEMG confidence interval endpoints are invalid.")
    }

    boot <- whole_sum_dwb_inference(
      fit = fit, abc = abc,
      bootstrap_reps = bootstrap_reps, confidence_level = confidence_level,
      interval_type = interval_type,
      finite_sample_correction = finite_sample_correction
    )

    est_ccemg[replication] <- beta_mg
    est_abc[replication] <- abc$beta_abc
    se_ccemg[replication] <- se_standard
    se_boot[replication] <- boot$se
    conditional_sigma_values[replication] <- boot$conditional_sigma
    estimator_variance_values[replication] <- boot$estimator_variance
    bootstrap_sigma_mc_values[replication] <- boot$sigma_boot_mc
    unit_score_variance_values[replication] <- boot$unit_score_variance
    bandwidth_values[replication] <- boot$bandwidth

    coverage_ccemg[replication] <- beta_true >= ci_standard[1L] && beta_true <= ci_standard[2L]
    coverage_boot[replication] <- beta_true >= boot$lower && beta_true <= boot$upper
    size_ccemg[replication] <- !coverage_ccemg[replication]
    size_boot[replication] <- !coverage_boot[replication]
    interval_length_ccemg[replication] <- diff(ci_standard)
    interval_length_boot[replication] <- boot$upper - boot$lower

    power_ccemg[replication, ] <- power_null_grid < ci_standard[1L] |
      power_null_grid > ci_standard[2L]
    power_boot[replication, ] <- power_null_grid < boot$lower |
      power_null_grid > boot$upper
  }

  closest_power_null <- which.min(abs(power_null_grid - power_null))
  summary <- data.frame(
    N = N, T = Tn, alpha = alpha,
    Bias_CCEMG = sqrt(N * Tn) * mean(est_ccemg - beta_true),
    Bias_BOOT = sqrt(N * Tn) * mean(est_abc - beta_true),
    RMSE_CCEMG = sqrt(mean((est_ccemg - beta_true)^2)),
    RMSE_BOOT = sqrt(mean((est_abc - beta_true)^2)),
    Coverage_CCEMG = mean(coverage_ccemg),
    Coverage_BOOT = mean(coverage_boot),
    AvgLen_CCEMG = mean(interval_length_ccemg),
    AvgLen_BOOT = mean(interval_length_boot),
    Size_CCEMG = mean(size_ccemg),
    Size_BOOT = mean(size_boot),
    Power_CCEMG = mean(power_ccemg[, closest_power_null]),
    Power_BOOT = mean(power_boot[, closest_power_null]),
    MeanSE_CCEMG = mean(se_ccemg),
    MeanSE_BOOT = mean(se_boot),
    MeanConditionalSigma_DWB = mean(conditional_sigma_values),
    MeanEstimatorVariance_BOOT = mean(estimator_variance_values),
    MeanBootstrapSigma_MC = mean(bootstrap_sigma_mc_values),
    MeanUnitInfluenceVariance = mean(unit_score_variance_values),
    MeanDWB_Bandwidth = mean(bandwidth_values),
    DWB_Dimension = "cross_section",
    Bootstrap_Reps = bootstrap_reps,
    Simulation_Reps = replications,
    Beta_True = beta_true,
    Reference_Power_Null = power_null_grid[closest_power_null]
  )

  power_curve <- data.frame(
    N = N, T = Tn, alpha = alpha, beta_true = beta_true,
    beta_null = power_null_grid,
    Power_CCEMG = colMeans(power_ccemg),
    Power_BOOT = colMeans(power_boot)
  )

  list(summary = summary, power_curve = power_curve)
}

build_case_table <- function() {
  data.frame(
    case_id = c(
      "case1_hetero_rank_satisfied",
      "case2_hetero_rank_not_satisfied",
      "case3_homo_rank_satisfied",
      "case4_homo_rank_not_satisfied"
    ),
    case_name = c(
      "Heterogeneous; rank satisfied",
      "Heterogeneous; rank not satisfied",
      "Homogeneous; rank satisfied",
      "Homogeneous; rank not satisfied"
    ),
    slope_mode = c("heterogeneous", "heterogeneous", "homogeneous", "homogeneous"),
    rank_satisfied = c(TRUE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

format_alpha_tag <- function(alpha, slope_mode) {
  if (slope_mode == "homogeneous" || is.na(alpha)) return("alpha_homo")
  paste0("alpha_", gsub("\\.", "p", format(alpha, trim = TRUE)))
}

write_metric_matrix <- function(data, column_name, destination) {
  n_values <- sort(unique(data$N))
  t_values <- sort(unique(data$T))
  values <- matrix(NA_real_, nrow = length(t_values), ncol = length(n_values),
                   dimnames = list(NULL, as.character(n_values)))
  for (row in seq_len(nrow(data))) {
    values[match(data$T[row], t_values), match(data$N[row], n_values)] <- data[[column_name]][row]
  }
  output <- data.frame(T = t_values, values, check.names = FALSE)
  utils::write.csv(output, destination, row.names = FALSE)
}

export_metric_tables <- function(results, output_dir, cases_df) {
  methods <- c("CCEMG", "BOOT")
  metrics <- c("AvgLen", "Coverage", "Bias", "RMSE", "Size", "Power", "MeanSE")

  for (index in seq_len(nrow(cases_df))) {
    case_definition <- cases_df[index, ]
    case_data <- results[results$case_id == case_definition$case_id, , drop = FALSE]
    alpha_values <- unique(case_data$alpha)

    for (alpha_value in alpha_values) {
      selected <- if (is.na(alpha_value)) is.na(case_data$alpha) else case_data$alpha == alpha_value
      alpha_data <- case_data[selected, , drop = FALSE]
      alpha_tag <- format_alpha_tag(alpha_value, case_definition$slope_mode)
      base_tag <- paste(case_definition$case_id, alpha_tag, sep = "_")

      for (metric in metrics) {
        for (method in methods) {
          column_name <- paste0(metric, "_", method)
          destination <- file.path(output_dir, paste0(base_tag, "_", column_name, ".csv"))
          write_metric_matrix(alpha_data, column_name, destination)
        }
      }
    }
  }
}

set.seed(read_int_env("SIM_SEED", 123))

N_values <- as.integer(read_num_vector_env("N_VALUES", c(50, 100, 200, 400), minimum = 2))
T_values <- as.integer(read_num_vector_env("T_VALUES", c(50, 100, 200, 400), minimum = 8))
alpha_values <- read_num_vector_env("ALPHA_VALUES", c(0, 0.5), minimum = 0)

replications <- read_int_env("SIM_REPS", 1000)
bootstrap_reps <- read_int_env("BOOT_REPS", 2000, minimum = 20)
beta_true <- read_num_env("BETA_TRUE", 1)
heterogeneity_sd <- read_num_env("HETEROGENEITY_SD", 1)
eps_sd <- read_num_env("EPS_SD", 1)
power_null <- read_num_env("POWER_NULL", 0.95)
power_half_width <- read_num_env("POWER_GRID_HALF_WIDTH", 0.5)
power_grid_points <- read_int_env("POWER_GRID_POINTS", 21, minimum = 3)
confidence_level <- read_num_env("CONFIDENCE_LEVEL", 0.95)
finite_sample_correction <- tolower(Sys.getenv(
  "DWB_FINITE_SAMPLE_CORRECTION", unset = "true"
)) %in% c("1", "true", "yes", "on")
interval_type <- tolower(Sys.getenv("BOOT_INTERVAL", unset = "asymmetric"))

if (confidence_level <= 0 || confidence_level >= 1) stop("CONFIDENCE_LEVEL must lie in (0, 1).")
if (!(interval_type %in% c("asymmetric", "basic", "percentile", "symmetric"))) {
  stop("BOOT_INTERVAL must be 'asymmetric' (or alias 'basic'/'percentile') or 'symmetric'.")
}

power_null_grid <- sort(unique(c(
  seq(beta_true - power_half_width, beta_true + power_half_width, length.out = power_grid_points),
  beta_true,
  power_null
)))

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
results_file <- Sys.getenv("RESULTS_FILE", unset = file.path(output_dir, "all_cases_results.csv"))
power_results_file <- Sys.getenv("POWER_RESULTS_FILE",
                                 unset = file.path(output_dir, "power_curve_results.csv"))

cases_df <- build_case_table()
results_list <- list()
power_list <- list()
row_index <- 1L

cat(sprintf("Monte Carlo replications: %d | bootstrap replications: %d | interval: %s\n",
            replications, bootstrap_reps, interval_type))

for (case_index in seq_len(nrow(cases_df))) {
  case_definition <- cases_df[case_index, ]
  alpha_grid <- if (case_definition$slope_mode == "heterogeneous") alpha_values else NA_real_

  for (alpha_value in alpha_grid) {
    alpha_label <- if (is.na(alpha_value)) "homogeneous" else format(alpha_value, trim = TRUE)
    cat(sprintf("Running %s | alpha = %s\n", case_definition$case_name, alpha_label))

    for (N in N_values) {
      for (Tn in T_values) {
        cat(sprintf("  N = %d, T = %d\n", N, Tn))
        cell <- run_single_cell(
          N = N, Tn = Tn, alpha = alpha_value,
          slope_mode = case_definition$slope_mode,
          rank_satisfied = case_definition$rank_satisfied,
          replications = replications, beta_true = beta_true,
          heterogeneity_sd = heterogeneity_sd, eps_sd = eps_sd,
          power_null = power_null, power_null_grid = power_null_grid,
          bootstrap_reps = bootstrap_reps, confidence_level = confidence_level,
          interval_type = interval_type,
          finite_sample_correction = finite_sample_correction
        )

        metadata <- data.frame(
          case_id = case_definition$case_id,
          case_name = case_definition$case_name,
          slope_mode = case_definition$slope_mode,
          rank_condition = if (case_definition$rank_satisfied) "satisfied" else "not_satisfied",
          stringsAsFactors = FALSE
        )

        results_list[[row_index]] <- cbind(metadata, cell$summary)
        power_list[[row_index]] <- cbind(metadata[rep(1L, nrow(cell$power_curve)), , drop = FALSE],
                                         cell$power_curve)
        row_index <- row_index + 1L
      }
    }
  }
}

results <- do.call(rbind, results_list)
power_results <- do.call(rbind, power_list)
rownames(results) <- NULL
rownames(power_results) <- NULL

utils::write.csv(results, results_file, row.names = FALSE)
utils::write.csv(power_results, power_results_file, row.names = FALSE)
export_metric_tables(results, output_dir, cases_df)

cat("Finished simulation and CSV exports.\n")
cat(sprintf("Summary results: %s\n", results_file))
cat(sprintf("Power-curve results: %s\n", power_results_file))
cat(sprintf("Output directory: %s\n", output_dir))
