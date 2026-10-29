### Plots for standard CCEMG and ABC whole-sum dependent wild bootstrap.

rm(list = ls())

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

library(ggplot2)

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "outputs")
results_file <- Sys.getenv("RESULTS_FILE", unset = file.path(output_dir, "all_cases_results.csv"))
power_results_file <- Sys.getenv("POWER_RESULTS_FILE",
                                 unset = file.path(output_dir, "power_curve_results.csv"))

if (!file.exists(results_file)) stop(sprintf("Summary results not found: %s", results_file))
if (!file.exists(power_results_file)) stop(sprintf("Power results not found: %s", power_results_file))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_file)
power_results <- utils::read.csv(power_results_file)

required_summary <- c("case_id", "case_name", "alpha", "N", "T",
                      "Coverage_CCEMG", "Coverage_BOOT",
                      "AvgLen_CCEMG", "AvgLen_BOOT")
required_power <- c("case_id", "case_name", "alpha", "N", "T",
                    "beta_true", "beta_null", "Power_CCEMG", "Power_BOOT")

missing_summary <- setdiff(required_summary, names(results))
missing_power <- setdiff(required_power, names(power_results))
if (length(missing_summary)) stop(sprintf("Missing summary columns: %s", paste(missing_summary, collapse = ", ")))
if (length(missing_power)) stop(sprintf("Missing power columns: %s", paste(missing_power, collapse = ", ")))

method_levels <- c("Standard CCEMG", "ABC + whole-sum DWB")
method_colors <- c("Standard CCEMG" = "#D55E00",
                   "ABC + whole-sum DWB" = "#0072B2")
confidence_level <- suppressWarnings(as.numeric(Sys.getenv("CONFIDENCE_LEVEL", unset = "0.95")))
if (!is.finite(confidence_level)) confidence_level <- 0.95

format_alpha_text <- function(alpha_value) {
  if (is.na(alpha_value)) "homogeneous" else format(alpha_value, trim = TRUE)
}

format_alpha_slug <- function(alpha_value) {
  if (is.na(alpha_value)) "alpha_homo" else paste0("alpha_", gsub("\\.", "p", format(alpha_value, trim = TRUE)))
}

make_method_long <- function(data, ccemg_column, bootstrap_column) {
  first <- data[, c("N", "T"), drop = FALSE]
  first$Method <- method_levels[1L]
  first$Value <- data[[ccemg_column]]

  second <- data[, c("N", "T"), drop = FALSE]
  second$Method <- method_levels[2L]
  second$Value <- data[[bootstrap_column]]

  combined <- rbind(first, second)
  combined$Method <- factor(combined$Method, levels = method_levels)
  combined$N <- factor(combined$N, levels = sort(unique(data$N)))
  combined
}

plot_one_case_alpha <- function(summary_data, curve_data, case_id, case_name, alpha_value) {
  alpha_text <- format_alpha_text(alpha_value)
  alpha_slug <- format_alpha_slug(alpha_value)
  t_values <- sort(unique(summary_data$T))
  n_count <- length(unique(summary_data$N))
  plot_width <- max(10, 3.1 * n_count)

  coverage_data <- make_method_long(summary_data, "Coverage_CCEMG", "Coverage_BOOT")
  coverage_plot <- ggplot(coverage_data, aes(x = T, y = Value, color = Method, group = Method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ N, nrow = 1,
               labeller = labeller(N = function(x) paste0("N = ", x))) +
    geom_hline(yintercept = confidence_level, linetype = "dashed", color = "#555555") +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = t_values) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = sprintf("%s | alpha = %s", case_name, alpha_text),
         x = "T", y = "95% confidence-interval coverage", color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA))

  coverage_file <- file.path(output_dir, paste0(case_id, "_", alpha_slug, "_coverage.png"))
  ggsave(coverage_file, coverage_plot, width = plot_width, height = 5.7,
         dpi = 300, bg = "white")

  length_data <- make_method_long(summary_data, "AvgLen_CCEMG", "AvgLen_BOOT")
  length_plot <- ggplot(length_data, aes(x = T, y = Value, color = Method, group = Method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ N, nrow = 1,
               labeller = labeller(N = function(x) paste0("N = ", x))) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = t_values) +
    labs(title = sprintf("%s | alpha = %s", case_name, alpha_text),
         x = "T", y = "Average 95% confidence-interval length", color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA))

  length_file <- file.path(output_dir, paste0(case_id, "_", alpha_slug, "_avglen.png"))
  ggsave(length_file, length_plot, width = plot_width, height = 5.7,
         dpi = 300, bg = "white")

  curve_standard <- curve_data[, c("N", "T", "beta_true", "beta_null"), drop = FALSE]
  curve_standard$Method <- method_levels[1L]
  curve_standard$Power <- curve_data$Power_CCEMG

  curve_bootstrap <- curve_data[, c("N", "T", "beta_true", "beta_null"), drop = FALSE]
  curve_bootstrap$Method <- method_levels[2L]
  curve_bootstrap$Power <- curve_data$Power_BOOT

  curve_long <- rbind(curve_standard, curve_bootstrap)
  curve_long$Method <- factor(curve_long$Method, levels = method_levels)
  curve_long$N <- factor(curve_long$N, levels = sort(unique(curve_long$N)))
  curve_long$T <- factor(curve_long$T, levels = sort(unique(curve_long$T)))

  power_plot <- ggplot(curve_long, aes(x = beta_null, y = Power,
                                      color = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    geom_vline(aes(xintercept = beta_true), linetype = "dashed",
               color = "#777777", linewidth = 0.4) +
    geom_hline(yintercept = 1 - confidence_level, linetype = "dotted",
               color = "#777777", linewidth = 0.4) +
    facet_grid(rows = vars(N), cols = vars(T),
               labeller = labeller(N = function(x) paste0("N = ", x),
                                   T = function(x) paste0("T = ", x))) +
    scale_color_manual(values = method_colors) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = sprintf("%s | alpha = %s", case_name, alpha_text),
         x = expression("Null value " * beta[0]),
         y = "Rejection probability", color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA))

  power_file <- file.path(output_dir, paste0(case_id, "_", alpha_slug, "_power_curves.png"))
  ggsave(power_file, power_plot,
         width = max(10, 2.8 * length(unique(curve_long$T))),
         height = max(7, 2.1 * length(unique(curve_long$N))),
         dpi = 300, bg = "white")

  cat(sprintf("Saved coverage, interval-length, and power plots for %s (%s).\n",
              case_id, alpha_slug))
}

groups <- unique(results[, c("case_id", "case_name", "alpha")])
for (index in seq_len(nrow(groups))) {
  case_id <- groups$case_id[index]
  alpha_value <- groups$alpha[index]
  summary_selected <- (results$case_id == case_id) &
    (if (is.na(alpha_value)) is.na(results$alpha) else results$alpha == alpha_value)
  curve_selected <- (power_results$case_id == case_id) &
    (if (is.na(alpha_value)) is.na(power_results$alpha) else power_results$alpha == alpha_value)

  plot_one_case_alpha(results[summary_selected, , drop = FALSE],
                      power_results[curve_selected, , drop = FALSE],
                      case_id, groups$case_name[index], alpha_value)
}

cat("Finished all plot exports.\n")
cat(sprintf("Output directory: %s\n", output_dir))
