#!/usr/bin/env python3
"""Executable reference for CCEMG versus ABC whole-sum DWB with l_N = 1.

The R scripts are the primary implementation. This companion verifies the
single-sum unit influence, independent-multiplier covariance, and confidence
interval calculations, and can generate reproducible outputs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import lfilter
from scipy.stats import norm


COLORS = {"Standard CCEMG": "#D55E00", "ABC + whole-sum DWB": "#0072B2"}


@dataclass
class CCEFit:
    betas: np.ndarray
    gammas: np.ndarray
    Gammas: np.ndarray
    residual_y: np.ndarray
    residual_x: np.ndarray
    q_values: np.ndarray
    y_bar: np.ndarray
    x_bar: np.ndarray


@dataclass
class ABCResult:
    beta_abc: float
    unit_corrected_slopes: np.ndarray
    unit_bias_adjustment: np.ndarray
    b0: np.ndarray
    b1: np.ndarray


def safe_inverse(values: np.ndarray | float, eps: float = 1e-10) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    protected = np.where(np.abs(array) < eps, np.sign(array + eps) * eps, array)
    return np.where(np.isfinite(protected), 1.0 / protected, 0.0)


def fit_cce_panel(y: np.ndarray, x: np.ndarray, rows: np.ndarray | None = None) -> CCEFit:
    if rows is not None:
        y, x = y[rows], x[rows]
    t_count = y.shape[0]
    if t_count < 4:
        raise ValueError("At least four periods are needed for the individual CCE regression.")

    y_bar, x_bar = y.mean(axis=1), x.mean(axis=1)
    controls = np.column_stack((np.ones(t_count), y_bar, x_bar))
    coef_x = np.linalg.lstsq(controls, x, rcond=1e-10)[0]
    coef_y = np.linalg.lstsq(controls, y, rcond=1e-10)[0]
    x_tilde = x - controls @ coef_x
    y_tilde = y - controls @ coef_y
    denominator = np.sum(x_tilde**2, axis=0)
    numerator = np.sum(x_tilde * y_tilde, axis=0)
    betas = np.divide(numerator, denominator, out=np.zeros_like(numerator),
                      where=denominator > 1e-10)
    gamma_controls = np.linalg.lstsq(controls, y - x * betas, rcond=1e-10)[0]
    residual_y = y_tilde - x_tilde * betas
    return CCEFit(betas, gamma_controls[1:3].T, coef_x[1:3].T,
                  residual_y, x_tilde, denominator / t_count, y_bar, x_bar)


def analytical_bias_correction(fit: CCEFit) -> ABCResult:
    """Feasible Equation (44) correction using CCE residual variance estimates."""
    t_count, n_count = fit.residual_x.shape
    c_bar = np.vstack((np.mean(fit.gammas + fit.Gammas * fit.betas[:, None], axis=0),
                       np.mean(fit.Gammas, axis=0)))
    cplus_t = np.linalg.pinv(c_bar).T
    inverse_denom = safe_inverse(fit.q_values)
    sig_v_hat = np.sum(fit.residual_x**2, axis=0) / max(1, t_count - 3)
    sig_eps_hat = np.sum(fit.residual_y**2, axis=0) / max(1, t_count - 4)

    b0 = np.empty(n_count)
    b1 = np.empty(n_count)
    for index in range(n_count):
        b0[index] = inverse_denom[index] * sig_v_hat[index] * (
            np.array([fit.betas[index], 1.0]) @ cplus_t @ fit.gammas[index]
        )
        b1[index] = inverse_denom[index] * (
            fit.Gammas[index] @ cplus_t @ np.array([sig_eps_hat[index], 1.0])
        )
    bias = np.sqrt(t_count / n_count) * (b0.mean() + b1.mean())
    # Equation (44): subtract the estimated leading bias.
    beta_abc = fit.betas.mean() - bias / np.sqrt(n_count * t_count)
    unit_bias_adjustment = -(b0 + b1) / n_count
    unit_corrected_slopes = fit.betas + unit_bias_adjustment
    if not np.isclose(unit_corrected_slopes.mean(), beta_abc, rtol=1e-10, atol=1e-10):
        raise AssertionError("Unit-level ABC contributions do not reproduce beta_ABC.")
    return ABCResult(float(beta_abc), unit_corrected_slopes, unit_bias_adjustment, b0, b1)


def bartlett_draws(unit_influence: np.ndarray, bandwidth: int, bootstrap_reps: int,
                   rng: np.random.Generator) -> np.ndarray:
    n_count = unit_influence.size
    coefficients = np.zeros(n_count + bandwidth - 1)
    for offset in range(bandwidth):
        coefficients[offset:offset + n_count] += unit_influence
    innovations = rng.normal(size=(coefficients.size, bootstrap_reps))
    return coefficients @ innovations / np.sqrt(n_count * bandwidth)


def bootstrap_interval(fit: CCEFit, abc: ABCResult,
                       bootstrap_reps: int, confidence: float, rng: np.random.Generator,
                       interval_type: str = "asymmetric",
                       finite_sample_correction: bool = True) -> dict[str, float]:
    t_count, n_count = fit.residual_x.shape
    unit_influence = np.sqrt(t_count) * (abc.unit_corrected_slopes - abc.beta_abc)
    if finite_sample_correction and n_count > 1:
        unit_influence *= np.sqrt(n_count / (n_count - 1))
    bandwidth = 1
    draws = bartlett_draws(unit_influence, bandwidth, bootstrap_reps, rng)
    unit_distance = np.abs(np.subtract.outer(np.arange(n_count), np.arange(n_count)))
    bartlett_kernel = np.maximum(1 - unit_distance / bandwidth, 0)
    conditional_sigma = unit_influence @ bartlett_kernel @ unit_influence / n_count
    estimator_variance = conditional_sigma / (n_count * t_count)
    corrected_centered = abc.unit_corrected_slopes - abc.beta_abc
    variance_denominator = (
        n_count * (n_count - 1)
        if finite_sample_correction else n_count**2
    )
    corrected_slope_variance = corrected_centered @ corrected_centered / variance_denominator
    if not np.isclose(estimator_variance, corrected_slope_variance,
                      rtol=1e-12, atol=1e-14):
        raise AssertionError("The l_N=1 DWB covariance identity failed.")
    # Asymmetric basic-bootstrap interval.  The bootstrap draws approximate
    # sqrt(NT)*(beta_ABC-beta_0), so their upper quantile is subtracted for the
    # lower endpoint and their lower quantile for the upper endpoint.
    if interval_type in {"asymmetric", "basic", "percentile"}:
        low, high = np.quantile(draws, [(1 - confidence) / 2, (1 + confidence) / 2])
        lower = abc.beta_abc - high / np.sqrt(n_count * t_count)
        upper = abc.beta_abc - low / np.sqrt(n_count * t_count)
    else:
        critical = np.quantile(np.abs(draws), confidence)
        lower = abc.beta_abc - critical / np.sqrt(n_count * t_count)
        upper = abc.beta_abc + critical / np.sqrt(n_count * t_count)
    if not (np.isfinite(lower) and np.isfinite(upper) and lower <= upper):
        raise AssertionError("ABC bootstrap confidence interval endpoints are invalid.")
    return {"lower": float(lower), "upper": float(upper),
            "se": float(np.sqrt(max(estimator_variance, 0.0))),
            "conditional_sigma": float(conditional_sigma),
            "estimator_variance": float(estimator_variance),
            "bootstrap_sigma_mc": float(np.var(draws, ddof=1)),
            "unit_score_variance": float(np.var(unit_influence, ddof=1)),
            "bandwidth": bandwidth}


def stationary_ar1(rho: float, length: int, rng: np.random.Generator,
                   burn_in: int = 150) -> np.ndarray:
    innovations = rng.normal(size=length + burn_in)
    return lfilter([1.0], [1.0, -rho], innovations)[burn_in:]


def generate_panel(n_count: int, t_count: int, alpha: float,
                   heterogeneous: bool, rank_satisfied: bool,
                   beta_true: float, heterogeneity_sd: float, eps_sd: float,
                   rng: np.random.Generator) -> tuple[np.ndarray, ...]:
    f1 = stationary_ar1(0.5, t_count, rng)
    f2 = stationary_ar1(0.5, t_count, rng)
    beta_i = (beta_true + heterogeneity_sd * rng.normal(size=n_count) / n_count**alpha
              if heterogeneous else np.full(n_count, beta_true))

    gamma1 = rng.normal(1.0, 0.2, n_count)
    gamma2 = rng.normal(0.0 if rank_satisfied else 1.0, 0.2, n_count)
    Gamma1 = rng.normal(0.5, 0.5, n_count)
    Gamma2 = rng.normal(0.5, 0.5, n_count)
    a_i = rng.normal(size=n_count)
    A_i = rng.normal(size=n_count)

    x = np.empty((t_count, n_count))
    y = np.empty((t_count, n_count))
    for index in range(n_count):
        v = stationary_ar1(rng.uniform(0.05, 0.95), t_count, rng)
        eps = eps_sd * stationary_ar1(rng.uniform(0.05, 0.95), t_count, rng)
        x[:, index] = A_i[index] + Gamma1[index] * f1 + Gamma2[index] * f2 + v
        y[:, index] = (
            a_i[index] + beta_i[index] * x[:, index]
            + gamma1[index] * f1 + gamma2[index] * f2 + eps
        )
    return y, x


def self_tests() -> None:
    package_directory = Path(__file__).resolve().parent
    for source_name in ("MCsimuforshrinkingheterogeneity.R",
                        "data_visualisation_shrinkingheterogeneity.R"):
        source = (package_directory / source_name).read_text(encoding="utf-8")
        stack: list[tuple[str, int]] = []
        quotation: str | None = None
        comment = False
        escaped = False
        line = 1
        matching = {")": "(", "]": "[", "}": "{"}
        for character in source:
            if character == "\n":
                line += 1
                comment = False
                if quotation is None:
                    continue
            if comment:
                continue
            if quotation is not None:
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == quotation:
                    quotation = None
                continue
            if character == "#":
                comment = True
            elif character in ('"', "'", "`"):
                quotation = character
            elif character in "([{":
                stack.append((character, line))
            elif character in ")]}":
                if not stack or stack[-1][0] != matching[character]:
                    raise AssertionError(f"Unmatched R delimiter {character!r} in {source_name}:{line}.")
                stack.pop()
        if quotation is not None or stack:
            raise AssertionError(f"Unclosed quote or delimiter in {source_name}: {stack}.")
        if any(marker in source for marker in ("JCK", "jackknife", "var_ccemg_proposed",
                                                "hybrid_bootstrap_inference",
                                                "estimate_omega_v_split")):
            raise AssertionError(f"Unexpected removed-estimator code remains in {source_name}.")
        if source_name == "MCsimuforshrinkingheterogeneity.R":
            required_fragments = (
                "beta_mg + z_critical * se_standard",
                "bandwidth <- 1L",
                "standard_error <- sqrt(max(conditional_sigma, 0) / (N * Tn))",
                "abc <- analytical_bias_correction(fit)",
            )
            missing_fragments = [
                fragment for fragment in required_fragments if fragment not in source
            ]
            if missing_fragments:
                raise AssertionError(
                    "R implementation is missing required corrected paths: "
                    + ", ".join(missing_fragments)
                )

    rng = np.random.default_rng(217)
    t_count, n_count = 60, 20
    y, x = generate_panel(n_count, t_count, 0.5, True, True,
                          1.0, 1.0, 1.0, rng)
    fit = fit_cce_panel(y, x)
    direct_design = np.column_stack((np.ones(t_count), x[:, 0], fit.y_bar, fit.x_bar))
    direct_beta = np.linalg.lstsq(direct_design, y[:, 0], rcond=1e-10)[0][1]
    if not np.isclose(fit.betas[0], direct_beta, rtol=1e-10, atol=1e-10):
        raise AssertionError("FWL individual CCE slope does not match direct OLS.")

    score = rng.normal(size=n_count)
    score -= score.mean()
    bandwidth = 5
    draws = bartlett_draws(score, bandwidth, 30000, rng)
    lag = np.abs(np.subtract.outer(np.arange(n_count), np.arange(n_count)))
    kernel = np.maximum(1 - lag / bandwidth, 0)
    theoretical_variance = score @ kernel @ score / n_count
    simulated_variance = np.var(draws, ddof=1)
    if not np.isclose(simulated_variance, theoretical_variance, rtol=0.05):
        raise AssertionError("Bartlett bootstrap draw variance does not match its kernel covariance.")

    abc = analytical_bias_correction(fit)
    if not np.isclose(abc.unit_corrected_slopes.mean(), abc.beta_abc, atol=1e-10):
        raise AssertionError("Corrected unit slopes do not average to the ABC estimator.")
    result = bootstrap_interval(fit, abc, 299, 0.95, rng)
    if not (np.isfinite(abc.beta_abc) and result["lower"] <= result["upper"] and result["se"] >= 0):
        raise AssertionError("The whole-sum DWB did not produce a valid interval.")

    # Check the requested quantile orientation independently of bootstrap RNG.
    toy_draws = np.array([-3.0, -1.0, 0.5, 2.0, 7.0])
    toy_low, toy_high = np.quantile(toy_draws, [0.025, 0.975])
    toy_center, toy_scale = 1.25, np.sqrt(40.0)
    toy_interval = (toy_center - toy_high / toy_scale,
                    toy_center - toy_low / toy_scale)
    if not toy_interval[0] < toy_interval[1]:
        raise AssertionError("Asymmetric bootstrap quantiles have the wrong endpoint orientation.")

    iid_result = bootstrap_interval(fit, abc, 10000, 0.95, rng)
    corrected_centered = abc.unit_corrected_slopes - abc.beta_abc
    exact_iid_sigma = t_count * np.sum(corrected_centered**2) / (n_count - 1)
    if not np.isclose(iid_result["conditional_sigma"], exact_iid_sigma,
                      rtol=1e-10, atol=1e-10):
        raise AssertionError("Bandwidth-one DWB does not equal the corrected Pesaran covariance.")
    exact_iid_se = np.sqrt(exact_iid_sigma / (n_count * t_count))
    if iid_result["bandwidth"] != 1 or not np.isclose(
        iid_result["se"], exact_iid_se, rtol=1e-12, atol=1e-12
    ):
        raise AssertionError("Reported ABC SE does not use the exact l_N=1 covariance.")

    print("PASS: R source delimiters are balanced and only the requested estimator paths remain.")
    print(f"PASS: FWL CCE slope agrees with direct OLS ({direct_beta:.6f}).")
    print(f"PASS: cross-sectional Bartlett DWB variance {simulated_variance:.6f}; "
          f"target {theoretical_variance:.6f}.")
    print(f"PASS: corrected unit slopes average to beta_ABC = {abc.beta_abc:.6f}.")
    print(f"PASS: bandwidth-one DWB conditional variance matches {exact_iid_sigma:.6f}.")
    print("PASS: analytical bias correction and whole-sum DWB interval are finite.")


def simulation_cell(case: dict[str, object], alpha: float, n_count: int, t_count: int,
                    args: argparse.Namespace, null_grid: np.ndarray,
                    rng: np.random.Generator) -> tuple[dict[str, object], pd.DataFrame]:
    est_standard, est_abc = [], []
    cover_standard, cover_boot = [], []
    length_standard, length_boot = [], []
    se_standard_values, se_boot_values = [], []
    conditional_sigma_values, estimator_variance_values = [], []
    bootstrap_sigma_mc_values, unit_score_variance_values, bandwidth_values = [], [], []
    rejection_standard = np.empty((args.reps, null_grid.size), dtype=bool)
    rejection_boot = np.empty_like(rejection_standard)
    z_critical = norm.ppf((1 + args.confidence) / 2)

    for replication in range(args.reps):
        y, x = generate_panel(
            n_count, t_count, alpha, bool(case["heterogeneous"]),
            bool(case["rank_satisfied"]), args.beta_true, args.heterogeneity_sd,
            args.eps_sd, rng
        )
        fit = fit_cce_panel(y, x)
        beta_mg = fit.betas.mean()
        abc = analytical_bias_correction(fit)
        se_standard = np.sqrt(np.sum((fit.betas - beta_mg)**2) / (n_count * (n_count - 1)))
        ci_standard = np.array([beta_mg - z_critical * se_standard,
                                beta_mg + z_critical * se_standard])
        boot = bootstrap_interval(
            fit, abc, args.bootstrap_reps, args.confidence, rng, args.interval,
            not args.no_finite_sample_correction
        )

        est_standard.append(beta_mg)
        est_abc.append(abc.beta_abc)
        cover_standard.append(ci_standard[0] <= args.beta_true <= ci_standard[1])
        cover_boot.append(boot["lower"] <= args.beta_true <= boot["upper"])
        length_standard.append(ci_standard[1] - ci_standard[0])
        length_boot.append(boot["upper"] - boot["lower"])
        se_standard_values.append(se_standard)
        se_boot_values.append(boot["se"])
        conditional_sigma_values.append(boot["conditional_sigma"])
        estimator_variance_values.append(boot["estimator_variance"])
        bootstrap_sigma_mc_values.append(boot["bootstrap_sigma_mc"])
        unit_score_variance_values.append(boot["unit_score_variance"])
        bandwidth_values.append(boot["bandwidth"])
        rejection_standard[replication] = (null_grid < ci_standard[0]) | (null_grid > ci_standard[1])
        rejection_boot[replication] = (null_grid < boot["lower"]) | (null_grid > boot["upper"])

    est_standard, est_abc = np.asarray(est_standard), np.asarray(est_abc)
    standard_power, boot_power = rejection_standard.mean(axis=0), rejection_boot.mean(axis=0)
    null_position = int(np.argmin(np.abs(null_grid - args.power_null)))
    alpha_value = alpha if bool(case["heterogeneous"]) else np.nan
    common = {"case_id": case["case_id"], "case_name": case["case_name"],
              "slope_mode": "heterogeneous" if bool(case["heterogeneous"]) else "homogeneous",
              "rank_condition": "satisfied" if bool(case["rank_satisfied"]) else "not_satisfied",
              "alpha": alpha_value, "N": n_count, "T": t_count}
    summary = {**common,
               "Bias_CCEMG": np.sqrt(n_count * t_count) * np.mean(est_standard - args.beta_true),
               "Bias_BOOT": np.sqrt(n_count * t_count) * np.mean(est_abc - args.beta_true),
               "RMSE_CCEMG": np.sqrt(np.mean((est_standard - args.beta_true)**2)),
               "RMSE_BOOT": np.sqrt(np.mean((est_abc - args.beta_true)**2)),
               "Coverage_CCEMG": np.mean(cover_standard), "Coverage_BOOT": np.mean(cover_boot),
               "AvgLen_CCEMG": np.mean(length_standard), "AvgLen_BOOT": np.mean(length_boot),
               "Size_CCEMG": 1 - np.mean(cover_standard), "Size_BOOT": 1 - np.mean(cover_boot),
               "Power_CCEMG": standard_power[null_position], "Power_BOOT": boot_power[null_position],
               "MeanSE_CCEMG": np.mean(se_standard_values), "MeanSE_BOOT": np.mean(se_boot_values),
               "MeanConditionalSigma_DWB": np.mean(conditional_sigma_values),
               "MeanEstimatorVariance_BOOT": np.mean(estimator_variance_values),
               "MeanBootstrapSigma_MC": np.mean(bootstrap_sigma_mc_values),
               "MeanUnitInfluenceVariance": np.mean(unit_score_variance_values),
               "MeanDWB_Bandwidth": np.mean(bandwidth_values), "Bootstrap_Reps": args.bootstrap_reps,
               "DWB_Dimension": "cross_section",
               "Simulation_Reps": args.reps, "Beta_True": args.beta_true,
               "Reference_Power_Null": null_grid[null_position]}
    curves = pd.DataFrame({**common, "beta_true": args.beta_true, "beta_null": null_grid,
                           "Power_CCEMG": standard_power, "Power_BOOT": boot_power})
    return summary, curves


def export_plots(summary: pd.DataFrame, power: pd.DataFrame,
                 output_directory: Path, confidence: float) -> None:
    groups = summary[["case_id", "case_name", "alpha"]].drop_duplicates()
    for _, group in groups.iterrows():
        alpha = group["alpha"]
        summary_mask = summary.case_id.eq(group.case_id) & (
            summary.alpha.isna() if pd.isna(alpha) else summary.alpha.eq(alpha)
        )
        power_mask = power.case_id.eq(group.case_id) & (
            power.alpha.isna() if pd.isna(alpha) else power.alpha.eq(alpha)
        )
        data, curve = summary.loc[summary_mask], power.loc[power_mask]
        alpha_text = "homogeneous" if pd.isna(alpha) else f"{alpha:g}"
        alpha_slug = "alpha_homo" if pd.isna(alpha) else f"alpha_{alpha:g}".replace(".", "p")
        n_values = np.sort(data.N.unique())
        t_values = np.sort(data["T"].unique())

        for metric, ylabel, filename in [
            ("Coverage", "95% confidence-interval coverage", "coverage"),
            ("AvgLen", "Average 95% confidence-interval length", "avglen"),
        ]:
            figure, axes = plt.subplots(1, len(n_values), figsize=(4.0 * len(n_values), 4.0),
                                        squeeze=False)
            for axis, n_count in zip(axes[0], n_values):
                panel = data.loc[data.N.eq(n_count)].sort_values("T")
                axis.plot(panel["T"].values,
                          panel[f"{metric}_CCEMG"], marker="o", color=COLORS["Standard CCEMG"],
                          label="Standard CCEMG")
                axis.plot(panel["T"].values, panel[f"{metric}_BOOT"], marker="s",
                          color=COLORS["ABC + whole-sum DWB"], label="ABC + whole-sum DWB")
                if metric == "Coverage":
                    axis.axhline(confidence, color="0.4", linestyle="--", linewidth=0.9)
                    axis.set_ylim(0, 1.03)
                axis.set_title(f"N = {n_count}")
                axis.set_xlabel("T")
                axis.grid(alpha=0.22)
            axes[0, 0].set_ylabel(ylabel)
            figure.suptitle(f"{group.case_name} | alpha = {alpha_text}")
            handles, labels = axes[0, 0].get_legend_handles_labels()
            figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
            figure.tight_layout(rect=(0, 0.08, 1, 0.92))
            figure.savefig(output_directory / f"{group.case_id}_{alpha_slug}_{filename}.png",
                           dpi=160, facecolor="white")
            plt.close(figure)

        figure, axes = plt.subplots(len(n_values), len(t_values),
                                    figsize=(3.2 * len(t_values), 2.55 * len(n_values)),
                                    squeeze=False, sharex=True, sharey=True)
        for row, n_count in enumerate(n_values):
            for col, t_count in enumerate(t_values):
                axis = axes[row, col]
                panel = curve.loc[curve["N"].eq(n_count) &
                                  curve["T"].eq(t_count)].sort_values("beta_null")
                axis.plot(panel.beta_null, panel.Power_CCEMG, color=COLORS["Standard CCEMG"],
                          label="Standard CCEMG")
                axis.plot(panel.beta_null, panel.Power_BOOT,
                          color=COLORS["ABC + whole-sum DWB"], label="ABC + whole-sum DWB")
                axis.axvline(float(panel.beta_true.iloc[0]), color="0.45", linestyle="--", linewidth=0.8)
                axis.axhline(1 - confidence, color="0.45", linestyle=":", linewidth=0.8)
                axis.set_ylim(0, 1.03)
                axis.set_title(f"N = {n_count}, T = {t_count}", fontsize=9)
                axis.grid(alpha=0.18)
                if row == len(n_values) - 1:
                    axis.set_xlabel("Null value")
                if col == 0:
                    axis.set_ylabel("Rejection probability")
        figure.suptitle(f"{group.case_name} | alpha = {alpha_text}")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
        figure.tight_layout(rect=(0, 0.08, 1, 0.93))
        figure.savefig(output_directory / f"{group.case_id}_{alpha_slug}_power_curves.png",
                       dpi=160, facecolor="white")
        plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true", help="Run numerical self-tests only.")
    parser.add_argument("--output-dir", default="example_outputs", help="Output directory.")
    parser.add_argument("--reps", type=int, default=1000, help="Monte Carlo replications per cell.")
    parser.add_argument("--bootstrap-reps", type=int, default=2000, help="Bootstrap draws per replication.")
    parser.add_argument("--n-values", default="50,100,200,400", help="Comma-separated N values.")
    parser.add_argument("--t-values", default="50,100,200,400", help="Comma-separated T values.")
    parser.add_argument("--alpha-values", default="0,0.5", help="Comma-separated heterogeneity exponents.")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--beta-true", type=float, default=1.0)
    parser.add_argument("--heterogeneity-sd", type=float, default=1.0)
    parser.add_argument("--eps-sd", type=float, default=1.0)
    parser.add_argument("--power-null", type=float, default=0.95)
    parser.add_argument("--power-grid-half-width", type=float, default=0.5)
    parser.add_argument("--power-grid-points", type=int, default=21)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--interval",
                        choices=("asymmetric", "basic", "percentile", "symmetric"),
                        default="asymmetric",
                        help=("Bootstrap CI: asymmetric/basic uses "
                              "[beta-q_(1-alpha/2), beta-q_(alpha/2)]/sqrt(NT); "
                              "'percentile' is retained as a legacy alias."))
    parser.add_argument("--no-finite-sample-correction", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    self_tests()
    if args.verify_only:
        return

    output_directory = Path(args.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    n_values = [int(value) for value in args.n_values.split(",")]
    t_values = [int(value) for value in args.t_values.split(",")]
    alpha_values = [float(value) for value in args.alpha_values.split(",")]
    if min(n_values) < 2:
        raise ValueError("All N values must be at least 2.")
    if min(t_values) < 8:
        raise ValueError("All T values must be at least 8.")
    if args.reps < 1:
        raise ValueError("--reps must be at least 1.")
    if args.bootstrap_reps < 20:
        raise ValueError("--bootstrap-reps must be at least 20.")
    if not 0 < args.confidence < 1:
        raise ValueError("--confidence must lie in (0, 1).")

    null_grid = np.unique(np.append(np.linspace(args.beta_true - args.power_grid_half_width,
                                                 args.beta_true + args.power_grid_half_width,
                                                 args.power_grid_points),
                                    [args.beta_true, args.power_null]))

    cases = [
        {"case_id": "case1_hetero_rank_satisfied", "case_name": "Heterogeneous; rank satisfied",
         "heterogeneous": True, "rank_satisfied": True},
        {"case_id": "case2_hetero_rank_not_satisfied", "case_name": "Heterogeneous; rank not satisfied",
         "heterogeneous": True, "rank_satisfied": False},
        {"case_id": "case3_homo_rank_satisfied", "case_name": "Homogeneous; rank satisfied",
         "heterogeneous": False, "rank_satisfied": True},
        {"case_id": "case4_homo_rank_not_satisfied", "case_name": "Homogeneous; rank not satisfied",
         "heterogeneous": False, "rank_satisfied": False},
    ]

    summaries, curves = [], []
    for case in cases:
        exponents = alpha_values if case["heterogeneous"] else [0.5]
        for alpha in exponents:
            for n_count in n_values:
                for t_count in t_values:
                    print(f"Simulating {case['case_id']}, alpha={alpha}, N={n_count}, T={t_count}")
                    summary, curve = simulation_cell(case, alpha, n_count, t_count,
                                                     args, null_grid, rng)
                    summaries.append(summary)
                    curves.append(curve)

    summary_results = pd.DataFrame(summaries)
    power_results = pd.concat(curves, ignore_index=True)
    summary_results.to_csv(output_directory / "all_cases_results.csv", index=False)
    power_results.to_csv(output_directory / "power_curve_results.csv", index=False)

    metrics = ("Coverage", "AvgLen", "Bias", "RMSE", "Size", "Power", "MeanSE")
    for (case_id, alpha), group in summary_results.groupby(["case_id", "alpha"], dropna=False):
        alpha_slug = "alpha_homo" if pd.isna(alpha) else f"alpha_{alpha:g}".replace(".", "p")
        for metric in metrics:
            for method in ("CCEMG", "BOOT"):
                table = group.pivot(index="T", columns="N", values=f"{metric}_{method}")
                table.to_csv(output_directory / f"{case_id}_{alpha_slug}_{metric}_{method}.csv")

    export_plots(summary_results, power_results, output_directory, args.confidence)
    print(f"Wrote {len(summary_results)} summary rows and {len(power_results)} power rows.")
    print(f"Example results and plots: {output_directory.resolve()}")


if __name__ == "__main__":
    main()
