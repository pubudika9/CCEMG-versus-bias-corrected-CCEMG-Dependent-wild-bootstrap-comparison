# Corrected CCEMG versus ABC whole-sum DWB simulation (`l_N = 1`)

This package compares:

1. standard CCEMG with Pesaran's cross-sectional mean-group variance estimator; and
2. analytically bias-corrected CCEMG with a whole-sum cross-sectional multiplier bootstrap fixed at `l_N = 1`.

The jackknife method is not included.

## Corrections implemented

- The CCEMG confidence interval is correctly formed as

  ```text
  [beta_MG - z_(1-alpha/2) se_MG,
   beta_MG + z_(1-alpha/2) se_MG].
  ```

  The earlier duplicated lower endpoint, which forced zero interval lengths, zero coverage and rejection probability one, has been removed.

- The analytical bias correction uses the subtraction sign in Equation (44):

  ```text
  beta_ABC = beta_MG - (b0_hat + b1_hat)/sqrt(NT).
  ```

- The corrected individual slopes use the same sign and exactly average to `beta_ABC`:

  ```text
  beta_ABC,i = beta_hat_i - (b0_i + b1_i)/N.
  ```

- The whole-sum DWB bandwidth is fixed at `l_N = 1`. No automatic bandwidth or arbitrary cross-sectional ordering is used.

- With the finite-sample correction,

  ```text
  Z_i^c = sqrt[N/(N-1)] sqrt(T) (beta_ABC,i - beta_ABC).
  ```

  The exact conditional covariance and estimator variance are

  ```text
  Sigma_DWB = N^(-1) sum_i Z_i^c Z_i^c,

  Var(beta_ABC) = Sigma_DWB/(NT)
                = [N(N-1)]^(-1) sum_i
                  (beta_ABC,i - beta_ABC)^2.
  ```

  The reported ABC standard error uses this exact covariance. The empirical variance of the multiplier draws is retained only as a numerical diagnostic.

- The two bootstrap tails are evaluated separately. The asymmetric basic interval is

  ```text
  [beta_ABC - q*_(1-alpha/2)/sqrt(NT),
   beta_ABC - q*_(alpha/2)/sqrt(NT)].
  ```

- The projection used in the correction is consistent with the intercept-augmented CCE regression.
- Nuisance variances in the ABC correction are estimated from CCE residuals; latent simulated errors are not supplied to the estimator.
- The DGP now includes the individual fixed effects `a_i` and `A_i`, uses unit-variance idiosyncratic innovations by default, and uses the Section 4.1 loading distributions:

  ```text
  Rank satisfied:
    gamma_i1 ~ N(1, 0.2^2), gamma_i2 ~ N(0, 0.2^2)
    Gamma_i1 ~ N(0.5, 0.5^2), Gamma_i2 ~ N(0.5, 0.5^2)

  Rank deficient:
    gamma_i1 ~ N(1, 0.2^2), gamma_i2 ~ N(1, 0.2^2)
    Gamma_i1 ~ N(0.5, 0.5^2), Gamma_i2 ~ N(0.5, 0.5^2).
  ```

## Default experiment

- `N,T in {50,100,200,400}`
- heterogeneous designs: `alpha in {0,0.5}`
- homogeneous rank-satisfied and rank-deficient designs
- 1,000 Monte Carlo replications per cell
- 2,000 multiplier draws per replication
- true mean slope `beta_0 = 1`
- 95% confidence intervals
- `l_N = 1`
- random seed 123

## Run the R implementation

Install the required packages:

```r
install.packages(c("MASS", "ggplot2"))
```

Run:

```bash
bash run_hetero_pipeline.sh
```

Quick smoke test:

```bash
SIM_REPS=50 BOOT_REPS=199 N_VALUES=50 T_VALUES=50 \
  OUTPUT_DIR=outputs_smoke bash run_hetero_pipeline.sh
```

## Run the Python implementation

The Python implementation contains matching estimation, simulation, output and plotting paths:

```bash
python3 python_reference_simulation.py --verify-only
```

Quick simulation:

```bash
python3 python_reference_simulation.py \
  --reps 50 --bootstrap-reps 199 \
  --n-values 50 --t-values 50 \
  --output-dir validation_outputs
```

Dependencies are `numpy`, `scipy`, `pandas`, and `matplotlib`.

## Configuration

| Variable | Default | Meaning |
|---|---:|---|
| `SIM_REPS` | `1000` | Monte Carlo replications per cell |
| `BOOT_REPS` | `2000` | Multiplier draws per replication |
| `SIM_SEED` | `123` | Random seed |
| `N_VALUES` | `50,100,200,400` | Cross-sectional sizes |
| `T_VALUES` | `50,100,200,400` | Time dimensions |
| `ALPHA_VALUES` | `0,0.5` | Heterogeneity exponents |
| `BETA_TRUE` | `1` | Population-average slope |
| `HETEROGENEITY_SD` | `1` | Slope-innovation standard deviation |
| `EPS_SD` | `1` | Idiosyncratic-error innovation scale |
| `CONFIDENCE_LEVEL` | `0.95` | Confidence level |
| `BOOT_INTERVAL` | `asymmetric` | Basic interval with separate lower/upper quantiles |
| `DWB_FINITE_SAMPLE_CORRECTION` | `true` | Apply `sqrt[N/(N-1)]` |
| `OUTPUT_DIR` | `outputs` | Output directory |

`l_N=1` is fixed in both implementations and is not a tunable configuration setting.

## Output

- `all_cases_results.csv`
- `power_curve_results.csv`
- separate performance tables for CCEMG and ABC
- exact ABC conditional covariance and estimator-variance diagnostics, plus the Monte Carlo multiplier-draw variance for comparison
- coverage plots
- average confidence-interval-length plots
- power-curve plots

The obsolete results from the earlier defective CCEMG confidence-interval code are deliberately not included.

## Validation

Both implementations contain checks for:

- correct CCEMG lower and upper endpoints;
- the Equation (44) subtraction sign;
- equality of the ABC estimator and the mean corrected individual slope;
- the exact `l_N=1` covariance identity;
- use of the exact covariance for the reported ABC standard error;
- ordered, finite confidence-interval endpoints.
