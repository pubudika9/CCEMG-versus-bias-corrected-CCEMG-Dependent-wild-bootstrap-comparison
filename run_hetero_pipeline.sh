#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

if ! command -v Rscript >/dev/null 2>&1; then
  printf 'Error: Rscript was not found. Install R, MASS, and ggplot2 before running this pipeline.\n' >&2
  exit 127
fi

SIM_REPS="${SIM_REPS:-1000}"
BOOT_REPS="${BOOT_REPS:-2000}"
SIM_SEED="${SIM_SEED:-123}"
N_VALUES="${N_VALUES:-50,100,200,400}"
T_VALUES="${T_VALUES:-50,100,200,400}"
ALPHA_VALUES="${ALPHA_VALUES:-0,0.5}"
BETA_TRUE="${BETA_TRUE:-1}"
HETEROGENEITY_SD="${HETEROGENEITY_SD:-1}"
EPS_SD="${EPS_SD:-1}"
CONFIDENCE_LEVEL="${CONFIDENCE_LEVEL:-0.95}"
BOOT_INTERVAL="${BOOT_INTERVAL:-asymmetric}"
DWB_FINITE_SAMPLE_CORRECTION="${DWB_FINITE_SAMPLE_CORRECTION:-true}"
OUTPUT_DIR="${OUTPUT_DIR:-outputs}"
RESULTS_FILE="${RESULTS_FILE:-${OUTPUT_DIR}/all_cases_results.csv}"
POWER_RESULTS_FILE="${POWER_RESULTS_FILE:-${OUTPUT_DIR}/power_curve_results.csv}"

export SIM_REPS BOOT_REPS SIM_SEED N_VALUES T_VALUES ALPHA_VALUES
export BETA_TRUE HETEROGENEITY_SD EPS_SD CONFIDENCE_LEVEL BOOT_INTERVAL
export DWB_FINITE_SAMPLE_CORRECTION OUTPUT_DIR RESULTS_FILE POWER_RESULTS_FILE

printf 'Running standard CCEMG versus ABC whole-sum DWB simulations with l_N=1...\n'
Rscript MCsimuforshrinkingheterogeneity.R

printf 'Generating coverage, interval-length, and power-curve plots...\n'
Rscript data_visualisation_shrinkingheterogeneity.R

printf 'Summary results: %s\n' "$RESULTS_FILE"
printf 'Power curves: %s\n' "$POWER_RESULTS_FILE"
printf 'All tables and plots: %s\n' "$OUTPUT_DIR"
