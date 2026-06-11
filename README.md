# Gaussian Invariant MCMC (GI-MALA) — Reproducibility Guide

This repository contains the code to reproduce all computational results in:

> **"Gaussian Invariant Markov Chain Monte Carlo"**  

## Software Requirements

### MATLAB
- **MATLAB R2020a or later** (tested on R2023b)
- Statistics and Machine Learning Toolbox
- No third-party packages; all utilities are in `aGrad/code/toolbox/`

### R
- **R 4.0 or later** (tested on R 4.3.2)
- Required packages: `MASS` (≥ 7.3), `ggplot2` (≥ 3.4.0), `coda` (≥ 0.19-4), `ks` (≥ 1.14)

---

## Quick Start — Reproduce Everything

The full pipeline requires two steps: run the MATLAB simulations, then run the R scripts.

### Step 1 — MATLAB (from the repository root)

```matlab
run_all_experiments
```

This single script runs all three MCMC experiments (binary classification, Cox process, GP regression) over 10 repeats, saves `.mat` results to `results/`, and produces the following tables in `diagrams/`:

| Tables produced | Description |
|---|---|
| `Table3.txt`, `Table4.txt`, `Table9.txt`, `Table10.txt`, `Table11.txt` | ESS comparison — binary classification (Heart, Australian, German, Pima, Ripley) |
| `Table5.txt` | ESS comparison — Log-Gaussian Cox process |
| `Table15.txt` | ESS comparison — GP regression |
| `Table7.txt`, `Table12.txt` | Variance reduction — binary classification |
| `Table8.txt` | Variance reduction — Cox process |


### Step 2 — R (from the repository root in R)

All figures and R-based tables are produced by self-contained R scripts. Run in any order:

```r
source("make_table2_logistic.R")    # Table 2  — ESS comparison, Bayesian logistic regression
source("make_table6_logistic.R")    # Table 6  — VR factors, Bayesian logistic regression
source("tmcmc.R")                   # Figure 1 + Table 13 — VR for tail probability
source("tmcmc_ess.R")               # Table 14 + Figure 3 — ESS on univariate Student-t
```

> Edit `data_dir` at the top of `make_table2_logistic.R` and `make_table6_logistic.R` to point to the folder containing the logistic-regression datasets before running.


---

## Complete Output Index

| Paper item | Produced by | Output file in `diagrams/` |
|---|---|---|
| **Table 2** | `make_table2_logistic.R` | `table2_logistic.tex` |
| **Table 3** | `make_results_binaryclassification` | `Table3.txt` (Heart) |
| **Table 4** | `make_results_binaryclassification` | `Table4.txt` (Australian) |
| **Table 5** | `make_results_cox` | `Table5_N64.txt`, `Table5_N32.txt` |
| **Table 6** | `make_table6_logistic.R` | `Table6.tex` |
| **Table 7** | `make_results_VR_binaryclassification` | `Table7.txt` |
| **Table 8** | `make_results_VR_girolami` | `Table8.txt` |
| **Table 9** | `make_results_binaryclassification` | `Table9.txt` (German) |
| **Table 10** | `make_results_binaryclassification` | `Table10.txt` (Pima) |
| **Table 11** | `make_results_binaryclassification` | `Table11.txt` (Ripley) |
| **Table 12** | `make_results_VR_binaryclassification` | `Table12.txt` |
| **Table 13** | `tmcmc.R` | `Table13.tex` |
| **Table 14** | `tmcmc_ess.R` | `Table14.tex` |
| **Table 15** | `make_results_regression` | `Table15.txt` |
| **Figure 1** | `tmcmc.R` | `Figure1.pdf` |
| **Figure 3** | `tmcmc_ess.R` | `Figure3.pdf` |

---

## Notes

- **Random seeds** are fixed at the top of `run_all_experiments.m` (`randn('seed',1234); rand('seed',1234)`) for exact reproducibility.
- **Parallel execution**: `run_all_experiments.m` uses `parfor` with up to 10 workers. Replace with `for` if running without the Parallel Computing Toolbox.
- **Robustness**: all `make_results_*.m` scripts skip any method whose `.mat` files are absent — no errors on partial runs.
- **Output directory**: `diagrams/` is created automatically if it does not exist.

