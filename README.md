# Gaussian Invariant MCMC (GI-MALA) — Reproducibility Guide

This repository contains the code to reproduce all computational results in:

> **"Gaussian Invariant Markov Chain Monte Carlo"**  
> Submitted to *Journal of the American Statistical Association*

---

## Software Requirements

### MATLAB
- **MATLAB R2020a or later** (tested on R2023b)
- Required toolboxes: Statistics and Machine Learning Toolbox, Parallel Computing Toolbox (for `parfor` acceleration; optional — replace `parfor` with `for` if unavailable)
- No additional third-party MATLAB packages are required; all utility functions are included in `aGrad/code/toolbox/`

### R
- **R 4.0 or later** (tested on R 4.3.2)
- Required packages (install via `install.packages(...)` before running):
  - `stats` (base, no install needed)
  - `MASS` (version ≥ 7.3)
  - `ggplot2` (version ≥ 3.4.0)
  - `coda` (version ≥ 0.19-4) — used for ESS diagnostics
  - `ks` (version ≥ 1.14) — used for kernel density estimation in `tmcmc_ess.R`

---

## Repository Structure

```
.
├── README.md                    ← this file
├── run_all_experiments.m        ← master wrapper (runs everything in order)
├── run_VR_experiments.m         ← VR wrapper (100 reps × multiple T values for Tables 7, 8, 12)
│
├── demos_binaryclassification.m ← Step 1a: run binary classification simulations
├── demos_logGaussianCox.m       ← Step 1b: run Log-Gaussian Cox process simulations
│
├── demBinaryClassification_Marg_fixedhypers_gimala.m    ← GI-MALA sampler demo (binary classif.)
├── demBinaryClassification_Marg_fixedhypers_gimala_VR.m ← GI-MALA + variance reduction demo
├── demBinaryClassification_Marg_fixedhypers_pMALA.m     ← preconditioned MALA demo (binary classif.)
├── demLogGaussianCoxGirolamiMarg_gimala.m               ← GI-MALA demo (Cox process)
├── demLogGaussianCoxGirolamiMarg_gimala_VR.m            ← GI-MALA + VR demo (Cox process)
├── demLogGaussianCoxGirolamiMarg_pMALA.m                ← pMALA demo (Cox process)
│
├── gpsampAuxMarg_fixedhypers_gimala.m      ← GI-MALA sampler implementation (core algorithm)
├── gpsampAuxMarg_fixedhypers_gimala_Gauss.m← GI-MALA sampler (Gaussian likelihood variant)
├── gpsampAuxMarg_fixedhypers_pMALA.m       ← preconditioned MALA sampler implementation
├── ergMean.m                               ← ergodic mean utility
│
├── tmcmc.R       ← Step 1c: variance reduction for tail probability (Section 5.2.1 / Figure 1; Appendix A.2 / Table 13)
├── tmcmc_ess.R   ← Step 1d: ESS comparison on univariate Student-t (Appendix A.2 / Table 14, Figure 3)
│
├── make_table2_logistic.R ← Step 1f: ESS comparison, Bayesian logistic regression (Table 2, Section 5.1.1)
├── make_table6_logistic.R ← Step 1f: VR factors, Bayesian logistic regression (Table 6, Section 5.2)
├── LogisticToy.R          ← exploratory logistic regression script (basis for the above)
│
├── RMHMC/
│   ├── demLogGaussianCoxGirolamiRHMC.m  ← RMHMC baseline demo (Cox process)
│   └── LGC_RMHMC_LV.m                  ← RMHMC sampler for Log-Gaussian Cox
│
├── aGrad/
│   └── code/
│       ├── demos_binaryclassification_fixedhypers.m  ← runner for fixed-hyperparameter baselines
│       ├── demos_binaryclassification_learnhypers.m  ← runner for learned-hyperparameter baselines
│       ├── demos_girolami.m                          ← runner for Cox process baselines
│       ├── demos_informlik.m                         ← runner for informative-likelihood experiment
│       ├── make_results_binaryclassification_fixedhypers.m ← Step 2a: generate figures/tables
│       ├── make_results_girolami.m                         ← Step 2b: generate figures/tables
│       ├── make_results_girolami_mMALA_RHMC.m              ← Step 2c: generate RHMC comparison
│       ├── make_results_VR_binaryclassification.m          ← Step 2d: Tables 7, 12 (VR factors)
│       ├── make_results_VR_girolami.m                      ← Step 2d: Table 8 (VR factors)
│       ├── toolbox/                    ← shared utility functions (kernels, likelihood, ESS, etc.)
│       └── data/                       ← benchmark datasets
│
└── results/
    ├── LogisticRegression_GP/   ← MCMC output files saved here (binary classification)
    └── Cox_regression/          ← MCMC output files saved here (Cox process)
```

---

## Step-by-Step Instructions to Reproduce All Results

> **Estimated total runtime:** ~2–4 hours on a modern multi-core machine with 10 parallel workers.  
> For a quick test, set `Repeats = 1` and `mcmcoptions.T = 100` in any script.

### Step 1 — Run MCMC simulations

All simulations can be run in one go by calling the master wrapper from MATLAB:

```matlab
run_all_experiments
```

Alternatively, run each experiment individually as described below.

---

#### Step 1a — Binary Classification (Sections 5.1.1 and 5.1.2)

From the repository root in MATLAB:

```matlab
addpath aGrad/code
addpath aGrad/code/toolbox
addpath aGrad/data

% Run all baseline methods (aGrad-z, aGrad-u, mGrad, Ellipt, pCN, pCNL)
demos_binaryclassification        % calls aGrad/code/demos_binaryclassification_fixedhypers.m

% Run the new GI-MALA variants introduced in this paper
demos_binaryclassification        % also calls demBinaryClassification_Marg_fixedhypers_gimala
                                  %              demBinaryClassification_Marg_fixedhypers_gimala_VR
                                  %              demBinaryClassification_Marg_fixedhypers_pMALA
```

Results are saved in `results/LogisticRegression_GP/` as `.mat` files, one per method × dataset × repeat.

---

#### Step 1b — Log-Gaussian Cox Process (Section 5.1.3)

```matlab
addpath aGrad/code
addpath aGrad/code/toolbox
addpath aGrad/data
addpath RMHMC

demos_logGaussianCox   % runs all Cox process methods including GI-MALA, pMALA, and RMHMC
```

Results are saved in `results/Cox_regression/` as `.mat` files.

---

#### Step 1e — Variance Reduction Experiments for Tables 7, 8, 12 (Section 5.2)

From the repository root in MATLAB:

```matlab
run_VR_experiments
```

This script runs `demBinaryClassification_Marg_fixedhypers_gimala_VR` for 100 repeats × 4 chain lengths (T ∈ {1,000; 10,000; 50,000; 200,000}) and `demLogGaussianCoxGirolamiMarg_gimala_VR` for 100 repeats × 2 chain lengths (T ∈ {1,000; 10,000}).

> **Runtime warning**: the full run is computationally intensive (100 reps × T=200,000 for 5 datasets). For a quick smoke-test, set `Repeats = 20` and `T_list_binary = [1000]` at the top of `run_VR_experiments.m` before running.

Results are saved as `<Dataset>_repeat<r>_Marg_fixedhypers_VR_<T>.mat` in the appropriate `results/` subdirectory.

> **Note on Table 6** (Bayesian logistic regression VR, Section 5.2): a VR demo for the fixed-coefficient logistic regression setting does not currently exist in this repository. Table 6 cannot be reproduced from the present code.

---

#### Step 1f — Bayesian Logistic Regression (Section 5.1.1 / **Table 2** and Section 5.2 / **Table 6**)

Open R (version ≥ 4.0) and run:

```r
source("make_table2_logistic.R")   # ESS comparison (Table 2, Section 5.1.1)
source("make_table6_logistic.R")   # Variance reduction factors (Table 6, Section 5.2)
```

Both scripts read data from the path set in `data_dir` at the top of each file (default: `/Users/angelos/Desktop/vr_pCN/data`). Edit this path if your data lives elsewhere. Required packages: `MASS`, `coda`.

| Script | Output file | Corresponds to |
|---|---|---|
| `make_table2_logistic.R` | `table2_logistic.tex` | **Table 2** (Section 5.1.1): ESS for GI-MALA vs. pMALA, Heart/Australian/German (10 repeats) |
| `make_table6_logistic.R` | `table6_logistic.tex` | **Table 6** (Section 5.2): VR factors, Heart/Australian, n ∈ {1k, 10k, 50k, 200k} (100 repeats) |

> **Runtime warning for `make_table6_logistic.R`**: 100 reps × T=200,000 is expensive. For a quick check, set `Repeats <- 20` and `T_list <- c(1000)` at the top of the script.

---

#### Step 1c — Variance Reduction for Tail Probability Estimation (Section 5.2.1 / **Figure 1**)

Open R (version ≥ 4.0) and run:

```r
source("tmcmc.R")
```

This script runs GI-MALA on a univariate Student-t target for combinations of degrees-of-freedom `nu` and threshold values `c`, and computes variance-reduction ratios. The output corresponds to **Figure 1** (Section 5.2.1) and **Table 13** (Appendix A.2, tabular version of Figure 1 for all b values).

---

#### Step 1d — ESS Comparison on Univariate Student-t (Appendix A.2 / **Table 14**, **Figure 3**)

```r
source("tmcmc_ess.R")
```

This script compares GI-MALA, MALA, and RWM on a univariate Student-t for varying degrees of freedom. Outputs are written to the folder specified by `base_out_dir` (edit this path at the top of the script before running):

| Output file | Corresponds to |
|---|---|
| `table_ess_t_selftuned.tex` | **Table 14** (Appendix A.2): ESS, acceptance rate, step size γ for GI-MALA, MALA, RWM |
| `density_vs_true_*.pdf` | **Figure 3** (Appendix A.2): estimated stationary densities vs. true Student-t |

---

### Step 2 — Generate Figures and Tables

After all MCMC runs have completed (Step 1a–1b), generate the paper figures and tables from MATLAB:

#### Step 2a — Binary Classification Figures and Tables

```matlab
cd aGrad/code
make_results_binaryclassification_fixedhypers
```

This script reads the `.mat` files from `results/LogisticRegression_GP/` and writes to `aGrad/diagrams/`:

| Output file pattern | Corresponds to |
|---|---|
| `Heart_table_fixedhypers.txt` | **Table 3** (Section 5.1.2): ESS comparison, Heart dataset (d=270) |
| `Australian_table_fixedhypers.txt` | **Table 4** (Section 5.1.2): ESS comparison, Australian dataset (d=690) |
| `German_table_fixedhypers.txt` | **Table 9** (Appendix A.1): ESS comparison, German dataset (d=1000) |
| `Pima_table_fixedhypers.txt` | **Table 10** (Appendix A.1): ESS comparison, Pima dataset (d=532) |
| `Ripley_table_fixedhypers.txt` | **Table 11** (Appendix A.1): ESS comparison, Ripley dataset (d=250) |
| `<Dataset>_*_logL_fixedhypers.{eps,pdf}` | **Figures 4, 5** (Appendix A.3): log-likelihood trace plots (Heart, Australian) |

---

#### Step 2b — Log-Gaussian Cox Process Figures and Tables

```matlab
cd aGrad/code
make_results_girolami
```

Writes to `aGrad/diagrams/`:

| Output file pattern | Corresponds to |
|---|---|
| `GaussianCox_table*.txt` | **Table 5** (Section 5.1.3): ESS comparison, Cox process (d=4096) |
| `GaussianCox_*_logL*.{eps,pdf}` | **Figure 6** (Appendix A.3): log-likelihood trace plots and boxplots |
| `GaussianCox_*_meanfield.{eps,pdf}` | posterior mean of the intensity field (Section 5.1.3) |
| `GaussianCox_*_varfield.{eps,pdf}` | posterior variance of the intensity field (Section 5.1.3) |
| `GaussianCox_*_expfield.{eps,pdf}` | posterior mean of exp(F) — inferred intensity (Section 5.1.3) |

---

#### Step 2d — Variance Reduction Tables (Section 5.2 / **Tables 7, 8, 12**)

After `run_VR_experiments` has completed (Step 1e), from `aGrad/code`:

```matlab
make_results_VR_binaryclassification
make_results_VR_girolami
```

Writes to `aGrad/diagrams/`:

| Output file | Corresponds to |
|---|---|
| `VR_GP_classification_table7.txt` | **Table 7** (Section 5.2): VR factors, GP binary classification (Heart d=270, Australian d=690) |
| `VR_GP_classification_table12.txt` | **Table 12** (Appendix A.1): VR factors, remaining GP classification datasets (German d=1000, Pima d=532, Ripley d=250) |
| `VR_GaussianCox_table8.txt` | **Table 8** (Section 5.2): VR factors, log-Gaussian Cox process (d=4096) |

---

#### Step 2c — RHMC Comparison (included in **Table 5**)

```matlab
cd aGrad/code
make_results_girolami_mMALA_RHMC
```

Processes the RMHMC results and adds them to the Cox process comparison (Table 5, Section 5.1.3).

---

## Notes

- **Random seeds** are fixed at the top of each demo script (`randn('seed',1); rand('seed',1)`) to ensure exact reproducibility.
- **Parallel execution**: `demos_binaryclassification.m` uses `parfor` with 10 workers. On machines with fewer cores, reduce the worker count or switch to a `for` loop — results will be identical but slower.
- **Output directory**: All figures are written to `aGrad/diagrams/`. Ensure this directory exists before running `make_results_*.m` scripts.
