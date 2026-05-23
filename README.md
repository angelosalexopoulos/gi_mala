# Gaussian Invariant MCMC (GI-MALA) — Reproducibility Guide

This repository contains the code to reproduce all computational results in:

> **"Gaussian Invariant Markov Chain Monte Carlo"**  
> Submitted to *Journal of the American Statistical Association*

---

## Software Requirements

### MATLAB
- **MATLAB R2020a or later** (tested on R2023b)
- Required toolboxes: Statistics and Machine Learning Toolbox
- No additional third-party MATLAB packages required; all utility functions are in `aGrad/code/toolbox/`

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
├── run_all_experiments.m        ← Step 1: runs all MCMC samplers for paper Tables 3–5, 9–11, 15
├── run_VR_experiments.m         ← Step 1 (VR): runs GI-MALA-VR for Tables 7, 8, 12 (heavy, ~days)
│
├── make_results_binaryclassification.m  ← Step 2: Tables 3,4,9,10,11; Figures 4,5
├── make_results_cox.m                   ← Step 2: Table 5; Figure 6
├── make_results_regression.m            ← Step 2: Table 15
├── make_results_VR_binaryclassification.m ← Step 2 (VR): Tables 7, 12
├── make_results_VR_girolami.m             ← Step 2 (VR): Table 8
│
├── demRegressInformLikelMarg_fixedhypers_gimala.m       ← GI-MALA demo (GP regression)
├── demRegressInformLikelMarg_fixedhypers_pMALA.m        ← pMALA demo (GP regression)
│
├── demBinaryClassification_Marg_fixedhypers_gimala.m    ← GI-MALA demo (binary classif.)
├── demBinaryClassification_Marg_fixedhypers_gimala_VR.m ← GI-MALA-VR demo (binary classif.)
├── demBinaryClassification_Marg_fixedhypers_pMALA.m     ← pMALA demo (binary classif.)
├── demLogGaussianCoxGirolamiMarg_gimala.m               ← GI-MALA demo (Cox process)
├── demLogGaussianCoxGirolamiMarg_gimala_VR.m            ← GI-MALA-VR demo (Cox process)
├── demLogGaussianCoxGirolamiMarg_pMALA.m                ← pMALA demo (Cox process)
│
├── gpsampAuxMarg_fixedhypers_gimala.m       ← GI-MALA sampler (core algorithm)
├── gpsampAuxMarg_fixedhypers_gimala_Gauss.m ← GI-MALA sampler (Gaussian likelihood variant)
├── gpsampAuxMarg_fixedhypers_pMALA.m        ← pMALA sampler implementation
├── ergMean.m                                ← ergodic mean utility
│
├── tmcmc.R                ← R: VR for tail probability (Figure 1 / Table 13)
├── tmcmc_ess.R            ← R: ESS on univariate Student-t (Table 14 / Figure 3)
├── LogisticToy.R          ← R: logistic-regression toy MCMC → Figure 3
├── make_table2_logistic.R ← R: ESS comparison, Bayesian logistic regression (Table 2)
├── make_table6_logistic.R ← R: VR factors, Bayesian logistic regression (Table 6)
├── make_figure4_logistic.R ← R: Figure 4 (GP binary classif. trace + boxplot)
├── make_figure5_regression.R ← R: Figure 5 (GP regression trace + boxplot)
├── make_figure6_cox.R      ← R: Figure 6 (Cox process trace + boxplot)
│
├── RMHMC/
│   ├── demLogGaussianCoxGirolamiRHMC.m  ← RHMC baseline demo (Cox process)
│   └── LGC_RMHMC_LV.m                  ← RHMC sampler for Log-Gaussian Cox
│
├── aGrad/
│   ├── code/
│   │   ├── demBinaryClassification_Marg_fixedhypers.m    ← mGrad baseline (binary classif.)
│   │   ├── demBinaryClassification_Ellipt_fixedhypers.m  ← Elliptical baseline (binary classif.)
│   │   ├── demBinaryClassification_pCN_fixedhypers.m     ← pCN baseline (binary classif.)
│   │   ├── demBinaryClassification_pCNL_fixedhypers.m    ← pCNL baseline (binary classif.)
│   │   ├── demLogGaussianCoxGirolamiMarg.m   ← mGrad baseline (Cox process)
│   │   ├── demLogGaussianCoxGirolamiEllipt.m ← Elliptical baseline (Cox process)
│   │   ├── demLogGaussianCoxGirolamipCN.m    ← pCN baseline (Cox process)
│   │   ├── demLogGaussianCoxGirolamipCNL.m   ← pCNL baseline (Cox process)
│   │   ├── toolbox/   ← shared utility functions (kernels, likelihood, ESS, etc.)
│   │   └── ...        ← legacy make_results scripts (not used for paper tables)
│   └── data/          ← benchmark datasets
│
├── results/
│   ├── LogisticRegression_GP/   ← Step 1 output: binary classification .mat files
│   ├── Cox_regression/          ← Step 1 output: Cox process .mat files
│   └── GPregression/            ← Step 1 output: GP regression .mat files
│
└── diagrams/                    ← Step 2 output: figures (.eps, .pdf) and LaTeX snippets (.txt)
```

---

## Full Reproduction Pipeline

```
Step 1 (MATLAB, from repo root)          Step 2 (MATLAB, from repo root)
────────────────────────────────         ────────────────────────────────────
run_all_experiments       ──►  results/**/*.mat
                                         │
                                         ▼
                                make_results_binaryclassification  ──►  diagrams/
                                make_results_cox                   ──►  diagrams/

run_VR_experiments        ──►  results/**/*_VR_*.mat
                                         │
                                         ▼
                                make_results_VR_binaryclassification  ──►  diagrams/
                                make_results_VR_girolami               ──►  diagrams/

R scripts (from repo root in R)  ──►  diagrams/
```

**All scripts are run from the repository root. No `cd` required.**

---

## Step 1 — Run MCMC Simulations

```matlab
run_all_experiments
```

Runs all methods for both experiments and saves results to `results/`. Set `Repeats = 1` at the
top of the script for a quick smoke-test; `Repeats = 10` for full paper results.

> **Runtime:** RHMC in Experiment 2 dominates (~2–6 hours per repeat). All other methods
> per repeat: ~1–3 hours total. For a quick check skip RHMC by commenting out
> `demLogGaussianCoxGirolamiRHMC(r)`.

### Experiment 1 — Binary Classification

Datasets: Australian, German, Heart, Pima, Ripley.

| Demo called | Method | Output folder | Output `.mat` filename |
|---|---|---|---|
| `demBinaryClassification_Marg_fixedhypers_pMALA` | pMALA | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_Marg_fixedhypers_pMALA.mat` |
| `demBinaryClassification_Marg_fixedhypers` | mGrad | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_Marg_fixedhypers.mat` |
| `demBinaryClassification_Ellipt_fixedhypers` | Ellipt | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_Ellipt_fixedhypers.mat` |
| `demBinaryClassification_pCN_fixedhypers` | pCN | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_pCN_fixedhypers.mat` |
| `demBinaryClassification_pCNL_fixedhypers` | pCNL | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_pCNL_fixedhypers.mat` |
| `demBinaryClassification_Marg_fixedhypers_gimala` | **GI-MALA** | `results/LogisticRegression_GP/` | `<Dataset>_repeat<r>_Marg_fixedhypers_gimala.mat` |

### Experiment 2 — Log-Gaussian Cox Process

| Demo called | Method | Output folder | Output `.mat` filename |
|---|---|---|---|
| `demLogGaussianCoxGirolamiMarg` | mGrad | `results/Cox_regression/` | `logGaussianCoxGirolami_Marg_repeat<r>.mat` |
| `demLogGaussianCoxGirolamiMarg_pMALA` | pMALA | `results/Cox_regression/` | `logGaussianCoxGirolami_Marg_repeat<r>_Marg_pMALA.mat` |
| `demLogGaussianCoxGirolamiEllipt` | Ellipt | `results/Cox_regression/` | `logGaussianCoxGirolami_Ellipt_repeat<r>.mat` |
| `demLogGaussianCoxGirolamipCN` | pCN | `results/Cox_regression/` | `logGaussianCoxGirolami_pCN_repeat<r>.mat` |
| `demLogGaussianCoxGirolamipCNL` | pCNL | `results/Cox_regression/` | `logGaussianCoxGirolami_pCNL_repeat<r>.mat` |
| `demLogGaussianCoxGirolamiRHMC` | RHMC | `results/Cox_regression/` | `logGaussianCoxGirolami_RHMC_repeat<r>_RHMC.mat` |
| `demLogGaussianCoxGirolamiMarg_gimala` | **GI-MALA** | `results/Cox_regression/` | `logGaussianCoxGirolami_Marg_repeat<r>_gimala.mat` |

### Experiment 3 — GP Regression with Informative Likelihood

Dataset: `regressinformlik_d1000.mat`, sigma2 = 0.1², d ≈ 1001.

| Demo called | Method | Output folder | Output `.mat` filename |
|---|---|---|---|
| `demRegressInformLikelMarg_fixedhypers_pMALA` | pMALA | `results/GPregression/` | `regression_repeat<r>_pMALA.mat` |
| `demRegressInformLikelMarg_fixedhypers` | mGrad | `results/GPregression/` | `regression_repeat<r>_Marg.mat` |
| `demRegressInformLikelEllipt_fixedhypers` | Ellipt | `results/GPregression/` | `regression_repeat<r>_Ellipt.mat` |
| `demRegressInformLikelpCN_fixedhypers` | pCN | `results/GPregression/` | `regression_repeat<r>_pCN.mat` |
| `demRegressInformLikelpCNL_fixedhypers` | pCNL | `results/GPregression/` | `regression_repeat<r>_pCNL.mat` |
| `demRegressInformLikelMarg_fixedhypers_gimala` | **GI-MALA** | `results/GPregression/` | `regression_repeat<r>_gimala.mat` |

### Variance-reduction experiments — run separately

```matlab
run_VR_experiments
```

Produces `*_VR_<T>.mat` files for T ∈ {1k, 10k, 50k, 200k} (binary) and T ∈ {1k, 10k} (Cox),
over 100 repeats each. Runtime: ~days. For a quick check set `Repeats = 20` and reduce `T_list`.

---

## Step 2 — Generate Figures and LaTeX Tables

All four scripts are run from the **repository root**. Output goes to `diagrams/`.
Any method whose `.mat` files are missing is silently skipped (no error).

### Binary Classification (Tables 3, 4, 9, 10, 11 / Figures 4, 5)

```matlab
make_results_binaryclassification
```

| Output file in `diagrams/` | Paper reference |
|---|---|
| `Table3.txt` | **Table 3** (Section 5.1.2) — Heart |
| `Table4.txt` | **Table 4** (Section 5.1.2) — Australian |
| `Table9.txt` | **Table 9** (Appendix A.1) — German |
| `Table10.txt` | **Table 10** (Appendix A.1) — Pima |
| `Table11.txt` | **Table 11** (Appendix A.1) — Ripley |
| `<Dataset>_<method>_logL_fixedhypers.pdf` | **Figures 4, 5** (Appendix A.3) |

Each table contains rows for: mGrad, pMALA, Ellipt, pCN, pCNL, **GI-MALA**.  
Columns: Time(s), step size δ, ESS (Min, Med, Max), Min ESS/s.

### Log-Gaussian Cox Process (Table 5 / Figure 6)

```matlab
make_results_cox
```

| Output file in `diagrams/` | Paper reference |
|---|---|
| `Table5_N64.txt` | **Table 5** (Section 5.1.3), N=64 |
| `Table5_N32.txt` | **Table 5** (Section 5.1.3), N=32 |
| `GaussianCox_<method>_logL1.pdf` | **Figure 6** (Appendix A.3) |

`table1` (N=64) contains: mGrad, pMALA, Ellipt, pCN, pCNL, RHMC, **GI-MALA**.  
`table2` (N=32) contains: mGrad, Ellipt, pCN, pCNL (pMALA, RHMC and GI-MALA run at N=64 only).

### GP Regression (Table 15)

```matlab
make_results_regression
```

| Output file in `diagrams/` | Paper reference |
|---|---|
| `Table15.txt` | **Table 15** — GP regression, mGrad / pMALA / Ellipt / pCN / pCNL / GI-MALA |

### Variance Reduction (Tables 7, 8, 12)

```matlab
make_results_VR_binaryclassification   % Tables 7, 12
make_results_VR_girolami               % Table 8
```

| Output file in `diagrams/` | Paper reference |
|---|---|
| `Table7.txt` | **Table 7** (Section 5.2) |
| `Table12.txt` | **Table 12** (Appendix A.1) |
| `Table8.txt` | **Table 8** (Section 5.2) |

---

## R-Based Experiments

Self-contained; do not depend on the MATLAB runs. Run from the repository root in R.

```r
source("make_table2_logistic.R")    # Table 2 (Section 5.1.1)
source("make_table6_logistic.R")    # Table 6 (Section 5.2)
source("tmcmc.R")                   # Figure 1 / Table 13 (Section 5.2.1 / Appendix A.2)
source("tmcmc_ess.R")               # Table 14 (Appendix A.2)
source("LogisticToy.R")             # Figure 3 (Section 5.1.1, logistic toy)
source("make_figure4_logistic.R")   # Figure 4 (GP binary classification)
source("make_figure5_regression.R") # Figure 5 (GP regression)
source("make_figure6_cox.R")        # Figure 6 (Log-Gaussian Cox process)
```

> `make_table*.R` and `LogisticToy.R` read data from the path set in `data_dir` at the top of each file. Edit before running.
> `make_figure4/5/6_*.R` read from `results/` (produced by Step 1) and require the `R.matlab` package.

| Script | Output file | Paper reference |
|---|---|---|
| `make_table2_logistic.R` | `table2_logistic.tex` | **Table 2** |
| `make_table6_logistic.R` | `table6_logistic.tex` | **Table 6** |
| `tmcmc.R` | `Figure1.pdf`, `Table13.tex` | **Figure 1 / Table 13** |
| `tmcmc_ess.R` | `Table14.tex` | **Table 14** |
| `LogisticToy.R` | `Figure3.pdf` | **Figure 3** |
| `make_figure4_logistic.R` | `Figure4.pdf` | **Figure 4** |
| `make_figure5_regression.R` | `Figure5.pdf` | **Figure 5** |
| `make_figure6_cox.R` | `Figure6.pdf` | **Figure 6** |

---

## Complete Output Index

| Paper item | Produced by | Output file in `diagrams/` |
|---|---|---|
| **Table 2** | `make_table2_logistic.R` | `table2_logistic.tex` |
| **Table 3** | `make_results_binaryclassification` | `Table3.txt` (Heart) |
| **Table 4** | `make_results_binaryclassification` | `Table4.txt` (Australian) |
| **Table 5** | `make_results_cox` | `Table5_N64.txt`, `Table5_N32.txt` |
| **Table 6** | `make_table6_logistic.R` | `table6_logistic.tex` |
| **Table 7** | `make_results_VR_binaryclassification` | `Table7.txt` |
| **Table 8** | `make_results_VR_girolami` | `Table8.txt` |
| **Table 9** | `make_results_binaryclassification` | `Table9.txt` (German) |
| **Table 10** | `make_results_binaryclassification` | `Table10.txt` (Pima) |
| **Table 11** | `make_results_binaryclassification` | `Table11.txt` (Ripley) |
| **Table 12** | `make_results_VR_binaryclassification` | `Table12.txt` |
| **Table 15** | `make_results_regression` | `Table15.txt` |
| **Table 13** | `tmcmc.R` | `Table13.tex` |
| **Table 14** | `tmcmc_ess.R` | `Table14.tex` (per GI-target subfolder) |
| **Figure 1** | `tmcmc.R` | `Figure1.pdf` |
| **Figure 3** | `LogisticToy.R` | `Figure3.pdf` |
| **Figure 4** | `make_figure4_logistic.R` | `Figure4.pdf` |
| **Figure 5** | `make_figure5_regression.R` | `Figure5.pdf` |
| **Figure 6** | `make_figure6_cox.R` | `Figure6.pdf` |

---

## Notes

- **Random seeds** are fixed at the top of each demo script (`randn('seed',1); rand('seed',1)`) for exact reproducibility.
- **Serial execution**: all experiments use plain `for` loops. `parfor` was removed because MATLAB parallel workers do not inherit the working directory reliably.
- **Robustness**: `make_results_*.m` scripts skip any method whose `.mat` files are absent — partial runs produce valid tables for completed methods.
- **Output directory**: all figures and LaTeX snippets go to `diagrams/` at the repository root, created automatically by `run_all_experiments.m`.
- **aGrad-z / aGrad-u**: these methods are not included in the paper tables and are not run by `run_all_experiments.m`. Their demo scripts remain in `aGrad/code/` for reference only.
