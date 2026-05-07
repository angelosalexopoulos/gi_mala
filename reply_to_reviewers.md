# Reply to Reviewers
## Manuscript: "Gaussian Invariant Markov Chain Monte Carlo"
## Journal of the American Statistical Association

We thank the reviewers for their careful reading and constructive feedback. Below we provide a point-by-point response to all comments. All changes are reflected in the revised manuscript and the accompanying code repository.

---

## Reviewer 3 (Reproducibility Review)

We thank the reviewer for the detailed and helpful reproducibility assessment. We have made substantial changes to the code repository to address each of the points raised.

---

**Comment 1:**
> While the authors appear to provide the necessary functions and scripts needed to reproduce the results, it is not immediately clear which exact display items in the manuscript will be generated. The authors should make more of an exact connection to display items in the README.

**Response:**
We have completely rewritten the `README.md` to provide an explicit, step-by-step mapping between scripts and every figure and table in the paper. Specifically:

- **Section "Step-by-Step Instructions"** now lists, for each step, which exact output files are produced and which paper display item (Table/Figure number) each file corresponds to.
- A summary table is provided for each results-generation script (`make_results_binaryclassification_fixedhypers.m`, `make_results_girolami.m`, `tmcmc.R`, `tmcmc_ess.R`), listing output file patterns alongside their corresponding display items.
- All figures and tables in the paper that include results are covered, including the main body (Tables 2–8, Figures 1–2), and supplementary material (Appendix A.1: Tables 9–12; Appendix A.2: Table 13, Table 14, Figure 3; Appendix A.3: Figures 4–6).
- For the variance reduction tables (Tables 6–8, Section 5.2), we have added `run_VR_experiments.m` and two new post-processing scripts (`make_results_VR_binaryclassification.m`, `make_results_VR_girolami.m`) that run 100 independent repeats for multiple chain lengths and compute the variance ratio ranges reported in the paper. We note that Table 6 (Bayesian logistic regression VR) cannot currently be reproduced from the provided code; a VR demo for the fixed-coefficient logistic regression setting is not included in this submission.

We also note that the code uses **MATLAB** (R2020a or later) and **R** (4.0 or later), not Python; no Python dependencies exist in this submission. We have added explicit version requirements for both languages and for all R packages used.

---

**Comment 2:**
> Part 3 of the ACC form should also provide the same instructions and connections to the display items in the manuscript.

**Response:**
We have updated Part 3 of the ACC form to mirror the step-by-step instructions and script-to-display-item mapping now present in the README. A reader following only the ACC form will be able to identify which script to run and which figure or table it produces, without needing to consult the README separately.

---

**Comment 3:**
> It would be helpful to have a wrapper function like `submit_experiment.sh` so that readers do not need to track individual scripts and parameter changes.

**Response:**
We have added `run_all_experiments.m`, a master MATLAB wrapper script that runs all MCMC experiments end-to-end in the correct order. Specifically, the script:

1. Sets all paths automatically.
2. Ensures output directories exist before writing results.
3. Runs the binary classification experiments (Experiment 1) — both the baseline methods and the new GI-MALA variants introduced in this paper.
4. Runs the Log-Gaussian Cox process experiments (Experiment 2).
5. Prints, upon completion, the exact commands needed to generate all figures and tables from the saved results.

Each section of the wrapper is preceded by a block comment that names the experiment, explains what it computes, identifies the output location, and states which paper display item(s) the outputs correspond to. We also direct users to `tmcmc.R` (Figure 1, Table 13) and `tmcmc_ess.R` (Table 14, Figure 3) for the R-based experiments in Appendix A.2.

---

**Comment 4:**
> It would also be helpful to have all notebooks and individual scripts more thoroughly commented, with each section describing what the code is doing methodologically.

**Response:**
We have substantially improved the inline documentation of the core algorithm file and all new demo scripts:

- **`gpsampAuxMarg_fixedhypers_gimala.m`** (the core GI-MALA sampler): added a full function-level docstring specifying all inputs, outputs, and their connection to the paper equations; added section-level comments at each key algorithmic step (preconditioner computation, proposal parameters, MH correction factor, step-size adaptation).
- **`demBinaryClassification_Marg_fixedhypers_gimala.m`**: added a function-level header describing inputs, outputs (file locations and variable names), and the paper display items the outputs feed into.
- **`demLogGaussianCoxGirolamiMarg_gimala.m`**: same treatment.
- **`tmcmc.R`**: added a file-level header block describing the purpose, target distribution, free parameters, output format, expected runtime, and which paper table it produces.
- **`tmcmc_ess.R`**: already contained a thorough header block; no changes needed.

---

**Comment 5:**
> The code uses Python, which incorporates libraries and other dependencies. These need to be listed in the ACC form with their version numbers.

**Response:**
We note that the code in this submission is implemented entirely in **MATLAB** and **R**; there is no Python code. We have added explicit software requirements to both the README and the ACC form:

**MATLAB:**
- MATLAB R2020a or later (tested on R2023b)
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (optional; used for `parfor` in `demos_binaryclassification.m`)

**R:**
- R 4.0 or later (tested on R 4.3.2)
- `MASS` ≥ 7.3
- `ggplot2` ≥ 3.4.0
- `coda` ≥ 0.19-4
- `ks` ≥ 1.14

All MATLAB utility functions are self-contained in `aGrad/code/toolbox/` and do not require any external MATLAB packages.
