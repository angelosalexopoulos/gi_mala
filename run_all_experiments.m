% run_all_experiments.m
%
% Master wrapper script: reproduces all core MCMC simulation results in the paper
% "Gaussian Invariant Markov Chain Monte Carlo" (JASA submission).
%
% USAGE:
%   Run this file from the repository root in MATLAB:
%       >> run_all_experiments
%
% OUTPUT:
%   All MCMC samples are saved as .mat files under results/:
%       results/LogisticRegression_GP/   -- binary classification experiments
%       results/Cox_regression/          -- Log-Gaussian Cox process experiments
%
%   After this script completes, generate figures and tables (from repo root):
%       >> make_results_binaryclassification   % Tables 3,4,9,10,11; Figures 4,5
%       >> make_results_cox                    % Table 5; Figure 6
%
%   For the variance-reduction experiments (Tables 7, 8, 12 / Section 5.2):
%       >> run_VR_experiments                      % 100 reps x multiple T values (~days)
%       >> make_results_VR_binaryclassification    % Tables 7, 12
%       >> make_results_VR_girolami                % Table 8
%
%   For the R-based experiments (from repo root in R):
%       source("make_table2_logistic.R")  # ESS comparison, Bayesian logistic regression (Table 2, Section 5.1.1)
%       source("make_table6_logistic.R")  # VR factors, Bayesian logistic regression (Table 6, Section 5.2)
%       source("tmcmc.R")                 # VR for tail probability estimation (Figure 1 / Table 13, Section 5.2.1 / Appendix A.2)
%       source("tmcmc_ess.R")             # ESS on univariate Student-t (Table 14 / Figure 3, Appendix A.2)
%
% RUNTIME: Set Repeats = 1 for a smoke-test (~hours); Repeats = 10 for full results (~days).
%          RHMC in Experiment 2 dominates runtime (several hours per repeat).
%
% REQUIREMENTS:
%   MATLAB R2020a or later.

% -------------------------------------------------------------------------
% Setup paths — all relative to repository root (run from root)
% -------------------------------------------------------------------------
addpath aGrad/code
addpath aGrad/code/toolbox
addpath aGrad/data
addpath RMHMC

% Ensure output directories exist
if ~exist('results', 'dir'),                      mkdir('results'); end
if ~exist('results/LogisticRegression_GP', 'dir'),mkdir('results/LogisticRegression_GP'); end
if ~exist('results/Cox_regression', 'dir'),       mkdir('results/Cox_regression'); end
if ~exist('results/GPregression', 'dir'),         mkdir('results/GPregression'); end
if ~exist('diagrams', 'dir'),                     mkdir('diagrams'); end

% -------------------------------------------------------------------------
% EXPERIMENT 1: Binary Classification (Section 5.1.2 / Appendix A.1)
%
% Runs 6 MCMC methods on 5 benchmark datasets (Australian, German, Heart,
% Pima, Ripley), each repeated 10 times with fixed kernel hyperparameters.
%
% Produces: results/LogisticRegression_GP/<Dataset>_repeat<r>_<Method>.mat
% Used in: Tables 3,4,9,10,11 (ESS comparison); Figures 4,5 (log-likelihood traces)
% -------------------------------------------------------------------------
fprintf('\n=== EXPERIMENT 1: Binary Classification ===\n');
fprintf('Running methods (mGrad, pMALA, Ellipt, pCN, pCNL, GI-MALA)...\n');

randn('seed', 1234);
rand('seed', 1234);
Repeats = 10;
for r = 1:Repeats
%parfor (r=1:Repeats,10)
    fprintf('  Repeat %d / %d\n', r, Repeats);
    demBinaryClassification_Marg_fixedhypers_pMALA(r);      % pMALA      → results/LogisticRegression_GP/
    demBinaryClassification_Marg_fixedhypers_gimala(r);     % GI-MALA → results/LogisticRegression_GP/
    demBinaryClassification_Marg_fixedhypers(r);            % mGrad      → results/LogisticRegression_GP/
    demBinaryClassification_Ellipt_fixedhypers(r);          % Ellipt     → results/LogisticRegression_GP/
    demBinaryClassification_pCN_fixedhypers(r);             % pCN        → results/LogisticRegression_GP/
    demBinaryClassification_pCNL_fixedhypers(r);            % pCNL       → results/LogisticRegression_GP/
end
fprintf('Experiment 1 complete.\n');

% -------------------------------------------------------------------------
% EXPERIMENT 2: Log-Gaussian Cox Process (Section 5.1.3)
%
% Runs all MCMC methods (including RMHMC baseline) on the 64x64 Gaussian
% Cox spatial intensity field from Girolami & Calderhead (2011).
%
% Produces: results/Cox_regression/logGaussianCoxGirolami_Marg_<Method>_repeat<r>.mat
% Used in: Table 5 (ESS comparison), Figure 6 (log-L traces and field images)
% -------------------------------------------------------------------------
fprintf('\n=== EXPERIMENT 2: Log-Gaussian Cox Process ===\n');
fprintf('Running methods (mGrad, pMALA, Ellipt, pCN, pCNL, RHMC, GI-MALA)...\n');

randn('seed', 1234);
rand('seed', 1234);
for r = 1:Repeats
%parfor (r=1:Repeats,10)    
    fprintf('mGrad  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamiMarg(r);            % mGrad      → results/Cox_regression/
    fprintf('pMALA  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamiMarg_pMALA(r);      % pMALA      → results/Cox_regression/
    fprintf('gi-mala  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamiMarg_gimala(r);     % GI-MALA → results/Cox_regression/
    fprintf('Ellipt Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamiEllipt(r);          % Ellipt     → results/Cox_regression/
    fprintf('pCN  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamipCN(r);             % pCN        → results/Cox_regression/
    fprintf('pCNL  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamipCNL(r);            % pCNL       → results/Cox_regression/
    fprintf('RHMC  Repeat %d / %d\n', r, Repeats);
    demLogGaussianCoxGirolamiRHMC(r);            % RHMC       → results/Cox_regression/
end
fprintf('Experiment 2 complete.\n');

% -------------------------------------------------------------------------
% EXPERIMENT 3: GP Regression with Informative Likelihood (Table 15)
%
% Runs 6 MCMC methods on a 1D regression problem with informative Gaussian
% likelihood (sigma2 = 0.1^2, d ~ 1001) from regressinformlik_d1000.mat.
%
% Produces: results/GPregression/regression_repeat<r>_<Method>.mat
% Used in: Table 15 (ESS comparison)
% -------------------------------------------------------------------------
fprintf('\n=== EXPERIMENT 3: GP Regression with Informative Likelihood ===\n');
fprintf('Running methods (mGrad, pMALA, Ellipt, pCN, pCNL, GI-MALA)...\n');

randn('seed', 1234);
rand('seed', 1234);
%for r = 1:Repeats
parfor (r=1:Repeats,10)    
    fprintf('  Repeat %d / %d\n', r, Repeats);
    demRegressInformLikelMarg_fixedhypers_pMALA(r);      % pMALA   → results/GPregression/
    demRegressInformLikelMarg_fixedhypers(r);            % mGrad   → results/GPregression/
    demRegressInformLikelEllipt_fixedhypers(r);          % Ellipt  → results/GPregression/
    demRegressInformLikelpCN_fixedhypers(r);             % pCN     → results/GPregression/
    demRegressInformLikelpCNL_fixedhypers(r);            % pCNL    → results/GPregression/
    demRegressInformLikelMarg_fixedhypers_gimala(r);     % GI-MALA → results/GPregression/
end
fprintf('Experiment 3 complete.\n');


make_results_binaryclassification   %% Tables 3,4,9,10,11; Figures 4,5
make_results_cox                    %% Table 5; Figure 6
make_results_regression             %% Table 15
run_VR_experiments
make_results_VR_binaryclassification   %% Tables 7, 12
make_results_VR_girolami               %% Table 8
