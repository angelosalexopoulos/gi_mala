% run_all_experiments.m
%
% Master wrapper script: reproduces all MCMC simulation results in the paper
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
%   After this script completes, generate figures and tables by running:
%       >> cd aGrad/code
%       >> make_results_binaryclassification_fixedhypers   % Section 5.X tables/figures
%       >> make_results_girolami                            % Section 5.X tables/figures
%       >> make_results_girolami_mMALA_RHMC                % RMHMC comparison figure
%
%   For the R-based experiments (Table 1 and Appendix A.3), run:
%       source("tmcmc.R")       # variance-reduction table (Table 1)
%       source("tmcmc_ess.R")   # ESS comparison (Tables A1-A3, Figures A1-A2)
%
% RUNTIME: approx. 2-4 hours on a modern multi-core machine (10 parallel workers).
%          Set Repeats = 1 for a quick smoke-test.
%
% REQUIREMENTS:
%   MATLAB R2020a or later; Parallel Computing Toolbox (optional, for parfor).

% -------------------------------------------------------------------------
% Setup paths
% -------------------------------------------------------------------------
addpath aGrad/code
addpath aGrad/code/toolbox
addpath aGrad/data
addpath RMHMC

% Ensure output directories exist
if ~exist('results/LogisticRegression_GP', 'dir')
    mkdir('results/LogisticRegression_GP');
end
if ~exist('results/Cox_regression', 'dir')
    mkdir('results/Cox_regression');
end
if ~exist('aGrad/diagrams', 'dir')
    mkdir('aGrad/diagrams');
end

% -------------------------------------------------------------------------
% EXPERIMENT 1: Binary Classification (Section 5.X)
%
% Runs 6 MCMC methods on 5 benchmark datasets (Australian, German, Heart,
% Pima, Ripley), each repeated 10 times with fixed kernel hyperparameters.
%
% Produces: results/LogisticRegression_GP/<Dataset>_repeat<r>_<Method>.mat
% Used in: Table X (ESS comparison) and Figure X (log-likelihood traces)
% -------------------------------------------------------------------------
fprintf('\n=== EXPERIMENT 1: Binary Classification ===\n');
fprintf('Running baseline methods (aGrad-z, aGrad-u, mGrad, Ellipt, pCN, pCNL)...\n');
demos_binaryclassification          % baseline methods from aGrad/code/

fprintf('Running GI-MALA variants introduced in this paper...\n');

randn('seed', 1);
rand('seed', 1);
Repeats = 10;
parfor (r = 1:Repeats, 10)
    % GI-MALA (proposed method, Section 3)
    demBinaryClassification_Marg_fixedhypers_gimala(r);
    % GI-MALA with variance reduction (Section 4)
    demBinaryClassification_Marg_fixedhypers_gimala_VR(r);
    % Preconditioned MALA baseline (for comparison)
    demBinaryClassification_Marg_fixedhypers_pMALA(r);
end
fprintf('Experiment 1 complete.\n');

% -------------------------------------------------------------------------
% EXPERIMENT 2: Log-Gaussian Cox Process (Section 5.X)
%
% Runs all MCMC methods (including RMHMC baseline) on the 64x64 Gaussisan
% Cox spatial intensity field from Girolami & Calderhead (2011).
%
% Produces: results/Cox_regression/logGaussianCoxGirolami_Marg_<Method>_repeat<r>.mat
% Used in: Table X (ESS comparison), Figure X (log-L traces), Figure X (field images)
% -------------------------------------------------------------------------
fprintf('\n=== EXPERIMENT 2: Log-Gaussian Cox Process ===\n');
fprintf('Running baseline methods (Ellipt, pCN, pCNL, mGrad, RHMC)...\n');
demos_logGaussianCox                % runs all methods including RMHMC

fprintf('Running GI-MALA variants introduced in this paper...\n');

randn('seed', 1);
rand('seed', 1);
for r = 1:10
    % GI-MALA (proposed method)
    demLogGaussianCoxGirolamiMarg_gimala(r);
    % GI-MALA with variance reduction
    demLogGaussianCoxGirolamiMarg_gimala_VR(r);
    % Preconditioned MALA baseline
    demLogGaussianCoxGirolamiMarg_pMALA(r);
end
fprintf('Experiment 2 complete.\n');

% -------------------------------------------------------------------------
% DONE — remind the user to generate figures next
% -------------------------------------------------------------------------
fprintf('\n=== All MCMC experiments complete ===\n');
fprintf('To generate figures and tables, run:\n');
fprintf('  cd aGrad/code\n');
fprintf('  make_results_binaryclassification_fixedhypers\n');
fprintf('  make_results_girolami\n');
fprintf('  make_results_girolami_mMALA_RHMC\n');
fprintf('For the R-based experiments (Table 1 / Appendix A.3):\n');
fprintf('  source("tmcmc.R")      %% Table 1\n');
fprintf('  source("tmcmc_ess.R")  %% Tables A1-A3, Figures A1-A2\n');
