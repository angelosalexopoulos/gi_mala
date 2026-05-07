% run_VR_experiments.m
%
% Runs all variance-reduction experiments needed for Tables 6-8 and Table 12.
%
% Calls:
%   demBinaryClassification_Marg_fixedhypers_gimala_VR(rep, T)
%     for T in {1000, 10000, 50000, 200000} and rep = 1..100
%     -> results/LogisticRegression_GP/   (Tables 7 and 12, Section 5.2 / Appendix A.1)
%
%   demLogGaussianCoxGirolamiMarg_gimala_VR(rep, T)
%     for T in {1000, 10000} and rep = 1..100
%     -> results/Cox_regression/          (Table 8, Section 5.2)
%
% NOTE on Table 6 (Bayesian logistic regression VR, Section 5.2):
%   No VR demo exists for the fixed-coefficient logistic regression setting.
%   Table 6 cannot be reproduced from the current code.
%
% RUNTIME WARNING:
%   Each GP binary classification run takes O(d^2) time per iteration.
%   100 reps x T=200000 x 5 datasets is computationally intensive (~days on
%   a single machine). For a quick verification, reduce Repeats to 20 and
%   T_list_binary to {1000} before running the full experiment.
%
% After all runs complete, generate tables with:
%   cd aGrad/code
%   make_results_VR_binaryclassification    % Tables 7, 12
%   make_results_VR_girolami                % Table 8

randn('seed', 42);
rand('seed',  42);

addpath aGrad/code;
addpath aGrad/code/toolbox;
addpath aGrad/data;

% Ensure output directories exist
if ~exist('results/LogisticRegression_GP', 'dir')
    mkdir('results/LogisticRegression_GP');
end
if ~exist('results/Cox_regression', 'dir')
    mkdir('results/Cox_regression');
end

Repeats        = 100;
T_list_binary  = [1000, 10000, 50000, 200000];
T_list_cox     = [1000, 10000];

% -------------------------------------------------------------------------
% Experiment VR-1: GP binary classification (Tables 7 and 12)
% -------------------------------------------------------------------------
fprintf('\n=== VR Experiment 1: GP Binary Classification ===\n');
fprintf('Running %d repeats x %d T values x 5 datasets\n', Repeats, length(T_list_binary));

for T = T_list_binary
    fprintf('\n  T = %d\n', T);
    parfor (r = 1:Repeats, 10)
        demBinaryClassification_Marg_fixedhypers_gimala_VR(r, T);
    end
end

% -------------------------------------------------------------------------
% Experiment VR-2: Log-Gaussian Cox process (Table 8)
% -------------------------------------------------------------------------
fprintf('\n=== VR Experiment 2: Log-Gaussian Cox Process ===\n');
fprintf('Running %d repeats x %d T values\n', Repeats, length(T_list_cox));

for T = T_list_cox
    fprintf('\n  T = %d\n', T);
    for r = 1:Repeats
        demLogGaussianCoxGirolamiMarg_gimala_VR(r, T);
    end
end

fprintf('\n=== All VR experiments complete. ===\n');
fprintf('Generate tables with:\n');
fprintf('  cd aGrad/code\n');
fprintf('  make_results_VR_binaryclassification   %% Tables 7, 12\n');
fprintf('  make_results_VR_girolami               %% Table 8\n');
