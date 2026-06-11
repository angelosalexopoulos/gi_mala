%patch_logL.m
%
% Patches existing .mat result files that were saved BEFORE the MATLAB demo
% scripts were updated to write a flat top-level 'LogL' variable.
%
% For each file the script:
%   1. Checks whether 'LogL' already exists at the top level (skip if yes).
%   2. Loads the summary struct, handles cell-array wrappers (Cox demos).
%   3. Extracts the post-burnin LogL: real(s.LogL(burnin+1 : end)).
%   4. Appends 'LogL' to the existing file with  save(..., '-append').
%
% Run once from the project root (gi-mala-new/) in MATLAB before executing
% the R figure scripts.
%
% Usage:
%   cd /path/to/gi-mala-new
%   patch_logL

fprintf('\n=== patch_logL: patching existing result files ===\n\n');

np_total = 0;   % patched
ns_total = 0;   % skipped (already has LogL)
nm_total = 0;   % missing file
nw_total = 0;   % warning / error

% =========================================================
% 1.  BINARY CLASSIFICATION  (LogisticRegression_GP/)
% =========================================================
datasets = {'Australian','German','Heart','Pima','Ripley'};
Reps     = 10;
burnin   = 5000;

logistic_methods = {
    '_pCN_fixedhypers',   'summarypCN';
    '_pCNL_fixedhypers',  'summarypCNL';
    '_Ellipt_fixedhypers','summaryEllipt';
    '_Marg_fixedhypers',  'summaryMarg';
};

fprintf('--- Binary Classification (LogisticRegression_GP/) ---\n');
for di = 1:numel(datasets)
    ds = datasets{di};
    for r = 1:Reps
        for mi = 1:size(logistic_methods,1)
            suffix = logistic_methods{mi,1};
            sname  = logistic_methods{mi,2};
            fname  = sprintf('results/LogisticRegression_GP/%s_repeat%d%s.mat', ...
                             ds, r, suffix);
            [np_total, ns_total, nm_total, nw_total] = ...
                patch_one(fname, sname, burnin, false, ...
                          np_total, ns_total, nm_total, nw_total);
        end
    end
end

% =========================================================
% 2.  GP REGRESSION  (GPregression/)
% =========================================================
burnin_reg = 10000;

regression_methods = {
    '_Marg',   'summaryMarg';
    '_Ellipt', 'summaryEllipt';
    '_pCN',    'summarypCN';
    '_pCNL',   'summarypCNL';
};

fprintf('--- GP Regression (GPregression/) ---\n');
for r = 1:Reps
    for mi = 1:size(regression_methods,1)
        suffix = regression_methods{mi,1};
        sname  = regression_methods{mi,2};
        fname  = sprintf('results/GPregression/regression_repeat%d%s.mat', r, suffix);
        [np_total, ns_total, nm_total, nw_total] = ...
            patch_one(fname, sname, burnin_reg, false, ...
                      np_total, ns_total, nm_total, nw_total);
    end
end

% =========================================================
% 3.  COX PROCESS  (Cox_regression/)
%     The summary structs are wrapped in a cell array {ds}.
% =========================================================
burnin_cox = 5000;
dataName_cox = 'logGaussianCoxGirolami';

cox_methods = {
    '_Marg',   'summaryMarg';    % mGrad
    '_Ellipt', 'summaryEllipt';
    '_pCNL',   'summarypCNL';
    '_pCN',    'summarypCN';
};

fprintf('--- Cox Process (Cox_regression/) ---\n');
for r = 1:Reps
    for mi = 1:size(cox_methods,1)
        infix  = cox_methods{mi,1};
        sname  = cox_methods{mi,2};
        % File naming: logGaussianCoxGirolami_{method}_repeat{r}.mat
        % (mGrad is stored without method tag in the file name)
        if strcmp(infix, '_Marg')
            fname = sprintf('results/Cox_regression/%s_Marg_repeat%d.mat', ...
                            dataName_cox, r);
        elseif strcmp(infix, '_Ellipt')
            fname = sprintf('results/Cox_regression/%s_Ellipt_repeat%d.mat', ...
                            dataName_cox, r);
        elseif strcmp(infix, '_pCNL')
            fname = sprintf('results/Cox_regression/%s_pCNL_repeat%d.mat', ...
                            dataName_cox, r);
        elseif strcmp(infix, '_pCN')
            fname = sprintf('results/Cox_regression/%s_pCN_repeat%d.mat', ...
                            dataName_cox, r);
        end
        % is_cell = true because Cox demos use summaryX{ds} cell arrays
        [np_total, ns_total, nm_total, nw_total] = ...
            patch_one(fname, sname, burnin_cox, true, ...
                      np_total, ns_total, nm_total, nw_total);
    end
end

% =========================================================
% Summary
% =========================================================
fprintf('\n=== Done ===\n');
fprintf('  Patched  : %d files\n', np_total);
fprintf('  Skipped  : %d files (LogL already present)\n', ns_total);
fprintf('  Missing  : %d files (not found)\n', nm_total);
fprintf('  Warnings : %d files (error during patching)\n', nw_total);

% =========================================================
% Helper function
% =========================================================
function [np, ns, nm, nw] = patch_one(fpath, sname, burnin, is_cell, np, ns, nm, nw)
%PATCH_ONE  Append flat post-burnin LogL to one .mat file.
%
%  fpath    - path to the .mat file
%  sname    - name of the summary struct inside the file
%  burnin   - number of burn-in samples to discard
%  is_cell  - true if the struct is wrapped in a cell array (Cox demos)

    if ~exist(fpath, 'file')
        fprintf('  [MISSING]  %s\n', fpath);
        nm = nm + 1;
        return
    end

    % Check if LogL already exists at the top level
    try
        w = whos('-file', fpath);
    catch ME
        fprintf('  [ERROR whos] %s : %s\n', fpath, ME.message);
        nw = nw + 1;
        return
    end

    names = {w.name};
    if ismember('LogL', names)
        % Already patched — skip silently
        ns = ns + 1;
        return
    end

    % Load the summary struct
    try
        data = load(fpath, sname);
    catch ME
        fprintf('  [ERROR load] %s : %s\n', fpath, ME.message);
        nw = nw + 1;
        return
    end

    if ~isfield(data, sname)
        fprintf('  [ERROR]  struct ''%s'' not found in %s\n', sname, fpath);
        nw = nw + 1;
        return
    end

    s = data.(sname);

    % Unwrap cell array (Cox process demos store summaryX{ds})
    if is_cell || iscell(s)
        s = s{1};
    end

    if ~isfield(s, 'LogL')
        fprintf('  [ERROR]  LogL field missing in struct %s in %s\n', sname, fpath);
        nw = nw + 1;
        return
    end

    % Extract post-burnin samples, strip imaginary artefacts
    try
        raw  = real(s.LogL(:));
        nraw = numel(raw);
        if nraw <= burnin
            fprintf('  [WARN]  %s: only %d samples (burnin=%d), using all\n', ...
                    fpath, nraw, burnin);
            LogL = raw;
        else
            LogL = raw(burnin+1:end);
        end
    catch ME
        fprintf('  [ERROR extract] %s : %s\n', fpath, ME.message);
        nw = nw + 1;
        return
    end

    % Append to file
    try
        save(fpath, 'LogL', '-append');
        fprintf('  [PATCHED]  %s  (%d samples)\n', fpath, numel(LogL));
        np = np + 1;
    catch ME
        fprintf('  [ERROR save] %s : %s\n', fpath, ME.message);
        nw = nw + 1;
    end
end
