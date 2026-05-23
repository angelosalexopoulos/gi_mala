function demRegressInformLikelpCNL_fixedhypers(rep)
% demRegressInformLikelpCNL_fixedhypers  (pCNL)
%
% pCNL sampler for GP regression with informative likelihood.
% Data: regressinformlik_d1000.mat, sigma2 = 0.1^2 (config i=2).
%
% Output: results/GPregression/regression_repeat<rep>_pCNL.mat
% Run from the repository root.

addpath aGrad/code/toolbox/;
addpath aGrad/data/;

% --- Data (sigma2 = 0.1^2, d ~ 1001) ---
i = 2;
load regressinformlik_d1000.mat;
X      = XX{i};
Y      = YY{i};
sigma2 = sigma2vals(i);

% --- Model ---
options = gpsampOptions('regression');
model   = gpsampCreate(Y, X, options);
model.Likelihood.logtheta   = log(sigma2);
model.constraints.kernHyper = 'fixed';
model.constraints.likHyper  = 'fixed';

model.K      = kernCompute(model.GP, model.X);
[model.U, model.Lambda, ~] = svd(model.K);
model.Lambda = diag(model.Lambda);
model.L      = diag(model.Lambda.^0.5) * (model.U');

% --- MCMC options ---
mcmcoptions.T          = 5000;
mcmcoptions.Burnin     = 10000;
mcmcoptions.StoreEvery = 1;
mcmcoptions.Langevin   = 1;

% --- Run ---
tic;
[model, samples, accRates] = gpsamppCN_fixedhypers(model, mcmcoptions);
elapsedTime = toc;

% --- Save ---
summarypCNL            = summaryStatistics(samples);
summarypCNL.elapsed    = elapsedTime;
summarypCNL.accRates   = accRates;
summarypCNL.delta      = model.delta;
summarypCNL.beta       = model.beta;
summarypCNL.eff_LogL   = mcmc_ess(samples.LogL(mcmcoptions.Burnin+1:end));
LogL = samples.LogL(mcmcoptions.Burnin+1:end);

save(['results/GPregression/regression_repeat' num2str(rep) '_pCNL.mat'], 'summarypCNL', 'LogL');
fprintf('Saved: results/GPregression/regression_repeat%d_pCNL.mat\n', rep);
