function demRegressInformLikelMarg_fixedhypers_gimala(rep)
% demRegressInformLikelMarg_fixedhypers_gimala
%
% GI-MALA sampler for GP regression with informative likelihood.
% Data: regressinformlik_d1000.mat, sigma2 = 0.1^2 (config i=2).
%
% Output: results/GPregression/regression_repeat<rep>_gimala.mat
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
model.delta                 = 1.0;
model.constraints.kernHyper = 'fixed';
model.constraints.likHyper  = 'fixed';

model.K      = kernCompute(model.GP, model.X);
[model.U, model.Lambda, ~] = svd(model.K);
model.Lambda = diag(model.Lambda);
model.sigma2 = sigma2;

% --- MCMC options ---
mcmcoptions.T          = 5000;
mcmcoptions.Burnin     = 10000;
mcmcoptions.StoreEvery = 1;
mcmcoptions.Langevin   = 1;
mcmcoptions.opt        = 0.574;

% --- Run ---
tic;
[model, samples, accRates] = gpsampAuxMarg_fixedhypers_gimala_Gauss(model, mcmcoptions);
elapsedTime = toc;

% --- Save ---
summaryMarg = summaryStatistics(samples);
ess         = summaryMarg.eff_F;
LogL        = summaryMarg.LogL;
eff_LogL    = mcmc_ess(samples.LogL);
delta       = model.delta;

save(['results/GPregression/regression_repeat' num2str(rep) '_gimala.mat'], ...
     'ess', 'eff_LogL', 'elapsedTime', 'delta', 'LogL', 'accRates');
fprintf('Saved: results/GPregression/regression_repeat%d_gimala.mat\n', rep);
