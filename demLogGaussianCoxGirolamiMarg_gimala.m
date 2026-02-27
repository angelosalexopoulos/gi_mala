function demLogGaussianCoxGirolamiMarg_gimala(rep)

addpath aGrad/code/toolbox/; 
addpath aGrad/data/;

dataName = 'logGaussianCoxGirolami_Marg';
storeRes = 1;

%for ds=[2 1]
ds=1;
%load the data
% Hyperparameters of model
load TestData64;
N = 64;
s = 1.91;
b = 1/33;
mu = log(126) - s/2;
m = 1/(N^2);
x_r = 0:1/(N - 1):1;
y_r = 0:1/(N - 1):1;
[xs, ys] = meshgrid(x_r, y_r);
Xt = [];
Xt(:, 1) = xs(:);
Xt(:, 2) = ys(:);

if ds== 2
% work in a subsampled field
N  = 64/ds;
s  = 1.91;
b  = 1/33;
mu = log(126) - s/2;
m  = 1/(N^2);
x_r = 0:1/(N - 1):1;
y_r = 0:1/(N - 1):1;
[xs, ys] = meshgrid(x_r, y_r);
Xt = [];
Xt(:, 1) = xs(:);
Xt(:, 2) = ys(:);
YY = reshape(Y, 64, 64);
XX = reshape(X, 64, 64);
YY = YY(1:2:64, 1:2:64) + YY(2:2:64, 1:2:64) + YY(1:2:64, 2:2:64) + YY(2:2:64, 2:2:64);
XX = XX(1:2:64, 1:2:64);
Y = YY(:);
X = XX(:);
end

means = zeros(1, size(X,1));
means_cv = zeros(1, size(X,1));



% model options
options = gpsampOptions('LogGaussianCox');
options.kern = 'sre';
% create the model
model = gpsampCreate(Y, Xt, options);
model.Likelihood.m = m;
model.Likelihood.logtheta = mu;
model.GP.logtheta = [log((32)*b) log(s)];
%model.GP.logtheta = [0 0];
model.constraints.kernHyper = 'fixed';
model.constraints.likHyper = 'fixed';
model.FF = X(:);
model.delta = 0.1/sqrt(size(X,1));

mcmcoptions.T = 5000;
mcmcoptions.Burnin = 5000;
mcmcoptions.StoreEvery = 1;
mcmcoptions.Langevin = 1;
mcmcoptions.opt =0.82;


% precompute the inverse covariacne matrix Q
model.K = kernCompute(model.GP, model.X);
[model.U, model.Lambda, tmp] = svd(model.K);
model.U = real(model.U);
model.Lambda = diag(model.Lambda);



   tic;
   [model samples accRates] = gpsampAuxMarg_fixedhypers_gimala(model, mcmcoptions);
   elapsedTime = toc;

 % compute statistics
   summaryMarg = summaryStatistics(samples);
   ess = summaryMarg.eff_F;
   LogL= summaryMarg.LogL;
   summaryMarg.elapsed = elapsedTime;
   summaryMarg.accRates = accRates;
   summaryMarg.delta = model.delta;
   eff_LogL = mcmc_ess(samples.LogL);
   save(['results/Cox_regression/' dataName '_repeat' num2str(rep) '_gimala.mat'], ...
       'elapsedTime','accRates','eff_LogL','ess','LogL');


end
