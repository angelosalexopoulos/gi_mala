function demLogGaussianCoxGirolamiRHMC(rep)

addpath aGrad/code/toolbox/; 

addpath RMHMC/;
dataName = 'logGaussianCoxGirolami_RHMC';


[samples, elapsedTime] = LGC_RMHMC_LV();

samples.F = real(samples.F);
 % compute statistics 
   summaryMarg = summaryStatistics(samples);
   ess = summaryMarg.eff_F;
   LogL= summaryMarg.LogL;
   summaryMarg.elapsed = elapsedTime;
   eff_LogL = mcmc_ess(real(samples.LogL));
   save(['results/Cox_regression/' dataName '_repeat' num2str(rep) '_RHMC.mat'], ...
       'elapsedTime','eff_LogL','ess','LogL');
