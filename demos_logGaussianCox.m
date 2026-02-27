
%
randn('seed', 1);
rand('seed', 1);

addpath aGrad/code;
addpath RMHMC/;

% how many times to repeat the experiments 
Repeats = 10;
%parfor (r=1:10,10)
for r=1:10
   demLogGaussianCoxGirolamiEllipt(r);
   demLogGaussianCoxGirolamipCN(r); 
   demLogGaussianCoxGirolamipCNL(r); 
   demLogGaussianCoxGirolamiMarg(r);
   demLogGaussianCoxGirolamiMarg_pMALA(r);
   demLogGaussianCoxGirolamiMarg_gimala(r);
   demLogGaussianCoxGirolamiRHMC(r);
end
