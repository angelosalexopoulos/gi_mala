% - Binary Classification examples
%
randn('seed', 1);
rand('seed', 1);

addpath aGrad/code;

% how many times to repeat the experiments 
Repeats = 10;
parfor (r=1:Repeats,10)
%for r=1:Repeats
   demBinaryClassification_Marg_fixedhypers(r);
   demBinaryClassification_Marg_fixedhypers_gimala(r);
   demBinaryClassification_Marg_fixedhypers_pMALA(r);
   demBinaryClassification_Ellipt_fixedhypers(r);
   demBinaryClassification_pCN_fixedhypers(r); 
   demBinaryClassification_pCNL_fixedhypers(r); 
end
