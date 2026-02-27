
This repository contains MATLAB demo scripts that run repeated MCMC experiments across:
- **models / experiments** (e.g., binary classification, regression with informative likelihood, etc.)
- **sampling methods** (e.g., preconditioned MALA, pCN variants, elliptical samplers)
- optionally **variance reduction (VR)** when indicated in the script name.

The typical workflow is:
1. Set the random seeds (for reproducibility).
2. Choose the number of repeats.
3. Run each method for each repeat `r` (often in parallel via `parfor`).

## 1) Binary classification experiment calls\

Example call list that executed within the file demos_binaryclassifcation.m file

Repeats = 10;
parfor \cf0 (r=1:Repeats,10)
 %for r=1:Repeats
   demBinaryClassification_Marg_fixedhypers(r);
   demBinaryClassification_Marg_fixedhypers_gimala(r);
   demBinaryClassification_Marg_fixedhypers_pMALA(r);
   demBinaryClassification_Ellipt_fixedhypers(r);
   demBinaryClassification_pCN_fixedhypers(r); 
   demBinaryClassification_pCNL_fixedhypers(r); 
end

The files tmcmc.R and tmcmc_ess.R contain R code that implements the experiment in Section 5.2.1 and in A.3.
