function [model, samples, accRates] = gpsampAuxMarg_fixedhypers_gimala(model, mcmcoptions)
% gpsampAuxMarg_fixedhypers_gimala  --  GI-MALA sampler (fixed kernel hyperparameters)
%
% Implements the Gaussian-Invariant MALA (GI-MALA) algorithm from Section 3
% of the paper, operating on the marginalised (whitened) GP latent space.
% The step-size delta is adapted during burn-in via a stochastic approximation
% scheme targeting the acceptance rate specified in mcmcoptions.opt.
%
% INPUTS
%   model         struct with fields:
%     .y            (n x 1) observations
%     .X            (n x D) input locations
%     .F            (1 x n) initial latent GP values
%     .K            (n x n) prior covariance matrix (precomputed)
%     .U, .Lambda   eigendecomposition of K: K = U * diag(Lambda) * U'
%     .delta        initial step size (scalar); will be adapted during burn-in
%     .GP           kernel hyperparameter struct (fixed, not updated)
%     .Likelihood   likelihood struct with field .type (e.g. 'Sigmoid', 'LogGaussianCox')
%     .constraints  struct; .kernHyper and .likHyper should be 'fixed'
%   mcmcoptions   struct with fields:
%     .T            number of post-burn-in iterations to store
%     .Burnin       number of burn-in iterations
%     .StoreEvery   store every k-th sample (thinning; set to 1 for no thinning)
%     .Langevin     flag (1 = use gradient; retained for interface compatibility)
%     .opt          target acceptance rate for step-size adaptation (e.g. 0.75)
%
% OUTPUTS
%   model         updated model struct (.F and .delta reflect final MCMC state)
%   samples       struct with fields:
%     .F            (T x n) stored latent samples (post burn-in)
%     .LogL         (T x 1) log-likelihood at each stored sample
%     .accprob      (T x 1) Metropolis acceptance probability at each step
%     .deltax       (T x 1) diagonal preconditioner value at each stored step
%   accRates      struct with field .F: mean acceptance rate (post burn-in)
%
% PAPER CONNECTION
%   The proposal distribution is defined in eq. (X) of the paper.
%   The MH correction factor (corrFactor) corresponds to eq. (X).
%   Step-size adaptation follows the scheme in Section X.X.



% ---- Initialise dimensions and storage arrays ----
BurnInIters = mcmcoptions.Burnin;
Iters = mcmcoptions.T;
StoreEvery = mcmcoptions.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F = zeros(num_stored, n);
samples.G = zeros(num_stored, n);
samples.PGhat = zeros(num_stored, n);
samples.LogL = zeros(num_stored,1);
samples.deltax = zeros(num_stored,1);
samples.deltax2 = zeros(num_stored,1);
samples.ff = zeros(1, num_stored);
samples.Fy = zeros(num_stored,n);
Y = model.y;
F = model.F;


%samples.LogL = zeros(1, BurnInIters + Iters);
samples.accprob = zeros(num_stored,1);
samples.ff = zeros(1, BurnInIters + Iters);

% ---- Evaluate initial log-likelihood and gradient ----
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:));

LambdaInv = 1./model.Lambda;

% Gradient and diagonal curvature of the log-likelihood at current F.
% der2F is the diagonal of the Hessian; delta_x = mean(-d^2 log p / dF^2)
% serves as a scalar preconditioner (see Section 3 of the paper).
gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
[derF, der2F] = gradloglikHandle(model.Likelihood, model.y, F);

% ---- Precompute GI-MALA proposal parameters ----
% delta_x: scalar preconditioner from mean negative curvature
delta_x     = mean(-(der2F));
delta_x_inv = 1/delta_x;

% Spectral quantities for the GI-MALA proposal mean and covariance
delta_xLambdaH = 1./(model.Lambda + delta_x_inv);
gx = (F + delta_x_inv*derF')*model.U;

% Mean shift for the proposal (scaled eigenbasis representation)
partOfMeanSampA  = gx'.*delta_xLambdaH.*model.Lambda*model.delta;
% Standard deviation for proposal noise (eigenspace)
partOfMeanSampB  = sqrt( (2*model.delta-model.delta*model.delta)/delta_x )*sqrt(model.Lambda.*delta_xLambdaH);
% Control-variate term (used in VR variant)
partOfMeanSampCV = partOfMeanSampA'*model.U'/model.delta;

% Log-determinant and quadratic form entering the MH log-ratio
partOfMeanMH = -0.5*sum(log(model.Lambda*delta_x+1.0));
partOfMeanMH = partOfMeanMH-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));


cnt = 0;

Langevin = 1;
range = 0.05;
opt = mcmcoptions.opt;

% adaption step size
epsilon = 0.05;

acceptHistF = zeros(1, BurnInIters + Iters);

% ---- Main MCMC loop (burn-in + sampling) ----
for it = 1:(BurnInIters + Iters)

    % --- Propose new state vector F (GI-MALA proposal, eq. X) ---
    Fnew = (1-model.delta)*F + (randn(1, n).*partOfMeanSampB' + partOfMeanSampA')*model.U';

    % Evaluate likelihood at the proposed state
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));

    % --- Compute MH acceptance ratio ---
    % Gradient and curvature at the proposed state (needed for reverse proposal)
    [derFnew, der2Fnew] = gradloglikHandle(model.Likelihood, model.y, Fnew);
    %der2Fnew = grad2loglikHandle(model.Likelihood, model.y, Fnew);

    delta_y = mean(-(der2Fnew));
    
    
    delta_y_inv = 1/delta_y;
    delta_yLambdaH = 1./(model.Lambda + delta_y_inv);
    delta_yLambdaH_Lambda = delta_yLambdaH.*model.Lambda;
    gy = (Fnew + delta_y_inv*derFnew')*model.U;
    partOfMeanSampAnew = (gy'.*delta_yLambdaH_Lambda)*model.delta;
    partOfMeanSampBnew = sqrt( (2*model.delta-model.delta*model.delta)/delta_y )*sqrt(delta_yLambdaH_Lambda);

    partOfMeanSampCVnew = partOfMeanSampAnew'*model.U'/model.delta;


    partOfMeanMHnew = -0.5*sum(log(model.Lambda*delta_y+1.0));
    partOfMeanMHnew = partOfMeanMHnew-0.5*(model.delta/(2-model.delta))*(gy*(delta_yLambdaH.*gy'));


    % Log proposal density q(Fnew | F) and q(F | Fnew) (quadratic terms)
    hxy = 0.5*(delta_x/(2*model.delta-model.delta^2))*sum( (Fnew-F -(model.delta/delta_x)*derF').^2 );
    hxy = hxy + partOfMeanMH;
    hyx = 0.5*(delta_y/(2*model.delta-model.delta^2))*sum( (F-Fnew -(model.delta/delta_y)*derFnew').^2 );
    hyx = hyx + partOfMeanMHnew;
    % corrFactor = log q(F|Fnew) - log q(Fnew|F); added to log-likelihood ratio in MH step
    corrFactor = hxy - hyx;

    [accept, uprob] = metropolisHastings(newLogLik + corrFactor, oldLogLik, 0, 0);

    acceptHistF(it) = accept;

    if accept == 1
        F = Fnew;
        derF = derFnew;
        der2F = der2Fnew;
        partOfMeanSampA = partOfMeanSampAnew;
        partOfMeanSampB = partOfMeanSampBnew;
        partOfMeanMH = partOfMeanMHnew;
        oldLogLik = newLogLik;
        delta_x = delta_y;
        partOfMeanSampCV = partOfMeanSampCVnew;
    end

    %if model.Likelihood.type ~= 'Gaussian'
    % --- Step-size adaptation (burn-in only) ---
    % Every 5 iterations, check the acceptance rate over the last 50 steps.
    % If it falls outside [opt-range, opt+range], rescale delta proportionally.
    % Adaptation is frozen after burn-in to ensure valid MCMC ergodicity.
    if mod(it, 5) == 0
        if (it >= 50)
            accRateF = mean(acceptHistF((it-49):it))*100;
            if (it <= BurnInIters)
                if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
                    model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;
                end
            end
        end
    end
    %end
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.F(cnt,:) = F;
        samples.G(cnt,:) = F/model.delta;
        samples.Fy(cnt,:) = Fnew;
        samples.accprob(cnt,1) = uprob;
        samples.LogL(cnt,1) = oldLogLik;
        samples.deltax(cnt,1) = -mean(der2F);
        samples.deltax2(cnt,1) = max(-der2F);
       
        samples.PGhat(cnt,:) = partOfMeanSampCV;

        %
    end

    %samples.LogL(it) = oldLogLik;
    samples.ff(it) = 0.5*(F*F');
    %
    %if mod(it,1000) == 0
    %it
    %end
end
%
%
model.F = F;
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100;
model.delta = model.delta;
