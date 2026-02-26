function [model, samples, accRates] = gpsampAuxMarg_fixedhypers_pMALA(model, mcmcoptions)
% Proposal:
%   q(y|x) = N(y | x + gamma A_x( grad g(x) - Sigma0^{-1} x ), 2*gamma*A_x )
% where A_x = (delta_x I + Sigma0^{-1})^{-1}, delta_x scalar depends on x.
% Sigma0 eigendecomp: Sigma0 = U * diag(Lambda) * U'

BurnInIters = mcmcoptions.Burnin;
Iters       = mcmcoptions.T;
StoreEvery  = mcmcoptions.StoreEvery;

[n, ~] = size(model.X);
num_stored = floor(Iters/StoreEvery);

samples.F       = zeros(num_stored, n);
samples.G       = zeros(num_stored, n);
samples.Fy      = zeros(num_stored, n);
samples.LogL    = zeros(num_stored, 1);
samples.accprob = zeros(num_stored, 1);
samples.ff      = zeros(1, BurnInIters + Iters);

Y = model.y;
F = model.F;                  % expect 1xn row; we'll enforce row below

U = model.U;
Lambda = model.Lambda(:)';    % 1xn row
LambdaInv = 1 ./ Lambda;

% step size gamma (you previously called it model.delta)
gamma = model.delta;

% Likelihood handles
loglikHandle     = str2func(['logL' model.Likelihood.type]);
gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);

% Ensure F is row
if iscolumn(F), F = F'; end

oldLogLik = sum(loglikHandle(model.Likelihood, Y, F(:)));

% initial grad/curvature and delta_x
[derF, der2F] = gradloglikHandle(model.Likelihood, model.y, F);
derF  = ensure_row(derF);
der2F = ensure_row(der2F);
delta_x = mean(-(der2F));   % scalar depending on x

cnt = 0;
acceptHist = zeros(1, BurnInIters + Iters);

% adaptation params (kept from your code)
range = 0.05;
opt   = 0.574;
epsilon = 0.05;

for it = 1:(BurnInIters + Iters)

    % --- current in eigenbasis ---
    xU     = F * U;          % 1xn
    gradxU = derF * U;       % 1xn

    % A_x eigenvalues: a_i = lambda_i / (1 + delta_x * lambda_i)
    ax = Lambda ./ (1 + delta_x .* Lambda);

    % mean in eigenbasis:
    % meanU = xU + gamma * A_x(grad - Sigma0^{-1}x)
    meanU = xU + gamma .* (ax .* (gradxU - xU .* LambdaInv));

    % sample yU ~ N(meanU, 2*gamma*ax)
    yU = meanU + randn(1, n) .* sqrt(2*gamma .* ax);

    % back to original coordinates
    Fnew = yU * U';

    % likelihood at proposal
    newLogLik = sum(loglikHandle(model.Likelihood, Y, Fnew(:)));

    % grad/curvature at proposal and delta_y
    [derFnew, der2Fnew] = gradloglikHandle(model.Likelihood, model.y, Fnew);
    derFnew  = ensure_row(derFnew);
    der2Fnew = ensure_row(der2Fnew);
    delta_y = mean(-(der2Fnew));

    % --- MH ratio via target + proposal terms ---
    % prior in eigenbasis: log N(0, Lambda)
    oldPrior = log_prior_eig(xU, Lambda);
    newPrior = log_prior_eig(yU, Lambda);

    % log q(y|x)
    logq_y_given_x = log_q_eig(xU, yU, gradxU, delta_x, gamma, Lambda, LambdaInv);

    % log q(x|y)
    gradyU = derFnew * U;
    logq_x_given_y = log_q_eig(yU, xU, gradyU, delta_y, gamma, Lambda, LambdaInv);

    logAlpha = (newLogLik - oldLogLik) + (newPrior - oldPrior) + (logq_x_given_y - logq_y_given_x);

    % accept/reject (avoid dependency on metropolisHastings.m)
    if log(rand) < logAlpha
        accept = 1;
        F = Fnew;
        derF = derFnew;
        delta_x = delta_y;
        oldLogLik = newLogLik;
    else
        accept = 0;
    end
    uprob = min(1, exp(min(0, logAlpha)));

    acceptHist(it) = accept;

    % adapt gamma during burn-in (same logic as yours)
    if mod(it,5) == 0 && it >= 50
        accRate = mean(acceptHist((it-49):it)) * 100;
        if it <= BurnInIters
            if (accRate > (100*(opt+range))) || (accRate < (100*(opt-range)))
                gamma = gamma + (epsilon*((accRate/100 - opt)/opt)) * gamma;
                model.delta = gamma;
            end
        end
    end

    % store after burn-in
    if (it > BurnInIters) && (mod(it,StoreEvery) == 0)
        cnt = cnt + 1;
        samples.F(cnt,:)       = F;
        samples.G(cnt,:)       = F / model.delta;
        samples.Fy(cnt,:)      = Fnew;
        samples.accprob(cnt,1) = uprob;
        samples.LogL(cnt,1)    = oldLogLik;
    end

    samples.ff(it) = 0.5*(F*F');
end

model.F = F;
accRates.F = mean(acceptHist(BurnInIters+1:end)) * 100;
end


