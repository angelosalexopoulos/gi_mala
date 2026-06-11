function [model samples accRates] = gpsampAuxMarg_fixedhypers_gimala_VR(model, mcmcoptions)
% gpsampAuxMarg_fixedhypers_gimala_VR
%
% GI-MALA sampler for latent Gaussian models.
% Variant used by the variance-reduction (VR) demo: stores the extra fields
% samples.Fy, samples.PGhat and samples.accprob that are needed to compute
% the control-variate estimator.  The standard (ESS-comparison) sampler
% gpsampAuxMarg_fixedhypers_gimala does NOT store these fields so that
% memory usage and timing are comparable with the other samplers.

BurnInIters = mcmcoptions.Burnin;
Iters       = mcmcoptions.T;
StoreEvery  = mcmcoptions.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F       = zeros(num_stored, n);
samples.Fy      = zeros(num_stored, n);
samples.PGhat   = zeros(num_stored, n);
samples.LogL    = zeros(num_stored, 1);
samples.accprob = zeros(num_stored, 1);
samples.ff      = zeros(1, BurnInIters + Iters);
Y = model.y;
F = model.F;

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:));

model.delta = 2/sqrt(n);
LambdaInv = 1./model.Lambda;



gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
[derF der2F] = gradloglikHandle(model.Likelihood, model.y, F);
derF  = derF(:);
der2F = der2F(:);

delta_x     = mean(-(der2F));
delta_x_inv = 1/delta_x;
delta_xLambdaH = 1./(model.Lambda + delta_x_inv);
gx = (F + delta_x_inv*derF')*model.U;

partOfMeanSampA  = gx'.*delta_xLambdaH.*model.Lambda*model.delta;
partOfMeanSampB  = sqrt( (2*model.delta-model.delta*model.delta)/delta_x )*sqrt(model.Lambda.*delta_xLambdaH);
partOfMeanSampCV = partOfMeanSampA'*model.U'/model.delta;

partOfMeanMH = -0.5*sum(log(model.Lambda*delta_x+1.0));
partOfMeanMH = partOfMeanMH-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));

cnt = 0;

range   = 0.05;
opt     = mcmcoptions.opt;
epsilon = 0.05;

acceptHistF = zeros(1, BurnInIters + Iters);

for it = 1:(BurnInIters + Iters)

    Fnew = (1-model.delta)*F + (randn(1, n).*partOfMeanSampB' + partOfMeanSampA')*model.U';

    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));

    [derFnew der2Fnew] = gradloglikHandle(model.Likelihood, model.y, Fnew);
    derFnew  = derFnew(:);
    der2Fnew = der2Fnew(:);

    delta_y     = mean(-(der2Fnew));
    delta_y_inv = 1/delta_y;
    delta_yLambdaH        = 1./(model.Lambda + delta_y_inv);
    delta_yLambdaH_Lambda = delta_yLambdaH.*model.Lambda;
    gy = (Fnew + delta_y_inv*derFnew')*model.U;
    partOfMeanSampAnew  = (gy'.*delta_yLambdaH_Lambda)*model.delta;
    partOfMeanSampBnew  = sqrt( (2*model.delta-model.delta*model.delta)/delta_y )*sqrt(delta_yLambdaH_Lambda);
    partOfMeanSampCVnew = partOfMeanSampAnew'*model.U'/model.delta;

    partOfMeanMHnew = -0.5*sum(log(model.Lambda*delta_y+1.0));
    partOfMeanMHnew = partOfMeanMHnew-0.5*(model.delta/(2-model.delta))*(gy*(delta_yLambdaH.*gy'));

    hxy = 0.5*(delta_x/(2*model.delta-model.delta^2))*sum( (Fnew-F -(model.delta/delta_x)*derF').^2 );
    hxy = hxy + partOfMeanMH;
    hyx = 0.5*(delta_y/(2*model.delta-model.delta^2))*sum( (F-Fnew -(model.delta/delta_y)*derFnew').^2 );
    hyx = hyx + partOfMeanMHnew;
    corrFactor = hxy - hyx;

    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0);

    acceptHistF(it) = accept;


    % keep samples after burn in
    if (it > BurnInIters) & (mod(it,StoreEvery) == 0)
        cnt = cnt + 1;
        samples.F(cnt,:)       = F;
        samples.Fy(cnt,:)      = Fnew;
        samples.PGhat(cnt,:)   = partOfMeanSampCV;
        samples.LogL(cnt,1)    = oldLogLik;
        samples.accprob(cnt,1) = uprob;
    end

    if accept == 1
        %accept
        F = Fnew;
        derF = derFnew;
        der2F = der2Fnew;
        gx = gy;
        partOfMeanSampA = partOfMeanSampAnew;
        partOfMeanSampB = partOfMeanSampBnew;
        partOfMeanMH = partOfMeanMHnew;
        oldLogLik = newLogLik;
        delta_x = delta_y;
        delta_x_inv = delta_y_inv;
        delta_xLambdaH=delta_yLambdaH;
        partOfMeanSampCV = partOfMeanSampCVnew;
    end

    % Adapt proposal during burnin
    if mod(it,5) == 0
        if (it >= 50)
            accRateF = mean(acceptHistF((it-49):it))*100;
            if (it <= BurnInIters)
                if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
                    model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;
                    model.delta = min(max(model.delta, 1e-6), 1.99);
                    partOfMeanSampA = gx'.*delta_xLambdaH.*model.Lambda*model.delta;
                    partOfMeanSampB = sqrt( (2*model.delta-model.delta*model.delta)/delta_x )*sqrt(model.Lambda.*delta_xLambdaH);%
                    partOfMeanMH = -0.5*sum(log(model.Lambda*delta_x+1.0));
                    partOfMeanMH = partOfMeanMH-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));
                    partOfMeanSampCV = partOfMeanSampA'*model.U'/model.delta;
                end
            end
        end
    end



    samples.ff(it) = 0.5*(F*F');
end

model.F    = F;
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100;
