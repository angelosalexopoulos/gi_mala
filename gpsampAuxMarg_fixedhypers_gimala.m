function [model samples accRates] = gpsampAuxMarg_fixedhypers_angelos(model, mcmcoptions)
%function [model samples accRates] = gpsampAuxMarg_fixedhypers(model, mcmcoptions)
%
% The pCN MALA for latent Gaussian models as described in overleaf notes
% pCN_Poisson
%
% Inputs:
%         -- model: the structure that contains the log likelihood and the latent
%                    Gaussian prior parameters such as the covariance matrix
%                    with its pre-computed eigendecomposition (see demos)
%         -- mcmcoptions: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model:
%         -- model: as above. The outputed model is updated to contain the
%                   state vector of the final MCMC iteration
%                   as well as the learned step size delta
%         -- samples: the structure that contrains the samples
%         -- accRates: acceptance rate
%
% 


BurnInIters = mcmcoptions.Burnin;
Iters = mcmcoptions.T;
StoreEvery = mcmcoptions.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F    = zeros(num_stored, n);
samples.LogL = zeros(num_stored, 1);
Y = model.y;
F = model.F;

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:));

model.delta = 2/sqrt(n); %understand this as the gamma parameter appearing in overleaf notes
%model.delta = 1.0;

gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
[derF der2F] = gradloglikHandle(model.Likelihood, model.y, F);

%derF  = derF(:);
%der2F = der2F(:);

delta_x = -mean(der2F);



delta_x_inv = 1/delta_x;
delta_xLambdaH = 1./(model.Lambda + delta_x_inv);
gx = (F + delta_x_inv*derF')*model.U;

partOfMeanSampA = gx'.*delta_xLambdaH.*model.Lambda*model.delta;
partOfMeanSampB = sqrt( (2*model.delta-model.delta*model.delta)/delta_x )*sqrt(model.Lambda.*delta_xLambdaH);%

partOfMeanMH1 = -0.5*sum(log(model.Lambda*delta_x+1.0));
partOfMeanMH2 = partOfMeanMH1-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));


cnt = 0;

Langevin = 1;
range = 0.05;
opt = mcmcoptions.opt;

% adaption step size
epsilon = 0.05;

acceptHistF = zeros(1, BurnInIters + Iters);

for it = 1:(BurnInIters + Iters)
    %
    % Propose new state vector F
    Fnew =(1-model.delta)*F + (randn(1, n).*partOfMeanSampB' +  partOfMeanSampA' )*model.U';

    % perform an evaluation of the likelihood p(Y | F)
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));

    % Metropolis-Hastings to accept-reject the proposal

    [derFnew der2Fnew] = gradloglikHandle(model.Likelihood, model.y, Fnew);
    %derFnew  = derFnew(:);
    %der2Fnew = der2Fnew(:);


    delta_y = -mean(der2Fnew);

    
    delta_y_inv = 1/delta_y;
    delta_yLambdaH = 1./(model.Lambda + delta_y_inv);
    delta_yLambdaH_Lambda = delta_yLambdaH.*model.Lambda;
    gy = (Fnew + delta_y_inv*derFnew')*model.U;
    partOfMeanSampAnew = (gy'.*delta_yLambdaH_Lambda)*model.delta;
    partOfMeanSampBnew = sqrt( (2*model.delta-model.delta*model.delta)/delta_y )*sqrt(delta_yLambdaH_Lambda);

    partOfMeanMHnew1 = -0.5*sum(log(model.Lambda*delta_y+1.0));
    partOfMeanMHnew2 = partOfMeanMHnew1-0.5*(model.delta/(2-model.delta))*(gy*(delta_yLambdaH.*gy'));


    hxy = 0.5*(delta_x/(2*model.delta-model.delta^2))*sum( (Fnew-F -(model.delta/delta_x)*derF').^2 );
    hxy = hxy +partOfMeanMH2;
    hyx = 0.5*(delta_y/(2*model.delta-model.delta^2))*sum( (F-Fnew -(model.delta/delta_y)*derFnew').^2 );
    hyx = hyx +partOfMeanMHnew2;
    corrFactor = hxy - hyx;

    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0);

    acceptHistF(it) = accept;

    if accept == 1
        %accept
        F = Fnew;
        derF = derFnew;
        der2F = der2Fnew;
        gx = gy;
        partOfMeanSampA = partOfMeanSampAnew;
        partOfMeanSampB = partOfMeanSampBnew;
        partOfMeanMH2 = partOfMeanMHnew2;
        partOfMeanMH1 = partOfMeanMHnew1;
        oldLogLik = newLogLik;
        delta_x = delta_y;
        delta_x_inv = delta_y_inv;
        delta_xLambdaH=delta_yLambdaH;
    end

    %if model.Likelihood.type ~= 'Gaussian'
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
                    %partOfMeanMH = -0.5*sum(log(model.Lambda*delta_x+1.0));
                    partOfMeanMH2 = partOfMeanMH1-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));
                end
            end
        end
    end
    %end
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.F(cnt,:)    = F;
        samples.LogL(cnt,1) = oldLogLik;
   

        %
    end
end
%
%
model.F = F;
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100;
model.delta = model.delta;
