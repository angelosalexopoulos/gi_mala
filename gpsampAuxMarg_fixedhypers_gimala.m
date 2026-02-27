function [model samples accRates] = gpsampAuxMarg_fixedhypers_gimala(model, mcmcoptions)


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

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:));


LambdaInv = 1./model.Lambda;



% Quantities that need to be updated when the state vector
% changes
gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
[derF,der2F] = gradloglikHandle(model.Likelihood, model.y, F);


delta_x = mean(-(der2F));


delta_x_inv = 1/delta_x;
delta_xLambdaH = 1./(model.Lambda + delta_x_inv);
gx = (F + delta_x_inv*derF')*model.U;

partOfMeanSampA = gx'.*delta_xLambdaH.*model.Lambda*model.delta;
partOfMeanSampB = sqrt( (2*model.delta-model.delta*model.delta)/delta_x )*sqrt(model.Lambda.*delta_xLambdaH);%
partOfMeanSampCV = partOfMeanSampA'*model.U'/model.delta;

partOfMeanMH = -0.5*sum(log(model.Lambda*delta_x+1.0));
partOfMeanMH = partOfMeanMH-0.5*(model.delta/(2-model.delta))*(gx*(delta_xLambdaH.*gx'));


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
    % new gradient
    [derFnew,der2Fnew] = gradloglikHandle(model.Likelihood, model.y, Fnew);
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


    hxy = 0.5*(delta_x/(2*model.delta-model.delta^2))*sum( (Fnew-F -(model.delta/delta_x)*derF').^2 );
    hxy = hxy +partOfMeanMH;
    hyx = 0.5*(delta_y/(2*model.delta-model.delta^2))*sum( (F-Fnew -(model.delta/delta_y)*derFnew').^2 );
    hyx = hyx +partOfMeanMHnew;
    corrFactor = hxy - hyx;

    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0);

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
    % Adapt proposal during burnin
    if mod(it,5) == 0
        if (it >= 50)
            accRateF = mean(acceptHistF((it-49):it))*100;
            if (it <= BurnInIters)
                if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
                    %
                    model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;
                    %
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
