function ESS = mcmc_ess(Samples)
% Samples: NumOfSamples x NumOfParameters (rows = iterations)

% If a vector is given, make it a column (NumOfSamples x 1)
if isvector(Samples)
    Samples = Samples(:);
end

% If it looks like parameters are in rows (few rows, many cols), transpose
% (common mistake when storing samples as NumParams x NumSamples)
if size(Samples,1) < size(Samples,2)
    Samples = Samples.';
end

[NumOfSamples, NumOfParameters] = size(Samples);

% If chain too short, ESS equals number of samples (no room for autocorr)
if NumOfSamples < 2
    ESS = NumOfSamples * ones(1, NumOfParameters);
    return;
end

% Choose a safe lag: must be >=1. Also don't use insanely large lags.
MaxLag = min(NumOfSamples - 1, 1000);  % you can change 1000 to e.g. floor(NumOfSamples/2)

ACs = zeros(MaxLag+1, NumOfParameters);
for i = 1:NumOfParameters
    %ACs(:,i) = autocorr(Samples(:,i), MaxLag);  % this syntax works in many versions
    ACs(:,i) = autocorr(Samples(:,i), 'NumLags', MaxLag);
    % If your MATLAB expects name/value, use:
    % ACs(:,i) = autocorr(Samples(:,i), 'NumLags', MaxLag);
end

M = floor(size(ACs,1)/2);
Gamma = zeros(M, NumOfParameters);
for j = 1:M
    Gamma(j,:) = ACs(2*j-1,:) + ACs(2*j,:);
end

for j = 2:M
    Gamma(j,:) = min(Gamma(j-1:j,:), [], 1);
end

MonoEst = zeros(1, NumOfParameters);
for i = 1:NumOfParameters
    PosGammas = find(Gamma(:,i) > 0);
    if isempty(PosGammas)
        MonoEst(i) = 1;
    else
        MonoEst(i) = -ACs(1,i) + 2*sum(Gamma(1:PosGammas(end), i));
        if MonoEst(i) < 1
            MonoEst(i) = 1;
        end
    end
end

ESS = NumOfSamples ./ MonoEst;
end