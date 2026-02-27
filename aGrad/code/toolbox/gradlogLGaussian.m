function [derF, der2F] = gradlogLGaussian(lik, Y, F)
%function derF = gardlogLGaussian(Y,F,sigma2)
%
%Description: Gradient of the Log likelihood wrt mean for the one dimensional 
%             Gaussian distribution   
%

% derivative for standard GP regression
derF = (Y(:) - F(:))/exp(lik.logtheta); 
der2F = -1.0/exp(lik.logtheta);