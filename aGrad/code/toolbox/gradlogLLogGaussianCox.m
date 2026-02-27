function [derF,der2F] = gradlogLLogGaussianCox(lik, Y, F)
%function derF = gradlogLLogGaussianCox(lik, Y, F)
%
%Description: Derivative of the Poisson GP log-likelihood useful for counts data 
%

forder = (lik.m*exp(lik.logtheta(1)))*exp(F(:));
derF = Y(:) - forder;         
der2F = -forder; 