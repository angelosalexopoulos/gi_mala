function [derF, der2F] = gradlogLSigmoid(lik, Y, F)
%function derF = gradlogLSigmoid(lik, Y, F)
%
%Description: Derivative of sigmoid GP log-likelihood useful for binary classification
%
expF = exp(-F(:));
% derivative for binary classification  
derF = 0.5*(Y(:)+1) - 1./(1+expF);  % derivative for binary classification  
der2F = -expF./( ( 1+expF  ).^2 );  % second derivative for binary classification  
