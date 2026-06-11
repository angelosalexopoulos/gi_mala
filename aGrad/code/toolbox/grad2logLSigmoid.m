function der2F = gradlogLSigmoid(lik, Y, F)
%function derF = gradlogLSigmoid(lik, Y, F)
%
%Description: Derivative of sigmoid GP log-likelihood useful for binary classification
%

% derivative for binary classification  
der2F = -(exp(-F(:)))./( ( 1+exp(-F(:))  ).^2 );  % second derivative for binary classification  
