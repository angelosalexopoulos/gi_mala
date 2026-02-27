function [ergmean] = ergMean(x)

ergmean = cumsum(x)./((1:length(x))');


