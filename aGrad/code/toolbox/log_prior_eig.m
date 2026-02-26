function lp = log_prior_eig(xU, Lambda)
    % log N(xU | 0, diag(Lambda)) up to additive constant
    lp = -0.5 * sum( (xU.^2) ./ Lambda );
end