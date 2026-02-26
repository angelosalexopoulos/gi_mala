function lq = log_q_eig(xU, yU, gradU, delta_x, gamma, Lambda, LambdaInv)
    % q(y|x) in eigenbasis:
    % meanU = xU + gamma * a_x .* (gradU - xU ./ Lambda)
    % cov   = 2*gamma * diag(a_x),   a_x = Lambda ./ (1 + delta_x * Lambda)
    ax = Lambda ./ (1 + delta_x .* Lambda);
    meanU = xU + gamma .* (ax .* (gradU - xU .* LambdaInv));
    varU  = 2*gamma .* ax;

    diff = yU - meanU;
    % full Gaussian log-density (constant cancels in differences but safe)
    lq = -0.5 * sum( log(2*pi .* varU) + (diff.^2) ./ varU );
end