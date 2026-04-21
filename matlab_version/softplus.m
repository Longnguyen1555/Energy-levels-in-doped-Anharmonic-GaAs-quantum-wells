function y = softplus(x)
% Stable log(1 + exp(x)).
    y = log1p(exp(-abs(x))) + max(x, 0.0);
end
