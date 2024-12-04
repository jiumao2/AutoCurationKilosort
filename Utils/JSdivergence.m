function JS_divergence = JSdivergence(mu1, sigma1, mu2, sigma2)

% Mixture distribution parameters
mu_M = (mu1 + mu2) / 2;
sigma_M = sqrt((sigma1^2 + sigma2^2) / 2);

% Compute KL divergences
KL_P_M = log(sigma_M / sigma1) + (sigma1^2 + (mu1 - mu_M)^2) / (2 * sigma_M^2) - 0.5;
KL_Q_M = log(sigma_M / sigma2) + (sigma2^2 + (mu2 - mu_M)^2) / (2 * sigma_M^2) - 0.5;

% Calculate JS divergence
JS_divergence = 0.5 * (KL_P_M + KL_Q_M);

end