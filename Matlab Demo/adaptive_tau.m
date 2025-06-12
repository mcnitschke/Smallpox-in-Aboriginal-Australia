function tau = adaptive_tau(S, P, remaining_time)
% Compute an adaptive time step (tau) based on system changes.
% Poisson -> variance = mean.

epsilon = 0.03; % Tolerance for changes

% Compute estimated mean and variance of changes
mean_change = P;
variance_change = P; 

% Compute tau candidates
tau1 = min(epsilon * (S ./ mean_change),[],'all');
tau2 = min((epsilon^2) * (S.^2 ./ variance_change),[],'all');

% Ensure finite values
tau1(isnan(tau1) | isinf(tau1)) = Inf;
tau2(isnan(tau2) | isinf(tau2)) = Inf;

% Select minimum tau
tau = min([tau1; tau2]);
tau = min(tau, remaining_time); 
tau = max(tau, 1e-6);

end
