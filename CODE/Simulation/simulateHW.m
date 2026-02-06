function xnext = simulateHW(g, a, sigma, xprev, tprev, tnext)
% SIMULATEHW   Simulates one step of the Hull–White short-rate process
%   g      : matrix of Gaussian random numbers (same size as xprev)
%   a      : mean reversion speed
%   sigma  : volatility of the Hull–White process (σ parameter)
%   xprev  : matrix/array of values x(tprev)
%   tprev  : datetime scalar, previous time point
%   tnext  : datetime scalar, next time point
%
% Returns xnext, a matrix of the same size as xprev/g, containing the
% simulated values x(tnext).

% Calculate the time increment τ in years using ACT/365 convention
ACT_365 = 3; 
tau = yearfrac(tprev, tnext, ACT_365);

% Compute the conditional mean term of the Ornstein-Uhlenbeck process:
%   mean_term = xprev * exp(-a * τ)
% This reflects the decay of the previous state xprev over the interval τ
mean_term = xprev .* exp(-a * tau);

% Compute the incremental variance of x over the interval [tprev, tnext]:
%   var_term = (σ^2 / (2a)) * (1 - exp(-2a * τ))
% This is the exact variance for the OU process over time τ
var_term  = sigma^2 * (1 - exp(-2 * a * tau)) / (2 * a);

% Compute the standard deviation from the variance
std_term  = sqrt(var_term);

% Finally, generate the next state xnext by adding Gaussian noise:
%   xnext = mean_term + std_term .* g'
% where g' are independent standard normals
xnext = mean_term + std_term .* g';
end
