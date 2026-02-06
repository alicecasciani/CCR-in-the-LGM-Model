function xnext = simulateLGM(g, a, sigma, xprev, tprev, tnext)
% SIMULATELGM   Simulates one step of the Linear Gaussian Markov (LGM) process
%   g      : matrix of Gaussian random numbers (same size as xprev)
%   a      : mean reversion speed
%   sigma  : volatility of the LGM process (σ parameter)
%   xprev  : matrix/array of values x(tprev)
%   tprev  : datetime scalar, previous time point
%   tnext  : datetime scalar, next time point
%
% Returns xnext, a matrix of the same size as xprev/g, containing the
% simulated values x(tnext).

% Calculate the time increment τ in years using ACT/365 convention
ACT_365 = 3; 
tau = yearfrac(tprev, tnext, ACT_365);

% Compute the conditional mean term for the Ornstein–Uhlenbeck process:
%   mean_term = xprev * exp(-a * τ)
% This represents the decay of xprev over the time interval τ
mean_term = xprev .* exp(-a * tau);

% Compute the incremental variance of x over [tprev, tnext]:
%   var_term = (σ^2 / (2a)) * (1 - exp(-2a * τ))
% This is the exact variance contribution for the OU process across τ
var_term  = sigma^2 * (1 - exp(-2 * a * tau)) / (2 * a);

% Compute the standard deviation from the variance
std_term  = sqrt(var_term);

% Generate the next state xnext by adding Gaussian noise:
%   xnext = mean_term + std_term .* g'
% where g' are independent standard normals
xnext = mean_term + std_term .* g';
end
