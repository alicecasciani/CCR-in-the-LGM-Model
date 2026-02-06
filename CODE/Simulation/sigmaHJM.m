function s = sigmaHJM(a, sigma, u, T)
% SIGMAHJM  Volatility function Hull–White model
%
%   INPUTS
%   a     : mean reversion speed
%   sigma : volatility level
%   u     : current time (in years)
%   T     : maturity time (in years)
%
%   OUTPUT
%   sigma

s = sigma^2 * (1 - exp(-a * (T - u))) / a;
end