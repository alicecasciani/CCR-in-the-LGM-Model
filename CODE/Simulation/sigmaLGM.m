function s = sigmaLGM(a, u, T)
% SIGMAHJM  Volatility function Linear Gaussian Markov (LGM) model
%   a     : mean reversion speed
%   u     : current time (in years)
%   T     : maturity time (in years)

s = (1 - exp(-a * (T - u))) / a;
end