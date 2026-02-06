function sigmaLGM = calibrateJamshidian(prices, expiry, tenors, dates, discounts, aLGM, time_grid)
% CALIBRATEJAMSHIDIAN   Bootstraps LGM volatilities via Jamshidian’s trick
%   sigmaLGM = calibrateJamshidian(prices, expiry, tenors, dates, discounts, aLGM, time_grid)
%
%   INPUTS:
%     prices      – 1×n vector of market‐observed swaption prices to match
%     expiry      – 1×n vector of swaption expiry dates (datetime) 
%     tenors      – 1×n vector of swaption tenors in years (integer)
%     dates       – vector of datetime points for the full discount curve
%     discounts   – vector of P(0, dates(k)) discount factors
%     aLGM        – mean‐reversion speed κ of the LGM model
%     time_grid   – vector of time grid points at which σ intervals change 
%
%   OUTPUT:
%     sigmaLGM    – 1×n vector of calibrated piecewise σ_i parameters for the LGM model
%
%   The function calibrates each σ_i sequentially:
%     1. The first swaption’s σ(1) is found via a root‐search on Jamshidian’s formula.
%     2. Each subsequent σ(i) is found by minimizing the squared pricing error of
%        a swaption pricing function JamshidianLGM that depends on all previously
%        calibrated σ(1:i−1) plus σ(i) itself.


n = numel(expiry);
sigmaLGM = zeros(1, n);

% --- 1) Calibrate the first swaption’s volatility via fzero
sigmaLGM(1) = calibFirstSwaption(prices(1), expiry(1), tenors(1), dates, discounts, aLGM);
% --- 2) Calibrate the remaining swaptions one by one
for i = 2:n
    sigmaLGM(i) = calibNthSwaption(prices(i), expiry(i), tenors(i), dates, discounts, aLGM, sigmaLGM(1:i-1), time_grid(1:i+1));
end
end