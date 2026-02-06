function sigma = calibFirstSwaption(price, expDate, tenor, dates, discounts, aLGM)
% CALIBFIRSTSWAPTION   Solve for σ(1) by equating model price to market price
%   sigma = calibFirstSwaption(price, expDate, tenor, dates, discounts, aLGM)
%
%   INPUTS:
%     price       – market price of the first swaption
%     expDate     – expiry date of the swaption (datetime)
%     tenor       – tenor in years (integer) of the underlying swap
%     dates       – vector of datetime points for the discount curve
%     discounts   – vector of P(0, dates(k)) discount factors
%     aLGM        – mean‐reversion speed of the LGM model
%
%   OUTPUT:
%     sigma       – calibrated volatility σ(1) for the first swaption
%
%   STEPS:
%     1. Construct the set of yearly payment dates for the underlying swap.
%     2. Compute the forward swap rate SR from the discount curve.
%     3. Define a root‐finding objective f0(s) = market_price − JamshidianModelPrice(s).
%     4. Solve f0(s) = 0 via fzero, starting from an initial guess of 0.01.

% Compute the swap’s end date by adding tenor years
end_date     = businessdayoffset(expDate + calyears(tenor));

% Build the vector of yearly payment dates: from expiry, stepping by 1 calendar year until end_date
dates_yearly = businessdayoffset(expDate:calyears(1):end_date);

% Compute the forward swap rate SR on [expDate, end_date]
SR           = FwdSwapRate(dates, discounts, dates_yearly);

% Define the objective: difference between market price and Jamshidian‐LGM model price with σ = s.
f0           = @(s) price - Jamshidian(dates, discounts, dates_yearly, SR, s, aLGM);

% Solve for σ using fzero, with initial guess of 0.01 (1%).
sigma = fzero(f0, 0.01);
end