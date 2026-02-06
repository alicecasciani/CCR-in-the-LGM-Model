function sigma = calibNthSwaption(price, expDate, tenor, dates, discounts, aLGM, prevSigmas, timeGrid)
% CALIBNTHSWAPTION   Optimize σ(i) to minimize squared pricing error
%   sigma = calibNthSwaption(price, expDate, tenor, dates, discounts, aLGM, prevSigmas, timeGrid)
%
%   INPUTS:
%     price       – market price of the i-th swaption
%     expDate     – expiry date of the swaption (datetime)
%     tenor       – tenor in years (integer) of the underlying swap
%     dates       – vector of datetime points for the discount curve
%     discounts   – vector of P(0, dates(k)) discount factors
%     aLGM        – mean‐reversion speed κ of the LGM model
%     prevSigmas  – 1×(i−1) vector of previously calibrated σ(1:i−1)
%     timeGrid    – vector of time grid points corresponding to interval boundaries
%
%   OUTPUT:
%     sigma       – calibrated volatility σ(i) for the i-th swaption
%
%   STEPS:
%     1. Build the swap‐payment schedule from expDate for tenor years.
%     2. Compute the forward swap rate SR.
%     3. Define the objective function as the squared difference between the market
%        price and the model price from JamshidianLGM(...).
%     4. Use fmincon to find σ(i) ≥ 0 minimizing this squared error, starting at 0.02

% 1) Compute the underlying swap’s end date and yearly payment dates
end_date     = businessdayoffset(expDate + calyears(tenor));
dates_yearly = businessdayoffset(expDate:calyears(1):end_date);

% 2) Compute the forward swap rate from the full discount curve
SR           = FwdSwapRate(dates, discounts, dates_yearly);

% 3) Build the squared‐error objective
obj = @(s) (price - JamshidianLGM(dates, discounts, dates_yearly, SR, prevSigmas, s, aLGM, timeGrid)).^2;
% 4) Set up fmincon to solve: minimize obj(s) subject to s ≥ 0
opts = optimoptions('fmincon', 'Display','off', 'Algorithm','interior-point', 'TolX',1e-8, 'TolFun',1e-8);
sigma = fmincon(obj, 0.02, [], [], [], [], 0, [], [], opts);
end