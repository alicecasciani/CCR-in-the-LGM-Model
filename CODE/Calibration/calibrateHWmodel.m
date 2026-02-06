function [aCalibrated, sigmaCalibrated, finalSSE, S, mkt_prices] = calibrateHWmodel(...
    dates, discounts, today, mkt_vols_raw, expiry, tenors, model)
% CALIBRATEHWMODEL   Calibrate Hull–White model parameters via Jamshidian
%   [aCalibrated, sigmaCalibrated, finalSSE, S, mkt_prices] = calibrateHWmodel(...
%        dates, discounts, today, mkt_vols_raw, expiry, tenors, model)
%
%   This function finds the optimal mean‐reversion speed a and volatility σ
%   of a one‐factor Hull–White model by minimizing the sum of squared errors
%   between market‐observed swaption prices (derived from Black or Bachelier vols)
%   and model‐computed swaption prices using Jamshidian’s decomposition.
%
%   INPUTS:
%     dates        – Vector of datetime points for the zero‐coupon curve
%     discounts    – Vector of P(0, dates(k)), the observed discount factors
%     today        – Valuation date as a datetime scalar
%     mkt_vols_raw – 1×n vector of market swaption volatilities (in decimal)
%     expiry       – 1×n vector of swaption expiry dates (datetime)
%     tenors       – 1×n vector of swaption tenors in years (integers)
%     model        – String: 'BLACK' or 'BACHELIER' to indicate pricing convention
%
%   OUTPUTS:
%     aCalibrated      – Calibrated HW mean‐reversion speed (a)
%     sigmaCalibrated  – Calibrated HW volatility parameter (σ)
%     finalSSE         – Final sum of squared pricing errors at optimum
%     S                – 1×n vector of forward swap rates for each swaption
%     mkt_prices       – 1×n vector of market‐implied swaption prices


% Number of swaptions to calibrate
n = length(mkt_vols_raw);
% Preallocate arrays
mkt_prices = zeros(1, n);
jamsh_prices = cell(1, n);
S = zeros(1, n);

% Loop over each swaption to compute S(i), mkt_prices(i), and jamsh_prices{i}
for i = 1:n
    % 1) Determine the underlying swap’s maturity date
    end_date = businessdayoffset(expiry(i) + calyears(tenors(i)));
    % 2) Construct the yearly payment dates of the swap
    dates_yearly = businessdayoffset(expiry(i):calyears(1):end_date);
    % 3) Compute the forward swap rate S(i) on [expiry, end_date]
    S(i) = FwdSwapRate(dates, discounts, dates_yearly);
    % 4) Convert the raw market volatility into a swaption price via Black or Bachelier formula
    mkt_prices(i) = price_swaption_model(dates, discounts, dates_yearly, today, mkt_vols_raw(i), S(i), model);
     % 5) Create an anonymous function jamsh_prices{i} that, for a given 
        %    x = [a; σ], returns the model price via Jamshidian’s method:
        %    Note: x(1) corresponds to a, and x(2) corresponds to σ.
    jamsh_prices{i} = @(x) Jamshidian(dates, discounts, dates_yearly, S(i), x(2), x(1));
end

% Define the objective function fObj(x):
    %   x is a 2×1 vector [a; σ]. 
    %   For each swaption j, jamsh_prices{j}(x) calls the model price given (a,σ).
    %   We compute the squared differences to the market prices and sum them all.
fObj = @(x) sum((mkt_prices - cellfun(@(F) F(x), jamsh_prices)).^2);


% Initial guess x0 = [a0; σ0]:
%   a0 = 0.1 (10% mean reversion speed)
%   σ0 = 0.02 if using Black, or σ0 = 0.08 if using Bachelier
x0 = [0.1; strcmpi(model, 'BLACK') * 0.02 + strcmpi(model, 'BACHELIER') * 0.08];
lb = [0; 0];

options = optimoptions('fmincon', ...
    'Display','off', ...        % niente output a video
    'Algorithm','interior-point', ...
    'TolX',1e-8, ...
    'TolFun',1e-8);

% Run fmincon to minimize fObj(x) subject to x ≥ [0; 0]
[xOpt, fval] = fmincon(fObj, x0, [], [], [], [], lb, [], [], options);

% Extract calibrated parameters and final sum of squared errors
aCalibrated = xOpt(1);
sigmaCalibrated = xOpt(2);
finalSSE = fval;

% Print summary of calibration results
fprintf('Calibration completed (%s model):\n', model);
fprintf('  a     = %.6f\n', aCalibrated);
fprintf('  sigma = %.6f\n', sigmaCalibrated);
fprintf('  SSE   = %.6e\n', finalSSE);
end
