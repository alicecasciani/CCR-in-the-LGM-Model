function [price] = Jamshidian(dates, discounts, dates_years, coupon, sigma, a)
% Jamshidian: Prices swaption using Jamshidian’s decomposition
    % in a one-factor Hull‐White (HJM) framework.
    %
    % Inputs:
    %   dates         : bootstrap dates
    %   discounts     : vector of known discount factors corresponding to “dates”
    %   dates_years   : vector of dates from t_alpha to t_omega 
    %   coupon        : annual coupon rate -> swaption strike
    %   sigma         : Hull‐White volatility parameter
    %   a             : Hull‐White mean‐reversion speed parameter
    %
    % Output:
    %   price         : the price via Jamshidian’s approach

% Set day‐count convention codes for later yearfrac calls:
ACT_360 = 2;  %  Actual/360 (for accrual between payment dates)
ACT_365 = 3;  %  Actual/365 (for yearfractions from today)

% Set the reference date
today = dates(1);
% Set t_alpha (swaption expiry)
t_alpha = dates_years(1);

% Compute the accrual fractions (“deltas”) between successive cashflow dates
% using Actual/360 day count.  
deltas = yearfrac(dates_years(1:end-1), dates_years(2:end), ACT_360);

% Compute yearfraction from “today” to each future date in dates_years using ACT/365
yf = yearfrac(today, dates_years, ACT_365);

% Build the vector of cash‐flow amounts (“coupons”) at each payment date:
%   coupon * accrual_delta for each coupon period, and 
%   the final period gets (1 + coupon*delta) to include redemption of principal.
coupons = coupon * deltas;
coupons(end) = 1 + coupons(end);   % add principal redemption in final cashflow

% Obtain discount factors at the cashflow times “dates_years”
disc = interpolation(discounts, dates, dates_years);
% Compute the forward‐discount factors for periods 1..end by dividing
% by the “t_alpha” discount (disc(1)).
fwd_disc = disc(2:end)/disc(1);

% Define the Hull‐White instantaneous volatility function σ_HJM(s, T):
%   σ_HJM(s,T) = (σ / a) * (1 − exp(−a*(T−s)))
sigmaHJM = @(s,T) sigma/a * (1 - exp(-a * (T - s)));

% Prepare a cell array of function handles fwd_stoc_disc{i}(x), each of which
% gives the stochastic forward‐discount factor when the short‐rate shift is x.
fwd_stoc_disc = cell(1,length(coupons));

for i = 1:length(coupons)
    % Precompute σ_HJM at (t_alpha,t_alpha)
    sigma_prev = sigmaHJM(yf(1),yf(1));    % often zero, since T−s = 0
    % Compute σ_HJM at (t_alpha, T_i) 
    sigma_next = sigmaHJM(yf(1),yf(i+1));

    % Define the integran
    integrand =  @(t) sigmaHJM(t,yf(1)).^2 - sigmaHJM(t,yf(i+1)).^2;
    % Numerically integrate from 0 to t_alpha
    intSigma = quadgk(integrand, 0, yf(1));

    % Create a function handle for the stochastic forward‐discount factor as a function of x
    fwd_stoc_disc{i} = @(x) fwd_disc(i) * exp(-x/sigma * (sigma_next-sigma_prev)-0.5*intSigma);
end

% We will solve price_CB(x) = 1 for x = x* (Jamshidian’s root).
price_CB = @(x) - 1;
for i = 1:length(coupons)
    price_CB = @(x) price_CB(x) + fwd_stoc_disc{i}(x)*coupons(i);
end
x_star = fzero(price_CB,0);

% Having found x*, we can compute the “strike” zero‐coupon bond prices K_i 
% for each maturity T_i that will replicate each coupon payment via a put.
K_i = zeros(1,length(coupons));
for i = 1:length(coupons)
    sigma_prev = sigmaHJM(yf(1),yf(1));
    sigma_next = sigmaHJM(yf(1),yf(i+1));
    integrand =  @(t) sigmaHJM(t,yf(1)).^2 - sigmaHJM(t,yf(i+1)).^2;
    intSigma = quadgk(integrand, 0, yf(1));
    % The strike K_i is the forward‐discount factor under shift x*
    K_i(i) = fwd_disc(i) * exp(-x_star/sigma * (sigma_next-sigma_prev) - 0.5*intSigma);
end


% Now compute the price of a put on a zero‐coupon bond expiring at T_i with strike K_i:
% pricePutZCB(K, today, t_alpha, T, sigma, a, discounts, dates)
%   – external function that prices the put on a ZCB under Hull‐White.
Puts  = zeros(1, length(K_i));
for i = 1:length(K_i)
    K   = K_i(i);                 % strike zero‐coupon price for maturity i
    T   = dates_years(i+1);       % maturity time (in years from today)
    Puts(i) = pricePutZCB(K, today, t_alpha, T, sigma, a, discounts, dates);
end

% Finally, the Jamshidian price is:
%   Price = Σ_{i=1}^N [ coupon_i * Put_i ]
price = sum(coupons.*Puts);

end
