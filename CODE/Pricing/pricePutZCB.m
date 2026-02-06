function price = pricePutZCB(strike, t0, start, expiry, sigma, a, discounts, dates)
% pricePutZCB Price of a European put option on a zero coupon bond under
% Hull-White Model
%   price = pricePutZCB(strike, t0, start, expiry, sigma, a, discounts, dates)
%
%   INPUTS:
%     strike    - Option strike rate
%     t0        - Valuation time
%     start     - Option start date
%     expiry    - Option expiry date
%     sigma     - Calibrated volatility
%     a         - Hull-White mean-reversion speed
%     discounts - Discount factors at 'dates'
%     dates     - Timeline for discount factors
%
%   OUTPUT:
%     price     - Put Option price

% Day count convention
ACT_365 = 3;

% Compute year fractions from valuation to bond maturity and option expiry
t1 = yearfrac(t0, start, ACT_365);   % time to bond maturity
t2 = yearfrac(t0, expiry, ACT_365);  % time to option expiry

% Discount factors at start and expiry
B1 = interpolation(discounts, dates, start);
B2 = interpolation(discounts, dates, expiry);

% Implied Black volatility for zero-coupon bond under HW
sigma_p = sigma*sqrt((1-exp(-2*a*t1))/(2*a))*(1-exp(-a*(t2-t1)))/a;
% Compute d-argument for Black formula
h = 1/sigma_p*log(B2/(B1*strike)) + sigma_p/2; 

% Put price under Black model
price = strike*B1*normcdf(-h+sigma_p) - B2*normcdf(-h);

end