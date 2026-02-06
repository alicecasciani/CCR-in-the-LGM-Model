function price = pricePutLGM(strike, t0, start, expiry, sigma, a, discounts, dates, time_grid)
% pricePutLGM Price of a European put option on a zero coupon bond under
% LGM model
%   price = pricePutLGM(strike, t0, start, expiry, sigma, a, discounts, dates, time_grid)
%
%   INPUTS:
%     strike    - Option strike rate
%     t0        - Valuation time
%     start     - Option start date
%     expiry    - Option expiry date
%     sigma     - Vector of calibrated volatilities along time grid
%     a         - Mean-reversion speed
%     discounts - Discount factors at 'dates'
%     dates     - Timeline for discount factors
%     time_grid - Grid of monitoring dates for sigma
%
%   OUTPUT:
%     price     - Option price

% Day count convention
ACT_365 = 3;

% Discount factors at start and expiry
B1 = interpolation(discounts, dates, start);
B2 = interpolation(discounts, dates, expiry);

% Compute integrated variance across intervals
sigma_square = zeros(1, length(sigma));
H=@(s,T) (1 - exp(-a * (T-s)))/a;
for i=1:length(sigma)
    integrand=@(s) sigma(i).^2.*(H(s,yearfrac(t0,expiry,ACT_365))-H(s,yearfrac(t0,time_grid(end),ACT_365))).^2;
    int = quadgk(integrand,yearfrac(t0,time_grid(i)),yearfrac(t0,time_grid(i+1),ACT_365));
    sigma_square(i) = int;
end
sigma_square = sum(sigma_square);

% Black-style d1, d2
d1 = log(B2/(strike*B1))./sqrt(sigma_square) + 0.5*sqrt(sigma_square);
d2 = log(B2/(strike*B1))./sqrt(sigma_square) - 0.5*sqrt(sigma_square);

% Option price under risk-neutral measure
price = -B2*normcdf(-d1)+B1*strike*normcdf(-d2);

end