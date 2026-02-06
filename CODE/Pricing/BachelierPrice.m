function price=BachelierPrice(dates, discounts, dates_yearly, today, mkt_vol, strike)
% BachelierPrice Computes payer swaption price under Bachelier model
%   price = BachelierPrice(dates, discounts, dates_yearly, today, mkt_vol, strike)
%
%   INPUTS:
%     dates         - Full timeline for discount factors
%     discounts     - Corresponding discount factors
%     dates_yearly  - Payment grid for the swap (including start and end)
%     today         - Valuation date
%     mkt_vol       - Bachelier volatility
%     strike        - Strike swap rate
%
%   OUTPUT:
%     price         - Payer swaption present value

% Interpolate discount factors for payment dates
B = interpolation(discounts, dates, dates_yearly(2:end)); 

% Day count conventions
EU_30_360 = 6; % for swap payments
ACT_365   = 3; % time-to-maturity measurement for Bachelier

% Year fractions between payment dates (30/360)
yf = yearfrac(dates_yearly(1:end-1),dates_yearly(2:end), EU_30_360);
% Forward swap rate
S = FwdSwapRate(dates, discounts, dates_yearly);
% Bachelier d-value
d=(S-strike)/(mkt_vol*sqrt(yearfrac(today,dates_yearly(1),ACT_365)));

% Price formula: annuity * normal integrals
price=dot(B,yf)*mkt_vol*sqrt(yearfrac(today,dates_yearly(1),EU_30_360))*(normcdf(d)*d+normpdf(d));
end