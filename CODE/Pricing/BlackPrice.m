function [price] = BlackPrice(dates, discounts, dates_yearly, today, mkt_vol, strike)
% BlackPrice Computes payer swaption price under Black model
%   price = BlackPrice(dates, discounts, dates_yearly, today, mkt_vol, strike)
%
%   INPUTS:
%     dates         - Full timeline for discount factors
%     discounts     - Corresponding discount factors
%     dates_yearly  - Payment grid for the swap (including start and end)
%     today         - Valuation date
%     mkt_vol       - Black volatility
%     strike        - Strike swap rate
%
%   OUTPUT:
%     price         - Payer swaption present value

% Discount factor to swap start
B = interpolation(discounts, dates, dates_yearly(1)); 

% Day count conventions
EU_30_360=6;
ACT_365=3;

% Calculate year fractions between consecutive payment dates
yf = yearfrac(dates_yearly(1:end-1),dates_yearly(2:end), EU_30_360);

% Interpolate discount factors for payment dates
discounts_inter = interpolation(discounts, dates, dates_yearly);
% Compute forward discount factors
fwd_disc = discounts_inter(2:end)./discounts_inter(1:end-1);


% Calculate Basis Point Value (BPV)
BPV = sum(yf.*fwd_disc);
% Forward Swap Rate
S = FwdSwapRate(dates, discounts, dates_yearly);

% d1 and d2 Black term
diff = yearfrac(today, dates_yearly(1), ACT_365);
d1 = log(S./strike)./(mkt_vol*sqrt(diff)) + 0.5*mkt_vol.*sqrt(diff); 
d2 = d1 - mkt_vol*sqrt(diff); 

% Black swaption price
price = B*BPV*(S*normcdf(d1)-strike*normcdf(d2)); 

end