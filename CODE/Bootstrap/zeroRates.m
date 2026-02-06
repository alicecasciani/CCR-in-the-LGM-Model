function zRates = zeroRates(dates, discounts)
% INPUT:
%   dates     : Complete set of dates used in the bootstrap
%   discounts : Vector of discount factors obtained from bootstrap
%
% OUTPUT:
%   zRates    : Vector of zero rates computed for each date

    % Set the settlement date
    t0 = datetime('28-Jun-2022');

    % Calculate the time intervals (in years) between the settlement date and each provided date
    % using day count convention 3 (ACT/365)
    deltas = yearfrac(t0, dates, 3);
      
    % Compute the zero-coupon rates using the formula:
    %     zRate = -ln(discount factor) / time interval,
    % which comes from the relation: discount factor = exp(-zRate * time interval)
    zRates = -log(discounts)./deltas;
end
