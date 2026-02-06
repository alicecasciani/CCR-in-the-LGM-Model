function [dates, discounts, zRates]=bootstrap(datesSet, ratesSet)
% This function bootstraps the discount factor curve
% from the available market data. We consider three different
% kind of instruments: deposits, futures and swaps.
% INPUTS
% datesSet: structure with the dates of the deposits - futures - swaps
% ratesSet: structure with the rates of the deposits - futures - swaps
% OUTPUT
% dates format datenum
% discounts: discount factors
% zero rates
[dates, discounts] = bootstrapDeposits(datesSet, ratesSet);
[dates, discounts] = bootstrapFutures(datesSet, ratesSet, dates, discounts);
[dates, discounts] = bootstrapSwaps(datesSet, ratesSet, dates, discounts);

% Compute Zero Rates
zRates = 100*zeroRates(dates, discounts);
% plot results
% discount factors
figure;
plot(dates, discounts,'b-*', 'LineWidth', 1.5);
grid on;
title('Discount Factors');
% zero-rates
figure;
plot(dates,zRates,'b-*', 'LineWidth', 1.5); 
grid on;
title('Zero Rates');

end     % function bootstrap