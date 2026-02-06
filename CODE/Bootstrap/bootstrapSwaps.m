function [dates, discounts] = bootstrapSwaps(datesSet, ratesSet, dates, discounts)
% Boostrapping the deposit rates to get the discount factors.
% INPUT
% datesSet: structure with the dates of the deposits
% ratesSet: structure with the rates of the deposits
% dates: previously defined dates. Deposits and Futures date
% discounts: previously defined discounts. Deposits and Futures discounts
% OUTPUT
% dates: dates of the deposits, format datenum
% discounts: discount factors of the deposits

% we extract the relevant dates and rates from the structs.
t0 = datesSet.settlement;
expiries = datesSet.swaps;
rates = ratesSet.swaps;
EU_30_360 = 6;
% To begin calibrating on interest rate swaps, we need to infer
% the discount factor at one year, as the complete set of swap rates 
% starts from 2 years. We compute the one year discount factor by
% linear interpolation. Recall that this discount factor is not supposed to
% be used in the swap curve. However, it is needed in the next points,
% for instance to compute the DV01 or price the bond. We could've 
oneyeardf = interpolation(discounts, dates, businessdayoffset(datetime('28-Jun-2023')));

[complete, rates] = swap_spline(expiries, rates);
yf_complete = yearfrac(t0, complete, EU_30_360);

deltas = yearfrac(t0, complete(1), EU_30_360);
deltas = [deltas, yf_complete(2:end)-yf_complete(1:end-1)];

discount_swaps = zeros(1,length(yf_complete));
discount_swaps(1) = oneyeardf;

for i=2:length(yf_complete)
    numerator = 1-rates(i)*sum(deltas(1:i-1).*discount_swaps(1:i-1));

    denominator = 1+deltas(i)*rates(i);
    discount_swaps(i) = numerator/denominator;
end

discounts = [discounts; discount_swaps(2:end)'];
dates = [datetime(dates,'ConvertFrom','datenum'); complete(2:end)'];
end