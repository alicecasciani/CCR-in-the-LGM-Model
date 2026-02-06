function [dates, discounts] = bootstrapFutures(datesSet, ratesSet, dates, discounts)
% Boostrapping the deposit rates to get the discount factors.
% INPUT
% datesSet: structure with the dates of the deposits
% ratesSet: structure with the rates of the deposits
% dates: previously defined dates. Deposits date
% discounts: previously defined discounts. Deposits discounts
% OUTPUT
% dates: dates of the deposits, format datenum
% discounts: discount factors of the deposits

% We consider only the first 7 futures, as they are the most liquid.
N_futures = 7;
% We extract the relevant dates and rates from the structs.
% recall that the first element of dates is the settlement date.
settlements = datesSet.futures(1:N_futures,1);
expiries = datesSet.futures(1:N_futures,2);
dates = [dates; expiries];
rates = ratesSet.futures(1:N_futures);

% We compute the forward discount factors
maturieties = yearfrac(settlements, expiries, 2);
fwd_discounts = 1./(1+maturieties.*rates);

% the first element of the dates vector is related to the 
% deposit. Same goes for rates. So we begin modifying the discounts vector
% from the second element, using the start = 1 variable.
start = 1;
for i=1:N_futures
    % The discount factor at the settlement date is computed by (linear) interpolation.
    settlementDF = interpolation(discounts, dates(1:start+i), settlements(i));
    discounts(start+i+1) = fwd_discounts(i)*settlementDF;
end

end     % function bootstrapFutures