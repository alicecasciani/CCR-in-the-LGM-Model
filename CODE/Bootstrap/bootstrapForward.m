function [dates, discounts] = bootstrapForward(datesSet, ratesSet, dates, discounts)
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
N_fra = 2;
% We extract the relevant dates and rates from the structs.
% recall that the first element of dates is the settlement date.
settlements = datesSet.fra(1:N_fra,1);
expiries = datesSet.fra(1:N_fra,2);
dates = [dates; expiries];
rates = ratesSet.fra(1:N_fra);
% We compute the forward discount factors
maturieties = yearfrac(settlements, expiries, 2);
fwd_discounts = 1./(1+maturieties.*rates);

% the first four elements of the dates vector are related to the 
% deposits. Same goes for rates. So we begin modifying the discounts vector
% from the fifth element, using the start = 4 variable.
start = 1;
for i=1:N_fra
    % The discount factor at the settlement date is computed by (linear) interpolation.
    settlementDF = interpolation(discounts, dates(1:start+i), settlements(i));
    discounts(start+i+1) = fwd_discounts(i)*settlementDF;
end

end 