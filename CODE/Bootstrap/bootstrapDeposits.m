function [dates,discounts] = bootstrapDeposits(datesSet, ratesSet)
% Boostrapping the deposit rates to get the discount factors.
% INPUT
% datesSet: structure with the dates of the deposits
% ratesSet: structure with the rates of the deposits
% OUTPUT
% dates: dates of the deposits, format datenum
% discounts: discount factors of the deposits

% The computations start from t0, the settlement date: 28-Jun-2022
t0 = datesSet.settlement;
% We only have one deposit. Later points of the curve
% are computed with more liquid assets.
N_depos = 1;

% We extract the relevant dates and rates from the structs. It's convenient to
% set the first date as t0, the settlement date. The corresponding
% discount factor will be 1.
dates = datesSet.depos;
dates = [t0; dates];
rates = ratesSet.depos;
rates = [1; rates];

% We compute the year fractions, with respect to the ACT/360 convention
ACT_360 = 2;
yf = yearfrac(t0, dates, ACT_360);
discounts = 1./(1+yf.*rates);

end     % function bootstrapDeposits