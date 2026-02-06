function [discounts] = interpolation(discounts, dates, target)
% interpolation: given a set of discount factors at known dates, compute
% the interpolated discount factor at a new target date via linear-zero-rate interpolation.
%
% Inputs:
%   discounts : vector of existing discount factors
%   dates     : corresponding vector of dates
%   target    : date at which we want the interpolated discount
%
% Output:
%   discounts : (re)assigned vector = the discount factors interpolated at
%   the target dates

% Day-count convention
ACT_365 = 3;

% If input dates are datetime arrays, convert them to MATLAB datenum
if isa(dates,'datetime')
    dates = datenum(dates);
end

if isa(target,'datetime')
    target = datenum(target);
end

% Reference date for all interpolations
today = datenum('28-Jun-2022');

% If the provided discounts vector has the same length as dates,
% assume the first element is a “t=0” discount (i.e. equals 1) and drop it.
% In that case, we only interpolate using the tail (from index 2 onward).
% Otherwise, assume discounts/dates already correspond one-to-one.
if numel(discounts)==numel(dates)
    dv = discounts(2:end);
    tv = dates(2:end);
else
    dv = discounts;
    tv = dates;
end

% Compute continuously compounded zero rates from today up to each known date tv:
%   zeroRate(i) = –ln(discount_i) / yearfrac(today, tv_i, ACT_365)
zeroRate = -log(dv) ./ yearfrac(today, tv, ACT_365);

% If there is only one zero rate available, that is our interpolated zero rate.
% Otherwise, linearly interpolate zeroRate vs. known dates (tv) to the target date.
if numel(zeroRate)<2
    interp_zr = zeroRate(end);
else
    % interp1 linearly interpolates the zero rate at “target”;
    % if target is outside range of tv, we hold the last zero rate constant
    interp_zr = interp1(tv, zeroRate, target, 'linear', zeroRate(end));
end


% Finally, convert the interpolated zero rate back to a discount factor:
%   discount(target) = exp(–z * yearfrac(today, target, ACT_365))
discounts = exp(-interp_zr .* yearfrac(today, target, ACT_365));

end
