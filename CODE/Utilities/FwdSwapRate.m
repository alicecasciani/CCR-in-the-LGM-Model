function swapRate = FwdSwapRate(dates, discounts, paymentDates)
% FwdSwapRate Computes the forward swap rate for an interest rate swap
%   swapRate = FwdSwapRate(dates, discounts, paymentDates) returns the fixed
%   rate that makes the net present value of the floating and fixed legs equal.
%
%   INPUTS:
%     dates        - Vector of dates corresponding to discount factors
%     discounts    - Vector of discount factors at each date
%     paymentDates - Vector of future payment dates for the swap
%
%   OUTPUT:
%     swapRate     - Fixed rate of the swap (annualized)

    % Define day count convention as European 30/360
    EU_30_360=6;

    % Calculate year fractions between consecutive payment dates
    yf = yearfrac(paymentDates(1:end-1),paymentDates(2:end), EU_30_360);

    % Interpolate discount factors for payment dates
    discounts_inter = interpolation(discounts,dates, paymentDates);
    % Compute the forward discounts
    fwd_disc = discounts_inter(2:end)./discounts_inter(1:end-1);
    
    % Calculate Basis Point Value (BPV)
    BPV = sum(yf.*fwd_disc);

    % Calculate swap rate as (1 - final discount factor) / BPV   
    swapRate = (1-discounts_inter(end)/discounts_inter(1))/BPV;
end