function results = paperApproachLGM(dates, discounts, today, S_black, S_bach, prices_black, prices_bach, expiry, tenors, H_fun)
% PAPERAPPROACHLGM   Calibrate LGM model using both Black and Bachelier approaches
%   results = paperApproachLGM(dates, discounts, today, ...
%                              S_black, S_bach, prices_black, prices_bach, ...
%                              expiry, tenors, H_fun)
%
%   This routine computes two separate calibrations:
%     1. Using Black‐implied swaption prices
%     2. Using Bachelier‐implied swaption prices
%   For each, it calls calibrateLGM(...) and then recovers σ_i from ζ_i.
%
%   INPUTS:
%     dates, discounts    – Full discount curve
%     today               – Valuation date (datetime)
%     S_black             – 1×N vector of forward swap rates for Black‐priced swaptions
%     S_bach              – 1×N vector of forward swap rates for Bachelier‐priced swaptions
%     prices_black        – 1×N vector of market swaption prices from Black vols
%     prices_bach         – 1×N vector of market swaption prices from Bachelier vols
%     expiry, tenors      – 1×N swaption expiries (datetime) and tenors (years)
%     H_fun               – Function handle H_fun(tau) returning H(τ)
%
%   OUTPUT:
%     results             – Struct with two fields: results.BLACK and results.BACHELIER
%       Each sub‐struct contains:
%         .zeta       – 1×N calibrated ζ_i
%         .y_star     – 1×N calibrated y*_i
%         .V          – 1×N model prices Vpay_opt
%         .sigma      – 1×N piecewise volatilities σ_i 


models = {'BLACK','BACHELIER'};
results = struct();
for m = 1:numel(models)
    modelType = models{m};
    % Call calibrateLGM to get ζ and y_star
    [zeta, y_star, V] = calibrateLGM(dates, discounts, today, ...
        selectMarketPrices(modelType, prices_black, prices_bach), expiry, tenors, ...
        selectRfix(modelType, S_black, S_bach), H_fun);
    % Recover piecewise σ from ζ
    sigma = computeSigmaFromZeta(zeta, today, expiry);
    % Store into the results struct
    results.(modelType).zeta    = zeta;
    results.(modelType).y_star = y_star;
    results.(modelType).V      = V;
    results.(modelType).sigma  = sigma;
end
end