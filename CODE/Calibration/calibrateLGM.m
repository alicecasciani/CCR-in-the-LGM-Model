function [zeta_all, y_star_all, V_all] = calibrateLGM(dates, discounts, today, mkt_prices, expiry, tenors, R_fix, H_fun)
% CALIBRATELGM   Calibrate the Linear Gaussian Markov (LGM) model parameters
%   [zeta_all, y_star_all, V_all] = calibrateLGM(dates, discounts, today, ...
%                                               mkt_prices, expiry, tenors, ...
%                                               R_fix, H_fun)
%
%   INPUTS:
%     dates       – Vector of datetime points for the zero‐coupon curve
%     discounts   – Vector of P(0, dates(k)), the observed discount factors
%     today       – Valuation date as a datetime scalar
%     mkt_prices  – 1×N vector of market‐observed swaption prices
%     expiry      – 1×N vector of swaption expiry dates (datetime)
%     tenors      – 1×N vector of swaption tenors in years (integers)
%     R_fix       – 1×N vector of fixed swap rates corresponding to each swaption
%     H_fun       – Function handle H_fun(tau) that returns H(tau) for each tau
%
%   OUTPUTS:
%     zeta_all    – 1×N vector of calibrated variance parameters ζ_e for each swaption
%     y_star_all  – 1×N vector of calibrated “threshold” y* for each swaption
%     V_all       – 1×N vector of optimal swaption prices computed with calibrated params


N = numel(mkt_prices);
% Preallocate outputs
y_star_all = zeros(1,N);
zeta_all    = zeros(1,N);
V_all       = zeros(1,N);

% Day‐count conventions
ACT_365 = 3; 
EU_30_360 = 6;
 
% To enforce monotonic (increasing) ζ_i across swaptions
z_prev = 0;
for i = 1:N
    % --- 1) Build the underlying swap’s payment schedule
    end_date      = businessdayoffset(expiry(i) + calyears(tenors(i)));
    dates_yearly  = businessdayoffset(expiry(i):calyears(1):end_date);

    % --- 2) Interpolate discount factors along those payment dates
    D        = interpolation(discounts, dates, dates_yearly);
    D0       = D(1);
    D        = D(2:end);

    % --- 3) Build the H vector at each payment date
    tau      = yearfrac(today, dates_yearly, ACT_365);
    H_all    = H_fun(tau);
    H0       = H_all(1);
    H        = H_all(2:end);

    % --- 4) Compute accrual fractions α_i for each annual period (30/360 EU)
    alpha    = yearfrac(dates_yearly(1:end-1), dates_yearly(2:end), EU_30_360); 

    % --- 5) Solve for (y_star, ζ_e) that satisfy both valuation equations
    sol = solve_swaption(alpha, R_fix(i), zeros(1,i), D, D0, H, H0, mkt_prices(i));
    y_star = sol(1);
    zeta   = sol(2);

    % --- 6) Enforce monotonicity of ζ across swaptions: ζ_i ≥ ζ_{i−1}
    zeta = max(zeta, z_prev);
    z_prev = zeta;

    % --- 7) Compute the swaption’s “optimal” model price with the calibrated parameters
    V_all(i) = Vpay_opt(y_star, zeta, alpha, R_fix(i), [], D, D0, H, H0);

    % --- 8) Store the calibrated values
    y_star_all(i) = y_star;
    zeta_all(i)    = zeta;
end
end
