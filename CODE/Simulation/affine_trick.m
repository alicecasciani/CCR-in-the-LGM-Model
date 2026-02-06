function [A, C] = affine_trick(t, pricing_grid, a, sigma, discount_factors, dates, discounts)
% AFFINETRICK   Precompute A(s, τ) and C(s, τ) for Hull–White ZCB prices
%   [A, C] = affine_trick(t, pricing_grid, a, sigma, discount_factors, dates, discounts)
%
%   Implements the affine decomposition for zero‐coupon bonds under Hull–White:
%     P(t, τ) = A(t, τ) * exp(−C(t, τ) * x_t),
%   where x_t is the OU short‐rate deviation at time t.
%
%   INPUTS:
%     t                – Valuation time s (datetime or datenum)
%     pricing_grid     – Vector of bond maturities τ_j (datetime array)
%     a                – HW mean‐reversion speed (κ)
%     sigma            – HW volatility parameter (σ)
%     discount_factors – P(0, τ_j) for each τ_j in pricing_grid
%     dates            – Dates corresponding to the full discount curve points
%     discounts        – Vector of P(0, dates(k)), the market zero‐coupon curve
%
%   OUTPUTS:
%     A, C             – Vectors A(t, τ_j) and C(t, τ_j) for each maturity τ_j

% Set day-count convention
ACT_365 = 3; 
% Preallocate output vectors
N = numel(pricing_grid);
A = zeros(N,1);
C = zeros(N,1);

% Compute the discount factor at t 
P0_t = interpolation(discounts, dates, t);

% Compute τ_t = yearfrac(today, t), where 'today' = dates(1) of the market data
today = dates(1);
tau_t = yearfrac(today, t, ACT_365);

% Loop over each maturity t_j in the pricing grid
for j = 1:N
    tj = pricing_grid(j);
    if tj > t
        % Compute τ = yearfrac(t, τ_j)
        tau = yearfrac(t, tj, ACT_365);
        % C(j) = (1 − exp(−a * τ)) / a
        C(j) = (1 - exp(-a * tau)) / a;
        % Forward discount factor = P(0, τ_j) / P(0, t)
        fwd_disc = discount_factors(j) / P0_t;

        % Compute the integral of σ^2 differences for A(j)
        %  Let T_tj = yearfrac(today, τ_j). Define integrand:
        %   integrand(u) = [σ_HJM(a, σ, u, T_tj)]^2 − [σ_HJM(a, σ, u, τ_t)]^2
        T_tj = yearfrac(today, tj, ACT_365);
        integrand = @(u) sigmaHJM(a, sigma, u, T_tj).^2 ...
                      - sigmaHJM(a, sigma, u, tau_t).^2;
        sigma_int = quadgk(integrand, 0, tau_t);
        
        % A(j) = fwd_disc * exp(−0.5 * sigma_int)
        A(j) = fwd_disc * exp(-0.5 * sigma_int);
    else
        % For τ_j ≤ t, set A(j) = 0 and C(j) = 0
        A(j) = 0;
        C(j) = 0;
    end
end
end
