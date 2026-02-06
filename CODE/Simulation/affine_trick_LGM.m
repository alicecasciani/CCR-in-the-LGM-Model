function [A, C] = affine_trick_LGM(t, pricing_grid, a, sigma, discount_factors, date_sigma, dates, discounts)
% AFFINETRICK_LGM    Precompute the functions A(t, t_j) and C(t, t_j) for
% Linear Gaussian Markov (LGM) model
%   t               : datetime scalar, valuation date (equals pricing_grid(1))
%   pricing_grid    : datetime vector [t = t_1, t_2, …, t_N] of maturity dates
%   a               : mean reversion speed of the LGM model
%   sigma           : vector of piecewise-constant volatilities (σ_i on each interval)
%   discount_factors: N×1 vector of zero-coupon bond prices P(0, t_j), aligned with pricing_grid
%   date_sigma      : vector of dates at which σ changes (boundaries for piecewise-constant σ)
%   dates           : vector of dates corresponding to 'discounts' (market curve)
%   discounts       : vector of discount factors P(0, dates)
%
% Outputs:
%   A : N×1 vector such that for each j,
%       B(t, t_j) = A(j) * exp(-C(j) * x(t))
%   C : N×1 vector of loadings on the state x(t) for each maturity t_j


% Set day-count convention
ACT_365 = 3; 
N = numel(pricing_grid);
% Preallocate output vectors
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
    if tj >t
        % 1) C_j = (1 - exp(-a*(t_j - t))) / a
        tau = yearfrac(t, tj, ACT_365);
        C(j) = (1 - exp(-a * tau)) / a;

        % 2) forward discount P(0,tj)/P(0,t)
        fwd_disc = discount_factors(j) / P0_t;

        % 3) Compute the integral σ_int needed for A(j)
        T_tj = yearfrac(today, tj, ACT_365);
        integrand = @(u) sigmaLGM(a, u, T_tj).^2 ...
                      - sigmaLGM(a, u, tau_t).^2;
        sigma_int=0;
        % Integrate piecewise over intervals where σ is constant, as defined by date_sigma
        for i=1:length(sigma)
                % If the valuation date t exceeds the upper boundary date_sigma(i+1),
                % integrate from date_sigma(i) to date_sigma(i+1)
                if t>date_sigma(i+1)
                    sigma_int = sigma_int + sigma(i)^2*quadgk(integrand, yearfrac(date_sigma(1), date_sigma(i),ACT_365), yearfrac(date_sigma(1),date_sigma(i+1),ACT_365));
                else
                % Otherwise, integrate only up to t, then break out of the loop
                    sigma_int = sigma_int + sigma(i)^2*quadgk(integrand, yearfrac(date_sigma(1), date_sigma(i),ACT_365), yearfrac(date_sigma(1),t,ACT_365));
                    break
                end

        end

        % 4) A_j
        A(j) = fwd_disc * exp(-0.5 * sigma_int);
    else
        % For t_j == t (j = 1) or any t_j < t, set A(j) = 0 and C(j) = 0
        A(j) = 0;
        C(j) = 0;
    end
end
end
