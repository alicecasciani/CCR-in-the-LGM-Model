function [price] = JamshidianLGM(dates, discounts, dates_years, coupon, sigma_prev, sigma, a, time_grid)
% JamshidianLGM: Prices a coupon‐bearing bond using Jamshidian’s decomposition
    % under a one‐factor LGM (Linear Gaussian Markov) model.
    %
    % Inputs:
    %   dates         : vector of bootstrap dates
    %   discounts     : vector of known discount factors corresponding to “dates”
    %   dates_years   : vector of dates from t_alpha to t_omega
    %   coupon        : annual coupon rate -> swaption strike
    %   sigma_prev    : vector of previously calibrated σ‐parameters at each earlier time‐step.
    %                   (Length = number of expiries before the current one.)
    %   sigma         : the current σ parameter to calibrate for the latest time‐step.
    %   a             : mean‐reversion speed.
    %   time_grid     : vector of datenum points = { today, each swaption expiry }, 
    %                   used to approximate integrals in the LGM volatility accumulate.
    %
    % Output:
    %   price         : price of the coupon‐bearing bond (clean price) by summing 
    %                   each coupon amount times the price of the corresponding zero‐coupon put.


%Day‐count convention
ACT_360 = 2;   % Actual/360 (used to compute accrual deltas between payment dates)
ACT_365 = 3;   % Actual/365 (used to compute year fractions from today to each date)

% Set the reference date
today = dates(1);
% Set t_alpha (swaption expiry)
t_alpha = dates_years(1);

% Compute the accrual fractions (“deltas”) between successive cashflow dates
% using Actual/360 day count.
deltas = yearfrac(dates_years(1:end-1), dates_years(2:end), ACT_360);

% Compute yearfraction from “today” to each future date in dates_years using ACT/365
yf = yearfrac(today, dates_years, ACT_365);


% Build the vector of cash‐flow amounts (“coupons”) at each payment date:
%   coupon * accrual_delta for each coupon period, and 
%   the final period gets (1 + coupon*delta) to include redemption of principal.
coupons = coupon * deltas;
coupons(end) = 1 + coupons(end); % add principal redemption in final cashflow

% Obtain discount factors at the cashflow times “dates_years”
disc = interpolation(discounts, dates, dates_years);
% Compute the forward‐discount factors for periods 1..end by dividing
% by the “t_alpha” discount (disc(1)).
fwd_disc = disc(2:end)/disc(1);

% In the Gaussian short‐rate model, β(s,T) = (1 − e^(−a·(T−s))) / a
beta = @(s,T) (1 - exp(-a * (T-s)))/a;

% Prepare arrays for stochastic forward‐discount function handles
fwd_stoc_disc = cell(1,length(coupons));

% Convert the “time_grid” dates into year‐fractions from today
time_grid_yf=yearfrac(today,time_grid,ACT_365);
% Also convert “dates_years” into year‐fractions from today
dates_years_yf=yearfrac(today,dates_years,ACT_365);

for i = 1:length(coupons)
    % Compute β at (t_alpha, t_alpha) and at (t_alpha, T_i)
    beta_prev = beta(yf(1),yf(1));
    beta_next = beta(yf(1),yf(i+1));

    % Build an array to accumulate the integral from time_grid points 1..end‐2
    %     int_prev(j) = contribution from subinterval j based on previously calibrated sigma_prev(j)
    int_prev=zeros(1,length(time_grid_yf)-2);
    for j=1:length(time_grid_yf)-2
    int_prev(j)=sigma_prev(j).^2/(2*a^3).*exp(-2*a.*dates_years_yf(i+1)).*(exp(2*a*time_grid_yf(j+1))+exp(2*a.*(-time_grid_yf(j+1)+time_grid_yf(j)+dates_years_yf(i+1))) ...
        -4.*exp(a*(-time_grid_yf(j+1)+time_grid_yf(j)+2*dates_years_yf(i+1)))+4.*exp(a*(time_grid_yf(j)+dates_years_yf(i+1)))-exp(2*a*time_grid_yf(j))-4*exp(a*(time_grid_yf(j+1)+dates_years_yf(i+1))) ...
        +3.*exp(2*a.*dates_years_yf(i+1)));
    end

    % Add the contribution from the final subinterval (using the current σ)
    int=int_prev+sigma.^2/(2*a^3).*exp(-2*a.*dates_years_yf(i+1)).*(exp(2*a*time_grid_yf(end))+exp(2*a.*(-time_grid_yf(end)+time_grid_yf(end-1)+dates_years_yf(i+1))) ...
        -4.*exp(a*(-time_grid_yf(end)+time_grid_yf(end-1)+2*dates_years_yf(i+1)))+4.*exp(a*(time_grid_yf(end-1)+dates_years_yf(i+1)))-exp(2*a*time_grid_yf(end-1))-4*exp(a*(time_grid_yf(end)+dates_years_yf(i+1)))...
        +3.*exp(2*a.*dates_years_yf(i+1)));

    % Construct the stochastic forward‐discount function handle
    fwd_stoc_disc{i} = @(x) fwd_disc(i) * exp(-x* (beta_next-beta_prev)-0.5*sum(int));
end

% We will solve price_CB(x) = 1 for x = x* (Jamshidian’s root).
price_CB = @(x) - 1;
for i = 1:length(coupons)
    price_CB = @(x) price_CB(x) + fwd_stoc_disc{i}(x)*coupons(i);
end
x_star = fzero(price_CB,0);

% Having found x*, we can compute the “strike” zero‐coupon bond prices K_i 
% for each maturity T_i that will replicate each coupon payment via a put.
K_i = zeros(1,length(coupons));
for i = 1:length(coupons)
    beta_prev = beta(yf(1),yf(1));
    beta_next = beta(yf(1),yf(i+1));
    int_prev=zeros(1,length(time_grid_yf)-2);
    for j=1:length(time_grid_yf)-2
    int_prev(j)=sigma_prev(j).^2/(2*a^3).*exp(-2*a.*dates_years_yf(i+1)).*(exp(2*a*time_grid_yf(j+1))+exp(2*a.*(-time_grid_yf(j+1)+time_grid_yf(j)+dates_years_yf(i+1))) ...
        -4.*exp(a*(-time_grid_yf(j+1)+time_grid_yf(j)+2*dates_years_yf(i+1)))+4.*exp(a*(time_grid_yf(j)+dates_years_yf(i+1)))-exp(2*a*time_grid_yf(j))-4*exp(a*(time_grid_yf(j+1)+dates_years_yf(i+1))) ...
        +3.*exp(2*a.*dates_years_yf(i+1)));
    end
    int=int_prev+sigma.^2/(2*a^3).*exp(-2*a.*dates_years_yf(i+1)).*(exp(2*a*time_grid_yf(end))+exp(2*a.*(-time_grid_yf(end)+time_grid_yf(end-1)+dates_years_yf(i+1))) ...
        -4.*exp(a*(-time_grid_yf(end)+time_grid_yf(end-1)+2*dates_years_yf(i+1)))+4.*exp(a*(time_grid_yf(end-1)+dates_years_yf(i+1)))-exp(2*a*time_grid_yf(end-1))-4*exp(a*(time_grid_yf(end)+dates_years_yf(i+1)))...
        +3.*exp(2*a.*dates_years_yf(i+1)));
    % The strike K_i is the forward‐discount factor under shift x*
    K_i(i) = fwd_disc(i) * exp(-x_star* (beta_next-beta_prev)-0.5*sum(int));
end

% Now compute the price of a put on a zero‐coupon bond expiring at T_i with strike K_i:
% pricePutLGM(K, today, t_alpha, T, sigma, a, discounts, dates)
%   – external function that prices the put on a ZCB under LGM.
Puts = zeros(1, length(K_i));
for i = 1:length(K_i)
    K   = K_i(i);
    T   = dates_years(i+1);
    j_alpha = find(time_grid == t_alpha, 1);
    Puts(i) = pricePutLGM(K, today, t_alpha, T, [sigma_prev sigma], a, discounts, dates, time_grid(1:j_alpha));
end

% Finally, the Jamshidian price is:
%   Price = Σ_{i=1}^N [ coupon_i * Put_i ]
price = sum(coupons.*Puts);

end