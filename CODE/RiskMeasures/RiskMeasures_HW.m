function [expected_exposure_nocoll,potential_future_exposure_nocoll, peak_pfe_nocoll, expected_positive_exposure_nocoll, expected_exposure_delta_HW,expected_positive_exposure_delta_HW,potential_future_exposure_delta_HW,peak_pfe_delta_HW, expected_exposure_epsilon_HW,expected_positive_exposure_epsilon_HW,potential_future_exposure_epsilon_HW,peak_pfe_epsilon_HW, EPE_tot, peakPFE_tot] = RiskMeasures_HW(Portfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach, flag)
% RISKMEASURES_HW   Compute credit exposure measures under Hull–White
%
%   [expected_exposure_nocoll, potential_future_exposure_nocoll, peak_pfe_nocoll, ...
%    expected_positive_exposure_nocoll, ...
%    expected_exposure_delta_HW, expected_positive_exposure_delta_HW, ...
%    potential_future_exposure_delta_HW, peak_pfe_delta_HW, ...
%    expected_exposure_epsilon_HW, expected_positive_exposure_epsilon_HW, ...
%    potential_future_exposure_epsilon_HW, peak_pfe_epsilon_HW, ...
%    EPE_tot, peakPFE_tot] = RiskMeasures_HW(...
%        Portfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach, flag)
%
%   Computes expected exposure (EE), expected positive exposure (EPE),
%   potential future exposure (PFE), and peak PFE for a portfolio of
%   interest‐rate swaps under a one‐factor Hull–White short‐rate model,
%   both without collateral and with a simple collateralization rule.
%
%   INPUTS:
%     Portfolio       – Table/struct containing portfolio details:
%                       .FixedLegRate, .Maturity, .FixedLeg (Receive/Pay),
%                       .NotionalAmount, .CounterpartyName
%     dates           – Vector of datetime points for the zero‐coupon curve
%     discounts       – Vector of P(0, dates(k)), the observed discount factors
%     simulations_num – Number of Monte Carlo scenarios
%     alpha           – Confidence level for PFE quantile (e.g., 0.95)
%     aBach           – Mean‐reversion speed κ of the Hull–White model
%     sigmaBach       – Volatility σ of the Hull–White model
%     flag            – Binary flag (if 1, re‐calibrate the fixed‐leg cashflow to par)
%
%   OUTPUTS:
%     expected_exposure_nocoll       – 1×T vector of EE without collateral
%     potential_future_exposure_nocoll – 1×T vector of PFE without collateral
%     peak_pfe_nocoll                – Scalar: maximum PFE over all time points (no collateral)
%     expected_positive_exposure_nocoll – Scalar: average EE (no collateral) excluding final time
%     expected_exposure_delta_HW     – 1×T EE profile for “Delta” counterparty under collateral
%     expected_positive_exposure_delta_HW – Scalar: average EE for Delta (collateral)
%     potential_future_exposure_delta_HW – 1×T PFE profile for Delta (collateral)
%     peak_pfe_delta_HW              – Scalar: peak PFE for Delta (collateral)
%     expected_exposure_epsilon_HW   – 1×T EE profile for “Epsilon” counterparty under collateral
%     expected_positive_exposure_epsilon_HW – Scalar: average EE for Epsilon (collateral)
%     potential_future_exposure_epsilon_HW – 1×T PFE profile for Epsilon (collateral)
%     peak_pfe_epsilon_HW            – Scalar: peak PFE for Epsilon (collateral)
%     EPE_tot                        – Scalar: combined portfolio EPE (collateral)
%     peakPFE_tot                    – Scalar: combined portfolio peak PFE (collateral)
%
%   The function proceeds in three main steps:
%     1. Simulate short‐rate scenarios x(t) at each quarterly time point.
%     2. For each swap in the portfolio, compute Monte Carlo MtM at each time for both
%        “Delta” and “Epsilon” counterparties, then aggregate to obtain EE and PFE
%        without collateral.
%     3. Compute collateralized MtM using a simple rule (variation‐margining every quarter),
%        then compute EE and PFE under collateral for each counterparty and overall.


% 1) Setup and common definitions
rng(42);                              % Fix random seed for reproducibility
today = dates(1);                    % Valuation date
EU_30_360 = 6;                       % Yearfrac code for 30/360 European
ACT_365   = 3;                       % Yearfrac code for ACT/365

% HW mean‐reversion speed
a = aBach;

% Extract portfolio fields into vectors
swapRate = Portfolio.FixedLegRate;
maturity = Portfolio.Maturity; 
fixedLeg = Portfolio.FixedLeg;
notional = 10^6*Portfolio.NotionalAmount;
counterparty = Portfolio.CounterpartyName;

% Convert swapType to +1 for Receive, −1 for Pay
swapType = zeros(size(fixedLeg));             % inizializza a 0
swapType(strcmpi(fixedLeg,'Receive')) =  1;   % imposta +1 per Receive
swapType(strcmpi(fixedLeg,'Pay'))     = -1;   % imposta -1 per Pay

% Build time grid for each quarter up to 10 years
fixed_leg_payments_dates = businessdayoffset(today:calmonths(3):today+calyears(10));
% Interpolate discount factors at those quarterly dates
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

simulations_grid = fixed_leg_payments_dates(2:end);

% Determine how many “Delta” and “Epsilon” counterparty swaps are in the portfolio
nEpsilon = sum(strcmp(counterparty,'Epsilon'));   % conta quante celle sono 'A'
nDelta = sum(strcmp(counterparty,'Delta'));

% Precompute fixed‐leg cash flows per unit notional (constant over time)
fixed_cash_flow = swapRate.*deltas;

% Generate all standard normal shocks for each simulation and quarter
z=randn(simulations_num, length(simulations_grid));

% Preallocate 3D arrays to store per‐scenario, per‐time MtM for each counterparty
r = simulations_num; c = length(simulations_grid);      % dimensione fissa di ciascuna X_i
THW_Epsilon = zeros(r, c, nEpsilon);
THW_Delta = zeros(r, c, nDelta);

% To store simulated x(t) paths for all simulations and times
xts = zeros(simulations_num, 40);

% Simulate short‐rate deviation x(t) under Hull–White
countDelta = 1; 
countEpsilon = 1; 
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    if idx == 1
        prev_sim_date = today;
        x_t           = zeros(1, simulations_num);
    else
        prev_sim_date = simulations_grid(idx-1);
    end
         % One‐step simulate OU‐process x(t) from prev_sim_date to sim_date
         x_t = simulateHW(z(:, idx), aBach, sigmaBach, x_t, prev_sim_date, sim_date);
         xts(:, idx) = x_t;   % Store all paths at this time
end

% 3) Compute MtM and exposure without collateral
expected_exposure_nocoll = 0;
potential_future_exposure_nocoll = 0;

for j=1:length(maturity)
% Preallocate MtM matrices for this swap, per counterparty type
MtMs_Epsilon=zeros(simulations_num, 40);
MtMs_Delta=zeros(simulations_num, 40);
% Reconstruct time grid from “today” to swap maturity j
simulations_grid =   businessdayoffset(today:calmonths(3):today+calmonths(maturity(j)*12));

    if flag == 1
        % Re‐compute par‐swap rate if flag==1 (i.e., calibrate to par on each path)
        discounts_timegrid = interpolation(discounts, dates, simulations_grid); 
        deltas = yearfrac(simulations_grid(1:end-1), simulations_grid(2:end), ACT_365); 
        BPV = sum(discounts_timegrid(2:end).*deltas);
        ParswapRate = (1-discounts_timegrid(end))/BPV;
        fixed_cash_flow = ParswapRate*deltas;
        disp((ParswapRate-swapRate(j))/ParswapRate);
    end
    
simulations_grid = simulations_grid(2:end); % drop “today” from the grid
Nstep           = numel(simulations_grid);
    for idx = 1:length(simulations_grid)
        sim_date = simulations_grid(idx);
        x_t=xts(:, idx)';
    
        % Check if we are before the swap’s final cashflow date
        if (maturity(j)>=1) && sim_date < businessdayoffset(today+ calyears(maturity(j)))
        % Compute A and C for pricing from sim_date to all remaining dates
        [A, C] = affine_trick(sim_date, simulations_grid, aBach, sigmaBach, discounts_timegrid, dates, discounts);
    
                if counterparty(j) == "Epsilon"
                    % Build the MtM for Epsilon at sim_date (per scenario)
                    if flag == 1
                         term1    = exp( log(fixed_cash_flow' .* A) - C * x_t ); 
                    else
                         term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t ); 
                    end
                    %term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t );     %alla funzione bisogna dare la riga giusta a fixed_cash_flow
                    sum_term = sum(term1, 1);           % 1×Nscen
                    terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                    adjust   = (idx < Nstep);   % 1 if we still have future payments
                    MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                    MtMs_Epsilon(:,idx) = MtM;
                elseif counterparty(j) == "Delta"
                    % Build the MtM for Delta at sim_date (per scenario)
                    if flag == 1
                        term1    = exp( log(fixed_cash_flow' .* A) - C * x_t ); 
                    else
                        term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t ); 
                    end
                    sum_term = sum(term1, 1);           % 1×Nscen
                    terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                    adjust   = (idx < Nstep);
                    MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                    MtMs_Delta(:,idx) = MtM;
                end
        elseif sim_date < businessdayoffset(today+ calmonths(maturity(j)*12))
            [A, C] = affine_trick(sim_date, simulations_grid, aBach, sigmaBach, discounts_timegrid, dates, discounts);
    
                if counterparty(j) == "Epsilon"
                    if flag == 1
                        term1    = exp( log(fixed_cash_flow' .* A) - C * x_t ); 
                    else
                         term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t ); 
                    end
                    sum_term = sum(term1, 1);           % 1×Nscen
                    terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                    adjust   = (idx < Nstep);
                    MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                    MtMs_Epsilon(:,idx) = MtM;
                elseif counterparty(j) == "Delta"
                    if flag == 1
                         term1    = exp( log(fixed_cash_flow' .* A) - C * x_t ); 
                    else
                         term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t ); 
                    end
                    sum_term = sum(term1, 1);           % 1×Nscen
                    terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                    adjust   = (idx < Nstep);
                    MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                    MtMs_Delta(:,idx) = MtM;
                end
        end
    end

    % Aggregate MtM for this swap into portfolio exposures
    if counterparty(j) == "Epsilon"
        THW_Epsilon(:,:,countEpsilon) = MtMs_Epsilon;
        countEpsilon = countEpsilon +1; 

        % Compute EE and PFE for Epsilon (no collateral) at each time:
        credit_exposure_nocoll = max(MtMs_Epsilon, 0);
        prob = mean(MtMs_Epsilon > 0) ;
        new=mean(credit_exposure_nocoll)./prob;
        new( isnan(new) ) = 0;
        expected_exposure_nocoll = expected_exposure_nocoll+new;
        potential_future_exposure_nocoll =potential_future_exposure_nocoll + quantile(MtMs_Epsilon, alpha);
    else
        THW_Delta(:,:,countDelta) = MtMs_Delta;
        countDelta=countDelta+1;
        credit_exposure_nocoll = max(MtMs_Delta, 0);
        prob = mean(MtMs_Delta > 0);
        new=mean(credit_exposure_nocoll)./prob;
        new( isnan(new) ) = 0;
        expected_exposure_nocoll = expected_exposure_nocoll+new;
        potential_future_exposure_nocoll =potential_future_exposure_nocoll + quantile(MtMs_Delta, alpha);
    end
end

% Final aggregate no‐collateral measures
peak_pfe_nocoll =max(potential_future_exposure_nocoll);
expected_positive_exposure_nocoll = mean(expected_exposure_nocoll(1:end-1));

% 4) Compute collateralized MtM and exposures
% Sum individual swap MtM to get total MtM for each counterparty (per scenario, per time)
MtM_Epsilon_HW = sum(THW_Epsilon, 3);
MtM_Delta_HW = sum(THW_Delta, 3);

% Initialize collateral trajectories
collateral_Epsilon_HW = zeros(simulations_num, 40);
collateral_Delta_HW = zeros(simulations_num, 40);

MtM_coll_Epsilon_HW = MtM_Epsilon_HW;
MtM_coll_Delta_HW = MtM_Delta_HW;

% Reconstruct discount factors at the quarterly grid for discounting collateral
simulations_grid = fixed_leg_payments_dates(2:end);
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end)./disc_simgrid(1:end-1);

% Counterparty Epsilon
for i=1:39
    % Collateral margin calls happen semiannual
    if mod(i,2) == 0
        stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts, fwd_disc);
        MtM_coll_Epsilon_HW(:, i) = MtM_coll_Epsilon_HW(:, i) + collateral_Epsilon_HW(:, i-1)./stoch_disc;
        collateral_Epsilon_HW(:, i) = - MtM_coll_Epsilon_HW(:, i);
        MtM_coll_Epsilon_HW(:, i) = 0; % After posting collateral, net MtM = 0
    elseif i==1
        % no previous collateral to adjust
        MtM_coll_Epsilon_HW(:, i) = MtM_coll_Epsilon_HW(:, i) ;
    else
        stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts,fwd_disc);
        MtM_coll_Epsilon_HW(:, i) = MtM_coll_Epsilon_HW(:, i) + collateral_Epsilon_HW(:, i-1)./stoch_disc;
    end
end

% Counterparty Delta
for i=1:39
    if mod(i,2) == 0
        stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts,fwd_disc);
        MtM_coll_Delta_HW(:, i) = MtM_coll_Delta_HW(:, i) + collateral_Delta_HW(:, i-1)./stoch_disc;
        collateral_Delta_HW(:, i) = - MtM_coll_Delta_HW(:, i);
        MtM_coll_Delta_HW(:, i) = 0; 
    elseif i==1
        MtM_coll_Delta_HW(:, i) = MtM_coll_Delta_HW(:, i) ;
    else
        stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts,fwd_disc);
        MtM_coll_Delta_HW(:, i) = MtM_coll_Delta_HW(:, i) + collateral_Delta_HW(:, i-1)./stoch_disc;
    end
end

% Compute exposures under collateral for each counterparty
credit_exposure_epsilon_HW = max(MtM_coll_Epsilon_HW, 0); 
prob_epsilon_HW = mean(MtM_coll_Epsilon_HW > 0);
expected_exposure_epsilon_HW = mean(credit_exposure_epsilon_HW)./prob_epsilon_HW;

potential_future_exposure_epsilon_HW = quantile(MtM_coll_Epsilon_HW, alpha);
peak_pfe_epsilon_HW = max(potential_future_exposure_epsilon_HW);
% Replace NaNs (if prob=0) with zero
expected_exposure_epsilon_HW( isnan(expected_exposure_epsilon_HW) ) = 0;
potential_future_exposure_epsilon_HW( isnan(potential_future_exposure_epsilon_HW )) = 0;
expected_positive_exposure_epsilon_HW = mean(expected_exposure_epsilon_HW(1:end-1));

nT = numel(simulations_grid);
vEPE     = expected_positive_exposure_epsilon_HW*ones(1,nT);
vPeakPFE = peak_pfe_epsilon_HW*ones(1,nT);

% Plot exposure profile for Epsilon under HW collateralization
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_epsilon_HW,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_epsilon_HW,'-r','LineWidth',1.5);
plot(simulations_grid, vEPE,      '--g', 'LineWidth',1.2);
plot(simulations_grid, vPeakPFE,  '--k', 'LineWidth',1.2);
hold off;

legend( ...
  'Expected Exposure Epsilon', ...
  'Potential Future Exposure Epsilon', ...
  'EPE', ...
  'Peak PFE', ...
  'Location','best' ...
);

xlabel('Time','FontSize',12);
ylabel('Exposure','FontSize',12);
title('Exposure Profile HW – Counterparty Epsilon','FontSize',14);

% Repeat for Delta under HW collateralization
credit_exposure_delta_HW = max(MtM_coll_Delta_HW, 0); 
prob_delta_HW = mean(MtM_coll_Delta_HW > 0); 
expected_exposure_delta_HW = mean(credit_exposure_delta_HW)./prob_delta_HW;

potential_future_exposure_delta_HW = quantile(MtM_coll_Delta_HW, alpha);
peak_pfe_delta_HW = max(potential_future_exposure_delta_HW);

expected_exposure_delta_HW( isnan(expected_exposure_delta_HW) ) = 0;
potential_future_exposure_delta_HW( isnan(potential_future_exposure_delta_HW) ) = 0;
expected_positive_exposure_delta_HW = mean(expected_exposure_delta_HW(1:end-1));

nT = numel(simulations_grid);
vEPE     = expected_positive_exposure_delta_HW * ones(1,nT);
vPeakPFE = peak_pfe_delta_HW             * ones(1,nT);

figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_delta_HW,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_delta_HW,'-r','LineWidth',1.5);
plot(simulations_grid, vEPE,      '--g', 'LineWidth',1.2);
plot(simulations_grid, vPeakPFE,  '--k', 'LineWidth',1.2);
hold off;

legend( ...
  'Expected Exposure Delta', ...
  'Potential Future Exposure Delta', ...
  'EPE', ...
  'Peak PFE', ...
  'Location','best' ...
);

xlabel('Time','FontSize',12);
ylabel('Exposure','FontSize',12);
title('Exposure Profile HW – Counterparty Delta','FontSize',14);

% 5) Plot no‐collateral profile for overall portfolio
vEPE     = expected_positive_exposure_nocoll * ones(1, nT);
vPeakPFE = peak_pfe_nocoll * ones(1, nT);
 
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_nocoll,           '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_nocoll,   '-r', 'LineWidth',1.5);
plot(simulations_grid, vEPE,                        '--g', 'LineWidth',1.2);
plot(simulations_grid, vPeakPFE,                    '--k', 'LineWidth',1.2);
hold off;
legend( ...
  'Expected Exposure', ...
  'Potential Future Exposure', ...
  'EPE', ...
  'Peak PFE', ...
  'Location','best' ...
);
xlabel('Time','FontSize',12);
ylabel('Exposure','FontSize',12);
title('Exposure Profile – Overall Portfolio','FontSize',14);

% 6) Aggregate collateralized EPE/PFE across counterparties
nT = numel(simulations_grid);
EPE_tot = mean(expected_exposure_epsilon_HW(1:end-1)+expected_exposure_delta_HW(1:end-1));
peakPFE_tot = max(potential_future_exposure_epsilon_HW+potential_future_exposure_delta_HW);
vEPE     = EPE_tot*ones(1,nT);
vPeakPFE = peakPFE_tot*ones(1,nT);
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_epsilon_HW+expected_exposure_delta_HW,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_epsilon_HW+potential_future_exposure_delta_HW,'-r','LineWidth',1.5);
plot(simulations_grid, vEPE,      '--g', 'LineWidth',1.2);
plot(simulations_grid, vPeakPFE,  '--k', 'LineWidth',1.2);
hold off;
legend( ...
  'Expected Exposure Epsilon', ...
  'Potential Future Exposure Epsilon', ...
  'EPE', ...
  'Peak PFE', ...
  'Location','best' ...
);
xlabel('Time','FontSize',12);
ylabel('Exposure','FontSize',12);
title('Exposure Profile – Overall Portfolio','FontSize',14);