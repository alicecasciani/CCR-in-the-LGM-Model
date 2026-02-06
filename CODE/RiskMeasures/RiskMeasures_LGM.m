function [expected_exposure_nocoll,potential_future_exposure_nocoll, peak_pfe_nocoll, expected_positive_exposure_nocoll , expected_exposure_delta,expected_positive_exposure_delta,potential_future_exposure_delta,peak_pfe_delta, expected_exposure_epsilon,expected_positive_exposure_epsilon,potential_future_exposure_epsilon,peak_pfe_epsilon, EPE_tot, peakPFE_tot] = RiskMeasures_LGM(Portfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach, flag)
% RISKMEASURES_LGM   Compute credit exposure measures under LGM model
%
%   [expected_exposure_nocoll, potential_future_exposure_nocoll, peak_pfe_nocoll, ...
%    expected_positive_exposure_nocoll, ...
%    expected_exposure_delta, expected_positive_exposure_delta, ...
%    potential_future_exposure_delta, peak_pfe_delta, ...
%    expected_exposure_epsilon, expected_positive_exposure_epsilon, ...
%    potential_future_exposure_epsilon, peak_pfe_epsilon, ...
%    EPE_tot, peakPFE_tot] = RiskMeasures_LGM(...
%        Portfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach, flag)
%
%   Computes EE, EPE, PFE, and peak PFE for a portfolio of interest‐rate swaps under
%   a Linear Gaussian Markov (LGM) model, both without collateral and with collateralization.
%   Uses piecewise‐constant volatilities calibrated via Jamshidian’s trick.
%
%   INPUTS:
%     Portfolio       – Table/struct with .FixedLegRate, .Maturity, .FixedLeg, 
%                       .NotionalAmount, .CounterpartyName
%     jamshResults    – Struct containing calibrated LGM parameters from Jamshidian:
%                       .BACHELIER.sigma (1×n piecewise vol vector)
%     dates           – Vector of datetime points for zero‐coupon curve
%     discounts       – Vector of P(0, dates(k)), observed discount factors
%     simulations_num – Number of Monte Carlo scenarios
%     alpha           – Confidence level for PFE (e.g., 0.95)
%     aLGM            – Mean‐reversion speed κ in LGM
%     sigmaBach       – Fallback volatility for last interval (if needed)
%     flag            – If 1, re‐calibrate fixed‐leg to par at each time
%
%   OUTPUTS:
%     (Same structure as RiskMeasures_HW, but under the LGM process)


rng(42);  % Seed for reproducibility
today = dates(1); % Reference date
% Day-count convention
EU_30_360 = 6;
ACT_365 = 3;
% LGM mean reversion
a = aLGM;  

% Build the time grid of piecewise volatility knot points
time_grid = [today, today + calmonths(3)];
time_grid = businessdayoffset([time_grid, today + calyears(1):calyears(1):today + calyears(10)]);
date_sigma = time_grid;  % Each interval [date_sigma(i), date_sigma(i+1)] has constant σ

% Extract portfolio details
swapRate = Portfolio.FixedLegRate;
maturity = Portfolio.Maturity;
fixedLeg = Portfolio.FixedLeg;
notional = 10^6 * Portfolio.NotionalAmount;
counterparty = Portfolio.CounterpartyName;

% Swap type: +1 Receive, -1 Pay
swapType = zeros(size(fixedLeg));
swapType(strcmpi(fixedLeg,'Receive')) =  1;
swapType(strcmpi(fixedLeg,'Pay'))     = -1;

% Piecewise σ_i from Jamshidian‐calibration under Bachelier
sigmaLGM_Bachelier = jamshResults.BACHELIER;

% Build quarterly payment grid up to 10 years, and precompute discount factors
fixed_leg_payments_dates = businessdayoffset(today:calmonths(3):today + calyears(10));
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

simulations_grid = fixed_leg_payments_dates(2:end);
Nstep = length(simulations_grid);
nEpsilon = sum(strcmp(counterparty,'Epsilon'));
nDelta = sum(strcmp(counterparty,'Delta'));

fixed_cash_flow = swapRate .* deltas;

% Standard normal shocks for LGM
z = randn(simulations_num, Nstep); 

% Preallocate arrays for MtM trajectories
r = simulations_num; c = Nstep;
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);
xtsLGM = zeros(simulations_num, Nstep);

countDelta = 1;
countEpsilon = 1;
simulations_grid = businessdayoffset(today:calmonths(3):today+calmonths(10*12));
simulations_grid = simulations_grid(2:end);

% Simulate LGM factor x(t) stepwise with piecewise σ
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    if idx == 1
        prev_sim_date = today;
        x_t = zeros(1, simulations_num);
    else
        prev_sim_date = simulations_grid(idx-1);
    end

    % Determine which σ to use on this interval
    if idx == 1
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(1), x_t, prev_sim_date, sim_date);
    elseif idx == 2 || idx == 3
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(2), x_t, prev_sim_date, sim_date);
    elseif idx < 40
        indice = floor(idx/4) + 2;
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(indice), x_t, prev_sim_date, sim_date);
    elseif idx == 40
        x_t = simulateLGM(z(:, idx), aLGM, sigmaBach, x_t, prev_sim_date, sim_date);
    end

    xtsLGM(:, idx) = x_t;
end

% Initialize no‐collateral aggregates
expected_exposure_nocoll = 0;
potential_future_exposure_nocoll = 0;

% 3) Compute MtM for each swap
for j = 1:length(maturity)
    MtMs_Epsilon = zeros(simulations_num, 40);
    MtMs_Delta = zeros(simulations_num, 40);
    simulations_grid = businessdayoffset(today:calmonths(3):today + calmonths(maturity(j)*12));
    
    if flag == 1
        discounts_timegrid = interpolation(discounts, dates, simulations_grid);
        deltas = yearfrac(simulations_grid(1:end-1), simulations_grid(2:end), ACT_365);
        BPV = sum(discounts_timegrid(2:end) .* deltas);
        ParswapRate = (1 - discounts_timegrid(end)) / BPV;
        fixed_cash_flow = ParswapRate * deltas;
    end

    simulations_grid = simulations_grid(2:end);
    Nstep = numel(simulations_grid);

    for idx = 1:length(simulations_grid)
        sim_date = simulations_grid(idx);
        x_t = xtsLGM(:, idx)';

        if (maturity(j) >= 1 && sim_date < businessdayoffset(today + calyears(maturity(j)))) || sim_date < businessdayoffset(today + calmonths(maturity(j)*12))
            
            [A, C] = affine_trick_LGM(sim_date, simulations_grid, aLGM, sigmaLGM_Bachelier, discounts_timegrid, date_sigma, dates, discounts);

            if counterparty(j) == "Epsilon"
                if flag == 1
                    term1 = exp(log(fixed_cash_flow' .* A) - C * x_t);
                else
                    term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                end
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Epsilon(:, idx) = MtM;
            elseif counterparty(j) == "Delta"
                if flag == 1
                    term1 = exp(log(fixed_cash_flow' .* A) - C * x_t);
                else
                    term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
                end
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
            end
        end
    end

    if counterparty(j) == "Epsilon"
        TLGM_Epsilon(:, :, countEpsilon) = MtMs_Epsilon;
        countEpsilon = countEpsilon + 1;  
        credit_exposure_nocoll = max(MtMs_Epsilon, 0);
        prob = mean(MtMs_Epsilon > 0) ;
        new=mean(credit_exposure_nocoll)./prob;
        new( isnan(new) ) = 0;
        expected_exposure_nocoll = expected_exposure_nocoll+new;
        potential_future_exposure_nocoll =potential_future_exposure_nocoll + quantile(MtMs_Epsilon, alpha);
        
    else
        TLGM_Delta(:, :, countDelta) = MtMs_Delta;
        countDelta = countDelta + 1;
        credit_exposure_nocoll = max(MtMs_Delta, 0);
        prob = mean(MtMs_Delta > 0);
        new=mean(credit_exposure_nocoll)./prob;
        new( isnan(new) ) = 0;
        expected_exposure_nocoll = expected_exposure_nocoll+new;
        potential_future_exposure_nocoll =potential_future_exposure_nocoll + quantile(MtMs_Delta, alpha);
        
    end
end
peak_pfe_nocoll =max(potential_future_exposure_nocoll);
expected_positive_exposure_nocoll = mean(expected_exposure_nocoll(1:end-1));

% 4) Compute collateralized MtM and exposures under LGM
% Aggregate MtM per counterparty across all swaps
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

collateral_Epsilon = zeros(simulations_num, 40);
collateral_Delta = zeros(simulations_num, 40);
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;
simulations_grid = fixed_leg_payments_dates(2:end);
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end) ./ disc_simgrid(1:end-1);
indice3m = 1;
indice1y = 4;

% Epsilon
for i = 1:39
    if mod(i,2) == 0
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1) ./ stoch_disc;
        collateral_Epsilon(:, i) = -MtM_coll_Epsilon(:, i);
        MtM_coll_Epsilon(:, i) = 0;
    elseif i == 1
        % No adjustment at first step
    else
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1) ./ stoch_disc;
    end
end

% Delta
for i = 1:39
    if mod(i,2) == 0
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1) ./ stoch_disc;
        collateral_Delta(:, i) = -MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0;
    elseif i == 1
        % No adjustment at first step
    else
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1) ./ stoch_disc;
    end
end

% Compute exposures under collateral for Epsilon
credit_exposure_epsilon = max(MtM_coll_Epsilon, 0); 
prob_epsilon = mean(MtM_coll_Epsilon > 0);
expected_exposure_epsilon = mean(credit_exposure_epsilon)./prob_epsilon;
potential_future_exposure_epsilon = quantile(MtM_coll_Epsilon, alpha);
peak_pfe_epsilon = max(potential_future_exposure_epsilon);

expected_exposure_epsilon( isnan(expected_exposure_epsilon) ) = 0;
potential_future_exposure_epsilon( isnan(potential_future_exposure_epsilon )) = 0;
expected_positive_exposure_epsilon = mean(expected_exposure_epsilon(1:end-1));

nT = numel(simulations_grid);
vEPE     = expected_positive_exposure_epsilon*ones(1,nT);
vPeakPFE = peak_pfe_epsilon*ones(1,nT);

% Plot collateralized profile for Epsilon
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_epsilon,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_epsilon,'-r','LineWidth',1.5);
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
title('Exposure Profile – Counterparty Epsilon','FontSize',14);

% Repeat for Delta under collateral
credit_exposure_delta = max(MtM_coll_Delta, 0); 
prob_delta = mean(MtM_coll_Delta > 0); 
expected_exposure_delta = mean(credit_exposure_delta)./prob_delta;
potential_future_exposure_delta = quantile(MtM_coll_Delta, alpha);
peak_pfe_delta = max(potential_future_exposure_delta);

expected_exposure_delta( isnan(expected_exposure_delta) ) = 0;
potential_future_exposure_delta( isnan(potential_future_exposure_delta) ) = 0;
expected_positive_exposure_delta = mean(expected_exposure_delta(1:end-1));

nT = numel(simulations_grid);
vEPE     = expected_positive_exposure_delta * ones(1,nT);
vPeakPFE = peak_pfe_delta              * ones(1,nT);

% Plot collateralized profile for Epsilon
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_delta,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_delta,'-r','LineWidth',1.5);
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
title('Exposure Profile – Counterparty Delta','FontSize',14);

% plot CRR no collateral
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

% Aggregate collateralized exposures for entire portfolio
nT = numel(simulations_grid);
EPE_tot = mean(expected_exposure_epsilon(1:end-1)+expected_exposure_delta(1:end-1));
peakPFE_tot = max(potential_future_exposure_epsilon+potential_future_exposure_delta);
vEPE     = EPE_tot*ones(1,nT);
vPeakPFE = peakPFE_tot*ones(1,nT);
figure;
hold on;
grid on;
plot(simulations_grid, expected_exposure_epsilon+expected_exposure_delta,      '-b', 'LineWidth',1.5);
plot(simulations_grid, potential_future_exposure_epsilon+potential_future_exposure_delta,'-r','LineWidth',1.5);
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