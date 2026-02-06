function [collateral_Epsilon, collateral_Delta, expected_exposure_delta,expected_positive_exposure_delta,potential_future_exposure_delta,peak_pfe_delta, expected_exposure_epsilon,expected_positive_exposure_epsilon,potential_future_exposure_epsilon,peak_pfe_epsilon] = Daily_LGM_Collateral(Portfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach)
% Set random seed for reproducibility
rng(42);
% Reference date
today=dates(1);
% Day-count conventions
EU_30_360=6;
ACT_365=3;
% Mean-reversion speed for LGM
a=aLGM;

% Extract swap parameters from the Portfolio
swapRate = Portfolio.FixedLegRate;
maturity = Portfolio.Maturity; 
fixedLeg = Portfolio.FixedLeg;
notional = 10^6*Portfolio.NotionalAmount; 
counterparty = Portfolio.CounterpartyName;

% Assign +1 for "Receive", -1 for "Pay"
swapType = zeros(size(fixedLeg));             % inizializza a 0
swapType(strcmpi(fixedLeg,'Receive')) =  1;   % imposta +1 per Receive
swapType(strcmpi(fixedLeg,'Pay'))     = -1;   % imposta -1 per Pay

% Build the piecewise‐constant volatility grid:
time_grid = [today, today + calmonths(3)];
time_grid = businessdayoffset([time_grid, today + calyears(1):calyears(1):today + calyears(10)]);
date_sigma=time_grid;

% Volatility vector from Jamshidian calibration under Bachelier
sigmaLGM_Bachelier = jamshResults.BACHELIER;

% Build a daily timeline from today up to 10 years ahead
fixed_leg_payments_dates = today:caldays(1):today+calyears(10);
% Interpolate discount factors at each day in that timeline
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);
% Compute year‐fractions (30/360 EUR) between consecutive days
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

% Simulation grid excludes t = 0 (days 2..end)
simulations_grid = fixed_leg_payments_dates(2:end);
N=length(simulations_grid);

% Count how many swaps belong to each counterparty ("Epsilon" or "Delta")
nEpsilon = sum(strcmp(counterparty,'Epsilon'));   % conta quante celle sono 'A'
nDelta = sum(strcmp(counterparty,'Delta'));

% Fixed cash‐flow per day per swap
fixed_cash_flow = swapRate.*deltas;

% Generate standard normal shocks for each scenario and each day
z=randn(simulations_num, N);

% Preallocate 3D arrays to store per‐counterparty MtM trajectories:
%  - TLGM_Epsilon: (scenarios × days × #Epsilon swaps)
%  - TLGM_Delta:   (scenarios × days × #Delta swaps)
r = simulations_num; c = N;      
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);
% Preallocate array for the simulated LGM factor x(t)
xtsLGM = zeros(simulations_num, N);

countDelta = 1; 
countEpsilon = 1; 
% 1) Simulate the LGM factor x(t) at daily frequency
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    if idx == 1
        prev_sim_date = today;
        x_t           = zeros(1, simulations_num);
    else
        prev_sim_date = simulations_grid(idx-1);
    end

    % Select the appropriate piecewise sigma for this interval:
    if idx <= 90
         x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(1), x_t, prev_sim_date, sim_date);
    elseif idx > 90 && idx <= 365
           x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(2), x_t, prev_sim_date, sim_date);
    elseif idx<N
        indice = floor(idx/365) + 2;
        if indice>11
            indice=11;
        end
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(indice), x_t, prev_sim_date, sim_date);
    elseif idx == N
        x_t = simulateLGM(z(:, idx), aLGM, sigmaBach, x_t, prev_sim_date, sim_date);
    end
    % Store x(t) for all scenarios
    xtsLGM(:, idx) = x_t;
end

% 2) For each swap, compute uncollateralized MtM trajectories
for j=1:length(maturity)
% Preallocate MtM arrays for this swap
MtMs_Epsilon=zeros(simulations_num, N);
MtMs_Delta=zeros(simulations_num, N);
% Build this swap’s own daily timeline up to its maturity
    if maturity(j) >= 1
        simulations_grid =   today:caldays(1):today+calyears(maturity(j));
    else 
        simulations_grid =  today:caldays(1):today+calmonths(maturity(j)*12);
    end
simulations_grid = simulations_grid(2:end);
Nstep           = numel(simulations_grid);
sim_date = simulations_grid(1);
% Compute affine coefficients A(t), C(t)
[A, C] = affine_trick_LGM(sim_date, simulations_grid, aLGM, sigmaLGM_Bachelier, discounts_timegrid, date_sigma, dates, discounts);

% Loop over each day in this swap’s timeline
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    x_t=xtsLGM(:, idx)';
    % Determine if this day is before the swap’s maturity payment
    if (maturity(j)>=1) && sim_date < businessdayoffset(today+ calyears(maturity(j)))
            % Store MtM for the correct counterparty
            if counterparty(j) == "Epsilon"
                term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t );     %alla funzione bisogna dare la riga giusta a fixed_cash_flow
                sum_term = sum(term1, 1);           % 1×Nscen
                terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                adjust   = (idx < Nstep);
                MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                MtMs_Epsilon(:,idx) = MtM;
                A=[0;A(1:end-1)];
                C=[0;C(1:end-1)];
            elseif counterparty(j) == "Delta"
                term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t );     %alla funzione bisogna dare la riga giusta a fixed_cash_flow
                sum_term = sum(term1, 1);           % 1×Nscen
                terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                adjust   = (idx < Nstep);
                MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                MtMs_Delta(:,idx) = MtM;
                A=[0;A(1:end-1)];
                C=[0;C(1:end-1)];
            end
    elseif sim_date < businessdayoffset(today+ calmonths(maturity(j)*12))
            if counterparty(j) == "Epsilon"
                term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t );     %alla funzione bisogna dare la riga giusta a fixed_cash_flow
                sum_term = sum(term1, 1);           % 1×Nscen
                terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                adjust   = (idx < Nstep);
                MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                MtMs_Epsilon(:,idx) = MtM;
                A=[0;A(1:end-1)];
                C=[0;C(1:end-1)];
            elseif counterparty(j) == "Delta"
                 term1    = exp( log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t );     %alla funzione bisogna dare la riga giusta a fixed_cash_flow
                sum_term = sum(term1, 1);           % 1×Nscen
                terminal = A(end) * exp(-C(end) .* x_t);  % 1×Nscen
                adjust   = (idx < Nstep);
                MtM =  notional(j)*swapType(j) * ( sum_term + terminal - adjust );
                MtMs_Delta(:,idx) = MtM;
                A=[0;A(1:end-1)];
                C=[0;C(1:end-1)];
            end
    end
end
% After looping all dates for swap j, assign to TLGM arrays
if counterparty(j) == "Epsilon"
    TLGM_Epsilon(:,:,countEpsilon) = MtMs_Epsilon;
    countEpsilon = countEpsilon +1; 
else  % "Delta"
    TLGM_Delta(:,:,countDelta) = MtMs_Delta;
    countDelta=countDelta+1;
end
end

% Sum across all swaps for each counterparty to get portfolio-level MtM
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

% 3) Initialize collateral arrays
collateral_Epsilon = zeros(simulations_num, N);
collateral_Delta = zeros(simulations_num, N);
% start with uncollateralized MtM
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;
% Recompute the simulation grid and forward discount ratios for daily steps
simulations_grid = fixed_leg_payments_dates(2:end);
disc_simgrid = interpolation(discounts, dates, simulations_grid); 
fwd_disc = disc_simgrid(2:end)./disc_simgrid(1:end-1);
% Indices for 3-month and 1-year intervals (90 and 365 days)
indice3m=90;
indice1y=365;

% 4) Compute collateralized MtM for Epsilon (weekly)
for i=1:N-1
    if mod(i,7) == 0
        % Every 7th day, post collateral:
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y);
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1)./stoch_disc;
        % Post new collateral equal to negative MtM
        collateral_Epsilon(:, i) = - MtM_coll_Epsilon(:, i);
        % Reset MtM to zero after posting collateral
        MtM_coll_Epsilon(:, i) = 0; 
    elseif i==1
         % On the first step, no prior collateral to adjust
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) ;
    else
        % On non‐collateral days, accrue previous collateral at stochastic discount
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y);
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1)./stoch_disc;
    end
end

% 5) Compute collateralized MtM for Delta (weekly)
for i=1:N-1
    if mod(i,7) == 0
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1)./stoch_disc;
        collateral_Delta(:, i) = - MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0; 
    elseif i==1
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) ;
    else
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1)./stoch_disc;
    end
end

% 6) Compute exposure metrics for Epsilon
credit_exposure_epsilon = max(MtM_coll_Epsilon, 0); 
prob_epsilon = mean(MtM_coll_Epsilon > 0);
expected_exposure_epsilon = mean(credit_exposure_epsilon)./prob_epsilon;
potential_future_exposure_epsilon = quantile(MtM_coll_Epsilon, alpha);
peak_pfe_epsilon = max(potential_future_exposure_epsilon);
% Replace any NaN (due to zero probability) with 0
expected_exposure_epsilon( isnan(expected_exposure_epsilon) ) = 0;
potential_future_exposure_epsilon( isnan(potential_future_exposure_epsilon )) = 0;
expected_positive_exposure_epsilon = mean(expected_exposure_epsilon(1:end-1));

% 7) Compute exposure metrics for Delta
credit_exposure_delta = max(MtM_coll_Delta, 0); 
prob_delta = mean(MtM_coll_Delta > 0);
expected_exposure_delta = mean(credit_exposure_delta)./prob_delta;
potential_future_exposure_delta = quantile(MtM_coll_Delta, alpha);
peak_pfe_delta = max(potential_future_exposure_delta);
% Replace NaNs with 0
expected_exposure_delta( isnan(expected_exposure_delta) ) = 0;
potential_future_exposure_delta( isnan(potential_future_exposure_delta )) = 0;
expected_positive_exposure_delta = mean(expected_exposure_delta(1:end-1));

% Return:
%   collateral_Epsilon, collateral_Delta,
%   expected_exposure_delta, expected_positive_exposure_delta, potential_future_exposure_delta, peak_pfe_delta,
%   expected_exposure_epsilon, expected_positive_exposure_epsilon, potential_future_exposure_epsilon, peak_pfe_epsilon
end