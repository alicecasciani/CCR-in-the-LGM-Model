% FINANCIAL ENGINEERING - AY 2024/2025
% FINAL PROJECT: Project 1B 
% AUTHORS: Leonardo Bacchi - Alice Sofia Casciani - Salvatore Molino
clear all
close all
clc

tic
addpath("./Bootstrap");
addpath("./Calibration");
addpath("./Jamshidian");
addpath("./Pricing");
addpath("./Simulation");
addpath("./Utilities");
addpath("./RiskMeasures");
formatData='dd/mm/yyyy'; 
filename = "mktData.xlsx";
[datesSet, ratesSet] = readExcelData( filename, formatData);

%% BOOTSTRAP 
[dates, discounts, zRates] = bootstrap(datesSet, ratesSet);

%% POINT 1 - Hull-White Calibration

today = dates(1);
expiry = businessdayoffset([today + calmonths(3), today + calyears(1:10)]);
tenors = [1,9,8,7,6,5,4,3,2,1,2];

mkt_vols = [184.95,109.06,103.63,98.20,94.73,92.11,88.76,87.10,86.09,84.18,81.96];

% Calibrazione con Black
[aBlack, sigmaBlack, sseBlack, S_black, mkt_prices_black] = calibrateHWmodel(...
    dates, discounts, today, mkt_vols, expiry, tenors, 'BLACK');

% Calibrazione con Bachelier
[aBach, sigmaBach, sseBach, S_bach, mkt_prices_bach] = calibrateHWmodel(...
    dates, discounts, today, mkt_vols, expiry, tenors, 'BACHELIER');


%% POINT 3 - PAPER APPROACH - LGM Calibration
aLGM       = 0.1;
H_fun      = @(t) (1 - exp(-aLGM .* t)) ./ aLGM;
lgmResults = paperApproachLGM(dates, discounts, today, S_black, S_bach, ...
    mkt_prices_black, mkt_prices_bach, expiry, tenors, H_fun);

disp('Sigma calibrated with Black - paper approach:');
disp(lgmResults.BLACK.sigma);

disp('Sigma calibrated with Bachelier - paper approach:');
disp(lgmResults.BACHELIER.sigma);

%% POINT 3 - JAMSHIDIAN APPROACH - LGM Calibration
aLGM = 0.1;
models = {'BLACK', 'BACHELIER'};
jamshResults = struct();
time_grid = [today, today + calmonths(3)];
time_grid = businessdayoffset([time_grid, today + calyears(1):calyears(1):today + calyears(10)]);

for m = 1:numel(models)
    modelType = models{m};
    prices = selectMarketPrices(modelType, mkt_prices_black, mkt_prices_bach);
    jamshResults.(modelType) = calibrateJamshidian(prices, expiry, tenors, dates, discounts, aLGM, time_grid);
    fprintf('Sigma calibrated with %s - Jamshidian approach:\n', modelType);
    disp(jamshResults.(modelType));
end

%% Plot sigma piece-wise constant 
edges  = [0 0.25 1 2 3 4 5 6 7 8 9 10];
sBlack     = jamshResults.BLACK;      
sBachelier = jamshResults.BACHELIER;  
xBlack = [];  yBlack = [];
xBach  = [];  yBach  = [];
for k = 1:numel(sBlack)
    xBlack = [xBlack  edges(k)  edges(k+1)  NaN];       
    yBlack = [yBlack  sBlack(k) sBlack(k) NaN];   
    xBach  = [xBach   edges(k)  edges(k+1)  NaN];
    yBach  = [yBach   sBachelier(k) sBachelier(k) NaN];
end

figure;  hold on;  grid on;  box on;
plot(xBlack, yBlack,  'LineWidth',2, 'DisplayName', 'Black');
plot(xBach,  yBach,   'LineWidth',2, 'DisplayName', 'Bachelier');
plot(xBlack, sigmaBlack*ones(1, length(xBlack)), 'LineWidth',2, 'LineStyle','--', 'DisplayName', 'HW Black')
plot(xBach, sigmaBach*ones(1, length(xBach)), 'LineWidth',2, 'LineStyle','--', 'DisplayName', 'HW Bachelier')
xticks([0.25 1:10]);
xlabel('T (years)');
ylabel('\sigma');
title('\sigma(t) piecewise constant – Jamshidian approach');
legend('Location','best');

%% POINT 4 
filename = 'portfolio.xlsx';
[alphaPortfolio, betaPortfolio, gammaPortfolio] = readIRSportfolios(filename);
simulations_num = 250000; 
alpha = 0.975; 

%% LGM - alpha ptf
% flag = 1 -> PAR swap rate 
% flag = 0 -> Excel swap rate

tic
[EE_nocoll_alpha_LGM,PFE_nocoll_alpha_LGM, peak_pfe_nocoll_alpha_LGM, EPE_nocoll_alpha_LGM, expected_exposure_delta_alpha_LGM_coll, expected_positive_exposure_delta_alpha_LGM_coll, potential_future_exposure_delta_alpha_LGM_coll, peak_pfe_delta_alpha_LGM_coll, ...
 expected_exposure_epsilon_alpha_LGM_coll, expected_positive_exposure_epsilon_alpha_LGM_coll, potential_future_exposure_epsilon_alpha_LGM_coll, peak_pfe_epsilon_alpha_LGM_coll,  EPE_LGM_alpha, peakPFE_LGM_alpha] = ...
    RiskMeasures_LGM(alphaPortfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach, 0);
toc

%% LGM - beta ptf
tic
[EE_nocoll_beta_LGM,PFE_nocoll_beta_LGM, peak_pfe_nocoll_beta_LGM, EPE_nocoll_beta_LGM,expected_exposure_delta_beta_LGM_coll, expected_positive_exposure_delta_beta_LGM_coll, potential_future_exposure_delta_beta_LGM_coll, peak_pfe_delta_beta_LGM_coll, ...
 expected_exposure_epsilon_beta_LGM_coll, expected_positive_exposure_epsilon_beta_LGM_coll, potential_future_exposure_epsilon_beta_LGM_coll, peak_pfe_epsilon_beta_LGM_coll,  EPE_LGM_beta, peakPFE_LGM_beta] = ...
    RiskMeasures_LGM(betaPortfolio,  jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach, 0);
toc

%% LGM - gamma ptf
tic
[EE_nocoll_gamma_LGM,PFE_nocoll_gamma_LGM, peak_pfe_nocoll_gamma_LGM, EPE_nocoll_gamma_LGM,expected_exposure_delta_gamma_LGM_coll, expected_positive_exposure_delta_gamma_LGM_coll, potential_future_exposure_delta_gamma_LGM_coll, peak_pfe_delta_gamma_LGM_coll, ...
 expected_exposure_epsilon_gamma_LGM_coll, expected_positive_exposure_epsilon_gamma_LGM_coll, potential_future_exposure_epsilon_gamma_LGM_coll, peak_pfe_epsilon_gamma_LGM_coll,  EPE_LGM_gamma, peakPFE_LGM_gamma] = ...
    RiskMeasures_LGM(gammaPortfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach, 0);
toc

%% HW - alpha ptf
% flag = 1 -> PAR swap rate 
% flag = 0 -> Excel swap rate
tic
[EE_nocoll_alpha_HW,PFE_nocoll_alpha_HW, peak_pfe_nocoll_alpha_HW, EPE_nocoll_alpha_HW, expected_exposure_delta_alpha_HW_coll, expected_positive_exposure_delta_alpha_HW_coll, potential_future_exposure_delta_alpha_HW_coll, peak_pfe_delta_alpha_HW_coll, ...
 expected_exposure_epsilon_alpha_HW_coll, expected_positive_exposure_epsilon_alpha_HW_coll, potential_future_exposure_epsilon_alpha_HW_coll, peak_pfe_epsilon_alpha_HW_coll, EPE_HW_alpha, peakPFE_HW_alpha] = ...
    RiskMeasures_HW(alphaPortfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach, 0);
toc

%% HW - beta ptf
tic
[EE_nocoll_beta_HW,PFE_nocoll_beta_HW, peak_pfe_nocoll_beta_HW, EPE_nocoll_beta_HW, expected_exposure_delta_beta_HW_coll, expected_positive_exposure_delta_beta_HW_coll, potential_future_exposure_delta_beta_HW_coll, peak_pfe_delta_beta_HW_coll, ...
 expected_exposure_epsilon_beta_HW_coll, expected_positive_exposure_epsilon_beta_HW_coll, potential_future_exposure_epsilon_beta_HW_coll, peak_pfe_epsilon_beta_HW_coll,  EPE_HW_beta, peakPFE_HW_beta] = ...
    RiskMeasures_HW(betaPortfolio,  dates, discounts, simulations_num, alpha, aBach, sigmaBach, 0);
toc

%% HW - gamma ptf
tic
[EE_nocoll_gamma_HW,PFE_nocoll_gamma_HW, peak_pfe_nocoll_gamma_HW, EPE_nocoll_gamma_HW,expected_exposure_delta_gamma_HW_coll, expected_positive_exposure_delta_gamma_HW_coll, potential_future_exposure_delta_gamma_HW_coll, peak_pfe_delta_gamma_HW_coll, ...
 expected_exposure_epsilon_gamma_HW_coll, expected_positive_exposure_epsilon_gamma_HW_coll, potential_future_exposure_epsilon_gamma_HW_coll, peak_pfe_epsilon_gamma_HW_coll,  EPE_HW_gamma, peakPFE_HW_gamma] = ...
    RiskMeasures_HW(gammaPortfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach, 0);
toc

%% POINT 4.c weekly time step - alpha ptf - LGM model
simulations_num = 2500; 
tic
[collateral_Delta_alpha, collateral_Epsilon_alpha] = Weekly_LGM_Collateral(alphaPortfolio, jamshResults, dates, discounts, simulations_num, aLGM, sigmaBach);
toc
%% POINT 4.c weekly time step - beta ptf - LGM model
[collateral_Delta_beta, collateral_Epsilon_beta] = Weekly_LGM_Collateral(betaPortfolio, jamshResults, dates, discounts, simulations_num, aLGM, sigmaBach);

%% POINT 4.c weekly time step - gamma ptf - LGM model
[collateral_Delta_gamma, collateral_Epsilon_gamma] = Weekly_LGM_Collateral(gammaPortfolio, jamshResults, dates, discounts, simulations_num, aLGM, sigmaBach);

%% POINT 4.c weekly time step - alpha ptf - HW model
[collateral_Delta_alpha_HW, collateral_Epsilon_alpha_HW] = Weekly_HW_Collateral(alphaPortfolio, dates, discounts, simulations_num, aBach, sigmaBach);

%% POINT 4.c weekly time step - beta ptf - HW model
[collateral_Delta_beta_HW, collateral_Epsilon_beta_HW] = Weekly_HW_Collateral(betaPortfolio, dates, discounts, simulations_num, aBach, sigmaBach);

%% POINT 4.c weekly time step - gamma ptf - HW model
[collateral_Delta_gamma_HW, collateral_Epsilon_gamma_HW] = Weekly_HW_Collateral(gammaPortfolio, dates, discounts, simulations_num, aBach, sigmaBach);

%% POINT 4.c daily time step - LGM - alpha ptf
% WARNING: FROM THIS POINT ON, EACH SECTION TAKES ABOUT 7 MINUTES
simulations_num = 100;
[collateral_epsilon_daily_alpha_LGM, collateral_delta_daily_alpha_LGM, expected_exposure_delta_daily_alpha_LGM,expected_positive_exposure_delta_daily_alpha_LGM,potential_future_exposure_delta_daily_alpha_LGM,peak_pfe_delta_daily_alpha_LGM, expected_exposure_epsilon_daily_alpha_LGM,expected_positive_exposure_epsilon_daily_alpha_LGM,potential_future_exposure_epsilon_daily_alpha_LGM,peak_pfe_epsilon_daily_alpha_LGM] = Daily_LGM_Collateral(alphaPortfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach);

%% POINT 4.c daily time step - LGM - beta ptf
[collateral_epsilon_daily_beta_LGM, collateral_delta_daily_beta_LGM, expected_exposure_delta_daily_beta_LGM,expected_positive_exposure_delta_daily_beta_LGM,potential_future_exposure_delta_daily_beta_LGM,peak_pfe_delta_daily_beta_LGM, expected_exposure_epsilon_daily_beta_LGM,expected_positive_exposure_epsilon_daily_beta_LGM,potential_future_exposure_epsilon_daily_beta_LGM,peak_pfe_epsilon_daily_beta_LGM] = Daily_LGM_Collateral(betaPortfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach);

%% POINT 4.c daily time step - LGM - gamma ptf
[collateral_epsilon_daily_gamma_LGM, collateral_delta_daily_gamma_LGM, expected_exposure_delta_daily_gamma_LGM,expected_positive_exposure_delta_daily_gamma_LGM,potential_future_exposure_delta_daily_gamma_LGM,peak_pfe_delta_daily_gamma_LGM, expected_exposure_epsilon_daily_gamma_LGM,expected_positive_exposure_epsilon_daily_gamma_LGM,potential_future_exposure_epsilon_daily_gamma_LGM,peak_pfe_epsilon_daily_gamma_LGM] = Daily_LGM_Collateral(gammaPortfolio, jamshResults, dates, discounts, simulations_num, alpha, aLGM, sigmaBach);

%% POINT 4.c daily time step - HW - alpha ptf
[collateral_epsilon_daily_alpha_HW, collateral_delta_daily_alpha_HW, expected_exposure_delta_daily_alpha_HW,expected_positive_exposure_delta_daily_alpha_HW,potential_future_exposure_delta_daily_alpha_HW,peak_pfe_delta_daily_alpha_HW, expected_exposure_epsilon_daily_alpha_HW,expected_positive_exposure_epsilon_daily_alpha_HW,potential_future_exposure_epsilon_daily_alpha_HW,peak_pfe_epsilon_daily_alpha_HW] = Daily_HW_Collateral(alphaPortfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach);

%% POINT 4.c daily time step - HW - beta ptf
[collateral_epsilon_daily_beta_HW, collateral_delta_daily_beta_HW, expected_exposure_delta_daily_beta_HW,expected_positive_exposure_delta_daily_beta_HW,potential_future_exposure_delta_daily_beta_HW,peak_pfe_delta_daily_beta_HW, expected_exposure_epsilon_daily_beta_HW,expected_positive_exposure_epsilon_daily_beta_HW,potential_future_exposure_epsilon_daily_beta_HW,peak_pfe_epsilon_daily_beta_HW] = Daily_HW_Collateral(betaPortfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach);

%% POINT 4.c daily time step - HW - gamma ptf
[collateral_epsilon_daily_gamma_HW, collateral_delta_daily_gamma_HW, expected_exposure_delta_daily_gamma_HW,expected_positive_exposure_delta_daily_gamma_HW,potential_future_exposure_delta_daily_gamma_HW,peak_pfe_delta_daily_gamma_HW, expected_exposure_epsilon_daily_gamma_HW,expected_positive_exposure_epsilon_daily_gamma_HW,potential_future_exposure_epsilon_daily_gamma_HW,peak_pfe_epsilon_daily_gamma_HW] = Daily_HW_Collateral(gammaPortfolio, dates, discounts, simulations_num, alpha, aBach, sigmaBach);
