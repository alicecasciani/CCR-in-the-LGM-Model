function [alphaPortfolio, betaPortfolio, gammaPortfolio] = readIRSportfolios(filename)
% readIRSportfolios Reads three IRS portfolios from an Excel workbook
%   [alphaP, betaP, gammaP] = readIRSportfolios(filename) imports portfolios
%   from sheets 'Alpha', 'Beta', and 'Gamma', returning them as tables.
%
%   INPUT:
%     filename - Path to the Excel file containing portfolio data
%
%   OUTPUTS:
%     alphaPortfolio - Table of IRS trades in the Alpha portfolio
%     betaPortfolio  - Table of IRS trades in the Beta portfolio
%     gammaPortfolio - Table of IRS trades in the Gamma portfolio
   
    % Store current warning state and suppress all warnings
    origWarnState = warning;          % capture the current warning configuration
    warning('off', 'all');            % turn off every warning
    cleanupObj = onCleanup(@() warning(origWarnState)); 
    % Read Alpha portfolio
    alphaData = readtable(filename, 'Sheet', 'Alpha', 'Range', 'A2:H102');
    alphaData.Properties.VariableNames = {
        'NotionalAmount', 
        'Maturity', 
        'FixedLegPaymentFrequency', 
        'FixedLeg', 
        'FixedLegRate', 
        'FloatingLegPaymentFrequency', 
        'FloatingLeg', 
        'CounterpartyName'
    };
    alphaPortfolio = alphaData;
    
    % Read Beta portfolio
    betaData = readtable(filename, 'Sheet', 'Beta', 'Range', 'A2:H102');
    betaData.Properties.VariableNames = {
        'NotionalAmount', 
        'Maturity', 
        'FixedLegPaymentFrequency', 
        'FixedLeg', 
        'FixedLegRate', 
        'FloatingLegPaymentFrequency', 
        'FloatingLeg', 
        'CounterpartyName'
    };
    betaPortfolio = betaData;
    
    % Read Gamma portfolio
    gammaData = readtable(filename, 'Sheet', 'Gamma', 'Range', 'A2:H102');
    gammaData.Properties.VariableNames = {
        'NotionalAmount', 
        'Maturity', 
        'FixedLegPaymentFrequency', 
        'FixedLeg', 
        'FixedLegRate', 
        'FloatingLegPaymentFrequency', 
        'FloatingLeg', 
        'CounterpartyName'
    };
    gammaPortfolio = gammaData;
    
end
