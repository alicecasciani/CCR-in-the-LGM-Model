function [dates, rates] = readExcelData( filename, formatData)
% Reads data from excel
%  It reads bid/ask prices and relevant dates
%  All input rates are in % units
%
% INPUTS:
%  filename: excel file name where data are stored
%  formatData: data format in Excel
% 
% OUTPUTS:
%  dates: struct with settlementDate, deposDates, futuresDates, swapDates
%  rates: struct with deposRates, futuresRates, swapRates

%% Dates from Excel

%Settlement date
[~, settlement] = xlsread(filename, 1, 'B2');
%Date conversion
dates.settlement = datenum(settlement, formatData);

%% Depos
%Dates relative to depos
date_depositi = addtodate(dates.settlement, 3, "month");
dates.depos = businessdayoffset(date_depositi);


%% Fra 

%Date relative to forward: expiry dates
dates_fra = [addtodate(dates.settlement, 1, "month"), addtodate(dates.settlement, 4, "month");
             addtodate(dates.settlement, 2, "month"), addtodate(dates.settlement, 5, "month");
             addtodate(dates.settlement, 3, "month"), addtodate(dates.settlement, 6, "month")];

dates.fra(:,1) = businessdayoffset(dates_fra(:, 1));
dates.fra(:,2) = businessdayoffset(dates_fra(:, 2));

%% Futures

%Dates relative to futures: calc start & end
[~, date_futures_read] = xlsread(filename, 1, 'Q6:R13');
numberFutures = size(date_futures_read,1);

dates.futures=ones(numberFutures,2);
dates.futures(:,1) = datenum(date_futures_read(:,1), formatData);
dates.futures(:,2) = datenum(date_futures_read(:,2), formatData);


%% Swaps

%Date relative to swaps: expiry dates
today = datetime(2022,6,28);
maturities  = [2 3 4 5 6 7 8 9 10 11 12 15 20 25 30 40 50];
expiryDates = today + calyears(maturities);
dates.swaps = businessdayoffset(expiryDates);


%% Rates from Excel (Bids & Asks)

%Depos
tassi_depositi = xlsread(filename, 1, 'J6');
rates.depos = tassi_depositi / 100;

%Fra
tassi_fra = xlsread(filename, 1, 'J7:J9');
rates.fra = tassi_fra / 100;

%Futures
tassi_futures = xlsread(filename, 1, 'J10:J17');
%Rates from futures
tassi_futures = 100 - tassi_futures;
rates.futures = tassi_futures / 100;

%Swaps
tassi_swaps = xlsread(filename, 1, 'J18:J34');
rates.swaps = tassi_swaps / 100;

end % readExcelData