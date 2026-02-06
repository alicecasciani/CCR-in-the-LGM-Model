function dates_bus = businessdayoffset(dates)
% INPUT: dates to be checked
% OUTPUT: business dates, using a Following convention

dates_bus=dates;

for i=1:length(dates)
    if weekday(dates(i)) == 1
        dates_bus(i) = dates(i)+1;
     elseif weekday(dates(i)) == 7
         dates_bus(i) = dates(i)+2;
     end
end

end
