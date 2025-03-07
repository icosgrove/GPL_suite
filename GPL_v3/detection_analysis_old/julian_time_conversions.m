function [margin, datep, datem] = julian_time_conversions(date_string, amount_of_time, time_units)
% This function identifies how much julian time needs to be added/subtratced from
% the original julian time value to get a certain margin of time. 
% Inputs: date_string: Input the given date in DD-MM-YYYY HH-MM-SS or
% whatever format as long as it is a string. Enter the amount of time you
% want for the margin, e.g. '5' for 5 seconds. Enter the units for that
% number e.g. 'S' for 5 seconds.  
% Supported time units are: 'Y', 'Mo', 'D', 'H', 'M', 'S' for year, month, 
% day, hour, minute, second. 
% The output 'margin' is the amount of julian time to be added/subtracted to your
% date. The output 'datep' is the julian time for your original date with
% differnence added, and 'datem' is the julian time for your original date
% with difference subtracted. 
% Note that there is significant error for large time spans (order of
% months/years) See code for that error. 


% Calculate how much julian time is for each time unit
sec_diff = datenum('01-Jan-2022 12:00:01') - datenum('01-Jan-2022 12:00:00'); % No error present on order of whole seconds (1s)
min_diff = 60*sec_diff; % No error present on order of whole seconds for 1 minute
hour_diff = 60*min_diff; % Error: +/- 1 second (1 Hour)
day_diff = 24*hour_diff; % Error: +/- 1 second (1 day)
month_diff = 30.437*day_diff; % Error: Multiple hours due to using average amount of days in a month.
year_diff = 365*day_diff; % Error: +/- 2.5 mins (1 year) 


% Switch statement to produce outputs
switch time_units
    case 'Y' % Year
        margin = amount_of_time*year_diff;
        datep = datenum(date_string) + margin;
        datem = datenum(date_string) - margin;
    case 'Mo' % Month
        margin = amount_of_time*month_diff;
        datep = datenum(date_string) + margin;
        datem = datenum(date_string) - margin;
    case 'D' % Day
        margin = amount_of_time*day_diff;
        datep = datenum(date_string) + margin;
        datem = datenum(date_string) - margin;
    case 'H' % Hour
        margin = amount_of_time*hour_diff;
        datep = datenum(date_string) + margin;
        datem = datenum(date_string) - margin;
    case 'M' % Minute
        margin = amount_of_time*min_diff;
        datep = datenum(date_string) + margin;
        datem = datenum(date_string) - margin;        
    case 'S' % Second
        margin = amount_of_time*sec_diff;
        datep = date_string + margin;
        datem = date_string - margin;
    otherwise
        margin = 'Wrong time unit format';
        datep = nan;
        datem = nan;
end
