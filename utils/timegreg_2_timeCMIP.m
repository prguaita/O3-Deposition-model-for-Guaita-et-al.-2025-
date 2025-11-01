function [start_time, end_time] = timegreg_2_timeCMIP(start_date, end_date, base_date, model_name)
%TIMEGREG_2_TIMECMIP Converts calendar dates to CMIP model time units
%
%   This function converts input start and end dates into CMIP time,
%   defined as the number of days elapsed since a given base date (usually
%   1850-01-01). Different climate models handle leap years and calendar days
%   differently: some include leap years, others use a fixed 365- or 360-day year.
%
%   INPUTS:
%     start_date  - datetime for the start of the period
%     end_date    - datetime for the end of the period
%     base_date   - datetime representing the CMIP base date (reference zero)
%     model_name  - string specifying the climate model name, affects time calculation
%
%   OUTPUTS:
%     start_time  - model time corresponding to start_date (days since base_date)
%     end_time    - model time corresponding to end_date (days since base_date)
%
%   Supported model names:
%     'MPI-ESM-1-0-HAM' : Includes leap years, uses standard calendar days
%     'UKESM1-0-LL'     : Uses a 360-day calendar, with 5 days removed per year
%     'GFDL-ESM4'       : Includes leap years, standard calendar days
%     'MIROC-ES2H'      : Includes leap years, standard calendar days
%
%   Example:
%     [s,e] = timegreg_2_timeCMIP(datetime(2000,1,1), datetime(2000,12,31), datetime(1850,1,1), 'MPI-ESM-1-0-HAM');

% Calculate number of leap days from base_date to start_date and end_date
if base_date.Year == start_date.Year
    n_leap_start = 0;
else
    n_leap_start = sum(leapyear(base_date.Year:1:start_date.Year));
end
n_leap_end = sum(leapyear(base_date.Year:1:(end_date.Year-1)));

switch model_name
    case 'MPI-ESM-1-0-HAM'
        % Standard calendar with leap years
        start_time = 1 + days(start_date - base_date);
        end_time   = 1 + days(end_date - base_date);
    case 'UKESM1-0-LL'
        % 360-day calendar: subtract leap days and 5 extra days per year
        start_time = 1 + days(start_date - base_date) - n_leap_start - 5*(start_date.Year - base_date.Year);
        end_time   = 1 + days(end_date - base_date) - n_leap_end   - 5*(end_date.Year   - base_date.Year);
    case 'GFDL-ESM4'
        % Standard calendar with leap years
        start_time = 1 + days(start_date - base_date) - n_leap_start;
        end_time   = 1 + days(end_date - base_date)   - n_leap_end;
    case 'MIROC-ES2H'
        % Standard calendar with leap years
        start_time = 1 + days(start_date - base_date);
        end_time   = 1 + days(end_date - base_date);
    otherwise
        error('Invalid model name: %s', model_name)
end

end
