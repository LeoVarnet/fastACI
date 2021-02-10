function [date_str, clock_now] = Get_date_and_time_str
% function [date_str, clock_now] = Get_date_and_time_str
%
% TODO: Fix the number of digits for each time unit (to see '02' instead of '2')
%
% Author: Alejandro Osses and the fastACI team
% Created on: 3/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clock_now = fix(clock);

zeros_str = {'','','','',''};

for i = 2:5 % check months, days, hh, mm
    if clock_now(i)<10
        zeros_str{i}='0';
    end
end

date_str = sprintf('%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f',zeros_str{1},clock_now(1), ...
    zeros_str{2},clock_now(2),zeros_str{3},clock_now(3),zeros_str{4},clock_now(4), ...
    zeros_str{5},clock_now(5)); % No seconds