function [date_str, clock_now] = Get_date_and_time_str
% function [date_str, clock_now] = Get_date_and_time_str
%
% TODO: Fix the number of digits for each time unit (to see '02' instead of '2')
%
% Author: Alejandro Osses and the fastACI team
% Created on: 3/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clock_now = fix(clock);
date_str  = sprintf('%.0f_%.0f_%.0f_%.0f_%.0f',clock_now(1:5)); % No seconds