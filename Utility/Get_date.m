function [p,numeric_y] = Get_date
% function [p,numeric_y] = Get_date
%
% 1. Description: 
%       Returns day, month and year in a struct (different fields). Time
%       is also returned in format hh:mm:ss
% 
%       p - output struct with date (string) and time
%
% 2. Example:
%       p = Get_date;
%       disp([p.dd '/' p.mm '/' p.yyyy ' at ' p.time]);
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium 2013
% Created in    : 2013
% Last update on: 26/11/2014 
% Last use on   : 26/11/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    type = 'time';
end

DateNum = fix(clock);

p.dd    = Num2str(DateNum(3),2);
p.mm    = Num2str(DateNum(2),2);
p.yyyy  = Num2str(DateNum(1),4);

if strcmp(type,'date')
    y = [p.dd,'-',p.mm,'-',p.yyyy,' (dd-mm-yy)'];
end

p.time          = [num2str(DateNum(4)),':',num2str(DateNum(5)),':',num2str(DateNum(6)),' (hh:mm:ss)'];
p.date2print    = [p.yyyy,'-',p.mm,'-',p.dd,'-at-' Num2str(DateNum(4)),'h-',Num2str(DateNum(5)),'m-',Num2str(DateNum(6)),'s'];
p.date4plots    = [p.yyyy,'-',p.mm,'-',p.dd,',' Num2str(DateNum(4)),':',Num2str(DateNum(5)),''];
p.date4files    = [p.yyyy p.mm p.dd];
numeric_y = DateNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end