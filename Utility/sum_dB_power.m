function [Y,ymean] = sum_dB_power(X,floor_dB)
% function [Y,ymean] = sum_dB_power(X,floor_dB)
%
% 1. Description:
%       Y [dB] (related to a pressure of y [Pa]) is obtained by combining 
%       Xn SPL levels (with a pressure of xn [Pa]) in such a way that y is:
%           y = sqrt(sum(x.^2))
% 
% 2. Stand-alone example:
%       X = [60 60];
%       [Y,ymean] = sum_dB_power(X); % expected result: 63 dB
% 
%       X = [60 64; 60 60];
%       [Y,ymean] = sum_dB_power(X); % expected result: 65.5 and 63.0 dB, Y     for first and second row
%                                    %                  62.4 and 60.0 dB, Ymean for first and second row
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: sum_dB_arit.m
%       Old name: sum_db.m
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 06/06/2014
% Last update on: 21/07/2016 
% Last use on   : 19/12/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    floor_dB = -eps;
end

[N,M] = size(X);

x2 = 10.^(X/10); % Anti-log squared

Y  =  10*log10( sum(x2,2) );

if nargout > 1
    
    idx = find(x2<floor_dB);
    x2(idx) = floor_dB;
    
    ymean =  10*log10(1/M * sum(x2,2) );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
