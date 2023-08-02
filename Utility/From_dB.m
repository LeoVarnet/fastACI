function gain = From_dB(gain_dB,div)
% function gain = From_dB(gain_dB,div)
%
% 1. Description:
%       From_dB: Convert decibels to voltage gain (if div = 20, default).
%       gain = From_dB(gain_dB)
%
% Author        : Brett Swanson
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2018
% Last update on: 30/07/2014 
% Last use on   : 14/11/2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    div = 20;
end

gain = 10 .^ (gain_dB / div);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
