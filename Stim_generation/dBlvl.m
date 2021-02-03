function [scaledS] = setdBlvl(S,outputlvl,reflvl)
% [scaledS] = dBlvl(S,outputlvl,reflvl)
% multiplies input sound S by a factor so that it will be played at level
% outputlvl (in dB SPL) with reflvl the level of a reference sound at -18
% dB full scale (default reflvl=72.6 dB SPL)
%
% Léo Varnet - 2020

if nargin<3
    reflvl = 72.6; % Calibration Sennheiser HD 650
end    
refrms = 0.088982893153102; % rms of a pure tone at -18 dB FS

scaledS = (S/rms(S))*refrms*10^((-reflvl+outputlvl)/20);

end

