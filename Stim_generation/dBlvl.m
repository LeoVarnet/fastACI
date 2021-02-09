function [S_cal,dBFS,gain_factor] = dBlvl(S,outputlvl,reflvl)
% function [S_cal,dBFS,gain_factor] = dBlvl(S,outputlvl,reflvl)
% multiplies input sound S by a factor so that it will be played at level
% outputlvl (in dB SPL) with reflvl the level of a reference sound at -18
% dB full scale (default reflvl=72.6 dB SPL)
%
% Leo Varnet - 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    reflvl = 72.6; % Calibration Sennheiser HD 650
    warning('%s: Default reference level %.1f dB is being used (calibration if Sennheiser HD 650 are used)',upper(mfilename))
end    
refrms = 0.088982893153102; % rms of a pure tone with -18 dB FS peak amplitude (-21 dB FS RMS)
lvl_rms_FS  = 20*log10(refrms);
lvl_peak_FS  = lvl_rms_FS+3; % valid for a pure tone only: here only for the record
dBFS        = reflvl-lvl_rms_FS;

factor1 = (1/rms(S)); % this makes that 'S' has an RMS of 1 (or full scale)
factor2 = 10^((-reflvl+outputlvl)/20); % gain factor according to how many dB 
                      % above or below the reference level are required
gain_factor = factor1*refrms*factor2;
S_cal = S*gain_factor;