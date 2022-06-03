function definput = arg_fastACI_simulation_detect(definput)
% ARG_FASTACI_SIMULATION_DETECT
%
%   #Author: Alejandro Osses (2021-)

definput.keyvals.maxtimelag_ms=0;
definput.keyvals.expvar_cal = -40; % dB, default SNR for calibration in speech experiments
definput.groups.maxtimelag_osses2021={'maxtimelag_ms',50};