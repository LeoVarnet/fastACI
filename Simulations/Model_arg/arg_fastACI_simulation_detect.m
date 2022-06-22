function definput = arg_fastACI_simulation_detect(definput)
% ARG_FASTACI_SIMULATION_DETECT
%
%   #Author: Alejandro Osses (2021-)

definput.keyvals.maxtimelag_ms=0;
definput.keyvals.expvar_cal = -40; % dB, default SNR for calibration in speech experiments
definput.flags.bias_each_session = {'bias_global'   ,'no_bias_each_session', ... % these two options are equivalent
                                    'no_bias_global','bias_each_session'}; % these two options are equivalent
definput.groups.maxtimelag_osses2021={'maxtimelag_ms',50};