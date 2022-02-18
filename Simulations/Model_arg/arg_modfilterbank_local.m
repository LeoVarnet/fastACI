function definput = arg_modfilterbank_local(definput)
% ARG_MODFILTERBANK_LOCAL
%
%   #License: GPL
%   #Author: Alejandro Osses (2020): Initial version
%   #Author: Clara Hollomey (2021): Adapted to AMT
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
%
% References: dau1997mapI verhey1999 jepsen2008cmh

definput.flags.mfb = {'mfb','no_mfb'};
definput.flags.modfilterlimit   = {'mfc_upper_limit','no_mfc_upper_limit'};
definput.flags.modfilter_150Hz_LP = {'LP_150_Hz','no_LP_150_Hz','LP_150_Hz_att'};
definput.flags.phase_insens = {'phase_insens_hilbert','no_phase_insens'};
% Attenuation factor applied to mod filters above 10 Hz (only applied if phase_insens_hilbert is on):
definput.flags.att_factor       = {'att_factor','no_att_factor'}; 

definput.keyvals.mfb_script = 'osses2022_modfilterbank';
definput.keyvals.mfc_upper_limit_max = 1000; % Hz, max limit
definput.keyvals.mflow  = 2; % Hz, modbank_fmin
definput.keyvals.mfhigh = 150; % Hz, modbank_fmax
definput.keyvals.Q_mfb =  1; % Q factor for the filters
definput.flags.modfilter_phase = {'phase_insens_hilbert', 'no_phase_insens'};
definput.keyvals.phase_insens_cut = 10; % Hz 
definput.keyvals.modbank_Nmod    = []; % number of filters, for overalpped 
                               % filters choose 'modbank_Nmod'

% definput.keyvals.mfc_upper_limit_max = 1000; % Hz, maximum upper limit

definput.groups.mfb_osses2022a = {'mfc_upper_limit',... % previuous name: mfb_osses2022_local
                            'no_LP_150_Hz', ...
                            'att_factor', ...
                            'mflow',4.0451, ... % if lowest -3 dB point is 2.5 Hz => mfc*(sqrt(5)/2)-mfc*1/2 = 2.5 => mfc = 4.0451 Hz
                            'mfb_script','osses2022_modfilterbank'}; 
                        
definput.groups.mfb_osses2021 = {'mfc_upper_limit',...
                            'LP_150_Hz_att', ...
                            'att_factor', ...
                            'Q_mfb',2, ...
                            'mfc_upper_limit_max', 1000, ...
                            'mfb_script','modfilterbank'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
