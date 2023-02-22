function definput = arg_vandorpschuitman2013a(definput)
%ARG_VANDORPSCHUITMAN2013A(definput)
%   #Author: Alejandro Osses

%% General options of the model:
definput.keyvals.order  = 1; % order of the combined outer- middle-ear filter

% %% For binaural processor:
% definput.keyvals.muASW = 2.00e-2;
% definput.keyvals.nuASW = 5.63e2;
% definput.keyvals.muLEV = 2.76e-2;
% definput.keyvals.nuLEV = 6.80e2;

%% For central processor:
definput.keyvals.Psimin = 0.34; % 'Optimisation' by running m20160714_fourth_report_MBBM_testing_vanDorp
                                % criterion used: minimise error (mean absolute difference).
% Psimin_dip determined from the ratio abs(-7.49e-3/1.33e-3), which are the
% Psimin and Psimin_dip as suggested by van Dorp:
definput.keyvals.Psimin_dip = -definput.keyvals.Psimin/5.63; 
definput.keyvals.Tmin   = 63.1e-3; % 63.1 [ms]
definput.keyvals.framelen = 5; % s
definput.keyvals.hopsize  = 1; % s

% definput.flags.multiple_one_second={'no_multiple_hopsize','multiple_hopsize'};
definput.flags.exclude_silence={'no_exclude','exclude'};

definput.groups.publ_vandorpschuitman2013 = {'no_exclude','framelen',[],'hopsize',[]}; % 'no_multiple_hopsize'
definput.groups.publ_osses2017 = {'exclude','framelen',5,'hopsize',1}; % 'multiple_hopsize'

end
