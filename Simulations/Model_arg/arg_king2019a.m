function definput=arg_king2019a(definput)
% ARG_KING2019A
%
%   #Author: Alejandro Osses

%% Internal noise
definput.flags.internalnoiseadd  = {'internalnoiseadd' ,'no_internalnoiseadd'};
definput.flags.internalnoisemult = {'internalnoisemult','no_internalnoisemult'};
definput.flags.internalnoisemem  = {'internalnoisemem' ,'no_internalnoisemem'};

% Additive noise:
definput.keyvals.intnoiseadd_std  = 425; % a.u., default in Andrew's codes

% Multiplicative noise:
definput.keyvals.intnoisemult_ratio_dB  = 1; % dB, King et al 2019

% Memory noise:
definput.keyvals.intnoisemem_std  = 50; %a.u., King et al 2019
definput.keyvals.intnoisemem_tau  = 1.2; % s, Wallaert et al 2017 and King et al 2019