function cfg_out = speechACI_varnet2015_cfg(cfg_in)
% function str_stim = speechACI_varnet2015_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.response_names = {'da', 'ga'}; 
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation
cfg.displayN       = 1; % 'oui'
cfg.sessionsN      = 100; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 1; % 'out';%
cfg.randorder      = 1;
 
cfg.startvar = 0;  % old name 'm_start'
cfg.expvar_description = 'SNR (dB)';

cfg.maxvar = 10;
 
% Staircase algorithm parameters
if cfg.adapt == 1
	cfg.start_stepsize     = 4;
    cfg.min_stepsize       = 1;
    cfg.adapt_stepsize     = 50/100;
else
    error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);