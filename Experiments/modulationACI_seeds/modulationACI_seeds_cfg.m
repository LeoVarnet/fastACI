function cfg_out = modulationACI_seeds_cfg(cfg_in)
% function str_stim = modulationACI_seeds_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.response_names = {'pure tone', 'modulated tone'}; 
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation

bDebug = 0;
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;

%cfg_game.end_sessions   = [500 1000 1500 2000 2500]; 
cfg.sessionsN      = 300; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 1; % 'out';%
cfg.randorder      = 1;

cfg.startvar = -8;  % old name 'm_start'
cfg.expvar_description = 'modulation depth (dB)';

cfg.maxvar = 0;

% Staircase algorithm parameters
if cfg.adapt == 1
	cfg.start_stepsize     = 4;
	cfg.min_stepsize       = 1;
    cfg.adapt_stepsize     = 90/100;
else
    error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);