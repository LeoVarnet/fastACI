function cfg_out = localisationILD_cfg(cfg_in)
% function str_stim = localisationILD_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.Language       = 'EN'; % or 'EN' or 'FR'
cfg.response_names = {'left', 'right'}; 
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation

bDebug             = 1;
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;

%cfg_game.end_sessions   = [500 1000 1500 2000 2500]; 
cfg.sessionsN      = 500; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 2; % 'out';%
cfg.randorder      = 1;

cfg.startvar = 5;  % 5 dB
cfg.expvar_description = 'Level difference (dB)';

cfg.maxvar = 20; % basically this is not set
cfg.minvar = -20;

cfg.adapt_stepsize = 50/100; % might be overwritten below...

% Staircase algorithm parameters
switch cfg.adapt 
    case {1,'transformed-up-down'}
        % cfg.step_resolution = 'linear';
        cfg.start_stepsize     = 1; % dB
        cfg.min_stepsize       = 0.2; % Hz, King et al. (2019, JASA)
        cfg.adapt_stepsize     = 0.9; % King et al. (2019, JASA) Implicit choice
        cfg.step_resolution = 'linear';
    case {2, 'weighted-up-down'}
        cfg.start_stepsize     = 1; % dB
        cfg.rule = [1 1]; 
        target_score = .707;
        cfg.step_down  = 1;
        cfg.step_up    = cfg.step_down*target_score/(1-target_score); % Kaernbach1991, Eq. 1
        cfg.min_stepsize = 1/cfg.step_up;
    otherwise
        error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);
