function cfg_out = modulationFM_cfg(cfg_in)
% function str_stim = modulationFM_cfg(cfg_in)
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
cfg.response_names = {'pure tone', 'modulated tone'}; 
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation

bDebug             = 1;
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;

%cfg_game.end_sessions   = [500 1000 1500 2000 2500]; 
cfg.sessionsN      = 500; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 1; % 'out';%
cfg.randorder      = 1;

cfg.startvar = 40;  % i.e., 0.04 times 1000 Hz, according to King et al. 2019
cfg.expvar_description = 'Frequency deviation (Hz)';

cfg.maxvar = 1000; % basically this is not set
cfg.minvar =    0;

% Staircase algorithm parameters
switch cfg.adapt 
    case {1,'transformed-up-down'}
        % cfg.step_resolution = 'linear';
        cfg.start_stepsize     = 1.58; % Hz
        cfg.min_stepsize       = 1.1172; % Hz, King et al. (2019, JASA)
        cfg.adapt_stepsize     = 1/(sqrt(sqrt(2))); % King et al. (2019, JASA) Implicit choice
        cfg.step_resolution = 'multiplicative';
    otherwise
        error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);
