function cfg_out = modulationACI_cfg(cfg_in)
% function str_stim = modulationACI_cfg(cfg_in)
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

cfg.startvar = -1;  % dB, as in Varnet and Lorenzi (2022, JASA)
cfg.expvar_description = 'modulation depth (dB)';

cfg.maxvar = 0;

% Staircase algorithm parameters
switch cfg.adapt 
    case {1,'transformed-up-down'}
        cfg.start_stepsize     = 2;  % dB, Varnet and Lorenzi (2022, JASA)
        cfg.min_stepsize       = .5; % dB, Varnet and Lorenzi (2022, JASA)
        cfg.adapt_stepsize     = 90/100; % Varnet and Lorenzi (2022, JASA)
        cfg.step_resolution = 'linear';
    otherwise
        error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);
