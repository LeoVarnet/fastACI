function cfg_out = speechACI_varnet2013_cfg(cfg_in)
% function str_stim = speechACI_varnet2013_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.filename_target  = {'Aba.wav','Ada.wav'};
cfg.response_names   = {'aba', 'ada'}; 
cfg.response_correct_target = [1,2]; % reponse correcte pour chaque signal (signaux dans l'ordre alphabetique)
                   % Alda is 'da', Alga is 'ga', Arda is 'da', Arga is 'ga'
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation

cfg.Language       = 'FR'; % or 'EN'
bDebug = 1; 
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;

cfg.sessionsN      = 400; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 2; warning('default is 1 - now we are only testing')% 
cfg.randorder      = 1;
 
cfg.startvar = 0;  % old name 'm_start'
cfg.expvar_description = 'SNR (dB)';

cfg.maxvar = 10;
 
% Staircase algorithm parameters
cfg.start_stepsize = 2;
cfg.adapt_stepsize = 50/100;
cfg.min_stepsize   = 1;
cfg.step_resolution = 'linear';

switch cfg.adapt
    case {1,'transformed-up-down'}
        cfg.rule = [1 2]; % [up down]-rule: [1 2] = 1-up 2-down   
        cfg.step_up    = 1;
        cfg.step_down  = 1;
        
    case {2, 'weighted-up-down'}
        cfg.rule = [1 1]; 
        target_score = .707;
        cfg.step_down  = 1;
        cfg.step_up    = cfg.step_down*target_score/(1-target_score);
        
    otherwise
        error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);
