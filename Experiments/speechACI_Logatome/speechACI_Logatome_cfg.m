function cfg_out = speechACI_Logatome_cfg(cfg_in)
% function str_stim = speechACI_Logatome_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.filename_target  = {'ap_pa.wav','at_ta.wav'};
cfg.response_names   = {'pa', 'ta'}; 
cfg.response_correct_target = [1,2];
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation
 
bDebug = 1; 
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;
 
cfg.sessionsN      = 400; % CAUTION: Overwritten in the case of simulation
cfg.adapt          = 1; % 'out';%
cfg_out.randorder      = 1;
  
cfg.startvar = 0;  % old name 'm_start'
cfg.expvar_description = 'SNR (dB)';
 
% cfg.maxvar = 10;
 
% Staircase algorithm parameters
if cfg.adapt == 1
	cfg.start_stepsize     = 4;
    cfg.min_stepsize       = 1;
    cfg.adapt_stepsize     = 50/100;
else
    error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);