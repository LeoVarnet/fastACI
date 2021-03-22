function cfg_inout = modulationACI_seeds_set(cfg_inout)
% function cfg_inout = modulationACI_seeds_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
% function cfg_out = il_set(cfg_in)
%
% Function comparable to *_set.m functions from AFC toolbox

cfg = [];

%%% Parameters to create the targets:
cfg.fs        = 48000; % Hz, sampling frequency
cfg.fm        = 4;
cfg.fc        = 1000;
cfg.stim_dur  = .75;   % s,  stimulus suration (s)
cfg.SNR       = -10;

% cfg.expvar_start = -8;  % old name 'm_start'
% cfg.expvar_description = 'modulation depth (dB)';

cfg.fadein_s  = 0.075; % CAUTION: Overwritten in the case of simulation
cfg.SPL       = 65;
cfg.dBFS      = 93.61; % Calibration with Headphones HD 650

cfg.N_presentation = 1500;  % number of stimuli / condition
cfg.N_target   = 2;     % Number of conditions
cfg.noise_type = 'white';
cfg.randorder  = 1;

cfg_inout = Merge_structs(cfg,cfg_inout);