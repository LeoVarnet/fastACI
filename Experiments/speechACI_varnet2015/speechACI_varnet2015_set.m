function cfg_inout = speechACI_varnet2015_set(cfg_inout)
% function cfg_out = speechACI_varnet2015_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
 
dir_main = '/home/alejandro/Documents/Databases/data/fastACI/speechACI/';
dir_speech = [dir_main 'speech-samples' filesep];
dir_noise  = [dir_main 'NoiseStim_S11'  filesep];

dBFS = 100;
% lvl_target = 65;

%%% Parameters to create the targets:
cfg.fs        = 48000; % Hz, sampling frequency
% cfg.fm        = 4;
% cfg.fc        = 1000;
% cfg.stim_dur  = 1.5;   % s,  stimulus suration (s)
% cfg.SNR       = -10;

% cfg.fadein_s  = 0.075; % CAUTION: Overwritten in the case of simulation
% cfg.SPL       = 65;
cfg.dBFS      = dBFS;

% Change the following names:
cfg.N_noise   = 2500;  % number of stimuli / condition
cfg.N_signal  = 4;     % Number of conditions
cfg.N         = cfg.N_noise*cfg.N_signal;

cfg_inout.dir_speech = dir_speech;
cfg_inout.dir_noise  = dir_noise;

% TODO: change this name:
cfg_inout.dir_stim  = dir_noise;
 
cfg_inout = Merge_structs(cfg,cfg_inout);
