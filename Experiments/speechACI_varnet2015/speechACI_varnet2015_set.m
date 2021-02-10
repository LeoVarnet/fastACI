function cfg_out = speechACI_varnet2015_set(cfg_in)
% function cfg_out = speechACI_varnet2015_set(cfg_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end
 
cfg_out = cfg_in; % copying input to output struct
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
cfg.N_noise   = 1500;  % number of stimuli / condition
cfg.N_signal  = 2;     % Number of conditions
cfg.N         = cfg.N_noise*cfg.N_signal;

cfg_out.dir_speech = dir_speech;
cfg_out.dir_noise  = dir_noise;

% TODO: change this name:
cfg_out.dir_stim  = dir_noise;
 
cfg_out = Merge_structs(cfg,cfg_out);
