function cfg_out = modulationACI_set(cfg_in)
% function str_stim = modulationACI_set(cfg_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end
% function cfg_out = il_set(cfg_in)
%
% Function comparable to *_set.m functions from AFC toolbox

cfg = [];
cfg_out = cfg_in; % copying input to output struct

%%% Parameters to create the targets:
cfg.fs        = 48000; % Hz, sampling frequency
cfg.fm        = 4;
cfg.fc        = 1000;
cfg.stim_dur  = 1.5;   % s,  stimulus suration (s)
cfg.SNR       = -10;
cfg.m_start   = -8; 
cfg.fadein_s  = 0.075; % CAUTION: Overwritten in the case of simulation
cfg.SPL       = 65;
cfg.dBFS      = 93.61; % Calibration with Headphones HD 650

cfg.N_noise         = 1500;  % number of stimuli / condition
cfg.N_signal        = 2;     % Number of conditions
cfg.noise_type      = 'white';

if ~isfield(cfg_in,'experiment')
    fname = strsplit(mfilename,'_');
    fname = fname{1};
    cfg_in.experiment = fname;
end

dir_stim = ['/home/alejandro/Documents/Databases/data/fastACI/' cfg_in.experiment filesep];
if ~isdir(dir_stim)
    error('%s: The specified stimulus directory does not exist, please indicate an existing directory...')
end
cfg.folder_name = ['NoiseStims-' cfg.noise_type]; % nom du dossier a creer contenant les noises

cfg_out.dir_stim = dir_stim;

cfg_out = Merge_structs(cfg,cfg_out);