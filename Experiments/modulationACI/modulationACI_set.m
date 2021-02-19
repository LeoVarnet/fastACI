function cfg_inout = modulationACI_set(cfg_inout)
% function cfg_inout = modulationACI_set(cfg_inout)
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
cfg.stim_dur  = 1.5;   % s,  stimulus suration (s)
cfg.SNR       = -10;

% cfg.expvar_start = -8;  % old name 'm_start'
% cfg.expvar_description = 'modulation depth (dB)';

cfg.fadein_s  = 0.075; % CAUTION: Overwritten in the case of simulation
cfg.SPL       = 65;
cfg.dBFS      = 93.61; % Calibration with Headphones HD 650

cfg.N_noise    = 1500;  % number of stimuli / condition
cfg.N_signal   = 2;     % Number of conditions
cfg.noise_type = 'white';

if ~isfield(cfg_inout,'experiment')
    fname = strsplit(mfilename,'_');
    fname = fname{1};
    cfg_inout.experiment = fname;
end

dir_stim = ['/home/alejandro/Documents/Databases/data/fastACI/' cfg_inout.experiment filesep];
if ~isdir(dir_stim)
    error('%s: The specified stimulus directory does not exist, please indicate an existing directory...')
end
cfg.folder_name = ['NoiseStims-' cfg.noise_type]; % nom du dossier a creer contenant les noises

cfg_inout.dir_stim = dir_stim;

cfg_inout = Merge_structs(cfg,cfg_inout);