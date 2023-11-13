function cfg_inout = speechACI_Audika_set(cfg_inout)
% function cfg_out = speechACI_Audika_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_vars

if nargin == 0
    cfg_inout = [];
end

dir_data_experiment = [fastACI_paths('dir_data') cfg_inout.experiment_full filesep];

% dir_logatome_src = '/media/alejandro/My Passport/Databases/data-Speech/french/Logatome/'; % Logatome-average-power-speaker-S46M_FR.mat
dir_target = [dir_data_experiment cfg_inout.Subject_ID filesep 'speech-samples' filesep];

%%% First checking point:
if ~isfield(cfg_inout,'Cond_extra_1')
    str2use  = ['Script2_Passation_EN(' cfg_inout.experiment '-abda, ' cfg_inout.Subject_ID ', Condition);'];
    error('%s.m: You need to specify the test speech samples (''abda'' for ab-ba vs ad-da; or ''apta'' for ap-pa vs at-ta), for instance:\n\t %s',upper(mfilename),str2use);
end

if ~isfield(cfg_inout,'Cond_extra_2')
    str2use  = ['Script2_Passation_EN(' cfg_inout.experiment_full '-abda-S41F, ' cfg_inout.Subject_ID ', Condition);'];
    error('%s.m: You need to specify the speaker (''S41F'' for abda; or ''S46M'' for apta), for instance:\n\t %s',upper(mfilename),str2use);
end

%%% Second checking point:
if ~isfield(cfg_inout,'Condition')
    str2use  = ['Script2_Passation_EN(' cfg_inout.experiment_full ', ' cfg_inout.Subject_ID ', Condition);'];
    str2use2 = '''white'' (default), ''pink'', or ''SSN''';
    error('%s.m: Please enter one of the experimental conditions. %s\n\tCondition can be %s',upper(mfilename),str2use,str2use2);
end

switch lower(cfg_inout.Condition) % lower case
    case 'white'
        dir_name_noise = 'NoiseStim-white';
        noise_type = 'white';
    case 'bumpv1p2_10db'
        dir_name_noise = 'NoiseStim-bumpv1p2_10dB';
        noise_type = 'bumpv1p2_10dB';
    case 'smpsv1p3'
        dir_name_noise = 'NoiseStim-sMPSv1p3';
        noise_type = 'sMPSv1p3';
    otherwise
        error('Noise type not recognised for this experiment...')
end
dir_noise  = [dir_data_experiment cfg_inout.Subject_ID filesep dir_name_noise filesep];

%%% Parameters to create the targets:
cfg.fs        = 16000; % Hz, sampling frequency
cfg.dur_ramp  = 75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
%%%
dBFS = 100; % full scale value to store the waveforms
cfg.dBFS = dBFS;
%%%

cfg.bRove_level = 1; % New option as of 16/04/2021
cfg.Rove_range  = 2.5; % plus/minus this value, changed from 4 to 2.5 dB on 26/05/2021

if isfield(global_vars,'Language')
    Language = global_vars.Language;
else
    Language = 'FR'; % French by default
end
cfg.Language = Language; % or 'EN'

if isfield(global_vars,'N_presentation')
    N_presentation = global_vars.N_presentation;
else
    if ~isfield(cfg_inout,'N_presentation')
        N_presentation = 2000; % before 14/04/2022, this value was 2500
    else
        N_presentation = cfg_inout.N_presentation;
    end
end

% if isfield(global_vars,'dBFS')
%     cfg_inout = Ensure_field(cfg_inout,'is_experiment',1);
%     if cfg_inout.is_experiment
%         error('AO on 9/02/2023: Validate this further')
%         cfg_inout.dBFS = global_vars.dBFS;
%         cfg_inout.dBFS_playback = cfg_inout.dBFS;
%     end
% end

cfg.N_presentation = 10%N_presentation; % number of stimuli / condition
cfg.N_target  = 2;     % Number of conditions
cfg.N         = cfg.N_target*cfg.N_presentation;

cfg_inout.dir_data_experiment = dir_data_experiment;

% Change the following names:
if isfield(cfg_inout,'dir_target')
    if ~exist(cfg_inout.dir_target,'dir')
        if iswindows
            warning('Leo: I had to re-enable this option here to have the possibility to use noises from other subjects, but I remember that this was an error source for you...')
            disp('Pausing for 10 s... (this is a temporal message)')
            pause(10)
        end
        cfg_inout.dir_target = dir_target;
    end
else
    cfg_inout.dir_target = dir_target;
end

if isfield(cfg_inout,'dir_noise')
    if ~exist(cfg_inout.dir_noise,'dir')
        cfg_inout.dir_noise  = dir_noise;
    end
else
    cfg_inout.dir_noise  = dir_noise;
end
cfg_inout.noise_type = noise_type;
 
cfg_inout = Merge_structs(cfg,cfg_inout);
