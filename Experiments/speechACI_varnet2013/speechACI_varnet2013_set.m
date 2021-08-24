function cfg_inout = speechACI_varnet2013_set(cfg_inout)
% function cfg_out = speechACI_varnet2013_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
 
if ismac % lab's computer
    dir_data_experiment = fastACI_paths('dir_data'); % '/Users/leovarnet/ownCloud/Data/Projet fastACI/';
    warning('Leo: please copy what you have in %s into %sspeechACI_varnet2013%s (remove this warning once you have already done that)',dir_data_experiment,dir_data_experiment,filesep);
    
elseif isunix % Alejandro's computer
    dir_data_experiment = [fastACI_paths('dir_data') 'speechACI_varnet2013' filesep];
    % dir_data_experiment = [fastACI_paths('dir_data') cfg_inout.experiment_full filesep];
    
elseif ispc % Leo's computer
    dir_data_experiment = fastACI_paths('dir_data'); % 'C:\Users\Léo\ownCloud\Data\Projet fastACI/';
    warning('Leo: please copy what you have in %s into %sspeechACI_varnet2013%s (remove this warning once you have already done that)',dir_data_experiment,dir_data_experiment,filesep);
end

dir_target = [dir_data_experiment cfg_inout.Subject_ID filesep 'speech-samples' filesep];
if ~isfield(cfg_inout,'Condition')
    str2use  = ['Script2_Passation_EN(' cfg_inout.experiment ', ' cfg_inout.Subject_ID ', Condition);'];
    str2use2 = ['''white'' (default), ''pink'', or ''SSN'''];
    error('%s.m: Please enter one of the experimental conditions. %s\n\tCondition can be %s',upper(mfilename),str2use,str2use2);
    % dir_name_noise = 'NoiseStim';
else
    switch lower(cfg_inout.Condition) % lower case
        case 'white'
            dir_name_noise = 'NoiseStim';
        case 'pink'
            dir_name_noise = 'NoiseStim-pink';
        case 'ssn'
            dir_name_noise = 'NoiseStim-SSN';
        otherwise
            error('%s: Condition not recognised',upper(mfilename));
    end
end
dir_noise  = [dir_data_experiment cfg_inout.Subject_ID filesep dir_name_noise filesep];

dBFS = 100;

%%% Parameters to create the targets:
cfg.fs        = 44100; % Hz, sampling frequency
cfg.dur_ramp  = 75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
cfg.dBFS      = dBFS;

% Change the following names:
cfg.N_presentation = 2500; % 5000;  % number of stimuli / condition
cfg.N_target  = 2;     % Number of conditions
cfg.N         = cfg.N_target*cfg.N_presentation;

cfg_inout.dir_data_experiment = dir_data_experiment;

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

% TODO: change this name:
% cfg_inout.dir_stim  = dir_noise;
 
cfg_inout = Merge_structs(cfg,cfg_inout);
