function publ_osses2023b_FA_1_sim(dir_zenodo)
% function publ_osses2023b_FA_1_sim(dir_zenodo)
%
% 1. Description:
%      Script used to generate the simulations presented in publ_osses2023b
%
% See also: publ_osses2021c_DAGA_1_sim.m 
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelname = 'king2019';
Conditions = {'white'};
experiment = 'toneinnoise_ahumada1975';

if nargin == 0
    bReplicate = 1; % a new cfg_crea (and new waveforms), new template
    dir_zenodo = '';
else
    if exist(dir_zenodo,'dir')
        bReplicate = 0;
    end
end
bReproduce = ~bReplicate;

if bReplicate
    % Nothing to do, everything will be generated from the scratch
    fname_template_suffix = '';
end
if bReproduce
    % The simulations will be 'repeated', exactly as in the FA paper:
    dir_data_zenodo = [dir_zenodo '02-Raw-data' filesep 'fastACI' filesep 'Publications' filesep ...
        'publ_osses2023b' filesep 'data_king2019' filesep];
    
    % We need to copy the cfg_crea file to the local directory
    fname = 'cfgcrea_2023_04_19_00_07_king2019_toneinnoise_ahumada1975_white.mat';
    file_src = [dir_data_zenodo '0-init' filesep fname];
    
    dir_subj = [fastACI_dir_data experiment filesep];
    if ~exist(dir_subj,'dir'); mkdir(dir_subj); end
    
    dir_subj = [dir_subj modelname filesep];
    if ~exist(dir_subj,'dir'); mkdir(dir_subj); end
    
    dir_results = [dir_subj 'Results' filesep];
    if ~exist(dir_results,'dir'); mkdir(dir_results); end
    
    file_dst = [dir_results fname];
    copyfile(file_src, file_dst);
    
    fname_template_suffix = 'white-2023-4-19'; %keyvals.fname_template_suffix
    %%% Now, copying the template:
    fname = ['template-king2019-toneinnoise_ahumada1975-trial-1-' fname_template_suffix '.mat'];
    file_src = [dir_data_zenodo '2-template_matching' filesep fname];
    file_dst = [dir_subj fname];
    copyfile(file_src, file_dst);
    
    %%% Copy the waveforms (alternatively, create the wav files from the seeds):
    dir_src = [dir_zenodo '01-Stimuli' filesep 'fastACI_data' filesep experiment filesep ...
        modelname filesep 'NoiseStims-white' filesep];
    dir_dst = [dir_subj 'NoiseStims-white' filesep];
    copyfile(dir_src,dir_dst);
end

%%%
dir_here = [fastACI_basepath 'Local' filesep];
fname_cfg =  [dir_here modelname '_cfg.m'];

runs = {'run-1'};
Show_cell(runs);

% fprintf('Enter a number from the list before, i.e., between 1 and %.0f\n',length(runs));
% idx = input('Enter the name of the simulation to be run (you need all these runs to reproduce the figure papers).: ');
idx = 1;
run_str = runs{idx};

%%% New way (maybe add an option to check whether the file is up to date):
% Update_cfg_file(file_src, file_dst);
file_dst = fname_cfg;
dir_src = [fastACI_basepath 'Simulations' filesep 'Stored_cfg' filesep];
switch run_str
    case {'run-1'}
        file_src = [dir_src 'osses2023b_FA_king2019.m'];
end
if exist(file_dst,'file')
    fprintf('The configuration file %s exists already, do you want to overwrite it?\n',file_dst);
    fprintf('\tPress any button to overwrite it and continue, press ctrl+c to abort\n');
    pause;
end
copyfile(file_src,file_dst);

%%% Getting the extra parameters that are not explicit in the cfg files:
p = il_get_model_config_Forum_Acusticum(run_str, modelname);

if isfield(p,'thres_for_bias')
    thres_for_bias = p.thres_for_bias;
else
    thres_for_bias = 0;
end
if isfield(p,'in_std')
    in_std = p.in_std;
else
    in_std = 0;
end

for i = 1:length(Conditions)
    noise_type = Conditions{i};
    
    p = Get_date;
    if isempty(fname_template_suffix)
        % Then, using the default template suffix
        fname_template_suffix = [noise_type '-' p.date4files]; % trick to always get a new template
    end
    flags_here = {'thres_for_bias',thres_for_bias,'in_std',in_std,'fname_template_suffix',fname_template_suffix};
    
    data_passation.i_current = 1; % idle numbers
    cfg_game.N = 3200; % idle numbers

    while data_passation.i_current < cfg_game.N
        [cfg_game,data_passation] = fastACI_experiment(experiment,modelname,noise_type,flags_here{:});
    end
    
    folder_src = cfg_game.dir_results;
    folder_new = [fileparts(cfg_game.dir_results(1:end-1)) filesep 'Results-' run_str filesep];
    if exist(folder_new,'dir')
        folder_new = [fileparts(cfg_game.dir_results(1:end-1)) filesep 'Results-' run_str '-retest-on-' p.date4files filesep];
        
        if exist(folder_new,'dir')
            error('Rename manually the ''Results'' folder...')
        end
    end
    movefile(folder_src,folder_new);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = il_get_model_config_Forum_Acusticum(run_str,modelname)

p = [];
p.modelname = modelname;
p.bStore_template = 0; % New option added on 17/09/2021, not relevant for
                       % the simulations here
switch run_str
    case 'run-1'
        p.thres_for_bias = 0;
        p.in_std = 0;
        
    otherwise
        error('Run condition not recognised')
end

if ~isfield(p,'thres_for_bias')
    p.thres_for_bias = nan;
end
