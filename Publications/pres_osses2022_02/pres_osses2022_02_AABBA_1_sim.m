function pres_osses2022_02_AABBA_1_sim
% function pres_osses2022_02_AABBA_1_sim
%
% See also: g20211201_analysis_pipeline_simulations_Alejandro
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

available_models = {'dau1997','king2019','relanoiborra2019','osses2021','osses2022a'};
Show_cell(available_models);

bInput = input('Choose the model by entering the number from the list above: ');
modelname = available_models{bInput};

Conditions = {'white'};
experiment = 'speechACI_Logatome-abda-S43M';
noise_type = Conditions{1};

%%%
dir_here = [fastACI_basepath 'Simulations' filesep];
fname_cfg =  [dir_here modelname '_cfg.m'];

%%% New way (maybe add an option to check whether the file is up to date):
% Update_cfg_file(file_src, file_dst);
file_dst = fname_cfg;
dir_src = [fastACI_basepath 'Simulations' filesep 'Stored_cfg' filesep];

pref = 'osses2022_02_AABBA_';
switch modelname
    case 'dau1997'
        file_src = [dir_src pref 'dau1997.m'];
        
    case 'king2019'
        file_src = [dir_src pref 'king2019.m'];
        
    case 'relanoiborra2019'
        file_src = [dir_src pref 'relanoiborra.m'];
        
    case 'osses2021'
        file_src = [dir_src pref 'osses2021.m'];
        
    case 'osses2022a'
        file_src = [dir_src pref 'osses2022a.m'];
end
if exist(file_dst,'file')
    fprintf('The configuration file %s exists already, do you want to overwrite it?\n',file_dst);
    fprintf('\tPress any button to overwrite it and continue, press ctrl+c to abort\n');
    pause;
end
copyfile(file_src,file_dst);

%%% Getting the extra parameters that are not explicit in the cfg files:
in_std = 0;

% fname_template_suffix = [noise_type '-10-rep-202201'];
fname_template_suffix = '-10-rep-202201';
flags_here = {'thres_for_bias',[],'in_std',in_std,'fname_template_suffix',fname_template_suffix};

bCalibration = 1;
if bCalibration
    %%% Runs the calibration first:
    [~,~,kv] = fastACI_experiment(experiment,modelname,noise_type,flags_here{:});
    Play_ready;

    thres_for_bias = kv.thres_for_bias;
    fprintf('%s: Simulations completed using thres_for_bias=%.4f, in_std=%.4f\n',upper(mfilename),thres_for_bias,in_std);
    %%%
end
flags_here = {'thres_for_bias',thres_for_bias,'in_std',in_std,'fname_template_suffix',fname_template_suffix};

cfg_game.N = 5000; % idle numbers, overwritten after the first pass of the while loop, below
data_passation.i_current = 1; % idle numbers
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










p = il_get_model_config_DAGA(run_str, modelname);

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
    fname_template_suffix = [noise_type '-' p.date4files]; % trick to always get a new template
    flags_here = {'thres_for_bias',thres_for_bias,'in_std',in_std,'fname_template_suffix',fname_template_suffix};
    
    data_passation.i_current = 1; % idle numbers
    cfg_game.N = 5000; % idle numbers

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
function p = il_get_model_config_DAGA(run_str,modelname)

switch modelname
    case 'osses2021'
        % Nothing to do
    otherwise
        error('This publication only used the model ''osses2021''')
end
p = [];
p.modelname = modelname;
p.bStore_template = 0; % New option added on 17/09/2021, not relevant for
                       % the simulations here
switch run_str
    case 'run-1'
        p.thres_for_bias = 0;
        p.in_std = 0;
        
    case 'run-3-m1p55'
        p.thres_for_bias = -1.55;
        p.in_std = 0;
        
    case 'run-3-p0p39'
        p.thres_for_bias = 0.39;
        p.in_std = 0;
        
    case 'run-3-p0p78'
        p.thres_for_bias = 0.78;
        p.in_std = 0;
        
    case 'run-4'
        p.in_std = 0;
        
    otherwise
        error('Run condition not recognised')
end

if ~isfield(p,'thres_for_bias')
    p.thres_for_bias = nan;
end
