function publ_osses2023b_Forum_Acusticum_1_sim
% function publ_osses2023b_Forum_Acusticum_1_sim
%
% 1. Description:
%      Script used to generate the simulations presented in publ_osses2023b
%
% See also: publ_osses2021c_DAGA_1_sim.m 
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelname = 'king2019';
Conditions = {'white'};
% experiment = 'toneinnoise_ahumada1975';
experiment = 'toneinnoise_ahumada1971';

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
        switch experiment
            case 'toneinnoise_ahumada1975'
                file_src = [dir_src 'king2019a_osses2023b_Forum_Acusticum.m'];
            case 'toneinnoise_ahumada1971'
                file_src = [dir_src 'king2019a_osses2023b_Forum_Acusticum_detlev10.m'];
        end
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
    fname_template_suffix = [noise_type '-' p.date4files]; % trick to always get a new template
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
