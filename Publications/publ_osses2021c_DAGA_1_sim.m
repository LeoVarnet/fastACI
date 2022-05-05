function publ_osses2021c_DAGA_1_sim
% function publ_osses2021c_DAGA_1_sim
%
% See also: g20210730_all_simulation_sessions
% Original name: g20210811_all_simulation_sessions
% See also: g20220414_simulations_Alejandro (new script for simulations)
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelname = 'osses2021'; % modelname = 'king2019';
Conditions = {'SSN'};%,'white','pink'};
experiment = 'speechACI_varnet2013'; % experiment = 'speechACI_Logatome-abda-S41F';

%%%
dir_here = [fastACI_basepath 'Simulations' filesep];
fname_cfg =  [dir_here modelname '_cfg.m'];

runs = {'run-1'; 'run-3-m1p55'; 'run-3-p0p39'; 'run-3-p0p78'; 'run-4'};
Show_cell(runs);

fprintf('Enter a number from the list before, i.e., between 1 and %.0f\n',length(runs));
idx = input('Enter the name of the simulation to be run (you need all these runs to reproduce the figure papers).: ');
run_str = runs{idx};

% %%% OLD WAY:
% templ_num = 10;
% p = model_cfg_osses2021c(run_str,modelname,templ_num);
% text_to_write = readfile_replace('model_cfg_replace.txt',p);
% 
% if exist(fname_cfg,'file')
%     fprintf('----------------------------------------------------------------------------\n')
%     fprintf('file %s exists, \npress any key to continue (will overwrite) or press ctrl+C to cancel \n',fname_cfg);
%     fprintf('----------------------------------------------------------------------------\n')
%     pause
% end
% 
% fid = fopen(fname_cfg, 'w');
% fwrite(fid, text_to_write);
% fclose(fid);
% %%%

%%% New way (maybe add an option to check whether the file is up to date):
% Update_cfg_file(file_src, file_dst);
file_dst = fname_cfg;
dir_src = [fastACI_basepath 'Simulations' filesep 'Stored_cfg' filesep];
switch run_str
    case {'run-1','run-3-m1p55','run-3-p0p39','run-3-p0p78'}
        file_src = [dir_src 'osses2021c_D1.m'];
    case 'run-4'
        file_src = [dir_src 'osses2021c_D2.m'];
end
if exist(file_dst,'file')
    fprintf('The configuration file %s exists already, do you want to overwrite it?\n',file_dst);
    fprintf('\tPress any button to overwrite it and continue, press ctrl+c to abort\n');
    pause;
end
copyfile(file_src,file_dst);

%%% Getting the extra parameters that are not explicit in the cfg files:
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
        movefile(folder_src,folder_new);
    end
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
        % p.decision_script = 'aci_detect';
        % p.template_script = 'model_template';
        % p.template_every_trial = 0;
        % p.templ_num = 10;
        % p.det_lev = -6;
        % p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0;
        p.in_std = 0;
        
    case 'run-3-m1p55'
        % p.decision_script = 'aci_detect';
        % p.template_script = 'model_template';
        % p.template_every_trial = 0;
        % p.templ_num = 10;
        % p.det_lev = -6;
        % p.type_decision = 'optimal_detector';
        p.thres_for_bias = -1.55;
        p.in_std = 0;
        
    case 'run-3-p0p39'
        % p.decision_script = 'aci_detect';
        % p.template_script = 'model_template';
        % p.template_every_trial = 0;
        % p.templ_num = 10;
        % p.det_lev = -6;
        % p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0.39;
        p.in_std = 0;
        
    case 'run-3-p0p78'
        % p.decision_script = 'aci_detect';
        % p.template_script = 'model_template';
        % p.template_every_trial = 0;
        % p.templ_num = 10;
        % p.det_lev = -6;
        % p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0.78;
        p.in_std = 0;
        
    case 'run-4'
        % p.decision_script = 'aci_detect';
        % p.template_script = 'model_template';
        % p.template_every_trial = 0;
        % p.templ_num = 10;
        % p.det_lev = -6;
        % p.type_decision = 'relanoiborra2019_decision';
        p.in_std = 0;
        
    otherwise
        error('Run condition not recognised')
end

if ~isfield(p,'thres_for_bias')
    p.thres_for_bias = nan;
end
