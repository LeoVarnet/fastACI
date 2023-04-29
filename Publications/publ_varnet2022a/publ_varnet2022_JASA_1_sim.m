function publ_varnet2022_JASA_1_sim
% function publ_varnet2022_JASA_1_sim
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelname = 'king2019'; 
Conditions = {'white'};
experiment = 'modulationACI'; 

%%%
dir_here = [fastACI_basepath 'Simulations' filesep];

str_in.config = 'varnet2022_XC-D';
Model_configuration_selector(modelname,str_in);

fastACI_experiment(experiment,modelname,experiment);

% disp('')
% 
% fname_cfg =  [dir_here modelname '_cfg.m'];
% 
% runs = {'run-1'; 'run-3-m1p55'; 'run-3-p0p39'; 'run-3-p0p78'; 'run-4'};
% Show_cell(runs);
% 
% fprintf('Enter a number from the list before, i.e., between 1 and %.0f\n',length(runs));
% idx = input('Enter the name of the simulation to be run (you need all these runs to reproduce the figure papers).: ');
% run_str = runs{idx};
% 
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
% 
% for i = 1:length(Conditions)
%     noise_type = Conditions{i};
%     
%     data_passation.i_current = 1; % idle numbers
%     cfg_game.N = 5000; % idle numbers
% 
%     while data_passation.i_current < cfg_game.N
%         [cfg_game,data_passation] = fastACI_experiment(experiment,modelname,noise_type);
%     end
%     
%     folder_src = cfg_game.dir_results;
%     folder_new = [fileparts(cfg_game.dir_results(1:end-1)) filesep 'Results-' run_str filesep];
%     if exist(folder_new,'dir')
%         folder_new = [fileparts(cfg_game.dir_results(1:end-1)) filesep 'Results-RETEST-' run_str filesep];
%         
%         if exist(folder_new,'dir')
%             error('Rename manually the ''Results'' folder...')
%         end
%         movefile(folder_src,folder_new);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function p = il_get_model_config_DAGA(run_str,modelname)
% 
% p = [];
% p.modelname = modelname;
% p.bStore_template = 0; % New option added on 17/09/2021, not relevant for
%                        % the simulations here
% switch run_str
%     case 'run-1'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.templ_num = 10;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0;
%         
%     case 'run-3-m1p55'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.templ_num = 10;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = -1.55;
%         
%     case 'run-3-p0p39'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.templ_num = 10;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0.39;
%         
%     case 'run-3-p0p78'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.templ_num = 10;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0.78;
%         
%     case 'run-4'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.templ_num = 10;
%         p.det_lev = -6;
%         p.type_decision = 'relanoiborra2019_decision';
%         
%     otherwise
%         error('Run condition not recognised')
% end
% 
% if ~isfield(p,'thres_for_bias')
%     p.thres_for_bias = nan;
% end
