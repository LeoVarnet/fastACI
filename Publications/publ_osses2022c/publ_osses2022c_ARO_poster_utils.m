function outs = publ_osses2022c_ARO_poster_utils(Subject_ID, noise_type, type_action)
% function outs = publ_osses2022c_ARO_poster_utils(Subject_ID, noise_type, type_action)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outs = [];

dir_target = [];
dir_noise  = [];

switch type_action
    case 'Get_filenames'
        experiment = 'speechACI_Logatome-abda-S43M';
        bGenerate_stimuli = 0;
        
        dir_res = [fastACI_dir_data filesep experiment filesep Subject_ID filesep 'Results' filesep];
        if ~exist(dir_res,'dir')
            bGenerate_stimuli = 1;    
        else
            file = Get_filenames(dir_res,['cfgcrea*' noise_type '.mat']);
            if isempty(file)
                bGenerate_stimuli = 1;
            end
        end
        
        if bGenerate_stimuli == 0 % i.e., if dir_res exists
            file = Get_filenames(dir_res,['cfgcrea*' noise_type '.mat']);
            if length(file) ~= 1
                error('More than one creafile was found on disk');
            end
            load([dir_res file{1}],'cfg_crea');

            cfg_crea = Check_cfg_crea_dirs(cfg_crea);

            dir_target = cfg_crea.dir_target;
            dir_noise  = cfg_crea.dir_noise;

            if ~exist(dir_target,'dir')
                error('Directory %s not found on disk',dir_target);
            end
            if ~exist(dir_noise,'dir')
                bGenerate_stimuli = 1;
                warning('Directory %s not found on disk',dir_noise);
            end

            switch Subject_ID
                case 'SLV'
                    switch noise_type
                        case 'white'
                        case 'sMPSv1p3'
                        case 'bumpv1p2_10dB'
                    end

                case 'SAO'
                    switch noise_type
                        case 'white'
                        case 'sMPSv1p3'
                        case 'bumpv1p2_10dB'
                    end
            end
            outs.dir_target = dir_target;
            outs.dir_noise  = dir_noise;
        end
        outs.bGenerate_stimuli = bGenerate_stimuli;
        % outs.fname_results = fname_results;
        
    case 'Get_flags'
        % N = 5000;
        % t_limits = [0 0.3425]; % s, excluding sound information after this timing
        % 
        % glmfct  = 'lassoglmslow'; 
        % TF_type = 'gammatone';
        % 
        % switch glmfct
        %     case 'lassoglmslow'
        %         N_lambda = 30;
        %         Lambdas = logspace(-4, -1, N_lambda);
        %         idx = find(Lambdas >= 10^-3);
        %         Lambdas = Lambdas(idx);
        % end
        % flags_for_input = {TF_type, ...
        %            glmfct, ... % 'dir_noise',dir_noise, 'dir_target',dir_target, ...
        %            'trialtype_analysis', 'total', ...
        %            'idx_trialselect', 1:N, ...
        %            'N_folds', 10, ... % 'dir_noise',cfg_game.dir_noise, 'dir_target',cfg_game.dir_target, ...
        %            'add_signal',0, ...
        %            'apply_SNR',0, ...
        %            'skip_if_on_disk',1, ...
        %            'expvar_after_reversal', 4, ...
        %            'lambda', Lambdas, ...
        %            'no_permutation', 'no_bias', ...
        %            't_limits', t_limits}; 
        % outs.flags_for_input = flags_for_input;
        % outs.glmfct = glmfct;
        % outs.TF_type = TF_type;
        
    case 'Model-calibrate'
        % % publ_osses2022c_utils('osses2022a','SSN','Model-calibrate')
        % 
        % % Completed on 2/02/2022 at 11:59...
        % model = 'osses2022a';
        % experiment = 'speechACI_varnet2013';
        % g20211207_calibrating_the_model(experiment,model);
        
    case 'Model-run'
        % % publ_osses2022c_utils('osses2022a','white','Model-run')
        % 
        % model = 'osses2022a';
        % if ~strcmp(Subject_ID,model)
        %     warning('Ignoring Subject_ID and using the model name...');
        % end
        % 
        % file_model_decision_config = [fastACI_basepath 'Interim_results' filesep 'osses2022a-optimal_detector-20220202-' noise_type '.mat'];
        % flags_here = {'file_model_decision_config',file_model_decision_config};
        % experiment = 'speechACI_varnet2013';
        % [cfg_game, data_passation] = fastACI_experiment(experiment,model,noise_type,flags_here{:});
        % disp('')
        
end

