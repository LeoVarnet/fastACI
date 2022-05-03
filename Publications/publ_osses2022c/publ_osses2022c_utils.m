function outs = publ_osses2022c_utils(Subject_ID, noise_type, type_action)
% function outs = publ_osses2022c_utils(Subject_ID, noise_type, type_action)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outs = [];

dir_target = [];
dir_noise  = [];

switch type_action
    case 'Get_filenames'
        [~,~,experiment] = publ_osses2021c_DAGA_0_checkdata({'osses2021c_S02'}); % IDLE
        switch Subject_ID
            case 'SLV'
                switch noise_type
                    case 'white'
                    %%% varnet2013_S02 (Sujet_Leo), loading folder, updating the target/noise folders:
                      % This is the old data of Leo, contrained to 5000 trials.    
                      [~,dir_subj] = publ_varnet2013_FrontHumNeurosci_0_checkdata;
                      dir_target = [dir_subj 'ListeSignal' filesep];
                      dir_noise  = [dir_subj 'ListeBruit' filesep];
                      fname_results = [dir_subj 'savegame_final.mat'];
                                            
                    case 'SSN'
                      % dir_target = outs.dir_target;
                      
                      dir_subj = [fastACI_dir_data experiment filesep 'osses2021c_S01' filesep];
                      dir_res = [dir_subj 'Results' filesep];
                      
                      fname_results = Get_filenames(dir_res,['savegame*' noise_type '.mat']);
                      if length(fname_results) ~= 1
                          error('More than one file was found...')
                      end
                      fname_results = [dir_res fname_results{1}];
                      dir_target = [dir_subj 'speech-samples' filesep];
                      dir_noise  = [dir_subj 'NoiseStim-SSN' filesep];
                      
                 end
            case 'SAO'
                switch noise_type
                    case 'white'
                      dir_subj = [fastACI_dir_data experiment filesep 'SAO-5000-trials' filesep];  
                      dir_res = [dir_subj 'Results' filesep];
                      filter_here = experiment; % no label for the white noise (avoiding to spot the SSN noise)
                      
                      fname_results = Get_filenames(dir_res,['savegame*' filter_here '.mat']);
                      if length(fname_results) ~= 1 
                          error('More than one file was found...')
                      end
                      fname_results = [dir_res fname_results{1}];
                      dir_target = [dir_subj 'speech-samples' filesep];
                      dir_noise  = [dir_subj 'NoiseStim' filesep];

                    case 'SSN'
                      subject_here = 'osses2021c_S02';
                      [bContinue,dir_subj] = publ_osses2021c_DAGA_0_checkdata({subject_here});
                      dir_subj = dir_subj{1};
                      
                      dir_res = [fastACI_dir_data experiment filesep subject_here filesep 'Results' filesep];
                      if bContinue == 1
                          fname_results = Get_filenames(dir_res,['savegame*' noise_type '.mat']);
                          if length(fname_results) ~= 1
                              error('More than one file was found...')
                          end
                          fname_results = [dir_res fname_results{1}];
                          dir_target = [dir_subj 'speech-samples' filesep];
                          dir_noise  = [dir_subj 'NoiseStim-SSN' filesep];
                      end
                end
        end
        outs.dir_target = dir_target;
        outs.dir_noise  = dir_noise;
        outs.fname_results = fname_results;
        
    case 'Get_flags'
        N = 5000;
        t_limits = [0 0.3425]; % s, excluding sound information after this timing

        glmfct  = 'l1glm'; % 'lassoglmslow'; 
        TF_type = 'gammatone';

        switch glmfct
            case 'l1glm' % 'lassoglmslow'
                N_lambda = 30;
                Lambdas = logspace(-4, -1, N_lambda);
                idx = find(Lambdas >= 10^-3);
                Lambdas = Lambdas(idx);
        end
        flags_for_input = {TF_type, ...
                   glmfct, ... % 'dir_noise',dir_noise, 'dir_target',dir_target, ...
                   'trialtype_analysis', 'total', ...
                   'idx_trialselect', 1:N, ...
                   'N_folds', 10, ... % 'dir_noise',cfg_game.dir_noise, 'dir_target',cfg_game.dir_target, ...
                   'add_signal',0, ...
                   'apply_SNR',0, ...
                   'skip_if_on_disk',1, ...
                   'expvar_after_reversal', 4, ...
                   'lambda', Lambdas, ...
                   'no_permutation', 'no_bias', ...
                   't_limits', t_limits, ...
                   'pyramid_script','imresize', ...
                   'pyramid_shape' ,0}; 
               
        outs.flags_for_input = flags_for_input;
        outs.glmfct = glmfct;
        outs.TF_type = TF_type;
        
    case 'Model-calibrate'
        % publ_osses2022c_utils('osses2022a','SSN','Model-calibrate')
        
        % Completed on 2/02/2022 at 11:59...
        model = 'osses2022a';
        experiment = 'speechACI_varnet2013';
        g20211207_calibrating_the_model(experiment,model);
        
    case 'Model-run'
        % publ_osses2022c_utils('osses2022a','white','Model-run')
        
        model = 'osses2022a';
        if ~strcmp(Subject_ID,model)
            warning('Ignoring Subject_ID and using the model name...');
        end
        
        file_model_decision_config = [fastACI_basepath 'Interim_results' filesep 'osses2022a-optimal_detector-20220202-' noise_type '.mat'];
        flags_here = {'file_model_decision_config',file_model_decision_config};
        experiment = 'speechACI_varnet2013';
        [cfg_game, data_passation] = fastACI_experiment(experiment,model,noise_type,flags_here{:});
        disp('')
        
end

