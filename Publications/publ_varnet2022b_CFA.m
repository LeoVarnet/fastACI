function [h,hname] = publ_varnet2022b_CFA(varargin)
% function [h,hname] = publ_varnet2022b_CFA(varargin)
%
% Generates the figures
%
% % To display Fig. 2 of Varnet, Lorenzi, Osses (2022, CFA) use :::
%     publ_varnet2022b_CFA('fig2');
%
% dir_data = fastACI_dir_data; % dir_data of this computer
% flags = {'dir_data',dir_data};
% publ_varnet2022b_CFA('fig2_mod22','SA',flags{:});
% publ_varnet2022b_CFA('fig2_mod22','SB',flags{:});
%
% publ_varnet2022b_CFA('fig2_abda13','SA',flags{:}); 
% publ_varnet2022b_CFA('fig2_abda13','SB',flags{:}); 
%
% publ_varnet2022b_CFA('fig2_abda21','SA',flags{:}); % DAGA
% publ_varnet2022b_CFA('fig2_abda21','SB',flags{:}); % DAGA
%
% publ_varnet2022b_CFA('fig2_abda22','SA',flags{:}); % ARO poster
% publ_varnet2022b_CFA('fig2_abda22','SB',flags{:}); % ARO poster
%
% Author: Alejandro Osses and Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_varnet2022b_CFA;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2_abda13','fig2_abda21','fig2_abda22','fig2_mod22'};
definput.flags.subject = {'SA','SB'};
definput.flags.publ = {'varnet2022b_CFA','varnet2022a_JASA'}; % default is varnet2022b_CFA
definput.keyvals.dir_out  = [];
definput.keyvals.dir_data = fastACI_paths('dir_data');
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

dir_data   = keyvals.dir_data;
dir_output = keyvals.dir_out;
dir_out_figs = [fastACI_paths('dir_output') 'varnet2022b_CFA' filesep];
if ~exist(dir_out_figs,'dir')
    mkdir(dir_out_figs);
end

%%%
flags_for_input = {'trialtype_analysis','total', ...'t1',...
    'N_folds', 10, ...
    'no_permutation', ...
    'no_bias', ...
    'no_plot', ...
    'pyramid_script','imresize', ...
    'pyramid_shape',0 ...
    };
DimCI = 'gammatone'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the participant ID:
if flags.do_fig2_abda13
    experiment = 'speechACI_varnet2013';
    exp_label = 'abda13-white';
    dir_expe = [dir_data experiment filesep];
    masker = 'white';
    
    if flags.do_SA
        participants = {'SLeo'}; 
    end
    if flags.do_SB
        if exist([dir_expe 'SAO' filesep],'dir')
            participants = {'SAO'}; % name after renaming
        else
            participants = {'SAO-5000-trials'}; % original name, as collected
        end        
    end
end
%%%
if flags.do_fig2_abda21
    experiment = 'speechACI_varnet2013';
    exp_label = 'abda21-SSN';
    dir_expe = [dir_data experiment filesep];
    masker = 'SSN';
    
    if flags.do_SA
        if exist([dir_expe 'S01' filesep],'dir')
            participants = {'S01'}; % as stored locally by Leo
        else
            participants = {'osses2021c_S01'}; % as downloaded from Zenodo
        end
    end
    if flags.do_SB
        if exist([dir_expe 'S02' filesep],'dir')
            participants = {'S02'}; % as stored locally by Leo
        else
            participants = {'osses2021c_S02'}; % as downloaded from Zenodo
        end
    end
end
%%%
if flags.do_fig2_abda22
    experiment = 'speechACI_Logatome-abda-S43M';
    exp_label = 'abda22-white';
    dir_expe = [dir_data experiment filesep];
    masker = 'white';
    
    if flags.do_SA
        participants = {'SLV'}; 
    end
    if flags.do_SB
        participants = {'SAO'};
    end
end
%%%
if flags.do_fig2_mod22
    experiment = 'modulationACI';
    exp_label = 'mod22';
    if flags.do_SA
        participants = {'S4'}; % S_LV2
    end
    if flags.do_SB
        participants = {'S1'}; % S_AO
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_subjects = length(participants);

for i_subject = 1:N_subjects
    % dir_local = [dir_data files{i_subject} filesep];    
    subj_here = participants{i_subject};
    dir_subj = [dir_data experiment filesep subj_here filesep];
    
    if flags.do_fig2_abda13 || flags.do_fig2_abda21 || flags.do_fig2_abda22
        cfg_game.experiment = experiment;
        cfg_game.dir_target = [dir_subj 'speech-samples' filesep];
    end
    
    if flags.do_fig2_abda13
        cfg_game.dir_noise  = [dir_subj 'NoiseStim'      filesep];
    end
    if flags.do_fig2_abda13 || flags.do_fig2_abda21 || flags.do_fig2_abda22
        dir_results = [dir_subj 'Results' filesep];

        files = Get_filenames(dir_results,'savegame_*.mat');
        if length(files) > 1
            Show_cell(files)
            fprintf('Looking for %s noises\n',masker);
            bInput = input(['Choose the one that you want to choose (expected=' masker '): ']);
            files = files(bInput);
        end
        fname_results = [dir_results files{1}];

        [cfg_game, data_passation, ListStim] = Convert_ACI_data_type(fname_results);
        if ~strcmp(cfg_game.Subject_ID,subj_here)
            cfg_game.Subject_ID = subj_here;
        end
    
        cfg_game = Check_cfg_crea_dirs(cfg_game,dir_data);
    end
    if flags.do_fig2_abda13
        if ~exist(cfg_game.dir_noise,'dir')
            switch cfg_game.Subject_ID
                case 'SAO-5000-trials'
                    cfg_game.dir_target = [dir_data experiment filesep subj_here filesep 'speech-samples' filesep];
                    cfg_game.dir_noise  = [dir_data experiment filesep subj_here filesep 'NoiseStim'      filesep];
            end
        end
    end
    if flags.do_fig2_abda13 || flags.do_fig2_abda21 || flags.do_fig2_abda22    
        dir_noise = cfg_game.dir_noise;
        dir_target = cfg_game.dir_target;

        flags_for_input{end+1} = 'dir_noise';
        flags_for_input{end+1} = dir_noise;
        flags_for_input{end+1} = 'dir_target';
        flags_for_input{end+1} = dir_target;

        %% compute ACI
        Data_matrix = [];

        %%% Put fname_results as an absolute path, it's easier:
        type = 'l1glm';% 'classic_revcorr';%
        switch type
            case 'l1glm' % 'lassoglmslow'
                N_lambda = 30;
                Lambdas = logspace(-4, -1, N_lambda);
                idx = find(Lambdas >= 10^-3);
                Lambdas = Lambdas(idx);

                flags_for_input{end+1} = 'lambda';
                flags_for_input{end+1} = Lambdas;
        end
        flags_for_input{end+1} = 'dir_out';
        flags_for_input{end+1} = dir_output;
        
        % type = 'lassoglm';
        [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, DimCI, ...
            type, flags_for_input{:},'Data_matrix',Data_matrix);
    end
    
    if flags.do_fig2_mod22
        data = publ_varnet2022a_utils(subj_here,'Get_CIt',flags,keyvals);
        f = data.f;
        t = data.t;
        
        results = data.results;
        cfg_ACI = data.cfg_ACI;
    end
    %% display ACI
    figure('Position', [100 100 400 300])
    affichage_tf(squeeze(results.ACI),'CI', 'cfg', cfg_ACI); hold on
    
    if flags.do_fig2_mod22
        XT  = get(gca,'XTick');
        XTL = round(100*t(XT))/100;
        set(gca,'XTickLabel',XTL);
        xlabel('Time (s)');
            
        YT = get(gca,'YTick');
        YTL = round(f(round(YT)));
        set(gca,'YTickLabel',YTL);
        ylabel('Frequency (Hz)');
            
        % if flags.do_SA
        %     title('SA (S4 in varnet2022a, S-LV)')
        % end
        % if flags.do_SB
        %     title('SB (S1 in varnet2022a, S-AO)')
        % end
    end
    if flags.do_SA
        subj_label = 'SA';
    end
    if flags.do_SB
        subj_label = 'SB';
    end
    title_here = [exp_label ', ' subj_label];
    title(title_here,'interpreter','none')
    
    % xlim([0 0.75])
end
h = [];
hname = [];

h(end+1) = gcf;
hname{end+1} = [exp_label '-' subj_label];

for i = 1:length(h)
    opts = [];
    opts.format = 'epsc';
    Saveas(h(i),[dir_out_figs hname{i}],opts);
    
    opts.format = 'png';
    Saveas(h(i),[dir_out_figs hname{i}],opts);
end
%%% End of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
