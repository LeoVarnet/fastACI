function data = exp_fastACI2021_JASA(varargin)
%EXP_FASTACI2021_JASA 
%
%   Usage: data = exp_fastACI2021_JASA(flags)
%
%   To display Fig. 3 of Osses et al., (2022) use :::
%
%     out = exp_osses2022('fig3');
%
%   To display Fig. 4 of Osses et al., (2022) use :::
%
%     out = exp_osses2022('fig4');

%   #Author: Alejandro Osses and Leo Varnet (2020-2021)

data = [];
 
h = []; % to store figure handles
hname = []; % to store figure names
close all

% definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig1'}; 

definput.flags.plot={'plot','no_plot'};
definput.keyvals.models=[];
 
[flags,keyvals]  = ltfatarghelper({},definput,varargin);
 
if flags.do_missingflag
    flags.do_fig1 = 1;
    disp('Plotting Fig. 1 (default)')
    
    flags.do_missingflag = 0;
end
disp('')
% %%%
% FontSize = 13; % FontSize for all figure axes
% fs = 100e3; % 100 kHz of sampling rate, this is arbitrary
% dBFS = 94; % i.e., amplitude 1 = 1 Pa = 94 dB SPL re 2x10^{-5} Pa
% 
% %%% for Figs 3, 4, 5, 6 (confirm for the first figures)
% models = {'dau1997','zilany2014','verhulst2015','verhulst2018','bruce2018','king2019','relanoiborra2019','osses2021'};
% Colours = {'b',il_rgb('Green'),'r',il_rgb('LightSkyBlue'),il_rgb('Maroon'),'m','k',il_rgb('mediumorchid')}; % ,'k'};
% Markers = {'o','s','<','>','v','^','p','d'};
% MarkersSize = [10 10 10 10 10 10 10 10];
% LineStyle = {'-','-','-.','-','--','-','-','-'};
% LineWidth = [2 2 2 2 2 2 2 2];
% 
% if ~isempty(keyvals.models),
% 	warning('Not all models selected'); 
% 	idxs=keyvals.models;
% 	Markers = Markers(idxs); 
% 	MarkersSize = MarkersSize(idxs); 
% 	Colours = Colours(idxs); 
% 	LineStyle = LineStyle(idxs); 
% 	LineWidth = LineWidth(idxs); 
% 	models = models(idxs);
% end
% 
% N_models = length(models);
% 
% figure_handle = []; % multiple figures will be generated
% figure_name   = [];

paths = fastACI_paths;
paths_results = [paths.dir_fastACI_sim 'Results' filesep];

%% ------ FIG 1 Osses and Varnet 2021 (in preparation) --------------------
if flags.do_fig1
    
    % fname_results = [paths_results '20210421-SAO-5000-trials_speech_varnet2013' filesep 'savegame_2021_04_22_14_46_SAO-5000-trials_speechACI_varnet2013.mat'];
    
    % From Slack (#proj-fastaci), from Leo on 20/05/2021 at 8:51:
    % ''ACIs in white noise (WN) and spech shape noise (SSN) for apa-ata 
    %   categorization, calculated with the basic revcorr or our new "lasso 
    %   on Gaussian pyramid"... nice, isn't it ? (5000 trials each, lines 
    %   show formants and f0 trajectory)''
    
    YL = [-40 4100];
    
    %%% SLeo: types 1, 2, 3, 5
    
    type = input('Enter 1 to process Logatome data; 3 to process varnet2013 data; 6 PH white: ');
    flags_extra = {};
    switch type
        case 1
            fprintf('SLeo, Logatome, white noise data\n')
            data_folder   = '20210428-SLeoWN-speechACI_Logatome';
            data_mat_file = 'savegame_2021_05_05_13_51_SLeoWN_speechACI_Logatome.mat';
            data_folder_full = [paths.dir_data 'speechACI_Logatome-apta-S46M' filesep 'SLeo' filesep];
            dir_noise  = [data_folder_full 'NoiseStim-white' filesep];
            dir_target = [data_folder_full 'speech-samples'  filesep];
            
            t_limits = [0.0 1]; % No limit at all
            f_limits = [1 10000];
            fname_results = [paths_results data_folder filesep data_mat_file];
            
            flags_extra = {'trialtype_analysis', 'incorrect'};
        case 2
            fprintf('SLeo, Logatome, SSN data\n')
            data_folder   = '20210421-SLeo_speechACI_Logatome';
            data_mat_file = 'savegame_2021_04_28_13_10_SLeo_speechACI_Logatome.mat';
            data_folder_full = [paths.dir_data 'speechACI_Logatome-apta-S46M'  filesep 'SLeo' filesep];
            dir_noise  = [data_folder_full 'NoiseStim-SSN'  filesep];
            dir_target = [data_folder_full 'speech-samples' filesep];
            
            t_limits = [0.0 1]; % No limit at all
            f_limits = [1 10000];
            fname_results = [paths_results data_folder filesep data_mat_file];
            
            flags_extra = {'trialtype_analysis', 'incorrect'};
            
        case 2.1
            fprintf('SLeo, Logatome, pink noise data\n')
            data_folder   = '20210531-SLeo_abdaS41Fpink';
            data_mat_file = 'savegame_2021_05_31_15_39_SLeo_speechACI_Logatome-abda-S41F_pink.mat';
            data_folder_full = [paths.dir_data 'speechACI_Logatome-abda-S41F'  filesep 'SLeo' filesep];
            dir_noise = [data_folder_full 'NoiseStim-pink'  filesep];
            dir_target = [data_folder_full 'speech-samples' filesep];
            
            t_limits = [0.0 1]; % No limit at all
            f_limits = [1 10000];
            fname_results = [paths_results data_folder filesep data_mat_file];
            
            type = type*10; % 21
            
            YL = [-40 8100];
            
            flags_extra = {'trialtype_analysis', 'incorrect'};
            
        case 3
            fprintf('SLeo, speechACI_varnet2013, white noise data\n')
            % '/home/alejandro/Documents/Databases/data/fastACI/data_varnet2013/'
            data_folder   = 'Sujet_Leo_S1';
            data_mat_file = 'savegame_final.mat';
            data_folder_full = [paths.dir_data 'data_varnet2013' filesep data_folder filesep];
            dir_noise  = [data_folder_full 'ListeBruit'  filesep];
            dir_target = [data_folder_full 'ListeSignal' filesep];
            
            f_limits = [1 10000];
            t_limits = [0.0 342.5e-3]; 
            fname_results = [data_folder_full filesep data_mat_file];
            
        case 4
            fprintf('SAOtest, speechACI_varnet2013, white noise data\n')
            data_folder   = 'SAO-5000-trials'; 
            data_mat_file = ['Results' filesep 'savegame_2021_04_22_14_46_SAO-5000-trials_speechACI_varnet2013.mat'];
            data_folder_full = [paths.dir_data 'speechACI_varnet2013' filesep data_folder filesep]; %  data/fastACI/speechACI_varnet2013/SAO-5000-trials/NoiseStim
            dir_noise  = [data_folder_full 'NoiseStim'      filesep];
            dir_target = [data_folder_full 'speech-samples' filesep];
            
            f_limits = [1 10000];
            t_limits = [0.0 342.5e-3]; 
            fname_results = [data_folder_full filesep data_mat_file];     
            
            % flags_extra = {'trialtype_analysis', 'incorrect'};
        
        case 4.1
            fprintf('SAOtest, speechACI_varnet2013, speech-shaped noise data\n')
            data_folder   = 'SAO-5000-trials'; 
            data_mat_file = ['Results' filesep 'savegame_2021_04_21_15_21_SAO-5000-trials_speechACI_varnet2013_SSN.mat'];
            data_folder_full = [paths.dir_data 'speechACI_varnet2013' filesep data_folder filesep]; %  data/fastACI/speechACI_varnet2013/SAO-5000-trials/NoiseStim
            dir_noise  = [data_folder_full 'NoiseStim-SSN'  filesep];
            dir_target = [data_folder_full 'speech-samples' filesep];
            
            f_limits = [1 10000];
            t_limits = [0.0 342.5e-3]; 
            fname_results = [data_folder_full filesep data_mat_file];     
            
            flags_extra = {'trialtype_analysis', 'incorrect'};
            
            type = type*10; % 21
        case 5
            fprintf('SLeo, speechACI_varnet2013, white noise data (new pilot data)\n')
            data_folder   = 'SLeo';
            % /home/alejandro/Documents/Databases/data/fastACI/speechACI_varnet2013/SLeo/Results/
            data_mat_file = ['Results' filesep 'savegame_2021_05_10_12_45_SLeo_speechACI_varnet2013_white.mat'];
            data_folder_full = [paths.dir_data 'speechACI_varnet2013' filesep data_folder filesep];
            dir_noise  = [data_folder_full 'NoiseStim'      filesep];
            dir_target = [data_folder_full 'speech-samples' filesep];
            
            f_limits = [1 10000];
            t_limits = [0.0 342.5e-3]; 
            fname_results = [data_folder_full filesep data_mat_file];
            
            % flags_extra = {'trialtype_analysis', 'incorrect'};
            
        case 6 
            fprintf('PH, speechACI_varnet2013, white noise data (new pilot data)\n')
            data_folder   = 'S_PH';
            data_mat_file = ['Results' filesep 'savegame_2021_06_16_15_15_S_PH_speechACI_Logatome-abda-S41F_white.mat'];
            data_folder_full = [paths.dir_data 'speechACI_Logatome-abda-S41F' filesep data_folder filesep];
            dir_noise  = [data_folder_full 'NoiseStim-white' filesep];
            dir_target = [data_folder_full 'speech-samples'  filesep];
            
            f_limits = [1 10000];
            t_limits = [0 1]; 
            fname_results = [data_folder_full filesep data_mat_file];
       
            %%% From l20210617_testScript4debug.m
            % % testScript4debug
            
            % flags = {'dir_out', fullpath,'no_plot','spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',.75... %
            %     'f_limits', [0 4050],...
            %     };
             
            % flags_extra = {'trialtype_analysis', 't1', 'expvar_limits',[-20, -5]};
            flags_extra = {'trialtype_analysis', 'incorrect'};
    end
   
    %%%
    bUse_global_dir_out = 0;
    if bUse_global_dir_out
        dir_out = [fastACI_paths('dir_output_fastACI2021_JASA') 'Data_proc' filesep];
        if ~exist(dir_out,'dir')
            mkdir(dir_out);
        end
    else
        dir_out = [data_folder_full filesep 'Results_ACI' filesep];
    end
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    %%%
    
    glm_functions = {   'classic_revcorr', ...
                        'glmfitqp', ...
                        'lasso'};
    prefix_title = {'A. ', 'B. ', 'C. '};                
    for i = 1:length(glm_functions)
        glmfct = glm_functions{i};
    
        DimCI = 'spect'; % 'tf'
    
        %%% Extra flags used by AO:
        % flags = {DimCI,'f_limits',f_limits,'t_limits',t_limits,'dir_out',dir_out, ...
        %     'dir_noise',dir_noise};
    
        switch glmfct
            case 'glmfitqp'
                f_limits = [0 4050];
                
                fg_ACI = {'dir_noise', dir_noise, 'dir_target', dir_target, 'dir_out', dir_out, 'no_plot', ...
                  'idx_trialselect',[], ... % 'idx_trialselect',1:2400, ...
                  'f_limits',f_limits,'t_limits',t_limits, ...
                  'spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',0, ... %'spect_NFFT',1024,'spect_Nwindow',1024,'spect_overlap',.5...
                }; 
            case {'lasso','classic_revcorr'}
                
                if type == 3 && strcmp(glmfct,'lasso')
                    f_limits = [1 12000]; warning('Temporal arrangement...')
                end
                if (type == 21 || type == 4 || type == 4.1 || type == 5 || type == 6) && strcmp(glmfct,'lasso')
                    f_limits = [1 12000]; warning('Temporal arrangement...')
                end
                
                fg_ACI = {'dir_noise', dir_noise, 'dir_target', dir_target, 'dir_out', dir_out, 'no_plot', ...
                  'idx_trialselect',[], ... % 'idx_trialselect',1:2400, ...
                  'f_limits',f_limits, ...
                  't_limits',t_limits, ...
                  'spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',.75 ... %'spect_NFFT',1024,'spect_Nwindow',1024,'spect_overlap',.5...
                };
        end

        flags_to_use = [fg_ACI flags_extra];
        [ACI,cfg_ACI,results] = fastACI_getACI(fname_results,DimCI,glmfct,flags_to_use{:});

        if isfield(results,'ACI_perm')
            bPermutation = 1;
            str_suf = ' + permutation test';
        else
            bPermutation = 0;
            str_suf = '';
        end
        
        if bPermutation == 0
            figure;
            affichage_tf(ACI,'CI','cfg',cfg_ACI); % ,'caxis',[-1 1]);
        end
        if bPermutation == 1
            %%% Plotting permutation test:
            idxs         = find(ACI(:)>results.ACI_perm_CI_high(:) |  ACI(:)< results.ACI_perm_CI_low(:));
            % idxs_exclude = find(ACI(:)<=results.ACI_perm_CI_high(:) & ACI(:)>=results.ACI_perm_CI_low(:));

            ACI_ex = zeros(size(ACI));
            ACI_ex(idxs) = ACI(idxs);
            figure;
            affichage_tf(ACI_ex, 'CI', 'cfg',cfg_ACI);  
        end
        
        ylim(YL)
        YT = 0:500:YL(2);
        set(gca,'YTick',YT);
        set(gca,'YTickLabel',YT);
        title([prefix_title{i} name2figname(glmfct) str_suf])
        
        h(end+1) = gcf;
        hname{end+1} = ['fig01-type-' num2str(type) '-' glmfct];
    end
    
    disp('')
    % figure; 
    % % h=pcolor(t_X, f_X, -sumReWeight(:,:,idxlambda)); set(h, 'EdgeColor', 'none'); xlabel('time (s)'); ylabel('freq (Hz)'); colorbar; title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); colorbar; caxis([-1 1]*max(abs(caxis)));
    % h=pcolor(cfg_ACI.t, cfg_ACI.f, ACI); 
    % set(h, 'EdgeColor', 'none'); 
    % xlabel('time (s)'); 
    % ylabel('freq (Hz)'); colorbar; 
    % title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); 
    % colorbar; % caxis([-1 1]*max(abs(caxis)));
        
end

dir_out_eps = paths.dir_output_fastACI2021_JASA_eps;
dir_out_all = [paths.dir_output_fastACI2021_JASA 'Figures-raw' filesep]; mkdir(dir_out_all);
for i = 1:length(h)
    opts = [];
    opts.format = 'epsc';
    Saveas(h(i),[dir_out_eps hname{i}],opts); % to be used in the paper
    Saveas(h(i),[dir_out_all hname{i}],opts);
    
    opts.format = 'fig';
    Saveas(h(i),[dir_out_all hname{i}],opts);
    
    opts.format = 'png';
    Saveas(h(i),[dir_out_all hname{i}],opts);
end