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

% data = [];
 
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
    
    dir_out = [fastACI_paths('dir_output_fastACI2021_JASA') 'Data_proc' filesep];
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    
    % fname_results = [paths_results '20210421-SAO-5000-trials_speech_varnet2013' filesep 'savegame_2021_04_22_14_46_SAO-5000-trials_speechACI_varnet2013.mat'];
    
    % From Slack (#proj-fastaci), from Leo on 20/05/2021 at 8:51:
    % ''ACIs in white noise (WN) and spech shape noise (SSN) for apa-ata 
    %   categorization, calculated with the basic revcorr or our new "lasso 
    %   on Gaussian pyramid"... nice, isn't it ? (5000 trials each, lines 
    %   show formants and f0 trajectory)''
    
    type = 1;
    switch type
        case 1
            fprintf('SLeo, Logatome, white noise data\n')
            data_folder   = '20210428-SLeoWN-speechACI_Logatome';
            data_mat_file = 'savegame_2021_05_05_13_51_SLeoWN_speechACI_Logatome.mat';
            dir_noise  = [paths.dir_data 'speechACI_Logatome-apta-S46M' filesep 'SLeo' filesep 'NoiseStim-white' filesep];
            dir_target = [paths.dir_data 'speechACI_Logatome-apta-S46M' filesep 'SLeo' filesep 'speech-samples'  filesep];
        case 2
            fprintf('SLeo, Logatome, SSN data\n')
            data_folder   = '20210421-SLeo_speechACI_Logatome';
            data_mat_file = 'savegame_2021_04_28_13_10_SLeo_speechACI_Logatome.mat';
            dir_noise = [paths.dir_data 'speechACI_Logatome-apta-S46M'  filesep 'SLeo' filesep 'NoiseStim-SSN'  filesep];
            dir_target = [paths.dir_data 'speechACI_Logatome-apta-S46M' filesep 'SLeo' filesep 'speech-samples' filesep];
    end
    
    fname_results = [paths_results data_folder filesep data_mat_file];
    
    glm_functions = {   'classic_revcorr', ...
                        'glmfitqp', ...
                        'lasso'};
    prefix_title = {'A. ', 'B. ', 'C. '};                
    for i = 1:length(glm_functions)
        glmfct = glm_functions{i};
    
        % t_limits = [0.0 1]; % No limit at all
        DimCI = 'spect'; % 'tf'
    
        %%% Extra flags used by AO:
        % flags = {DimCI,'f_limits',f_limits,'t_limits',t_limits,'dir_out',dir_out, ...
        %     'dir_noise',dir_noise};
    
        switch glmfct
            case 'glmfitqp'
                f_limits = [0 4050];
                
                flags = {'dir_noise', dir_noise, 'dir_target', dir_target, 'dir_out', dir_out, 'no_plot', ...
                  'idx_trialselect',[], ... % 'idx_trialselect',1:2400, ...
                  'f_limits',f_limits, ...
                  'spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',0, ... %'spect_NFFT',1024,'spect_Nwindow',1024,'spect_overlap',.5...
                }; 
            case {'lasso','classic_revcorr'}
                flags = {'dir_noise', dir_noise, 'dir_target', dir_target, 'dir_out', dir_out, 'no_plot', ...
                  'idx_trialselect',[], ... % 'idx_trialselect',1:2400, ...
                  'spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',.75 ... %'spect_NFFT',1024,'spect_Nwindow',1024,'spect_overlap',.5...
                };
        end

        [ACI,cfg_ACI,results] = Script4_Calcul_ACI(fname_results,DimCI,glmfct,flags{:});

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