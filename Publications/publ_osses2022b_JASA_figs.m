function data = publ_osses2022b_JASA_figs(varargin)
% function data = publ_osses2022b_JASA_figs(varargin)
%
% Generates the figures
%
% % To display Fig. 1 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig1');
%
% % To display Fig. 2 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig2');
%
% % To display Fig. 4 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig4');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc

if nargin == 0
    help publ_osses2022b_JASA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig1','fig2','fig4'};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';
% 
% experiment = 'speechACI_varnet2013';
% dir_exp = [dir_fastACI_results experiment filesep];

if flags.do_fig1 || flags.do_fig2
    dir_data = fastACI_paths('dir_data');
    dir_subj = [dir_data 'speechACI_Logatome-abda-S43M' filesep 'S01' filesep];
    
    dir_noise = [dir_subj 'NoiseStim-white' filesep];
    
    basef = 8000;
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5, ...
        'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};
    
    CLim = [-8 -4.5];
    
% if mod(i_experiment,2) == 1
%         
%     else
%         flags_extra = {'colorbar','yes'};
%     end
%     set(gca,'YTickLabel',[]);
%     switch experiments{i_experiment}{1}
%         case 'modulationACI'
%             ACI_to_plot = ACI; %-squeeze(results.ACI);
%             warning('Temporal')
%         otherwise
%             ACI_to_plot = ACI; % squeeze(results.ACI);
%     end
%     outs = affichage_tf(ACI_to_plot,'CI', 'cfg', cfg_ACI, 'NfrequencyTicks', 8, flags_extra{:}); h    
end

if flags.do_fig1 % Speech-in-noise representation
    
    flags_extra = {'NfrequencyTicks',8,'colorbar','no'};
    
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    [aba, fs] = audioread([dir_subj 'speech-samples' filesep 'S43M_ab_ba.wav']);
    
    [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
    
    N_noises = 1000; % arbitrary number
    % SNR = -14; % dB, arbitrary value
    % gain_SNR = 10^(SNR/20);
    % aba_at_SNR = gain_SNR * aba;
        
    fname_noises = Get_filenames(dir_noise);
    fname_noises = fname_noises(1:N_noises); % names of the first N_noises
    
    G_Naba = nan([size(G_aba) N_noises]); % memory allocation
    G_N    = nan(size(G_Naba));
    
    for i_noise = 1:N_noises
        [noise, fs] = audioread([dir_noise fname_noises{i_noise}]);

        noisy_aba = noise+aba;
        [G_Naba(:,:,i_noise), fc, t, outs] = Gammatone_proc(noisy_aba, fs, flags_gamma{:});
        [G_N(:,:,i_noise), fc, t, outs]    = Gammatone_proc(noise    , fs, flags_gamma{:});
    end
    
    %%%
    G_Naba_avg = mean(G_Naba,3);
    G_N_avg    = mean(G_N,3);
    G_N_avg_std = G_N_avg + 2.2*std(G_N,[],3);

    % G_Naba_avg_thres = G_Naba_avg;
    %G_Naba_avg_thres(G_Naba_mean_thres<G_N_meanstd) = nan;

    idx_realisation = round(.733*N_noises);
    if N_noises ~= 1000
        warning('Less noise representations being used to derive the masking effects...');
    end
    G_Naba_thres = G_Naba(:,:,idx_realisation);
    %G_Naba_thres(G_Naba_thres<G_N_meanstd) = nan;

    figure('Position',[100 100 1000 250]); 
    il_tiledlayout(1,5,'TileSpacing','Compact'); % 'tight'); % none'); % 'Compact');
    
    x1 = -0.05; y1 = 1.05;
    x2 =  0.55; y2 = 0.92;
    x3 =  0.95; y3 = 0.85;
    
    il_nexttile(1); 
    affichage_tf(log(G_aba)', 'pow', t, fc, flags_extra{:}); 
    caxis(CLim); 
    title('Clean /aba/');
    text(x1,y1,'A','Units','Normalized','FontSize',12,'FontWeight','Bold')
    colorbar off;
    
    %nexttile; affichage_tf(log(G_ada)', 'pow', t, fc); caxis([-8 -4.5]); title('/ada/ target');colorbar off;ylabel('');set(gca,'YTickLabels',[]);xlabel('');set(gca,'XTickLabels',[]);%
    il_nexttile(2); 
    affichage_tf(log(G_Naba(:,:,1))', 'pow', t, fc, flags_extra{:}); 
    caxis(CLim); 
    title('Noisy /aba/','FontSize',8);
    text(x1,y1,'B','Units','Normalized','FontSize',12,'FontWeight','Bold');
    text(x2,y2,'N=1','Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
    % colorbar off;
    ylabel('');
    set(gca,'YTickLabels',[]);

    il_nexttile(3); 
    affichage_tf(log(G_Naba_avg)', 'pow', t, fc, flags_extra{:}); 
    caxis(CLim); 
    title('/aba/ and EM','FontSize',8);
    text(x1,y1,'C','Units','Normalized','FontSize',12,'FontWeight','Bold');
    text(x2,y2,sprintf('N=%.0f',N_noises),'Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
    text(x3,y3,'(floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');
    colorbar off;
    ylabel('');
    set(gca,'YTickLabels',[]);%
    
    il_nexttile(4); 
    affichage_tf(log(G_Naba_avg)', 'pow', t, fc, flags_extra{:}); 
    caxis(CLim); 
    title('/aba/ and EM+MM','FontSize',8);
    text(x1,y1,'D','Units','Normalized','FontSize',12,'FontWeight','Bold');
    text(x2,y2,sprintf('N=%.0f',N_noises),'Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
    text(x3,y3,'(floor+mod.floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');
    
    colorbar off;
    ylabel('');
    set(gca,'YTickLabels',[]);%
    
    il_nexttile(5); 
    affichage_tf(log(G_Naba_thres)', 'pow', t, fc, flags_extra{:}); 
    caxis(CLim); 
    title('    /aba/ and EM+MM+IM','FontSize',8);
    text(x1,y1,'E','Units','Normalized','FontSize',12,'FontWeight','Bold');
    text(x2,y2,'N=1','Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
    text(x3,y3,'(floor+mod.floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');
    
    ylabel('');
    set(gca,'YTickLabels',[]);%

    listImage = findobj('Type','Image');
    set(listImage(2),'AlphaData',1-0.5*(G_Naba_avg<G_N_avg_std)');   % before last panel
    set(listImage(1),'AlphaData',1-0.5*(G_Naba_thres<G_N_avg_std)'); % last panel 
    
    h(end+1) = gcf;
    hname{end+1} = 'fig1-schematic-masking';
end

if flags.do_fig2
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    [aba, fs] = audioread([dir_subj 'speech-samples' filesep 'S43M_ab_ba.wav']);
    [ada, fs] = audioread([dir_subj 'speech-samples' filesep 'S43M_ad_da.wav']);
    
    %%% Plot targets + 1 example of noise
    [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
    [G_ada, fc, t, outs] = Gammatone_proc(ada, fs, flags_gamma{:});

    warning('Conceptually wrong: The spectrograms are subjected to a natural logarithm...')
    disp('Pausing for 5 s')
    pause(5)
    
    figure('Position',[100 100 500 250]); 
    il_tiledlayout(1,2,'TileSpacing','Compact');
    il_nexttile(1); 
    affichage_tf(log(G_aba)', 'pow', t, fc); 
    caxis(CLim); 
    title('/aba/ target');
    colorbar off;
    
    il_nexttile(2); 
    affichage_tf(log(G_ada)', 'pow', t, fc); 
    caxis(CLim); 
    title('/ada/ target');
    ylabel('');
    set(gca,'YTickLabels',[]); %xlabel('')

    h(end+1) = gcf;
    hname{end+1} = 'fig2-spec-targets';
    
    disp('')    
end

if flags.do_fig4
    %%% Creating the frequency matrix:
    N_freq = 64;
    fhigh = 8000;
    fhigh_erb = freqtoaud(fhigh);
    delta_f = 0.5; % ERB
    flow_erb = fhigh_erb-delta_f*(N_freq-1);

    cfg_ACI = [];
    f = audtofreq(flow_erb:delta_f:fhigh_erb);
    cfg_ACI.f = f(:);

    %%% Creating the time matrix:
    binwidth = 10e-3; % 10 ms
    dur = .86;
    t = (0+binwidth:binwidth:dur);
    cfg_ACI.t = t;

    %%% Creating a B matrix, with arbitrary bumps located at indexes 'idxs':
    w = zeros(length(f),length(t));
    w = w(:);
    idxs = [10,150,560,1000,1300,2100,2193,2300,2599]; % Arbitrary indexes
    w(idxs) = 1;

    cfg_ACI.lasso_Nlevelmin = 2;
    cfg_ACI.lasso_Nlevelmax = 5;
    cfg_ACI.lasso_Nlevel    = 5;
    % The following is 'a prori' knowledge: the size of each dimension (freq, time)
    %    for levels 1 (omited), 2, 3, 4, and 5:
    cfg_ACI.lasso_Pyra_size = [0 0; 32 64; 16 32; 8 16; 4 8];

    %%% Convertion 'back' into the time-frequency domain:
    idxlambda = 1; % idle value
    keyvals.pyramid_script = 'imresize'; % old default value
    [~, cfg_ACI, sumReWeight] = Convert_lasso_B2ACI(w, cfg_ACI, idxlambda, keyvals);

    [~,~,sizeZ] = size(sumReWeight);
    switch sizeZ
        case 1 % then sumReWeight is bidimensional (because one lambda was simulated)
            % dim_f = 1; dim_t = 2;
            dim_lambda = [];
        otherwise % then sumReWeight is tridimensional (because more than one lambda was simulated):
            % dim_f = 2; dim_t = 3;
            dim_lambda = 1;
    end

    %% 5. Plot
    figure;
    if dim_lambda == 1
        affichage_tf( squeeze(sumReWeight(2,:,:)), 'CI', cfg_ACI.t, cfg_ACI.f);
    else
        affichage_tf( sumReWeight, 'CI', cfg_ACI.t, cfg_ACI.f);
    end
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');

    hold on;

    Pos = get(gcf,'Position');
    Pos(3:4) = [560   300];
    set(gcf,'Position',Pos);
    % set(gca,'YTick',[]);
    % set(gca,'XTick',[]);

    colorbar off

    there = [0.26,0.4];
    fhere = 5000*[1 1];
    plot(there,fhere,'b-','LineWidth',2);

    text(0.24,0.58,'ex. wide kernel','Units','Normalized','Color','b','FontWeight','bold');

    there = [0.61,0.65];
    fhere = 1450*[1 1];
    plot(there,fhere,'b-','LineWidth',2);

    text(0.57,0.14,'ex. narrow kernel','Units','Normalized','Color','b','FontWeight','bold');
    
    h(end+1) = gcf;
    hname{end+1} = 'fig4-Gaussbasis';    
end

data.h = h;
data.hname = hname;

% if flags.do_fig1a || flags.do_fig1b || flags.do_fig2 || flags.do_fig3a || flags.do_fig3b
%     
%     N_plots = length(folders)*length(noise_types);
%     thres = nan(13,N_plots);
% 
%     bPlot_ACI_norm = 1;
%     
%     for k = 1:length(noise_types)
%         noise_type = noise_types{k};
% 
%         filt = ['savegame*' noise_type '*.mat'];
% 
%         for i = 1:length(folders)
% 
%             data_folder_full = [dir_exp model filesep];
%             dir_where = [data_folder_full folders{i} filesep];
% 
%             files = Get_filenames(dir_where,filt);
% 
%             fname_results = [dir_where files{1}];
%             [cfg_game, data_passation] = Convert_ACI_data_type(fname_results);
%             
%             N_sessions = length(data_passation.resume_trial);
%             for j = 1:N_sessions
% 
%                 idxi = data_passation.resume_trial(j);
%                 if idxi == 0
%                     idxi = 1;
%                 end
%                 if j < N_sessions
%                     idxf = data_passation.resume_trial(j+1)-1;
%                 else
%                     idxf = cfg_game.N;
%                 end
% 
%                 thres(j,count) = prctile(data_passation.expvar(idxi:idxf),50);
%                 correct_score(j,count) = 100*sum(data_passation.is_correct(idxi:idxf))/(idxf-idxi+1);
%                 idx = find(data_passation.n_responses(idxi:idxf)==1);
%                 response_is_one(j,count) = 100*length(idx)/(idxf-idxi+1);
%                 idx = find(data_passation.n_responses(idxi:idxf)==2);
%                 response_is_two(j,count) = 100*length(idx)/(idxf-idxi+1);
%                 fprintf('\thres=%.2f dB, tidxi=%.0f, idxf=%.0f\n',thres(j,count),idxi,idxf);
% 
%             end
% 
%             Me(count) = prctile(thres(:,count),50);
%             errL(count) = Me(count) - prctile(thres(:,count),25);
%             errU(count) = prctile(thres(:,count),75) - Me(count);
% 
%             is_correct_all(count,:) = data_passation.is_correct;
% 
%             %%% Plotting the thresholds per session and global:
%             x_var = count;
%             if count == 1
%                 h(1) = figure; 
%             end
%             set(0, 'CurrentFigure', h(1));
% 
%             errorbar(x_var, Me(count),errL(count),errU(count)); hold on;
% 
%             plot(x_var*ones(size(thres(:,count))),thres(:,count),'bo');
% 
%             hname{1} = 'thres';
%             xlim([0.5 N_plots+.5]);
%             XT = 1:N_plots;
%             set(gca,'XTick',XT);
% 
%             count = count+1;
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             glmfct = 'lasso';
%             DimCI = 'gammatone';
%             add_signal = 0; % '1' to add signal in the ACI assessment, '0' to use noise alone trials
% 
%             if isempty(keyvals.dir_out)
%                 dir_out_ACI = [dir_where 'Results_ACI' filesep];
%             else
%                 if cfg_game.is_simulation
%                     % Creating a new subfolder for each model run. This is 
%                     %   important because all model runs will produce the same
%                     %   ACI name.
%                     dir_out_ACI = [keyvals.dir_out model '-' folders{i} filesep];
%                     if ~exist(dir_out_ACI,'dir')
%                         mkdir(dir_out_ACI);
%                     end
%                 else
%                     dir_out_ACI = keyvals.dir_out;
%                 end
%             end
%             if ~exist(dir_out_ACI,'dir')
%                 mkdir(dir_out_ACI);
%             end
%             dir_noise  = cfg_game.dir_noise;
%             dirname4waveforms = Get_subjectname_from_dirname(dir_noise);
%             
%             if ~exist(dir_noise,'dir')
%                 switch dirname4waveforms
%                     case {'SLeo','SLeoVarnet2013'} % this folder does not exist anymore:
%                         dir_new = 'osses2021c_S01';
%                     case 'SAO-5000-trials'
%                         dir_new = 'osses2021c_S02';
%                 end
%                 idx = strfind(dir_noise,dirname4waveforms);
%                 idx = idx + length(dirname4waveforms);
%                 dir_noise = [fastACI_paths('dir_data') experiment filesep dir_new filesep dir_noise(idx+1:end-1) filesep];
%             end
%             
%             if isfield(cfg_game,'dir_speech')
%                 cfg_game.dir_target = cfg_game.dir_speech;
%             end
%             dir_target = cfg_game.dir_target; 
% 
%             
%             if ~exist(dir_noise,'dir')
%                 idx = strfind(dir_noise,filesep);
%                 dir_noise = [data_folder_full dir_noise(idx(end-1)+1:end-1) filesep];
%                 cfg_game.dir_noise = dir_noise;
%             end
%             
%             if ~exist(dir_target,'dir')
%                 idx = strfind(dir_target,filesep);
%                 dir_target = [data_folder_full dir_target(idx(end-1)+1:end-1) filesep];
%                 cfg_game.dir_target = dir_target;
%             end
% 
%             %%% bCalculate:
%             f_limits = [1 10000];
%             t_limits = [0 1]; 
%             %%% end reading folders
% 
%             Data_matrix = [];
%             fname_ACI = [];
% 
%             trial_select = 5000; 
%             ACI_all    = [];
%             ACI_incorr = [];
% 
%             for ii = 1:length(trial_select)
%                 switch trial_select(ii)
%                     case 5000
%                         idx_trialselect = [];
%                     otherwise
%                         idx_trialselect = 1:trial_select(ii);
%                 end
%                 fg_ACI = {'dir_noise', dir_noise, 'dir_target', dir_target, ...
%                   'dir_out', dir_out_ACI, 'no_plot', ...
%                   'idx_trialselect', idx_trialselect, ...
%                   'f_limits',f_limits, ...
%                   't_limits',t_limits, ... % 'spect_NFFT',512,'spect_Nwindow',512,'spect_overlap',.75 ... %'spect_NFFT',1024,'spect_Nwindow',1024,'spect_overlap',.5...
%                   'skip_if_on_disk',1, ...
%                   'add_signal',add_signal, ...
%                   'N_perm',20, ...
%                   'pyramid_script','imresize' ...
%                 };
% 
%                 if isempty(Data_matrix)
%                     flags_to_use = fg_ACI;
%                     [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results,DimCI,glmfct,flags_to_use{:});
%                     flags_Data_matrix = {'Data_matrix',Data_matrix};
%                 else
%                     flags_to_use = [fg_ACI flags_Data_matrix];
%                     [ACI,cfg_ACI,results] = fastACI_getACI(fname_results,DimCI,glmfct,flags_to_use{:});
%                 end
%                 correct_text = '';
%                 if bPlot_ACI_norm == 1
%                     ACI = ACI / max(max(abs(ACI)));
%                 end
%                 ACI_all(:,:,ii) = ACI;
%                 htmp = figure;
%                 figure(htmp); set(0, 'CurrentFigure', htmp)
%                 [h(end+1),hname{end+1}] = il_plot_the_ACI(ACI, cfg_ACI, cfg_game);
%                 hname{end} = [hname{end} '-' folders{i}];
%                 
%                 xlim([0.05 0.55]); % warning('Temporal here')
%             end
% 
%         end
%     end
% 
%     disp('The response ''1'' was preferred in the following percentage of the times:')
%     bias = [prctile(response_is_one,75); median(response_is_one); prctile(response_is_one,25)];
% end
% 
% if flags.do_fig4
%     Pos34 = [700 250];
%     
%     for k = 1:length(model)
%         noise_type = 'SSN';
% 
%         filt = ['savegame*' noise_type '*.mat'];
% 
%         switch model{k}
%             case 'osses2021c_S01'
%                 folders = folders1;
%             case 'osses2021c_S02'
%                 folders = folders2;
%             case 'osses2021'
%                 folders = folders3;
%             otherwise
%                 error('Maybe a mistake?')
%         end
% 
%         for i = 1:length(folders)
% 
%             data_folder_full = [dir_exp model{k} filesep];
%             dir_where = [data_folder_full folders{i} filesep];
% 
%             files = Get_filenames(dir_where,filt);
% 
%             fname_results = [dir_where files{1}];
%             cfg_game = [];
%             data_passation = [];
%             load(fname_results,'cfg_game','data_passation');
% 
%             N_sessions = length(data_passation.resume_trial);
%             for j = 1:N_sessions
% 
%                 idxi = data_passation.resume_trial(j);
%                 if idxi == 0
%                     idxi = 1;
%                 end
%                 if j < N_sessions
%                     idxf = data_passation.resume_trial(j+1)-1;
%                 else
%                     idxf = cfg_game.N;
%                 end
% 
%                 thres(j,count) = prctile(data_passation.expvar(idxi:idxf),50);
%                 correct_score(j,count) = 100*sum(data_passation.is_correct(idxi:idxf))/(idxf-idxi+1);
%                 idx = find(data_passation.n_responses(idxi:idxf)==1);
%                 response_is_one(j,count) = 100*length(idx)/(idxf-idxi+1);
%                 idx = find(data_passation.n_responses(idxi:idxf)==2);
%                 response_is_two(j,count) = 100*length(idx)/(idxf-idxi+1);
%                 fprintf('\thres=%.2f dB, tidxi=%.0f, idxf=%.0f\n',thres(j,count),idxi,idxf);
% 
%             end
% 
%             Me(count) = prctile(thres(:,count),50);
%             errL(count) = Me(count) - prctile(thres(:,count),25);
%             errU(count) = prctile(thres(:,count),75) - Me(count);
% 
%             is_correct_all(count,:) = data_passation.is_correct;
% 
%             %%% Plotting the thresholds per session and global:
%             x_off = .1;
% 
%             x_var = count;
% 
%             if count == 1
%                 h(1) = figure; grid on
%             end
%             set(0, 'CurrentFigure', h(1));
% 
%             plot((x_var-x_off)*ones(size(thres(:,count))),thres(:,count),'o','Color',[.7 .7 .7],'MarkerFaceColor',[.7 .7 .7]); hold on
%             errorbar(x_var, Me(count),errL(count),errU(count),'ks-','LineWidth',2,'MarkerFaceColor','k'); 
% 
%             xlim([0.5 N_plots+.5]);
%             XT = 1:N_plots;
%             set(gca,'XTick',XT);
% 
%             if count == 1
%                 ylim([-22 2]); grid on
% 
%                 Pos = get(gcf,'Position');
%                 Pos(3:4) = Pos34;
%                 set(gcf,'Position',Pos);
%                 YT = -20:2:0;
%                 set(gca,'YTick',YT);
%                 ylabel({'Discrimination threshold';'SNR (dB)'})
%                 %title('(a)')
%                 xlabel('Listener ID')
%                 hname{1} = 'fig4-a-thres';
%             end
% 
%             Meb(count) = prctile(response_is_one(:,count),50);
%             errLb(count) = Meb(count) - prctile(response_is_one(:,count),25);
%             errUb(count) = prctile(response_is_one(:,count),75) - Meb(count);
% 
%             if count == 1
%                 h(2) = figure; 
%             end
%             set(0, 'CurrentFigure', h(2));
% 
%             if count == 1
%                 plot([0 N_plots+1],[50 50],'r--','LineWidth',2); grid on; hold on
% 
%                 Pos = get(gcf,'Position');
%                 Pos(3:4) = Pos34;
%                 set(gcf,'Position',Pos);
%                 YT = 10:10:90;
%                 set(gca,'YTick',YT);
%                 for kk = 1:length(YT)
%                     if kk == 1
%                         YTL{kk} = [' ' num2str(YT(kk))];
%                     else
%                         YTL{kk} = num2str(YT(kk));
%                     end
%                 end
%                 ylim([0 100]);
%                 set(gca,'YTickLabel',YTL);
% 
%                 ylabel({'Response bias';'towards /aba/ (%)'})
%                 % title('(b)')
%                 xlabel('Listener ID')
%                 hname{2} = 'fig4-b-bias';
%             end
% 
%             plot((x_var-x_off)*ones(size(response_is_one(:,count))),response_is_one(:,count),'o','Color',[.7 .7 .7],'MarkerFaceColor',[.7 .7 .7]); hold on
%             errorbar(x_var, Meb(count),errLb(count),errUb(count),'ks-','LineWidth',2,'MarkerFaceColor','k'); hold on;
% 
% 
%             xlim([0.5 N_plots+.5]);
%             XT = 1:N_plots;
%             set(gca,'XTick',XT);
% 
%             count = count+1;
% 
%         end
%     end
%     lab2use = {'S01','S02','-1.55 MU','0 MU','0.39 MU','0.78 MU','Dec2'};
%     for i = 1:length(h)
%         set(0, 'CurrentFigure', h(i));
%         set(gca,'XTickLabels',lab2use);
%         % xtickangle(-90);
% 
%         YL = get(gca,'YLim');
%         plot(2.5*[1 1],YL,'k-');
%         plot(6.5*[1 1],YL,'k-');
% 
%         switch i
%             case 1
%                 lab_here = '(a)';
%             case 2
%                 lab_here = '(b)';
%         end
%         text(0.02,0.95,lab_here,'Units','Normalized','FontWeight','Bold')
%         text(0.03,0.88,'Exp.','Units','Normalized','FontWeight','Bold')
%         text(0.33,0.88,'Decision 1','Units','Normalized','FontWeight','Bold')
%         text(0.88,0.88,'Dec. 2','Units','Normalized','FontWeight','Bold')
%     end
%     disp('The response ''1'' was preferred in the following percentage of the times:')
%     bias = [prctile(response_is_one,75); median(response_is_one); prctile(response_is_one,25)];
% end
% 
% bSave = 1;
% if nargout ~= 0
%     data.h = h;
%     data.hname = hname;
%     bSave = 0;
% end
% 
% if bSave
%     dir_out = fastACI_paths('dir_output');
%     if ~exist(dir_out,'dir')
%         dir_out = [pwd filesep 'outputs' filesep];
%         mkdir(dir_out);
%     end
% 
%     for i = 1:length(h)
%         opts = [];
% 
%         opts.format = 'epsc';
%         Saveas(h(i),[dir_out hname{i}],opts);
% 
%         opts.format = 'png';
%         Saveas(h(i),[dir_out hname{i}],opts);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [fig_handle,fig_name] = il_plot_the_ACI(ACI, cfg_ACI, cfg_game)
% 
% dir_where = cfg_game.dir_target;
% 
% h_outs = affichage_tf(ACI,'CI','cfg',cfg_ACI); hold on;
% 
% %%%
% par_formants.timestep = 0.01; % positive timestep 0.01
% par_formants.nformants = 5; % positive nformants 5
% par_formants.maxformant = 5500; % positive maxformant 5500
% par_formants.windowlength = 0.025; % positive windowlength 0.025
% par_formants.dynamicrange = 30; % positive dynamic range 20
% 
% par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
% par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
% par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)
% 
% par_formants.I_min = 50; % 60, arbitrary value
% try
%     outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
%     if isempty(outs_from_Praat.filesf0)
%         error('Files not found...')
%     end
% catch me
%     fprintf('%s: Praat not found on disk, loading pre-stored Praat files\n',upper(mfilename));
%     
%     dir_where_stored = [fastACI_basepath 'Stimuli' filesep 'varnet2013' filesep]; 
%     outs_from_Praat = Get_all_metrics_from_stored_Praat(dir_where_stored,par_formants);
% end
% 
% idx = find(outs_from_Praat.t_f0{1} < 0.089 | (outs_from_Praat.t_f0{1} > 0.164 & outs_from_Praat.t_f0{1} < .280) | outs_from_Praat.t_f0{1} > 0.404);
% outs_from_Praat.f0{1}(idx,:) = nan;
% outs_from_Praat.f0{2}(idx,:) = nan;
% outs_from_Praat.F{1}(idx,:) = nan;
% outs_from_Praat.F{2}(idx,:) = nan;
% 
% %%%
% Style  = {'-','--'};
% Colour = {[0.5 0.5 0.5],'k'};
% LW     = 1; 
% outs_from_Praat = affichage_tf_add_Praat_metrics(dir_where, cfg_ACI,outs_from_Praat,Style,Colour,LW);
% 
% %%% Plotting
% idx = strfind(cfg_ACI.fnameACI,filesep);
% fname_ACI_short = cfg_ACI.fnameACI(idx(end)+1:end);
% 
% title(name2figname(fname_ACI_short));
% fig_handle = gcf;
% fig_name = name2figname(fname_ACI_short);
% 
% Pos = get(fig_handle,'Position');
% Pos(3) = 400;
% set(fig_handle,'Position',Pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_tiledlayout(N,M,TileSpacing,TileSpacing_option)

if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout','file'); % tiledlayout.p
if bExist
    tiledlayout(N,M,TileSpacing,TileSpacing_option);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_nexttile(N)

bExist = exist('nexttile','file'); % nexttile
if bExist
    if nargin == 0
        nexttile;
    else
        nexttile(N);
    end
else
    if N == 1
        close;
    end
    figure;
end