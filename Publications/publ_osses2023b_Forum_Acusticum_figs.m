function data = publ_osses2023b_Forum_Acusticum_figs(varargin)
% function publ_osses2023b_Forum_Acusticum_figs(varargin)
%
% Generates the figures
%
% % To display Fig. 2 of Osses and Varnet, (2023, Forum Acusticum) use :::
%     publ_osses2023b_Forum_Acusticum_figs('fig2a');
%     publ_osses2023b_Forum_Acusticum_figs('fig2b');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc

if nargin == 0
    help publ_osses2023b_Forum_Acusticum_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2a','fig2b', ...
    'fig3', ... % supra-threshold
    'fig5'}; % ACI
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';

experiment = 'toneinnoise_ahumada1975';
dir_exp = [dir_fastACI_results experiment filesep];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2a || flags.do_fig2b  || flags.do_fig3   
    dir_noise = [dir_exp 'king2019' filesep 'NoiseStims-white' filesep]; 
    files = Get_filenames(dir_noise,'*.wav');
    
    var = load([dir_exp 'king2019/Results-run-1/savegame_2023_04_19_01_08_king2019_toneinnoise_ahumada1975_white.mat']);
    cfg_game = var.cfg_game;
    tone = cfg_game.signal;
    dBFS = cfg_game.dBFS;
    fs = cfg_game.fs;
    
    if flags.do_fig2a || flags.do_fig3
        N = 1000;
    end
    if flags.do_fig2b 
        N = 100;
        warning('temporal')
    end
    Ni = 1600;
    for i = 1:N
        [noises_TN(:,i),fs] = audioread([dir_noise files{i}]);
        
        noises_N(:,i) = audioread([dir_noise files{i+Ni}]);
    end
    if flags.do_fig3
        thres = median(var.data_passation.expvar);
        noises_TN_thres = noises_TN + From_dB(thres)*repmat(tone,[1 N]);
    end
    
    noises_TN = noises_TN + repmat(tone,[1 N]);
end

if flags.do_fig2a    
    t = (1:length(tone))/fs;
    toffset = 12e-3;
    idx = find(t>=.2+toffset & t<=.3-toffset);
            
    %%% Upsampling:
    % noises_TN = resample(noises_TN,2*fs,fs);
    % noises_N  = resample(noises_N ,2*fs,fs);
    % fs = 2*fs;
    
    % [x, x_dB, f]  = freqfft2(tone(idx),K,fs,win_type,dBFS);
    win_type = 'hanning';
    % win_type = 'rectangular';
    for i = 1:N
        K = fs/2; % length(idx); % fs/2;
        % if i == 1
        %     
        % end
        [x, y_dB_TN(:,i),f] = freqfft2(noises_TN(idx,i),K,fs,win_type,dBFS);
        [x, y_dB_N(:,i)]    = freqfft2(noises_N(idx,i) ,K,fs,win_type,dBFS);
    end
    
    y_TN = median(y_dB_TN,2);
    y_N  = median(y_dB_N,2);
    
    f_ERB = freqtoaud(f);
    
    figure;
    
    plot(f_ERB,y_TN,'b-','LineWidth',2); hold on; grid on;
    plot(f_ERB,y_N,'r--','LineWidth',2);
    % plot(f_ERB,x_dB,'m-','LineWidth',2);
    ylabel('Average amplitude (dB)');
    xlabel('Frequency (Hz)');
    XTL = [125 250 500 1000 2000];
    set(gca,'XTick',freqtoaud(XTL));
    set(gca,'XTickLabel',XTL);
    xlim(freqtoaud([100 2400]))
    
    Pos = get(gcf,'Position');
    Pos(4) = 300;
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2a-FFT';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2b
    basef = 500;
	for i = 1:N
        env_TN(:,i) = il_get_hilbert(noises_TN(:,i),fs,basef);
        env_N(:,i)  = il_get_hilbert(noises_N(:,i),fs,basef);
    end
    
    factor_to_Pa = 2; % from dBFS of 100 to dBFS of 94 dB SPL
    
    t = (1:size(env_TN,1))/fs;
    figure;
    plot(t,factor_to_Pa*env_TN,'Color',rgb('Gray')); hold on;
    plot(t,factor_to_Pa*env_N,'Color',rgb('Gray')); hold on;
    plot(t,factor_to_Pa*mean(env_TN,2),'b-','LineWidth',2); hold on; grid on;
    plot(t,factor_to_Pa*mean(env_N,2),'r--','LineWidth',2);
    
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    
    Pos = get(gcf,'Position');
    Pos(4) = 200;
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2b-Hilbert';
end

if flags.do_fig3    
    dur_sil_bef = round(100e-3*fs);
    sil_bef = zeros([dur_sil_bef 1]);
    t_sil_offset = dur_sil_bef/fs;
    
    dur_sil_after = round(1300e-3*fs);
    sil_after = zeros([dur_sil_after 1]);
    
    basef = 500;
    out_TN     = zeros([size(noises_TN,1)+dur_sil_bef+dur_sil_after 1]);
    out_TN_thres = zeros(size(out_TN));
    out_TN_mod = zeros(size(out_TN));
    out_TN_mod_thres = zeros(size(out_TN));
     
    out_N = zeros(size(out_TN));
    out_N_mod = zeros(size(out_TN));
    
    idx_mfc = 1:5; % 3:5;
    for i = 1:N
        [out,fc,mfc] = king2019([sil_bef; noises_TN(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'no_mfb','subfs',fs);
        out = squeeze(out);
        out_TN = out_TN + out;
        
        [out,fc,mfc] = king2019([sil_bef; noises_TN_thres(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'no_mfb','subfs',fs);
        out = squeeze(out);
        out_TN_thres = out_TN_thres + out;
        
        [out,fc,mfc] = king2019([sil_bef; noises_N(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'no_mfb','subfs',fs);
        out = squeeze(out);
        out_N = out_N + out;
        
        [out,fc,mfc] = king2019([sil_bef; noises_TN(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'subfs',fs);
        out = squeeze(out);
        
        out_TN_mod = out_TN_mod + sum( out(:,idx_mfc) ,2);
        
        %%%
        [out,fc,mfc] = king2019([sil_bef; noises_TN_thres(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'subfs',fs);
        out = squeeze(out);
        
        out_TN_mod_thres = out_TN_mod_thres + sum( out(:,idx_mfc) ,2);
        %%%
        [out,fc,mfc] = king2019([sil_bef; noises_N(:,i); sil_after],fs,'basef',basef,'flow',basef,'fhigh',basef,'subfs',fs);
        out = squeeze(out);
        
        out_N_mod = out_N_mod + sum( out(:,idx_mfc) ,2);
    end
    out_TN = out_TN/N;
    out_TN_thres = out_TN_thres/N;
    out_N = out_N/N;
    out_TN_mod = out_TN_mod/N;
    out_TN_mod_thres = out_TN_mod_thres/N;
    out_N_mod = out_N_mod/N;
    
    t = (1:size(out_TN,1))/fs - t_sil_offset;
    
    figure;
    plot(t,out_TN_thres,'m-'); hold on;
    plot(t,out_TN,'b-'); hold on;
    plot(t,out_N,'r-');
    
    factor = 10^4;
    if factor ~= 1
        factor_str = sprintf(' x 10^{%.0f}',-log10(factor));
    else
        factor_str = '';
    end
    figure;
    Pos = get(gcf,'Position');
    Pos(4) = 300;
    set(gcf,'Position',Pos);
    
    subplot(5,1,[1 2 3]);
    pl1 = plot(t,factor*out_TN_mod_thres,'m-','LineWidth',2); hold on; grid on
    pl2 = plot(t,factor*out_TN_mod,'b-'); hold on;
    pl3 = plot(t,factor*out_N_mod,'r-');
    % ylabel(['Amplitude (a.u.' factor_str ')']);
    title('A. Average TN and N trials','FontSize',10);
    ha =gca;
    set(gca,'XTickLabels','');
    
    legend([pl2 pl1 pl3],{'T+N, SNR=0 dB (supra)','T+N, SNR=-18 dB','N alone'});
    
    subplot(5,1,[4 5]);
    plot(t,factor*(out_TN_mod_thres-out_N_mod),'m-','LineWidth',2); hold on; grid on
    plot(t,factor*(out_TN_mod-out_N_mod),'b-'); hold on;
    
    % ylabel(['Amplitude (a.u.' factor_str ')']);
    title('B. Difference between TN and N','FontSize',10);
    ha(end+1) = gca;
    linkaxes(ha,'x');
    xlim([-.05 1.205]);
    xlabel('Time (s)');
    
    han=axes(gcf,'visible','off'); 
    %han.Title.Visible='on';
    %han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,['Amplitude (a.u.' factor_str ')']);
    
    set(gca,'YTick',[-2:2]);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig3-model-rep';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5

    N_lambda = 30;
    Lambdas = logspace(-4, -1, N_lambda);
    idx = find(Lambdas >= 2*10^-3);
    Lambdas = Lambdas(idx);

    glmfct = 'l1glm';
    TF_type = 'gammatone';
    flags_in = {    'trialtype_analysis','total', ...
                    'N_folds',10, ...
                    'no_permutation', ...
                    'bias', ...
                    'plot', ...
                    'pyramid_script','imresize', ...
                    'pyramid_shape',0, ...
                    'lambda',Lambdas, ...
                    'bwmul',0.43};
	fname_results = [dir_exp 'king2019' filesep 'Results-run-1' filesep 'savegame_2023_04_19_01_08_king2019_toneinnoise_ahumada1975_white.mat'];
    [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results,TF_type,glmfct,flags_in{:});
    
    hname{end+1} = 'fig5-ACI';
    h(end+1) = gcf;
    
    title('')
    Pos = get(gcf,'Position');
    Pos(4) = 300;
    set(gcf,'Position',Pos);
end

% if flags.do_fig3b
%     %%% Only the final fits:
%     model = 'osses2021'; 
%     noise_types = {'SSN'}; 
%     folders = {'Results-run-4'};
% end
% 
% if flags.do_fig4
%     bProceed = publ_osses2021c_DAGA_0_checkdata; % Checking if the experimental data is on disk
%     
%     model{1} = 'osses2021c_S01';  folders1 = {'Results'}; 
%     model{2} = 'osses2021c_S02';  folders2 = {'Results'}; 
%     model{3} = 'osses2021';       folders3 = {'Results-run-3-m1p55','Results-run-1', ...
%         'Results-run-3-p0p39','Results-run-3-p0p78','Results-run-4'}; 
%     
%     N_plots = (length(folders1)+length(folders2)+length(folders3));
%     thres = nan(15,N_plots);
% end
% 
% if flags.do_fig1a || flags.do_fig1b || flags.do_fig4
%     if bProceed == 0
%         error('Please follow the instructions to download the experimental data before you can successfully run this script again...')
%     end
% end
% 
% count = 1;
% 
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

bSave = 1;
if nargout ~= 0
    data.h = h;
    data.hname = hname;
    bSave = 0;
end

if bSave
    dir_out = fastACI_paths('dir_output');
    if ~exist(dir_out,'dir')
        dir_out = [pwd filesep 'outputs' filesep];
        mkdir(dir_out);
    end

    for i = 1:length(h)
        opts = [];

        opts.format = 'epsc';
        Saveas(h(i),[dir_out hname{i}],opts);

        opts.format = 'png';
        Saveas(h(i),[dir_out hname{i}],opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_get_hilbert(insig,fs,basef);

cutofffreq=30;

outsig = auditoryfilterbank(insig,fs,'basef',basef,'flow',basef,'fhigh',basef);
outsig = ihcenvelope(outsig,fs,'ihc_king2019');

% Extra LPF filter at 30 Hz:

ihc_filter_order = 2;
[b, a] = butter(1, cutofffreq*2/fs);
for ii=1:ihc_filter_order
    outsig = filter(b,a, outsig);
end
metric = outsig;
description = 'BB waveform + HWR + LPF at 30 Hz';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig_handle,fig_name] = il_plot_the_ACI(ACI, cfg_ACI, cfg_game)

dir_where = cfg_game.dir_target;

h_outs = affichage_tf(ACI,'CI','cfg',cfg_ACI); hold on;

%%%
par_formants.timestep = 0.01; % positive timestep 0.01
par_formants.nformants = 5; % positive nformants 5
par_formants.maxformant = 5500; % positive maxformant 5500
par_formants.windowlength = 0.025; % positive windowlength 0.025
par_formants.dynamicrange = 30; % positive dynamic range 20

par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

par_formants.I_min = 50; % 60, arbitrary value
try
    outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
    if isempty(outs_from_Praat.filesf0)
        error('Files not found...')
    end
catch me
    fprintf('%s: Praat not found on disk, loading pre-stored Praat files\n',upper(mfilename));
    
    dir_where_stored = [fastACI_basepath 'Stimuli' filesep 'varnet2013' filesep]; 
    outs_from_Praat = Get_all_metrics_from_stored_Praat(dir_where_stored,par_formants);
end

idx = find(outs_from_Praat.t_f0{1} < 0.089 | (outs_from_Praat.t_f0{1} > 0.164 & outs_from_Praat.t_f0{1} < .280) | outs_from_Praat.t_f0{1} > 0.404);
outs_from_Praat.f0{1}(idx,:) = nan;
outs_from_Praat.f0{2}(idx,:) = nan;
outs_from_Praat.F{1}(idx,:) = nan;
outs_from_Praat.F{2}(idx,:) = nan;

%%%
Style  = {'-','--'};
Colour = {[0.5 0.5 0.5],'k'};
LW     = 1; 
outs_from_Praat = affichage_tf_add_Praat_metrics(dir_where, cfg_ACI,outs_from_Praat,Style,Colour,LW);

%%% Plotting
idx = strfind(cfg_ACI.fnameACI,filesep);
fname_ACI_short = cfg_ACI.fnameACI(idx(end)+1:end);

title(name2figname(fname_ACI_short));
fig_handle = gcf;
fig_name = name2figname(fname_ACI_short);

Pos = get(fig_handle,'Position');
Pos(3) = 400;
set(fig_handle,'Position',Pos);
