function data = publ_osses2023b_FA_figs(varargin)
% function publ_osses2023b_FA_figs(varargin)
%
% Generates the figures
%
% % To display Fig. 2 of Osses and Varnet, (2023, Forum Acusticum) use :::
%     publ_osses2023b_FA_figs('fig2a');
%     publ_osses2023b_FA_figs('fig2b');
%     publ_osses2023b_FA_figs('fig3');
%     publ_osses2023b_FA_figs('fig4');
%     publ_osses2023b_FA_figs('fig5');
%     publ_osses2023b_FA_figs('fig6');
%
% % If you have the data from Zenodo stored locally, you can use:
%     publ_osses2023b_FA_figs('fig2a','zenodo','dir_zenodo',dir_zenodo);
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_osses2023b_FA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2a','fig2b', ...
    'fig3', ... % Ahumada's data
    'fig4', ... % supra-threshold, model
    'fig5', ...
    'fig6'}; % ACI

definput.flags.local={'local','zenodo'};
definput.keyvals.dir_zenodo=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

bZenodo = flags.do_zenodo;
bLocal = ~bZenodo;

dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';

experiment = 'toneinnoise_ahumada1975';
 
if bLocal
    dir_exp = [dir_fastACI_results experiment filesep];
    dir_noise = [dir_exp 'king2019' filesep 'NoiseStims-white' filesep]; % common parameter
    dir_savegame = [dir_exp 'king2019' filesep 'Results-run-1' filesep];
end
if bZenodo
    dir_zenodo = keyvals.dir_zenodo;
    dir_noise = [dir_zenodo '01-Stimuli' filesep 'fastACI_data' filesep ...
         experiment filesep 'king2019' filesep 'NoiseStims-white' filesep];
    dir_savegame = [dir_zenodo '02-Raw-data' filesep 'fastACI' filesep 'Publications' filesep ...
         'publ_osses2023b' filesep 'data_king2019' filesep '1-experimental_results' filesep];
end
var = load([dir_savegame 'savegame_2023_04_19_01_08_king2019_toneinnoise_ahumada1975_white.mat']);
cfg_game       = var.cfg_game;
data_passation = var.data_passation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2a || flags.do_fig2b  || flags.do_fig4   
    files = Get_filenames(dir_noise,'*.wav');
    
    tone = cfg_game.signal;
    dBFS = cfg_game.dBFS;
    fs = cfg_game.fs;
    
    N = 1000;
    
    Ni = 1600;
    for i = 1:N
        [noises_TN(:,i),fs] = audioread([dir_noise files{i}]);
        
        noises_N(:,i) = audioread([dir_noise files{i+Ni}]);
    end
    
    % Obtaining the 'experimental' threshold:
    thres = median(data_passation.expvar);
    
    noises_TN_thres = noises_TN + From_dB(thres)*repmat(tone,[1 N]);
    noises_TN       = noises_TN +                repmat(tone,[1 N]); % SNR=0 dB
end

if flags.do_fig2a    
    t = (1:length(tone))/fs;
    toffset = 12e-3;
    idx = find(t>=.2+toffset & t<=.3-toffset);
            
    win_type = 'hanning';
    % win_type = 'rectangular';
    for i = 1:N
        K = fs/2; % length(idx); % fs/2;
        [x, y_dB_TN(:,i),f]   = freqfft2(noises_TN(idx,i),K,fs,win_type,dBFS);
        [x, y_dB_TN_m18(:,i)] = freqfft2(noises_TN_thres(idx,i),K,fs,win_type,dBFS);
        [x, y_dB_N(:,i)]      = freqfft2(noises_N(idx,i) ,K,fs,win_type,dBFS);
    end
    
    y_TN = median(y_dB_TN,2);
    y_TN_m18 = median(y_dB_TN_m18,2); % -18.2
    y_N  = median(y_dB_N,2);
    
    f_ERB = freqtoaud(f);
    
    figure;
    plot(f_ERB,y_TN_m18,'m-','LineWidth',2); hold on; grid on;
    plot(f_ERB,y_TN,'b-','LineWidth',2); 
    plot(f_ERB,y_N,'r--','LineWidth',2);
    % plot(f_ERB,x_dB,'m-','LineWidth',2);
    ylabel('Average amplitude (dB)');
    xlabel('Frequency (Hz)');
    XTL = [125 250 500 1000 2000];
    set(gca,'XTick',freqtoaud(XTL));
    set(gca,'XTickLabel',XTL);
    xlim(freqtoaud([100 2400]))
    
    %%%
    [max_val, idx_max] = max(y_TN_m18);
    xvar = [f_ERB(idx_max) f_ERB(idx_max)+2.5]; % extending up to +2 ERB_N
    yvar = max_val*[1 1];
    plot(xvar,yvar,'m--','LineWidth',2);
    text(xvar(1)+2,max_val,sprintf('tone at\nSNR=%.1f dB',thres),'Color','m','FontSize',12);
    
    
    [max_val, idx_max] = max(y_TN);
    xvar = [f_ERB(idx_max) f_ERB(idx_max)+2.5]; % extending up to +2 ERB_N
    yvar = max_val*[1 1];
    plot(xvar,yvar,'b--','LineWidth',2);
    text(xvar(1)+2,max_val,sprintf('tone at\nSNR=0 dB'),'Color','b','FontSize',12);
    %%%
    
    Pos = get(gcf,'Position');
    Pos(4) = 250;
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2a-FFT';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2b
    basef = 500;
	for i = 1:N
        env_TN(:,i) = il_get_hilbert(noises_TN(:,i),fs,basef);
        env_TN_m18(:,i) = il_get_hilbert(noises_TN_thres(:,i),fs,basef);
        env_N(:,i)  = il_get_hilbert(noises_N(:,i),fs,basef);
    end
    
    factor_to_Pa = 2; % from dBFS of 100 to dBFS of 94 dB SPL
    
    t = (1:size(env_TN,1))/fs;
    figure;
    plot(t,factor_to_Pa*env_TN,'Color',rgb('Gray')); hold on;
    plot(t,factor_to_Pa*env_N,'Color',rgb('Gray')); hold on;
    plot(t,factor_to_Pa*mean(env_TN_m18,2),'m-','LineWidth',2); hold on; grid on;
    plot(t,factor_to_Pa*mean(env_TN,2),'b-','LineWidth',2); 
    plot(t,factor_to_Pa*mean(env_N,2),'r--','LineWidth',2);
    
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    YL = get(gca,'YLim');
    yoff = .003;
    ylim([YL(1)-yoff YL(2)+yoff])
    Pos = get(gcf,'Position');
    Pos(4) = 150;
    set(gcf,'Position',Pos);
    
    YT = 0:.005:YL(2);
    set(gca,'YTick',YT);
    YTL = [];
    for i=1:length(YT)
        if mod(i,2)==1
            YTL{i} = num2str(YT(i));
        else
            YTL{i} = '';
        end
    end
    set(gca,'YTickLabel',YTL);
    h(end+1) = gcf;
    hname{end+1} = 'fig2b-Hilbert';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig3
    % Weights in pixels:
    Weights_RM = [  18  -9  7 -8 18; ... % 600
                    0    0  0  0  0; ... % 550
                   -14   0 67 22 15; ... % 500
                   -20   9  0  0  0; ... % 450 Hz
                     0  14  0  0  0];    % 400
                 
    Weights_KL = [   0   0 -16  9 13; ... % 600
                     9   0 -11  0 10; ... % 550
                     0   9  57 18 25; ... % 500
                     0   0 -12  0  0; ... % 450 Hz
                    -4  -3   0  0 -12];   % 400
                
    val_max = max([Weights_KL(:); Weights_RM(:)]);
    val_max = 1.2*val_max;
    
    Weights_RM = Weights_RM/val_max;
    Weights_KL = Weights_KL/val_max;
    
    YT  = [0 1 2 3 4];
    YTL = [400 450 500 550 600];
    
    offsety = repmat([4; 3; 2; 1; 0],1,5);
    
    ti = [0:.1:.4];
    figure;
    tiledlayout(1,2,'tilespacing','compact');
    nexttile(1);
    for i = 1:size(Weights_RM,1)
        y_var = Weights_RM(i,:)+offsety(i,:);
        plot(ti,y_var,'b-'); hold on; grid on;
        idx = find(y_var~=offsety(i,1));
        if ~isempty(idx)
            plot(ti(idx),y_var(idx),'bo','MarkerFaceColor','b'); hold on; grid on;
        end
    end
    ylim([-.5 4.5]);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YTL);
    set(gca,'XTick',ti);
    title('Participant RM');
    
    ylabel('Frequency (Hz)');
    
    %%%
    nexttile(2);
    for i = 1:size(Weights_KL,1)
        y_var = Weights_KL(i,:)+offsety(i,:);
        plot(ti,y_var,'b-'); hold on; grid on;
        idx = find(y_var~=offsety(i,1));
        if ~isempty(idx)
            plot(ti(idx),y_var(idx),'bo','MarkerFaceColor','b'); hold on; grid on;
        end
    end
    ylim([-.5 4.5]);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel','');
    set(gca,'XTick',ti);
    title('Participant KL');
            
    han=axes(gcf,'visible','off'); 
    han.XLabel.Visible='on';
    
    xlabel('Starting time of the segment (s)');
    
    Pos = get(gcf,'Position');
    Pos(3:4) = [650 250];
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig3-Ahumada-data';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig4    
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
    
    % figure;
    % plot(t,out_TN_thres,'m-'); hold on;
    % plot(t,out_TN,'b-'); hold on;
    % plot(t,out_N,'r-');
    % 
    % factor = 10^4;
    % if factor ~= 1
    %     factor_str = sprintf(' x 10^{%.0f}',-log10(factor));
    % else
    %     factor_str = '';
    % end
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
    hname{end+1} = 'fig4-model-rep';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5
    expvar     = data_passation.expvar;
    is_correct = data_passation.is_correct;
    
    thres = median(expvar);
    
    var_min = floor(min(expvar)-.5);
    var_max = max(expvar);
    var_step = 2; % dB
    
    SNR_edges = var_min:var_step:var_max;
    
    % idx = find(expvar<=SNR_edges(1)); % below first edge there is nothing
    for i_edge = 1:length(SNR_edges)-1
        idx = find(expvar>=SNR_edges(i_edge) & expvar<SNR_edges(i_edge+1));
        N_bin(i_edge) = length(idx);
        correct_bin(i_edge) = sum(is_correct(idx))/N_bin(i_edge);
        SNR_bin(i_edge) = mean(SNR_edges(i_edge:i_edge+1));
    end
    
    idx = find(expvar>=SNR_edges(i_edge+1));
    N_bin(i_edge+1) = length(idx);
    correct_bin(i_edge+1) = sum(is_correct(idx))/N_bin(i_edge+1);
    SNR_bin(i_edge+1) = SNR_edges(i_edge+1)+var_step/2;
    
    SNR_edges(end+1) = SNR_edges(end)+var_step;
    % correct_bin(end+1) = correct_bin(end);
    
    figure;
    bar(SNR_bin,100*correct_bin,'BarWidth',1); grid on; hold on; % BarWidth in relative units
    xlim([SNR_edges(1)-2 SNR_edges(end)]);
    
    XT = SNR_edges(1):SNR_edges(end);
    for i = 1:length(XT)
        if mod(i,2)==1
            XTL{i} = '';
        else
            XTL{i} = num2str(XT(i));
        end
    end
    XT = [SNR_edges(1)-1 XT];
    set(gca,'XTick',XT);
    
    XTL = ['avg',XTL];
    
    set(gca,'XTickLabel',XTL);
    
    Pos = get(gcf,'Position');
    Pos(4) = 230; % 350;
    set(gcf,'Position',Pos);
    ylabel('Percentage correct (%)');
    xlabel('expvar, SNR (dB)');
    
    YL = get(gca,'YLim');
    plot(thres*[1 1],YL,'r--');
    
    N_session = 400;
    idx2use = 51:400; % visual inspection
    
    expvar_buf = reshape(expvar,[N_session, length(expvar)/N_session]);
    is_correct_buf = reshape(is_correct,[N_session, length(is_correct)/N_session]);
    score_session = 100*mean(is_correct_buf(idx2use,:));
    % x_var = (SNR_edges(end)+1)*ones(size(score_session));
    % plot(x_var,score_session,'bo-','LineWidth',2);
    Me = mean(score_session);
    errL = Me-prctile(score_session,5);
    errU = prctile(score_session,95)-Me;
    x_var = SNR_edges(1)-var_step/2;
    errorbar(x_var,Me,errL,errU,'ro-','LineWidth',2,'MarkerFaceColor','w');
    
    XL = get(gca,'XLim');
    plot(XL,50*[1 1],'k--','LineWidth',2);
    text(x_var-.5,50+4,sprintf('chance'),'FontWeight','Bold','FontSize',9);
    h(end+1) = gcf;
    hname{end+1} = 'fig5-performance';
    ylim([37 103])
    set(gca,'YTick',40:10:100);
    disp('')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig6
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
                    'lambda',Lambdas}; % , ...
                    % 'bwmul',0.43};
    if bZenodo
        dir_out_ACI = [keyvals.dir_zenodo filesep '03-Post-proc-data' filesep];
        flags_in(end+1:end+2) = {'dir_out',dir_out_ACI};
    end
	fname_results = [dir_savegame 'savegame_2023_04_19_01_08_king2019_toneinnoise_ahumada1975_white.mat'];
    [ACI,cfg_ACI,results,Data_matrix,extra_outs] = fastACI_getACI(fname_results,TF_type,glmfct,flags_in{:});
    
    hname{end+1} = 'fig6-ACI';
    h(end+1) = gcf;
    
    title('')
    Pos = get(gcf,'Position');
    Pos(4) = 300;
    set(gcf,'Position',Pos);
    
    tcolourbar = extra_outs.out_affichage.tcolorbar;
    set(tcolourbar,'TickLabels','');
    
    colourbar_map = 'DG_jet';
    my_map = Get_colourmap_rgb(colourbar_map);
    col_max = my_map(1,:);
    col_min = my_map(end,:);
    text(.95,1.05,cfg_ACI.response_names{1},'Units','Normalized','FontWeight','Bold','Color',col_min);
    text(.95,-.05,cfg_ACI.response_names{2},'Units','Normalized','FontWeight','Bold','Color',col_max);

    set(gca,'XTick',.1:.1:.4);
end

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
