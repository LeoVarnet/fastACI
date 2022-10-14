function data = publ_osses2022b_JASA_figs(varargin)
% function data = publ_osses2022b_JASA_figs(varargin)
%
% Generates the figures
%
% % To display Fig. 1 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig1'); % EM, MM, and IM for /aba/, /ada/
%
% % To display Fig. 2 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig2'); % T-F representations /aba/, /ada/
%
% % To display Fig. 4 of Osses and Varnet, (2022, JASA) use :::
%     publ_osses2022b_JASA_figs('fig4'); % illustration of the Gaussian basis 
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all, clc

if nargin == 0
    help publ_osses2022b_JASA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig1_old', ... % Schematic of EM, IM, and MM
    'fig1', ...
    'fig2a', ... % spectrogram white, bump, MPS noises
    'fig2b', ... % Other metrics to characterise white, bump, and MPS noises
    'fig3', ...  % Illustration of the Gaussian basis elements
    'fig4', ...
    'fig5_simu', ... % old name: fig4_simu
    'fig5', 'fig5_stats', ... % dprime and statistics
    'fig6', ... % ACIs participants, L1109
    'fig7', ...
    'fig8', ... % Cross predictions, L1608
    'fig8b_old', ... % Discarded significance test
    'fig9', ... % Cross predictions between noise
    'fig10', ...
    'fig8b','fig11_old', ... % Cross predictions 'between participant' for simulations
    'fig12', ... % Cross predictions between noise for simulations
    'fig1_suppl', ... % Audiograms
    'fig2_suppl', ... % Intellitest, L3209 (on 29/09/2022)
    'fig3_suppl', ...
    'fig3b_suppl', ...
    'fig5_suppl'}; % different lambdas
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';
% 
% experiment = 'speechACI_varnet2013';
% dir_exp = [dir_fastACI_results experiment filesep];

% Common variables:
noise_str         = {'white','bumpv1p2_10dB','sMPSv1p3'};
noise_types_label = {'white','bump','MPS'};
masker_colours    = {[0,0,1],[1,0,0],rgb('Green')}; % Leo's colours

experiment        = 'speechACI_Logatome-abda-S43M';
N_maskers = length(noise_str);

do_fig4_stats_N10 = 0; % by default no analysis of N=10 data
do_fig5_stats_N10 = 0; % by default no analysis of N=10 data
% End common variables
    
if flags.do_fig1_old || flags.do_fig1 || flags.do_fig2a || flags.do_fig2b
    dir_data = fastACI_paths('dir_data');
    dir_subj = [dir_data experiment filesep 'S01' filesep];
    
    dir_noise = [dir_subj 'NoiseStim-white' filesep];
    
    basef = 8000;
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5, ...
        'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};
    
    CLim = [-8 -4.5];
    
    flags_extra = {'NfrequencyTicks',8,'colorbar','no'};
end

if flags.do_fig1_old % Speech-in-noise representation
    
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    [aba, fs] = audioread([dir_subj 'speech-samples' filesep 'S43M_ab_ba.wav']);
    [ada, fs] = audioread([dir_subj 'speech-samples' filesep 'S43M_ad_da.wav']);
    
    [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
    G_ada                = Gammatone_proc(ada, fs, flags_gamma{:});
    
    N_noises = 1000; % arbitrary number
    % SNR = -14; % dB, arbitrary value
    % gain_SNR = 10^(SNR/20);
    % aba_at_SNR = gain_SNR * aba;
        
    fname_noises = Get_filenames(dir_noise);
    fname_noises = fname_noises(1:N_noises); % names of the first N_noises
    
    G_Naba = nan([size(G_aba) N_noises]); % memory allocation
    G_Nada = nan(size(G_Naba));
    G_N    = nan(size(G_Naba));
    
    for i_noise = 1:N_noises
        [noise, fs] = audioread([dir_noise fname_noises{i_noise}]);

        noisy_aba = noise+aba;
        noisy_ada = noise+ada;
        [G_Naba(:,:,i_noise), fc, t, outs] = Gammatone_proc(noisy_aba, fs, flags_gamma{:});
        G_Nada(:,:,i_noise)                = Gammatone_proc(noisy_ada, fs, flags_gamma{:});
        G_N(:,:,i_noise)                   = Gammatone_proc(noise    , fs, flags_gamma{:});
    end
    
    %%%
    G_Naba_avg = mean(G_Naba,3);
    G_Nada_avg = mean(G_Nada,3);
    G_N_avg    = mean(G_N,3);
    G_N_avg_std = G_N_avg + 2.2*std(G_N,[],3);
    
    % G_Naba_avg_thres = G_Naba_avg;
    %G_Naba_avg_thres(G_Naba_mean_thres<G_N_meanstd) = nan;

    idx_realisation = round(.733*N_noises);
    if N_noises ~= 1000
        warning('Less noise representations being used to derive the masking effects...');
    end
    G_Naba_thres = G_Naba(:,:,idx_realisation);
    G_Nada_thres = G_Nada(:,:,idx_realisation);
    %G_Naba_thres(G_Naba_thres<G_N_meanstd) = nan;

    N_targets = 2; % default is 1
    if N_targets ~= 1
        correction = .9; % slightly smaller
    else
        correction = 1;
    end
    figure('Position',[100 100 1000 250*N_targets*correction]); 
    il_tiledlayout(N_targets,5,'TileSpacing','Compact'); % 'tight'); % none'); % 'Compact');
    
    x1 = -0.05; y1 = 1.05;
    x2 =  0.55; y2 = 0.92;
    x3 =  0.95; y3 = 0.85;
    
    txt_targets = {'/aba/','/ada/'};
    for i = 1:N_targets
        
        txt_target = txt_targets{i};
        
        switch i 
            case 1
                G = G_aba;
            case 2
                G = G_ada;
        end
            
        G_dB = il_To_dB(G');
        il_nexttile(1+(i-1)*5); 
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); 
        caxis(CLim); 
        title(['Clean ' txt_target]);
        if i == 1
            text(x1,y1,'A','Units','Normalized','FontSize',12,'FontWeight','Bold');
            set(gca,'XTickLabels',[]);
            set(gca,'XLabel',[]);
        end
        % colorbar off;

        switch i 
            case 1
                G = G_Naba;
            case 2
                G = G_Nada;
        end
        %nexttile; affichage_tf(log(G_ada)', 'pow', t, fc); caxis([-8 -4.5]); title('/ada/ target');colorbar off;ylabel('');set(gca,'YTickLabels',[]);xlabel('');set(gca,'XTickLabels',[]);%
        il_nexttile(2+(i-1)*5); 
        G_dB = il_To_dB(G(:,:,1)');
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); 
        caxis(CLim); 
        title(['Noisy ' txt_target],'FontSize',8);
        if i == 1
            text(x1,y1,'B','Units','Normalized','FontSize',12,'FontWeight','Bold');
            set(gca,'XTickLabels',[]);
            set(gca,'XLabel',[]);
        end
        text(x2,y2,'N=1','Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
        % colorbar off;
        ylabel('');
        set(gca,'YTickLabels',[]);

        switch i 
            case 1
                G_avg = G_Naba_avg;
            case 2
                G_avg = G_Nada_avg;
        end
        il_nexttile(3+(i-1)*5); 
        G_dB = il_To_dB(G_avg');
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); 
        caxis(CLim); 
        title([txt_target ' and EM'],'FontSize',8);
        if i == 1
            text(x1,y1,'C','Units','Normalized','FontSize',12,'FontWeight','Bold');
            set(gca,'XTickLabels',[]);
            set(gca,'XLabel',[]);
        end
        text(x2,y2,sprintf('N=%.0f',N_noises),'Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
        text(x3,y3,'(floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');
        % colorbar off;
        ylabel('');
        set(gca,'YTickLabels',[]);%

        il_nexttile(4+(i-1)*5); 
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); % G_dB was not refreshed, but the 'mod threshold' will be drawn later
        caxis(CLim); 
        title([txt_target ' and EM+MM'],'FontSize',8);
        if i == 1
            text(x1,y1,'D','Units','Normalized','FontSize',12,'FontWeight','Bold');
            set(gca,'XTickLabels',[]);
            set(gca,'XLabel',[]);
        end
        text(x2,y2,sprintf('N=%.0f',N_noises),'Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
        text(x3,y3,'(floor+mod.floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');

        % colorbar off;
        ylabel('');
        set(gca,'YTickLabels',[]);%

        switch i 
            case 1
                G = G_Naba_thres;
            case 2
                G = G_Nada_thres;
        end
        il_nexttile(5+(i-1)*5);
        G_dB = il_To_dB(G');
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); 
        caxis(CLim); 
        title(['    ' txt_target ' and EM+MM+IM'],'FontSize',8);
        if i == 1
            text(x1,y1,'E','Units','Normalized','FontSize',12,'FontWeight','Bold');
            set(gca,'XTickLabels',[]);
            set(gca,'XLabel',[]);
        end
        text(x2,y2,'N=1','Units','Normalized','FontSize',10,'FontWeight','Bold','Color','white');
        text(x3,y3,'(floor+mod.floor)','Units','Normalized','FontSize',7,'FontWeight','Bold','Color','white','HorizontalAlignment','right');

        ylabel('');
        set(gca,'YTickLabels',[]);%

        listImage = findobj('Type','Image');
        set(listImage(2),'AlphaData',1-0.5*(G_avg<G_N_avg_std)');   % before last panel
        set(listImage(1),'AlphaData',1-0.5*(G    <G_N_avg_std)'); % last panel 
    end
    
    h(end+1) = gcf;
    hname{end+1} = 'fig1-schematic-masking';
end

if flags.do_fig1
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    fname_aba = [dir_subj 'speech-samples' filesep 'S43M_ab_ba.wav'];
    fname_ada = [dir_subj 'speech-samples' filesep 'S43M_ad_da.wav'];
    [aba, fs] = audioread(fname_aba);
    [ada, fs] = audioread(fname_ada);
    
    %%% Plot targets + 1 example of noise
    [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
    [G_ada, fc, t, outs] = Gammatone_proc(ada, fs, flags_gamma{:});

    warning('Conceptually wrong: The spectrograms are subjected to a natural logarithm...')
    disp('Pausing for 5 s')
    pause(5)
    
    figure('Position',[100 100 500 280]); 
    il_tiledlayout(1,2,'TileSpacing','tight');
    il_nexttile(1); 
    G_dB = il_To_dB(G_aba');
    affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); hold on;
    XTL = 0:.2:.8;
    XT  = XTL; XT(1) = min(t)/2;
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',XTL);
    caxis(CLim); 
    title('/aba/ target');
        
    cfg_in = [];
    cfg_in.f = fc;
    cfg_in.t = t;
    il_add_formants(fname_aba,cfg_in);
    
    % flags_extra{end} = 'on';
    il_nexttile(2); 
    G_dB = il_To_dB(G_ada');
    affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); hold on;
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',XTL);
    
    caxis(CLim); 
    title('/ada/ target');
    ylabel('');
    set(gca,'YTickLabels',[]); %xlabel('')

    % fname1 = [cfg_ACI.dir_target files{1}];
	% fname2 = [cfg_ACI.dir_target files{2}];
           
    il_add_formants(fname_ada,cfg_in);
    
    text(0.52,0.08,'f_0','Color','white','Units','Normalized');
    text(0.52,0.35,'F_1','Color','white','Units','Normalized');
    text(0.52,0.48,'F_2','Color','white','Units','Normalized');
    text(0.52,0.62,'F_3','Color','white','Units','Normalized');
    text(0.52,0.75,'F_4','Color','white','Units','Normalized');
    
    h(end+1) = gcf;
    hname{end+1} = 'fig1-spec-targets';
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2a
    
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5/4, ...
        'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};
    fname = 'Noise_00101.wav';
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    fname_no1 = [dir_subj 'NoiseStim-white'         filesep fname];
    fname_no2 = [dir_subj 'NoiseStim-bumpv1p2_10dB' filesep fname];
    fname_no3 = [dir_subj 'NoiseStim-sMPSv1p3'      filesep fname];
    
    [insig(:,1),fs] = audioread(fname_no1);
    insig(:,2) = audioread(fname_no2);
    insig(:,3) = audioread(fname_no3);
    
    figure('Position',[100 100 500 280]); 
    il_tiledlayout(1,3,'TileSpacing','tight');
    
    for i = 1:3
        %%% Plot targets + 1 example of noise
        [G, fc, t, outs] = Gammatone_proc(insig(:,i), fs, flags_gamma{:});
     
        dBFS = 100;
        floor2use = 0;
        G_dB(:,:,i) = max(20*log10(G)'+dBFS,floor2use); % il_To_dB(G');
    end
    
    dB_max = max(max(max(G_dB)));
    % G_dB(find(G_dB> dB_max))= dB_max;
    % G_dB(:,end,1) = 0;
    % G_dB(:,end,2) = 0;
    % G_dB(:,end,3) = 0;
    idx_t = find(t>=0.1 & t<=t(end)-100e-3);
    idx_f = find(round(fc>200));
    for i = 1:3
        if i ~= 3
            flags_extra = {'NfrequencyTicks',8,'colorbar','no'};
        else
            flags_extra = {'NfrequencyTicks',8,'colorbar','yes'};
        end
        
        il_nexttile(i); 
        affichage_tf(G_dB(idx_f,idx_t,i), 'pow', t(idx_t), fc(idx_f), flags_extra{:}); hold on;
        caxis([0 dB_max])
        
        if i ~= 1
            set(gca,'YTickLabel',[]);
            ylabel('');
        end
        title(noise_types_label{i})
    end
    disp('')
    h(end+1) = gcf;
    hname{end+1} = 'fig2a-spec-noises';   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2b 
	N_sounds = 1000; % arbitrary choice
        
    % Band levels are always computed
    do_kohlrausch2021_env = 1; % flags.do_fig2c; % Modulation spectrum
	do_V                  = 1; % flags.do_fig2d; % V metric
    fmod_xlim = [0 60];
    fc_ylim = [30 70];
    V_lim = [-6.8 .8];
    %%%
    
    figure('Position',[100 100 650 650]); % Pos(4) = 840 
    il_tiledlayout(3,3,'TileSpacing','tight');
    
    % i replaced by i_noise
	for i_noise = 1:length(noise_str)
        dir_where = [dir_subj 'NoiseStim-' noise_str{i_noise} filesep];
        suff = noise_types_label{i_noise};
        
        dBFS = 100;
        files1 = Get_filenames(dir_where,'*.wav');
        files1 = files1(1:N_sounds);
        lvls = [];
        
        for j = 1:N_sounds
            file = files1{j};
            [insig,fs] = audioread([dir_where file]);

            if do_kohlrausch2021_env
                [env_dB_full(:,j),xx,env_extra] = Get_envelope_metric(insig,fs, ...
                    'kohlrausch2021_env_noDC'); suff_DC = ''; fmod_ylim = [10 50];
                % [env_dB_full(:,j),xx,env_extra] = Get_envelope_metric(insig,fs, ...
                %     'kohlrausch2021_env'); suff_DC = '-with-DC'; fmod_ylim = [10 70];
            end
 
            [outsig1,fc] = auditoryfilterbank(insig,fs);
            if j == 1
                t = (1:size(outsig1,1))/fs;
            end
            lvls(j,:) = rmsdb(outsig1) + dBFS; % These are the band levels
 
            for i_fc = 1:length(fc)
                if do_V
                    [V1(j,i_fc),description,yenv1] = Get_envelope_metric(outsig1(:,i_fc),fs,'V');
                end
            end
            V_overall(j) = Get_envelope_metric(insig,fs,'V');
            
            if mod(j,50) == 1
                fprintf('\tProcessing sound %.0f of %.0f\n',j,N_sounds);
            end
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Band levels:
        L1me = prctile(lvls,50);
        L1perL = prctile(lvls, 5);
        L1perU = prctile(lvls,95);
        % errLL1 = L1me - prctile(lvls,5);
        % errUL1 = prctile(lvls,95)-L1me;

        il_nexttile(i_noise)
        semilogx(fc,L1perU,'Color',rgb('Gray')); hold on;
        plot(fc,L1perL,'Color',rgb('Gray'));
        plot(fc,L1me,'o-','Color',rgb('Maroon'),'MarkerFaceColor',rgb('Maroon')); grid on
        % errorbar(fc,L1me,errLL1,errUL1,'-','Color',Colours{i});
               
        if i_noise == 1
            data.L1_freq = fc;
        end
        data.L1me(:,i_noise) = L1me(:);
        data.LperL(:,i_noise) = L1perL(:);
        data.LperU(:,i_noise) = L1perU(:);
        
        if i_noise == 1
            ylabel('Band level (dB)')
        else
            set(gca,'YTickLabel',[]);
        end
        if i_noise == 2
            xlabel('Frequency (Hz)')
        end

        XT = [125 250 500 1000 2000 4000 8000];
        for i_freq=1:length(XT)
            if XT(i_freq)<1000
                XTL{i_freq} = num2str(XT(i_freq));
            else
                XTL{i_freq} = [num2str(XT(i_freq)/1000) 'k'];
            end
        end
        set(gca,'XTick',XT);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTick',fc_ylim(1)+5:5:fc_ylim(2)-5);
        
        ylim(fc_ylim);

        if i_noise == 1
            text(0,1.05,'B.','Units','Normalize','FontWeight','Bold','FontSize',13);
        end
        title(noise_types_label{i_noise});

        if do_V
            %%% Metric V
            V1me = prctile(V1,50);
            V1perL = prctile(V1, 5);
            V1perU = prctile(V1,95);
            % errL1 = V1me - prctile(V1,5);
            % errU1 = prctile(V1,95)-V1me;

            il_nexttile(3+i_noise)
            semilogx(fc,V1perU,'Color',rgb('Gray'),'LineWidth',2); hold on
            plot(fc,V1perL,'Color',rgb('Gray'),'LineWidth',2);
            plot(fc,V1me,'o-','Color',rgb('Maroon'),'MarkerFaceColor',rgb('Maroon')); hold on; grid on
            % errorbar(fc,V1me,errL1,errU1,'-','Color',Colours{i});

            if i_noise == 1
                ylabel('V (dB)')
            else
                set(gca,'YTickLabel',[]);
            end
            if i_noise == 2
                xlabel('Frequency (Hz)')
            end

            set(gca,'XTick',XT);
            set(gca,'XTickLabel',XTL);

            ylim(V_lim);
            
            if i_noise == 1
                text(0,1.05,'C.','Units','Normalize','FontWeight','Bold','FontSize',13);
            end
            % title(' ')
            % title(noise_types_label{i_noise});
            
            data.V1me(:,i_noise) = V1me(:);
            data.V1perL(:,i_noise) = V1perL(:);
            data.V1perU(:,i_noise) = V1perU(:);
            if i_noise == 1
                data.V1_freq = fc;
            end
        end

        if do_kohlrausch2021_env
            idx_env = find(env_extra.f_env<=fmod_xlim(2));
            
            f_env  = env_extra.f_env(idx_env);
            % extra.fs_env = env_extra.fs_env;
            env_dB_U = prctile(env_dB_full(idx_env,:),95,2);
            env_dB   = prctile(env_dB_full(idx_env,:),50,2);
            env_dB_L = prctile(env_dB_full(idx_env,:), 5,2);

            il_nexttile(6+i_noise)
            plot(f_env,env_dB_U,'-','Color',rgb('Gray'),'LineWidth',2); hold on;
            plot(f_env,env_dB_L,'-','Color',rgb('Gray'),'LineWidth',2);
            plot(f_env,env_dB  ,'-','Color',rgb('Maroon'),'LineWidth',2); grid on
            
            if i_noise == 1
                data.f_env = f_env;
            end
            data.env_dB_U(:,i_noise) = env_dB_U(:);
            data.env_dB_L(:,i_noise) = env_dB_L(:);
            data.env_dB(:,i_noise)   = env_dB(:);
            
            if i_noise == 1
                ylabel('Envelope spectrum (dB)')
            else
                set(gca,'YTickLabel',[]);
            end
            if i_noise==2
                xlabel('Modulation frequency (Hz)')
            end

            if i_noise == 1
                text(0,1.05,'D.','Units','Normalize','FontWeight','Bold','FontSize',13);
            end
            % title(' ')
            % title(noise_types_label{i_noise});
            
            deltaf = 10;
            XT = deltaf:deltaf:fmod_xlim(2)-deltaf;
            set(gca,'XTick',XT);
            set(gca,'YTick',fmod_ylim(1)+5:5:fmod_ylim(2)-5);
            xlim(fmod_xlim);
            ylim(fmod_ylim);

            % h(end+1) = gcf;
            % hname{end+1} = sprintf('Env-N-%.0f%s',N_sounds,suff);
        end
 
        data.V_overall(:,i_noise) = V_overall(:);
    end
    h(end+1) = gcf;
    hname{end+1} = sprintf('fig2b-analysis-noises-N-%.0f%s',N_sounds,suff_DC);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig3
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
    hname{end+1} = 'fig3-Gaussbasis';    
end

if flags.do_fig4 || flags.do_fig5 || flags.do_fig5_stats || flags.do_fig6 || ...
   flags.do_fig7 || flags.do_fig8 || flags.do_fig9       || flags.do_fig10 || ...
   flags.do_fig1_suppl || flags.do_fig2_suppl || flags.do_fig3_suppl || flags.do_fig3b_suppl || flags.do_fig8b_old

    %  g20220720_behavstats_from_l20220325
    %%% Analysis for 12 participants:
    Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'}; 
    N_subjects = length(Subjects);
end

if flags.do_fig5_simu || flags.do_fig12
    Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'}; 
    N_subjects = length(Subjects);
end

if flags.do_fig4 || flags.do_fig5_simu || flags.do_fig5 || flags.do_fig5_stats
    %%% If analysis for the first 10 participants: 
    % % Subjects = {'S01','S02','S03','S04','S05','S07','S08','S09','S10','S11'};
    
    % % dir_experiment = [fastACI_dir_data experiment filesep]; % this is the default directory
    % dir_experiment = ['/home/alejandro/Desktop/fastACI_today/fastACI_data/' experiment filesep]; % /S01/Results/Results_final

    Markers = {'o','s','^'};
    LW = [1 1 2];

    bLeo = 0;
    bAlejandro = ~bLeo;
    bExclude = 1;

    for i_subject = 1:N_subjects
        subject = Subjects{i_subject};
        
        if flags.do_fig4 || flags.do_fig5 || flags.do_fig5_stats % Human listeners:
            dir_subject = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep 'data_' Subjects{i_subject} filesep '1-experimental_results' filesep];
        end
        if flags.do_fig5_simu
            model = 'osses2022a_debug_0715';
            dir_subject = [fastACI_dir_data 'speechACI_Logatome-abda-S43M' filesep ...
                model filesep 'Results-' Subjects{i_subject} '-v1' filesep];
        end
        % dir_subject = [dir_experiment subject filesep];

        % switch subject
        %     case 'osses2021' % if model
        %         YL_fig1_2 = [-30 -10];
        %     otherwise % if a human participant
        %         YL_fig1_2 = [-20 0];
        % end

        for i_masker = 1:N_maskers
            masker      = noise_str{i_masker};
            dir_results = dir_subject;
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            if length(fname_results)~= 1
                error('Multiple savegame files were found, but it should only be one...')
            end
            fname_results = fname_results{end};

            cfg_game = [];
            data_passation = [];
            load([dir_results fname_results]);
            %subj = Get_subjectname_from_dirname(cfg_game.dir_noise);

            n_responses = data_passation.n_responses;
            N_trials    = length(n_responses); % completed number of trials
            % n_signal    = data_passation.n_targets(1:N_trials);
            SNR         = data_passation.expvar;

            resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
            N_sessions  = length(resume_trial)-1;

            % %%% Preparing for Fig. 1 and 2:
            % data_passation_tmp = data_passation;
            % 
            % 
            % bPlot_in_Script3 = 0;
            % % with option 'histogram-group' the histograms are defined on a
            % % common scale from -20 dB to -5 dB
            % data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation_tmp,'histogram-group',bPlot_in_Script3);

            % Fig. 2: Reversal analysis ---------------------------------------

            if N_sessions > 10
                idx_odd = find(mod(resume_trial(1:end-1)-1,400)~=0);
                fprintf('Participant %s has %.0f shorter sessions\n',subject,length(idx_odd));
                fprintf('\tI am going to merge that session with the previous one...\n',subject,length(idx_odd));
                for i_odd = 1:length(idx_odd)
                    fprintf('\t(session %.0f started at trial %.0f)\n',idx_odd(i_odd),resume_trial(idx_odd));
                end
                resume_trial(idx_odd) = [];
                N_sessions = length(resume_trial)-1;

                disp('')
            end
            
            idxs2use_all = [];
            for i_session = 1:N_sessions

                %%% Leo's way: Analysis only using the reversals: -------------
                % sessions are redefined as blocks of 400 trials
                idxi = 1+400*(i_session-1);
                idxf = 400*(i_session);
                try
                    [rev, idx] = Get_mAFC_reversals(SNR(idxi:idxf));
                catch
                    disp('')
                end

                % Then uses median:
                medianSNR_old(i_subject, i_masker, i_session) = median(rev);
                percSNR_L_old(i_subject, i_masker, i_session) = prctile(rev,5);
                percSNR_U_old(i_subject, i_masker, i_session) = prctile(rev,95);
                iscorr = data_passation.is_correct(idx);
                perc_corr_old(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);

                %%% Alejandro's way:
                idxi = resume_trial(i_session);
                idxf = resume_trial(i_session+1)-1;

                idxi_rev = idx(4); % reversal number 4 (including it) and later
                idxi100 = idxf-100; % reversal number 4 (including it) and later
                
                idxs2use = idxi+idxi_rev:idxf;
                idxs2use_all = [idxs2use_all idxs2use]; % collating the idxs of the kept trials
                SNR_here = SNR(idxs2use);
                medianSNR(i_subject, i_masker, i_session) = median(SNR_here);
                percSNR_L(i_subject, i_masker, i_session) = prctile(SNR_here,5);
                percSNR_U(i_subject, i_masker, i_session) = prctile(SNR_here,95);

                iscorr = data_passation.is_correct(idxs2use);
                perc_corr(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
                
                iscorr = data_passation.is_correct(idxi:idxf);
                perc_corr_all(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
                
                iscorr = data_passation.is_correct(idxi_rev:idxf);
                perc_corr_rev(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);

                iscorr = data_passation.is_correct(idxi100:idxf);
                perc_corr_100(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
                disp('')
            end

            %%% Preparing for Fig. 1 and 2:
            data_passation_tmp = data_passation;
            data_passation_tmp.expvar = data_passation_tmp.expvar(idxs2use_all);
            data_passation_tmp.n_responses = data_passation_tmp.n_responses(idxs2use_all);
            data_passation_tmp.n_targets = data_passation_tmp.n_targets(idxs2use_all);
            bPlot_in_Script3 = 0;
            % with option 'histogram-group' the histograms are defined on a
            % common scale from -20 dB to -5 dB
            data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation_tmp,'histogram-group',bPlot_in_Script3);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fig. 2: Performance as a function of SNR
            bin_centres = data_hist.bin_centres; 
            P_H(i_subject, i_masker, :)       = data_hist.H ./(data_hist.H +data_hist.M );
            P_FA(i_subject, i_masker, :)      = data_hist.FA./(data_hist.CR+data_hist.FA);
            dprime(i_subject, i_masker, :)    =   norminv(P_H(i_subject, i_masker, :)) - norminv(P_FA(i_subject, i_masker, :));    % Eq.  9 from Harvey2004
            criterion(i_subject, i_masker, :) = -(norminv(P_H(i_subject, i_masker, :)) + norminv(P_FA(i_subject, i_masker, :)))/2; % Eq. 12 from Harvey2004

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N_hist(i_subject, i_masker, :) = data_hist.N_hist;
            perc_correct(i_subject, i_masker, :) = 100*data_hist.prop_correct;
            perc_bias(i_subject, i_masker, :)    = 100*data_hist.bias_if_1;    

            if flags.do_fig4 || flags.do_fig5_simu
                fprintf('Median percentage correct, subj=%.0f, noise=%.0f: all rev=%.1f; after rev 4=%.1f; last 100t=%.1f\n', ...
                    i_subject,i_masker, ...
                    100*median(perc_corr_all(i_subject, i_masker, :)), ...
                    100*median(perc_corr_rev(i_subject, i_masker, :)), ...
                    100*median(perc_corr_100(i_subject, i_masker, :)) );
            end
        end

        if bLeo 
            medianSNR = medianSNR_old;
            % percSNR_L_old(i_subject, i_masker, i_session) = prctile(rev,5);
            % percSNR_U_old(i_subject, i_masker, i_session) = prctile(rev,95);
            % iscorr = data_passation.is_correct(idx);
            % perc_corr_old(
        end
    end
    
    if flags.do_fig4 || flags.do_fig5_simu
        if bAlejandro
            % Summary for the text:
            for i_subject = 1:N_subjects
                perc_here = squeeze(perc_corr_all(i_subject, :, :));
                Me_all(i_subject) = median(perc_here(:));
                Mi_all(i_subject) = min(perc_here(:));
                Ma_all(i_subject) = max(perc_here(:));
            end
            [Mi,idx_mi] = min(Me_all);
            [Ma,idx_ma] = max(Me_all);

            fprintf('...percentage correct varied between %.1f (S%.0f) and %.1f (S%.0f)\n', ...
                100*Mi,idx_mi,100*Ma,idx_ma);
            fprintf('...session by session scores between %.1f and %.1f\n', ...
                100*min(Mi_all),100*max(Ma_all));
            disp('')
        end
    end

    % % Fig. 0: percentage correct as a function of condition
    % figure
    % global_percentcorrect = mean(percentcorrect,3);
    % for i_masker = 1:N_maskers
    %     Me = squeeze(mean(global_percentcorrect(:,i_masker)));
    %     Sem = squeeze(std(global_percentcorrect(:,i_masker)))/sqrt(size(medianSNR,1));
    %     errorbar(i_masker,Me,Sem,'o','Color',maskercolors{i_masker},'MarkerFaceColor',maskercolors{i_masker});
    %     hold on
    %     for i_subject = 1:N_subjects
    %         plot(i_masker-0.1,squeeze(global_percentcorrect(i_subject,i_masker,:)),'.','Color',[maskercolors{i_masker},0.1]);
    %     end
    % end
    % xlabel('Masker condition')
    % ylabel('Correct response rate (%)')
    % title('group analysis of Correct response rate')
end

if flags.do_fig4 || flags.do_fig5_simu
    % Fig. 1: Group plot of SNR thresholds as a function of trial number
    figure('Position',[100 100 700 420]); % (hrev)
    % x_var = 400:400:4000;%resume_trial(2:end);
    x_var = 1:length(resume_trial)-1;
    x_var_label = resume_trial(1:end-1);

    hplt = [];
    for i_masker = 1:N_maskers
        medianSNR_here = squeeze(medianSNR(:,i_masker,:));
        Me  = mean(medianSNR_here);
        Sem = std(medianSNR_here)/sqrt(N_subjects);

        offx = 0.15*(i_masker-2);
        hplt(end+1) = errorbar(x_var+offx,Me,Sem,'-','Marker',Markers{i_masker}, ...
            'Color',masker_colours{i_masker},'MarkerFaceColor',masker_colours{i_masker},'LineWidth',LW(i_masker));
        hold on; grid on
        if i_masker == 1
            xlim([0 max(x_var)+1])
        end
    end
    for i_subject = 1:N_subjects
        meanSNR_participant(i_subject,:) = mean(squeeze(medianSNR(i_subject,:,:)));
        meanSNR_overall(i_subject) = mean(meanSNR_participant(i_subject,:),2);
    end
    [ma,ma_idx] = min(meanSNR_overall); % minimum is best performance
    [mi,mi_idx] = max(meanSNR_overall); % maximum is worst performance
    for i_subject = 1:N_subjects
        if i_subject == mi_idx || i_subject == ma_idx
            LW_here = 2;
            LS = '-'; % linestyle
        else
            LW_here = 1;
            LS = '--';
        end
        colour_grey = rgb('Gray');
        plot(x_var,meanSNR_participant(i_subject,:),LS,'Color',colour_grey,'LineWidth',LW_here);

        if i_subject == mi_idx
            text(0.85,0.75,Subjects{mi_idx},'Units','Normalized','FontWeight','Bold','Color',colour_grey);
        end
        if i_subject == ma_idx
            text(0.85,0.25,Subjects{ma_idx},'Units','Normalized','FontWeight','Bold','Color',colour_grey);
        end
    end

    set(gca,'XTick',x_var);
    set(gca,'XTickLabel',x_var_label);
    text(1700,-10,'S10','Color',[0.75,0.75,0.75])
    xlabel('Starting trial of the block #'); % xlim([0 4100])
    ylabel('SNR threshold (dB)')
    legend(hplt, noise_types_label)
    if bLeo
        title('Leo')
    end
    
    h(end+1) = gcf;
    hname{end+1} = 'fig4-SNRs';
    
    %%%
    % Run mixed ANOVA on medianSNR
    % % repeated measure anova
    % tdata = array2table(medianSNR(:,:));
    % [F_maskers,F_sessions] = meshgrid(1:N_maskers,1:N_sessions);
    % F_maskers = Maskers(F_maskers); F_maskers = F_maskers(:);
    % F_sessions = [cellfun(@num2str,num2cell(F_sessions),'UniformOutput',false)]; F_sessions = F_sessions(:);
    % factorNames = {'Masker','Session'};
    % within = table(F_maskers, F_sessions, 'VariableNames', factorNames);
    % % fit the repeated measures model
    % rm = fitrm(tdata,'Var1-Var30~1','WithinDesign',within);
    % [ranovatblb] = ranova(rm, 'WithinModel','Masker*Session');

    % Run differential mixed ANOVA on medianSNR
    [F_subjects,F_maskers,F_sessions] = ndgrid(1:N_subjects,1:N_maskers,1:N_sessions);
    [p,tbl,stats] = anovan(medianSNR(:),{F_maskers(:),F_sessions(:),F_subjects(:)},'varnames',{'maskers','sessions','subjects'},'random',3,'continuous',2,'display','off','model','linear');%
    fprintf(['\n*Results of mixed ANOVA on SNR thresholds with factors Masker and session (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})

    %%% Leo's analysis:
    % *Results of mixed ANOVA on SNR thresholds with factors Masker and session (and random factor subject)*
    % maskers: F(2,345) = 17.67, p = 0.000
    % sessions: F(1,345) = 35.50, p = 0.000

    %%% Alejandro's analysis (no data exclusion yet)
    % maskers: F(2,345) = 16.16, p = 0.000 % slightly smaller effect of masker
    % sessions: F(1,345) = 37.44, p = 0.000 % slightly larger effect of session

    %%% Latest results on 9/09/2022, after SNR exclusion:
    % *Results of mixed ANOVA on SNR thresholds with factors Masker and session (and random factor subject)*
    % maskers: F(2,345) = 15.87, p = 0.000
    % sessions: F(1,345) = 36.63, p = 0.000
    
    multcompare(stats,'display','off') % post-hoc analysis
    % Factors being compared (1-2; 1-3; 2-3)
    % 1.0000    2.0000   -0.7001   -0.4834   -0.2668    0.0000
    % 1.0000    3.0000   -0.2907   -0.0740    0.1426    0.7025
    % 2.0000    3.0000    0.1928    0.4094    0.6260    0.0000
    
    %stats.varnames {'maskers' } {'sessions'} {'subjects'}
    %
    % The fourth column shows the difference between the estimated group means. 
    % The third and fifth columns show the lower and upper limits for 95% CIs for the true mean difference. 
    % The sixth column contains the p-value for a hypothesis test that the corresponding mean difference is equal to zero. 
    % If All p-values are very small, indicates that the the dependent variable differs across all three factors.

    do_fig4_stats_N10 = input('Show analysis for 10 participants? (1=yes; 0=no): ');
    if do_fig4_stats_N10
        % Run differential mixed ANOVA on medianSNR
        idx_N = 1:N_subjects;
        idx_N([6 12]) = [];
    
        [F_subjects,F_maskers,F_sessions] = ndgrid(idx_N,1:N_maskers,1:N_sessions);
        SNR_here = medianSNR(idx_N,:,:);
        [p,tbl,stats] = anovan(SNR_here(:),{F_maskers(:),F_sessions(:),F_subjects(:)},'varnames',{'maskers','sessions','subjects'},'random',3,'continuous',2,'display','off','model','linear');%
        fprintf(['\n*Results of mixed ANOVA on SNR thresholds with factors Masker and session (and random factor subject)*\n'])
        fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
        fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})
        
        % Results on 4/10/2022:
        % *Results of mixed ANOVA on SNR thresholds with factors Masker and session (and random factor subject)*
        % maskers: F(2,287) = 15.82, p = 0.000
        % sessions: F(1,287) = 33.06, p = 0.000
        
        multcompare(stats,'display','off') % post-hoc analysis
        % 1.0000    2.0000   -0.7812   -0.5384   -0.2956    0.0000
        % 1.0000    3.0000   -0.3192   -0.0764    0.1664    0.7411
        % 2.0000    3.0000    0.2192    0.4620    0.7048    0.0000
    end
end % end do_fig4

if flags.do_fig5 
    % Fig. 2: Performance as a function of SNR
    P_H_resp2 = P_H; % Script3 assumes that the target sound is option '2'
    resp2_label = [data_hist.H_label '-trials (H)'];
    P_H_resp1 = 1-P_FA; % 'Hits' when the participant response was '1'
    resp1_label = [data_hist.CR_label '-trials (CR)'];

    figure('Position',[100 100 700 600]);
    il_tiledlayout(2,2,'TileSpacing','tight');
    
    XL = [-18.5 -8.5];
    for i_masker = 1:N_maskers
        
        il_nexttile(1)
        
        x_var = bin_centres;
        Score = 100*squeeze(P_H_resp1(:,i_masker,:));
        Me  = mean(Score);
        Sem = std(Score)/sqrt(size(Score,1));
        
        offx = 0.15*(i_masker-1);
        errorbar(x_var+offx,Me,Sem,'o-','Color',masker_colours{i_masker},'MarkerFaceColor',masker_colours{i_masker});
        hold on; grid on
        
        if i_masker == N_maskers
            % title('PC');
            % xlabel('SNR (dB)');
            ylabel('Percentage correct');
            ylim(100*[0.35 1.05]); 
            xlim(XL)
            % legend(resp2_label,'Location','SouthEast')
            set(gca,'XTickLabel','');
        end
        txt_options = {'Unit','Normalized','FontWeight','Bold','FontSize',14};
        text(0.05,0.92,'A. /aba/ trials (H)',txt_options{:});
        
        il_nexttile(2)
        
        Score = 100*squeeze(P_H_resp2(:,i_masker,:));
        Me = mean(Score);
        Sem = std(Score)/sqrt(size(Score,1));
        offx = 0.15*(i_masker-1);
        errorbar(x_var+offx,Me,Sem,'o-','Color',masker_colours{i_masker},'MarkerFaceColor', 'w');%maskercolors{i_masker});
        hold on; grid on
        
        if i_masker == N_maskers
            % title('PC');
            % xlabel('SNR (dB)');
            ylabel('Percentage correct');
            ylim(100*[0.35 1.05]); 
            xlim(XL)
            % legend(resp1_label,'Location','SouthEast')
            set(gca,'XTickLabel','');
        end
        text(0.05,0.92,'B. /ada/ trials (CR)',txt_options{:});
    end
    
    il_nexttile(3)
        
    for i_masker = 1:N_maskers
        x_var = bin_centres;
        Me = squeeze(mean(dprime(:,i_masker,:)));
        Sem = squeeze(std(dprime(:,i_masker,:)))/sqrt(size(dprime,1));
        offx = 0.15*(i_masker-1);
        errorbar(x_var+offx,Me,Sem,'o-','Color',masker_colours{i_masker},'MarkerFaceColor',masker_colours{i_masker});
        hold on; grid on;
        if i_masker == N_maskers
            % title('d prime');
            xlabel('SNR (dB)');
            ylabel('d''');
            xlim(XL)
            legend(noise_types_label, 'interpreter', 'none','Location','SouthEast')
        end
        
        ylim([-0.1 2.3]);
        set(gca,'YTick',0:0.2:2.2);
        
        text(0.05,0.92,'C. d''',txt_options{:});
    end

    il_nexttile(4)
        
    for i_masker = 1:N_maskers
        x_var = bin_centres;
        Me = squeeze(mean(criterion(:,i_masker,:)));
        Sem = squeeze(std(criterion(:,i_masker,:)))/sqrt(size(criterion,1));
        offx = 0.15*(i_masker-1);
        errorbar(x_var+offx,Me,Sem,'o-','Color',masker_colours{i_masker},'MarkerFaceColor',masker_colours{i_masker});
        hold on; grid on;
        %plot(bin_centres, criterion ,'-o','Color',colour_here); hold on, grid on;
        if i_masker == N_maskers
            xlim(XL)
            % title('criterion');
            xlabel('SNR (dB)');
            ylabel('c');
            
            ylim([-0.17 0.23]);
            set(gca,'YTick',-0.15:0.05:0.20);
        
            text(0.05,0.92,'D. Criterion c',txt_options{:});
        end
    end
    
    h(end+1) = gcf;
    hname{end+1} = 'fig5-Rates'; 
end

if flags.do_fig5_stats % Originally by Leo
    % Run mixed models on dprime & criterion
    
    % First we need to select a subset of SNRs where all 3 conditions are defined
    idx_SNR2analyse = 5:9; SNR2analyze = bin_centres(idx_SNR2analyse);
    dprime2analyse = dprime(:,:,idx_SNR2analyse);
    criterion2analyse = criterion(:,:,idx_SNR2analyse);
   
    [F_subjects,F_maskers,F_SNR] = ndgrid(1:N_subjects,1:N_maskers,SNR2analyze);
    [p,tbl,stats] = anovan(dprime2analyse(:),{F_maskers(:),F_SNR(:),F_subjects(:)},'varnames',{'maskers','SNR','subjects'},'random',3,'continuous',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on dprime with factors Masker and SNR (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})
    % multcompare(stats,'display','off')
    
    [p,tbl,stats] = anovan(criterion2analyse(:),{F_maskers(:),F_SNR(:),F_subjects(:)},'varnames',{'maskers','SNR','subjects'},'random',3,'continuous',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on criterion with factors Masker and SNR (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})
    %multcompare(stsats,'display','off')
    % Fig. 3: Performance as a function of SNR
    
    %%% Results on 9/09/2022:
    % *Results of Mixed ANOVA on dprime with factors Masker and SNR (and random factor subject)*
    % maskers: F(2,165) = 10.39, p = 0.000
    % SNR: F(1,165) = 1017.65, p = 0.000
    % 
    % *Results of Mixed ANOVA on criterion with factors Masker and SNR (and random factor subject)*
    % maskers: F(2,165) = 1.75, p = 0.178
    % SNR: F(1,165) = 9.77, p = 0.002   
    
    do_fig5_stats_N10 = input('Show analysis for 10 participants? (1=yes; 0=no): ');
end

if do_fig5_stats_N10 % requires that flags.do_fig5 = 1
    % Repeat the mixed ANOVAs on dprime & criterion, if requested
    
    idx_N = 1:N_subjects;
    idx_N([6 12]) = [];
    
    % First we need to select a subset of SNRs where all 3 conditions are defined
    dprime2analyse = dprime(idx_N,:,idx_SNR2analyse);
    criterion2analyse = criterion(idx_N,:,idx_SNR2analyse);
   
    fprintf('\tAnalysis using N=10\n');
    [F_subjects,F_maskers,F_SNR] = ndgrid(idx_N,1:N_maskers,SNR2analyze);
    [p,tbl,stats] = anovan(dprime2analyse(:),{F_maskers(:),F_SNR(:),F_subjects(:)},'varnames',{'maskers','SNR','subjects'},'random',3,'continuous',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on dprime with factors Masker and SNR (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})
    % multcompare(stats,'display','off')
    
    [p,tbl,stats] = anovan(criterion2analyse(:),{F_maskers(:),F_SNR(:),F_subjects(:)},'varnames',{'maskers','SNR','subjects'},'random',3,'continuous',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on criterion with factors Masker and SNR (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{5,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{5,3}, tbl{3,6}, tbl{3,7})
    %multcompare(stsats,'display','off')
    % Fig. 3: Performance as a function of SNR
    
    %%% Results on 4/10/2022:
    % *Results of Mixed ANOVA on dprime with factors Masker and SNR (and random factor subject)*
    % maskers: F(2,137) = 10.24, p = 0.000
    % SNR: F(1,137) = 788.09, p = 0.000
    % 
    % *Results of Mixed ANOVA on criterion with factors Masker and SNR (and random factor subject)*
    % maskers: F(2,137) = 1.41, p = 0.249
    % SNR: F(1,137) = 13.02, p = 0.000   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 6
% if flags.do_fig10 || flags.do_fig8b || flags.do_fig11_old || flags.do_fig12
%     % All simulation results were stored in one output folder
%     if isunix
%         dir_results = '/home/alejandro/Documents/Databases/data/fastACI_data_z_tmp/20220715_sim_Q1_osses2022a/';
%     else
%         error('Leo put your folder here...')
%         dir_results = [];
%     end
%     Subjects_ID = {'Results-S01-v1','Results-S02-v1','Results-S03-v1','Results-S04-v1', ...
%                    'Results-S05-v1','Results-S06-v1','Results-S07-v1','Results-S08-v1', ...
%                    'Results-S09-v1','Results-S10-v1','Results-S11-v1','Results-S12-v1'};
%     model_str = 'osses2022a';
% end

if flags.do_fig6 || flags.do_fig7 || flags.do_fig8 || flags.do_fig9 || ...
   flags.do_fig3_suppl || flags.do_fig8b_old
    glmfct = 'l1glm';
    bIs_Alejandro = isunix;
    bIs_Leo = ~bIs_Alejandro;
    
    flags_base = {'gammatone', ... % gt
                  glmfct, ... % l1glm
                  'no_bias', ... % nob
                  'expvar_after_reversal',4, ... %'idx_trialselect', 1:4000 ... % tr4000
                  'pyramid_script','imgaussfilt', ...
                  'pyramid_shape',-1 ...
                  }; % the extra fields

    if bIs_Alejandro
        % dir_where = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc/';
        dir_where = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc_Leo/';
            % Run-S01/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
            % Run-S02/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
            % ...
            % Run-S12/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat

        % Location of the savegames:
                   % /home/alejandro/Desktop/fastACI_today/fastACI_data/speechACI_Logatome-abda-S43M/S01/Results
        % dir_res = {'/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Publications/publ_osses2022b/',['1-experimental_results' filesep]};
        % dir_wav = '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/';
    end
    % for ii = 1:N_subjects
    %     fname_crosspred{ii} = [dir_where 'ACI-' subject '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
    %     if exist(fname_crosspred{ii},'file')
    %         bInclude(ii) = 1;
    %     end
    % end
end

if flags.do_fig6
     
    N_noises = length(noise_str);
     
    Pos4 = [800 800];
    N_subjs2plot = [6 6];
     
    for i_repeat = 1:2
        N_subj2plot = N_subjs2plot(i_repeat);
        figure('Position',[100 100 600 Pos4(i_repeat)]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
        il_tiledlayout(N_subj2plot,N_noises,'TileSpacing','tight');
        for i_subj = 1:N_subj2plot    
            dir_subj = dir_where;
            for i_noise = 1:N_noises
                ACI_fname = [dir_subj 'ACI-' Subjects{i_subj} '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4.mat'];
 
                ACI = [];
                cfg_ACI = [];
                results = [];
                info_toolbox = [];
                load(ACI_fname);
 
                bBox = 0;
                switch Subjects{i_subj}
                    case 'S10'
                        if i_noise == 1
                            bBox = 1;
                        end
                    case 'S11'
                        if i_noise == 1
                            bBox = 1;
                        end
                    otherwise
                        % Nothing to do
                end
                
                if i_noise == 3
                    bColourbar = 'yes';
                else
                    bColourbar = 'no';
                end
                flags_opts = {'NfrequencyTicks',5,'colorbar',bColourbar}; 
 
                il_nexttile(N_noises*(i_subj-1)+i_noise);
                affichage_tf(ACI,'CInorm', cfg_ACI.t, cfg_ACI.f,flags_opts{:}); hold on
                if i_subj ~= N_subj2plot
                    xlabel('')
                    set(gca,'XTickLabel','');
                end
                if i_noise ~= 1
                    ylabel('');
                    set(gca,'YTickLabel','');
                end
                
                %%% Adding a box to the null ACIs
                if bBox
                    YL_box = get(gca,'YLim');
                    XL_box = get(gca,'XLim');
                    xcoor = mean(XL_box);
                    ycoor = mean(YL_box);
                    offx = 0.98*abs(diff(XL_box)/2);
                    offy = 0.98*abs(diff(YL_box)/2);
                    il_plot_square_XY(xcoor,ycoor,'m',offx,offy,'--',2);
                end
                
                if i_subj == 1
                    title(noise_types_label{i_noise});
                    switch i_noise
                        case 1
                            if i_repeat == 1
                                text2show = 'A.';
                            else
                                text2show = 'D.';
                            end
                        case 2
                            if i_repeat == 1
                                text2show = 'B.';
                            else
                                text2show = 'E.';
                            end
                        case 3
                            if i_repeat == 1
                                text2show = 'C.';
                            else
                                text2show = 'F.';
                            end
                    end
                    text(0,1.15,text2show,'FontWeight','Bold','Units','Normalized','FontSize',14);
                end
                if i_noise == 3
                    text2show = Subjects{i_subj};
                    text(.7,.90,text2show,'FontWeight','Bold','Units','Normalized','FontSize',12);
                end
                
                disp('')
            end
        end
        Subjects(1:N_subj2plot) = []; % remove the first six labels, so that it's easier to repeat this code
         
        h(end+1) = gcf;
        hname{end+1} = ['fig06-ACI-set' num2str(i_repeat)];
    end
    disp('')
   
end % end do_fig6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7 || flags.do_fig8 || flags.do_fig9 || flags.do_fig3_suppl || ...
   flags.do_fig8b_old
    % Loading cross predictions
    
    dir_res = {[fastACI_paths('dir_fastACI') 'Publications' filesep 'publ_osses2022b' filesep],['1-experimental_results' filesep]};
           
    trialtype_analysis = 'total';

    
    for i_masker = 1:N_maskers
        masker  = noise_str{i_masker};

        for i_subject = 1:N_subjects
         
            subject = Subjects{i_subject};
            
            %%%
            if flags.do_fig7
                if i_subject == 1 % Only done once, at the beginning
                    bInclude = zeros(1,N_subjects);
                    
                    cross_suffix = '';
                    fname_crosspred = [];
                    for ii = 1:N_subjects
                        fname_crosspred{ii} = [dir_where 'ACI-' subject '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                        if exist(fname_crosspred{ii},'file')
                            bInclude(ii) = 1;
                        end
                    end
                end
            end
            %%%
            if flags.do_fig8 || flags.do_fig9 || flags.do_fig3_suppl || ...
                    flags.do_fig3b_suppl || flags.do_fig8b_old
                % 2.1. Checking whether you previously stored already the cross-prediction data:
                if flags.do_fig8 || flags.do_fig3_suppl || flags.do_fig8b_old
                    if i_subject == 1 
                        bInclude = zeros(1,N_subjects);
                        cross_suffix = '';
                        fname_crosspred = [];
                        for ii = 1:N_subjects
                            fname_crosspred{ii} = [dir_where 'ACI-' Subjects{ii} '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                            if exist(fname_crosspred{ii},'file')
                                bInclude(ii) = 1;
                            end    
                        end
                    end
                else
                    bInclude = zeros(1,N_subjects);
                    
                    cross_suffix = 'noise';
                    for ii = 1:N_maskers
                        % fname_crosspred{ii} = [dir_subject 'Results' filesep 'Results_ACI' filesep 'ACI-' subject '-speechACI_Logatome-' noise_str{ii} '-nob-gt-l1glm-rev4.mat'];
                        fname_crosspred{ii} = [dir_where 'ACI-' subject '-speechACI_Logatome-' noise_str{ii} '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                        if ii == i_masker
                            if exist(fname_crosspred{ii},'file')
                                bInclude(i_subject) = 1;
                            end
                        end
                    end
                end
            end
            
            if bInclude(i_subject)
                % loading data
                dir_results = [dir_res{1} 'data_' subject filesep dir_res{2}];
                files = Get_filenames(dir_results,['savegame_*' masker '.mat']);
                fname_results = [dir_results files{1}];

                [cfg_game,data_passation] = Convert_ACI_data_type(fname_results);

                dir_out_here = dir_where; %[dir_where subject filesep];

                cfg_game = Check_cfg_crea_dirs(cfg_game); % probably gives an error if the folders aren't visible

                flags_for_input = {...
                    'trialtype_analysis', trialtype_analysis, ...
                    'N_folds', 10, ...
                    'dir_noise',cfg_game.dir_noise, ...
                    'dir_target',cfg_game.dir_target, ...
                    'add_signal',0, ...
                    'apply_SNR',0, ...
                    'skip_if_on_disk',1, ...%'force_dataload', ...
                    'no_permutation', ...
                    'no_bias', ...
                    'no_plot'...
                    'dir_out',dir_out_here ...
                    };

                %%% 1. Obtaining the ACI of the current subject:
                [col1,cfg_ACI_here,res,Data_matrix_here] = fastACI_getACI(fname_results, ...
                    flags_base{:},flags_for_input{:}); % ,'force_dataload');
                % crosspred = [];
                
                % 2.1. Checking whether you previously stored already the cross-prediction data:
                [crosspred, fname_cross] = Check_if_crosspred(cfg_ACI_here.fnameACI,fname_crosspred,cross_suffix);

                bForceComplete = 0;

                if isempty(crosspred) || bForceComplete
                    % 2.2. If there is no cross-prediction in your computer, then
                    %      it obtains it
                    flags = {'ACI_crosspred',fname_crosspred,'Data_matrix',Data_matrix_here}; % the extra fields
                    [col1,cfg_ACI_here,res,col4] = fastACI_getACI(fname_results, flags_base{:},flags{:},flags_for_input{:});
                else
                    % (Part of 2.1): places 'crosspred' as a field of res. (local results)
                    res.crosspred = crosspred;
                end

                %%% 3. It places the cross-prediction reults in a cell array:
                ACI{i_subject,i_masker}     = col1;
                cfg_ACI{i_subject,i_masker} = cfg_ACI_here; 
                results{i_subject,i_masker} = res;
                % idxlambda(i_subject) = res.idxlambda;

                [cross,outs_cross] = Read_crosspred(fname_cross,res.idxlambda);
                PC_matrix(i_subject,:,i_masker)  = outs_cross.PA_mean_re_chance;
                Dev_matrix(i_subject,:,i_masker) = outs_cross.Dev_mean;
                Dev_matrix_plus_SEM(i_subject,:,i_masker) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
                Dev_matrix_min_SEM(i_subject,:,i_masker)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
                
                if flags.do_fig7
                    res_here = results{i_subject,i_masker};
                    % PC and DEV for incorrect trial
                    [PC_per_fold_inc,Dev_per_fold_inc] = il_tbtpred_inc(cfg_ACI_here,res_here,data_passation);

                    %%% Save some data for group-level analysis
                    idxlambda = res_here.idxlambda;

                    num1 = res_here.FitInfo.Dev_test./res_here.FitInfo.CV.TestSize;
                    num2 = num1(end,:); % referential null ACI 
                    G_Devtrial_opt(:,i_subject,i_masker) = num1(idxlambda,:) - num2;

                    num1 = res_here.FitInfo.PC_test;
                    num2 = num1(end,:); % referential null ACI 
                    G_PAtrial_opt(:,i_subject,i_masker) = num1(idxlambda,:) - num2;

                    G_OptDEV(i_subject,i_masker) = mean(G_Devtrial_opt(:,i_subject,i_masker));
                    G_OptPA(i_subject,i_masker)  = mean(G_PAtrial_opt( :,i_subject,i_masker));
                    
                    %%% Incorrect:
                    num1 = Dev_per_fold_inc;
                    num2 = num1(end,:); % referential null ACI 
                    G_Devtrial_opt_inc(:,i_subject,i_masker) = num1(idxlambda,:) - num2;

                    num1 = PC_per_fold_inc;
                    num2 = num1(end,:); % referential null ACI 
                    G_PAtrial_opt_inc(:,i_subject,i_masker) = num1(idxlambda,:) - num2;

                    G_OptDEV_inc(i_subject,i_masker) = mean(G_Devtrial_opt_inc(:,i_subject,i_masker));
                    G_OptPA_inc(i_subject,i_masker)  = mean(G_PAtrial_opt_inc(:,i_subject,i_masker));
                end % if flags.do_fig7
                
            end % end if bInclude
        end % end if i_subject
    end % end i_masker
    
end

if flags.do_fig7
    % Taken from l20220325_crossprediction_newversion.m and 
    %   g20220203_crossprediction.m
    
    %% function l20220502_analysis_pipeline_ACIstats
    figure('Position',[100 100 700*2-200 550]);
    il_tiledlayout(2,2,'TileSpacing','tight');
    
    crit_sig = [0 0 1.3 2.39]; % significance criterion
    
    XL = [0.5 3.5]; 
    il_nexttile(1);
    plot(XL,crit_sig(1)*[1 1],'k--'); hold on; grid on
    
    il_nexttile(2);
    plot(XL,crit_sig(2)*[1 1],'k--'); hold on; grid on
    
    il_nexttile(3);
    plot(XL,crit_sig(3)*[1 1],'k--'); hold on; grid on
    
    il_nexttile(4);
    plot(XL,crit_sig(4)*[1 1],'k--'); hold on; grid on
    
    for i_masker = 1:N_maskers
        for i_subject = 1:N_subjects

            if i_subject == 9 || i_subject == 10
                disp('')
            end
            offx = 0.05*(i_subject - ((N_subjects+2)/2-.5)); %  0.05*i_subject;
            % % Find optimal lambda
            % [Devtrial_opt,lambda_opt] = min(mean(G_Devtrial(:,:,i_subject,i_masker),2),[],1);

            % y2plot  = G_Devtrial_opt(lambda_opt,:,i_subject,i_masker);
            % y2plot2 = 100*G_PA_opt(lambda_opt,:,i_subject,i_masker);
            y2plot  = G_Devtrial_opt(:,i_subject,i_masker);
            y2plot2 = 100*G_PAtrial_opt(:,i_subject,i_masker);
            
            il_nexttile(1);
            Me = mean(y2plot);
            EB = 1.64*sem(y2plot);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(1) && min(Me_SEM) <= crit_sig(1) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx, Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
            
            il_nexttile(3);
            Me = mean(y2plot2);
            EB = 1.64*sem(y2plot2);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(3) && min(Me_SEM) <= crit_sig(3) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            elseif diff(Me_SEM) == 0
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx, Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
            
            %%% Incorrect trials:
            y2plot  = G_Devtrial_opt_inc(:,i_subject,i_masker);
            y2plot2 = 100*G_PAtrial_opt_inc(:,i_subject,i_masker);

            il_nexttile(2);
            Me = mean(y2plot);
            EB = 1.64*sem(y2plot);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(2) && min(Me_SEM) <= crit_sig(2) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx, Me , EB , 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

            il_nexttile(4);
            Me = mean(y2plot2);
            EB = 1.64*sem(y2plot2);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(4) && min(Me_SEM) <= crit_sig(4) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            elseif diff(Me_SEM) == 0
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx,Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on            
        end

        y2plot  = G_OptDEV(:,i_masker);
        % subplot(2,1,1)
        il_nexttile(1);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        y2plot  = 100*G_OptPA(:,i_masker);
        % subplot(2,1,2)
        il_nexttile(3);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        %%% Incorrect:
        y2plot  = G_OptDEV_inc(:,i_masker);
        il_nexttile(2);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        y2plot  = 100*G_OptPA_inc(:,i_masker);
        il_nexttile(4);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
    end

    XLabel = 'Masker';
    % subplot(2,1,1)
    il_nexttile(1);
    title('All trials')
    grid on
    xlim(XL);
    % xlabel(XLabel); 
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', []); 
    ylabel('Deviance / trial benefit (adim)')
    xcoor = 0.01;
    ycoor = 0.94;
    txt_opts = {'Units','Normalized','FontWeight','Bold','FontSize',14};
    text(xcoor,ycoor,'A',txt_opts{:});
    
    il_nexttile(2);
    title('Incorrect trials only')
    grid on
    xlim(XL);
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', []); 
    set(gca, 'YTickLabels',[]);
    text(xcoor,ycoor,'B',txt_opts{:});
    
    il_nexttile(3);
    grid on
    xlim(XL)
    ylim([-1 23])
    xlabel(XLabel); 
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', noise_types_label); 
    ylabel('Percent accuracy benefit (%)')
    text(xcoor,ycoor,'C',txt_opts{:});
    
    il_nexttile(4);
    
    xlim(XL)
    ylim([-1 23])
    xlabel(XLabel); 
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', noise_types_label); 
    set(gca, 'YTickLabels',[]);
    % ylabel('Percent accuracy benefit (%)')
    text(xcoor,ycoor,'D',txt_opts{:});
    
    h(end+1) = gcf;
	hname{end+1} = ['fig7-metrics-benefit'];
        
    disp('')
    
    % Test on Optimal DEV
    [F_subjects,F_maskers] = ndgrid(1:N_subjects,1:N_maskers);
    [p,tbl,stats] = anovan(G_OptDEV(:),{F_maskers(:),F_subjects(:)},'varnames',{'maskers','subjects'},'random',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on optimal DEV with factors Masker (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{4,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{4,3}, tbl{3,6}, tbl{3,7})
    %multcompare(stats,'display','off')
    [p,tbl,stats] = anovan(G_OptPA(:),{F_maskers(:),F_subjects(:)},'varnames',{'maskers','subjects'},'random',2,'display','off');
    fprintf(['\n*Results of Mixed ANOVA on optimal PA with factors Masker (and random factor subject)*\n'])
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{2,1}, tbl{2,3}, tbl{4,3}, tbl{2,6}, tbl{2,7})
    fprintf('%s: F(%d,%d) = %.2f, p = %.3f\n', tbl{3,1}, tbl{3,3}, tbl{4,3}, tbl{3,6}, tbl{3,7})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flags.do_fig8 || flags.do_fig3_suppl || flags.do_fig8b_old
    % Code by Leo, taken from l20220325_crossprediction_newversion.m
 
    XTL = [];
    XT = 1:N_subjects;
    for ii = 1:N_subjects
        txt2use = 'S';
        if ii < 10 
            txt2use = [txt2use '0'];
        end
        txt2use = [txt2use num2str(ii)];
        YTL{ii} = txt2use;
        if mod(ii,2)==0
            XTL{ii} = txt2use;
        else
            XTL{ii} = '';
        end
    end
    
    % figure('Position',[100 100 800 350]);
    figure('Position',[100 100 720 415]);
    il_tiledlayout(1,3,'TileSpacing','tight'); % 'tight');

    idx_sig = [0 0 0]; % idx of significant ACIs
    for i_column = 1:3
        switch mod(i_column,3)
            case 1 % white
                Colour_here = masker_colours{1};
            case 2 % bump
                Colour_here = masker_colours{2};
            case 0 % MPS
                Colour_here = masker_colours{3};
        end
        Dev_matrix_here = Dev_matrix(:,:,i_column);
        Dev_matrix_plus = Dev_matrix_plus_SEM(:,:,i_column);
        % Dev_matrix_min  = Dev_matrix_min_SEM(:,:,i_column);
        % Dev_diag(:,i_column) = diag(Dev_matrix_here);
        Dev_diag(:,i_column) = diag(Dev_matrix_plus);
        Dev_matrix_nodiag = Dev_matrix_here-repmat(Dev_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        if flags.do_fig3_suppl
            il_nexttile(i_column)
            % imagesc(1:N_subjects,1:N_subjects,Dev_matrix(:,:,i_column)); hold on
            imagesc(1:N_subjects,1:N_subjects,Dev_matrix_plus); hold on
            % caxis([-10 0])
            colormap('gray')
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            if i_column ~= 1
                set(gca,'YTickLabel',[]);
            end
            if i_column == 3
                c = colorbar;
            end
            %%% Plotting significant (or non-significant)
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0
                    % idx = find(Dev_matrix_nodiag(:,ii)<0); % criterion based on the diagonal
                    idx = find(Dev_matrix_plus(:,ii)<0); % criterion based on Dev + 1.64 SEM < 0
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                            if ii ~= jj % only counting cross predictions
                                idx_sig(i_column) = idx_sig(i_column)+1;
                            end
                        end
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                il_plot_square(ii,ii,Colour_here,offxy); % square of the diagonal
            end
            %%%
        end
        PC_matrix_here = PC_matrix(:,:,i_column);
        PC_diag(:,i_column) = diag(PC_matrix_here);
        PC_matrix_nodiag = PC_matrix_here-repmat(PC_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8
            il_nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([0 10])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            if i_column ~= 1
                set(gca,'YTickLabel',[]);
            end
            if i_column == 3
                c = colorbar;
            end
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus(:,ii)<0);
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                            if ii ~= jj % only counting cross predictions
                                idx_sig(i_column) = idx_sig(i_column)+1;
                            end
                        end
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
        end
    end
    fprintf('Delta PAexp values ranged between %.1f and %.1f\n', ...
        min(min(min(PC_matrix))),max(max(max(PC_matrix))));
    fprintf('ACIs leading to significant predictions using the CVD criterion: %.1f, %.1f, %.1f (white, bump, MPS)\n', ...
        idx_sig);
    FS_here = 10;

    nexttile(1);
    ylabel('Data from');
    set(gca,'YTickLabel',YTL);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'A.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    il_nexttile(2);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    text(0,1.05,'B.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    il_nexttile(3);
    text(0,1.05,'C.',opts_txt{:});
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    if flags.do_fig8
        % c.Label.String = 'Percentage accuracy benefit (%)';
        set(c,'Limits',[0 10.5]);
        fig_name = 'fig8a-crosspred_participant';
    end
    if flags.do_fig3_suppl
        % c.Label.String = 'Deviance / trial benefit (adim)';
        % Lc = get(c,'Limits');
        % set(c,'Limits',[Lc(1) 0])
        fig_name = 'suppl-fig3a-crosspred_participant';
    end
    set(gca,'FontSize',FS_here);

    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);

    for i_masker = 1:N_maskers
        nexttile(i_masker);
        ht = title(noise_types_label{i_masker});
        set(ht,'FontSize',FS_here);
        xlabel('ACI from');
    end
    h(end+1) = gcf;
    hname{end+1} = fig_name;

    %%%
    if flags.do_fig8b_old
        error('Uncomment below, removed on 27/09/2022')
        % PC_rem = [];
        % for idx_noise = 1:3
        %     PC_rem{idx_noise} = PC_matrix(:,:,idx_noise)-repmat(diag(PC_matrix(:,:,idx_noise))',N_subjects,1);
        %     idx = find(PC_rem{idx_noise}>0); 
        %     length(idx)
        % end
        % 
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %% Replicating the figure, but showing averages of cross predictions (between subjects)
        % figure('Position',[100 100 700*2-200 550]);
        % il_tiledlayout(2,2,'TileSpacing','tight');
        % 
        % XL = [0.5 3.5]; 
        % 
        % il_nexttile(1);
        % plot(XL,0*[1 1],'k--'); hold on; grid on
        % 
        % il_nexttile(2);
        % plot(XL,0*[1 1],'k--'); hold on; grid on
        % 
        % il_nexttile(3);
        % plot(XL,0*[1 1],'k--'); hold on; grid on
        % 
        % il_nexttile(4);
        % plot(XL,0*[1 1],'k--'); hold on; grid on
        % 
        % for i_masker = 1:N_maskers
        %     for i_subject = 1:N_subjects
        % 
        %         if i_subject == 9 || i_subject == 10
        %             disp('')
        %         end
        %         offx = 0.05*(i_subject - ((N_subjects+2)/2-.5)); %  0.05*i_subject;
        % 
        %         idx_here = 1:N_subjects;
        %         idx_here(i_subject) = []; % removing the diagonal
        %         y2plot  = nan([N_subjects 1]); 
        %         y2plot2 = PC_rem{i_masker}(idx_here,i_subject);
        % 
        %         MCol = 'w';
        %         il_nexttile(1);
        %         Me = nan; % y2plot;
        %         % EB = 1.64*sem(y2plot);
        %         % MCol = 'w';
        %         % Me_SEM = [Me+EB Me-EB];
        %         plot(i_masker+offx, Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
        % 
        %         il_nexttile(3);
        %         Me = median(y2plot2);
        %         % EB = 1.64*sem(y2plot2);
        %         errU = max(y2plot2)-Me;
        %         errL = Me-min(y2plot2);
        %         % MCol = 'w';
        %         % Me_SEM = [Me+EB Me-EB];
        %         errorbar(i_masker+offx, Me, errL, errU, 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
        % 
        %         %%% Incorrect trials:
        %         y2plot  = nan([N_subjects 1]); 
        %         y2plot2 = nan([N_subjects 1]); 
        % 
        %         il_nexttile(2);
        %         Me = mean(y2plot);
        %         % EB = 1.64*sem(y2plot);
        %         % MCol = 'w';
        %         % Me_SEM = [Me+EB Me-EB];
        %         plot(i_masker+offx, Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
        % 
        %         il_nexttile(4);
        %         Me = mean(y2plot2);
        %         % EB = 1.64*sem(y2plot2);
        %         % MCol = 'w';
        %         % Me_SEM = [Me+EB Me-EB];
        %         plot(i_masker+offx,Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on            
        %     end
        % end
        % 
        % XLabel = 'Masker';
        % % subplot(2,1,1)
        % il_nexttile(1);
        % title('All trials')
        % grid on
        % xlim(XL);
        % % xlabel(XLabel); 
        % set(gca, 'XTick', 1:N_maskers); 
        % set(gca, 'XTickLabels', []); 
        % ylabel('Deviance / trial benefit (adim)')
        % xcoor = 0.01;
        % ycoor = 0.94;
        % txt_opts = {'Units','Normalized','FontWeight','Bold','FontSize',14};
        % text(xcoor,ycoor,'NaN',txt_opts{:});
        % set(gca,'FontSize',FS_here);
        % 
        % il_nexttile(2);
        % title('Incorrect trials only')
        % grid on
        % xlim(XL);
        % set(gca, 'XTick', 1:N_maskers); 
        % set(gca, 'XTickLabels', []); 
        % set(gca, 'YTickLabels',[]);
        % text(xcoor,ycoor,'NaN',txt_opts{:});
        % set(gca,'FontSize',FS_here);
        % 
        % il_nexttile(3);
        % grid on
        % xlim(XL)
        % % ylim([-1 6])
        % xlabel(XLabel); 
        % set(gca, 'XTick', 1:N_maskers); 
        % set(gca, 'XTickLabels', noise_types_label); 
        % ylabel('Percent accuracy benefit (%)')
        % text(xcoor,ycoor,'D',txt_opts{:});
        % set(gca,'FontSize',FS_here);
        % 
        % il_nexttile(4);
        % 
        % xlim(XL)
        % % ylim([-1 6])
        % xlabel(XLabel); 
        % set(gca, 'XTick', 1:N_maskers); 
        % set(gca, 'XTickLabels', noise_types_label); 
        % set(gca, 'YTickLabels',[]);
        % % ylabel('Percent accuracy benefit (%)')
        % text(xcoor,ycoor,'NaN',txt_opts{:});
        % set(gca,'FontSize',FS_here);
        % 
        % h(end+1) = gcf;
        % hname{end+1} = ['fig8b-metrics-benefit-cross'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig9
    for j = 1:N_subjects
        for k = 1:N_maskers
            for i = 1:N_maskers 
                idxlambda = results{j,i}.idxlambda;

                PC_test  = results{j,i}.crosspred(k).PC_test;
                yvar = PC_test(idxlambda,:)-PC_test(end,:);
                PC_matrix(j,k,i) = 100*mean(yvar); % mean converted to a percentage

                Dev_test = results{j,i}.crosspred(k).Dev_test; 
                yvar = Dev_test(idxlambda,:)-Dev_test(end,:);
                Dev_matrix(j,k,i) = mean(yvar); % mean 
                Dev_matrix_plus_SEM(j,k,i) = mean(yvar)+1.64*sem(yvar);
                Dev_matrix_min_SEM(j,k,i)  = mean(yvar)-1.64*sem(yvar);
            end
        end
    end
        
    % Build the exact matrix to plot
    for j = 1:N_subjects
        for k = 1:N_maskers
            PC_pool_here{j}(k,:) = transpose(squeeze(PC_matrix(j,k,:)));
            Dev_pool_here{j}(k,:) = transpose(squeeze(Dev_matrix(j,k,:)));
        end
    end
    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 1:4
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool = PC;
    PC_pool_diag_removed = PC_rem_tmp;
    Dev_pool = Dev;

    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 5:8
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool  = [PC_pool; PC];
    PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
    Dev_pool = [Dev_pool; Dev];

    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 9:12
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool = [PC_pool; PC];
    PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
    Dev_pool = [Dev_pool; Dev];
    
    % figure;
    figure('Position',[100 100 1000 420]); % 560 420]);
    il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');

    il_nexttile(1);
    x_var = 1:size(Dev_pool,2);
    y_var = 1:size(Dev_pool,1);
    imagesc(x_var,y_var,Dev_pool); hold on
    colormap('gray')
    caxis([-10 0])

    c = colorbar;
    c.Label.String = 'Deviance / trial benefit (adim)';

    for ii = 1:4
        text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 5:8
        text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 9:12
        text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    %%% Drawing a separator between participants:
    plot(3.5*[1 1],[0 9.5],'Color','b');
    plot(3.5*[1 1],[0 9.5],'Color','b');

    plot(6.5*[1 1],[0 9.5],'Color','b');
    plot(6.5*[1 1],[0 9.5],'Color','b');

    plot(9.5*[1 1],[0 9.5],'Color','b');
    plot(9.5*[1 1],[0 9.5],'Color','b');

    plot([0 12.5],3.5*[1 1],'Color','b');
    plot([0 12.5],6.5*[1 1],'Color','b');
    plot([0 12.5],9.5*[1 1],'Color','b');
    set(gca,'FontSize',10);
    %%%

    disp('')
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'D.','Units','Normalized','FontSize',14,'FontWeight','Bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % figure;
    il_nexttile(2);
    x_var = 1:size(PC_pool,2);
    y_var = 1:size(PC_pool,1);
    imagesc(x_var,y_var,PC_pool); hold on;
    colormap('gray')
    caxis([0 10])

    %%% Drawing a separator between participants:
    plot(3.5*[1 1],[0 9.5],'Color','b');
    plot(3.5*[1 1],[0 9.5],'Color','b');

    plot(6.5*[1 1],[0 9.5],'Color','b');
    plot(6.5*[1 1],[0 9.5],'Color','b');

    plot(9.5*[1 1],[0 9.5],'Color','b');
    plot(9.5*[1 1],[0 9.5],'Color','b');

    plot([0 12.5],3.5*[1 1],'Color','b');
    plot([0 12.5],6.5*[1 1],'Color','b');
    plot([0 12.5],9.5*[1 1],'Color','b');
    set(gca,'FontSize',10);
    %%%

    PC_pool_diag_removed = transpose(PC_pool_diag_removed);
    offxy = 0.5;
    offxy = offxy*.9;
    for jj = 1:size(PC_pool_diag_removed,1)
        idx_here = find(PC_pool_diag_removed(jj,:)>0);
        for ii = 1:length(idx_here)
            il_plot_square(jj,idx_here(ii),'m',offxy,'--');
        end
    end

    c = colorbar;
    c.Label.String = 'Percentage accuracy benefit (%)';

    for ii = 1:4
        text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 5:8
        text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 9:12
        text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
    end

    disp('')
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

    h(end+1) = gcf;
    hname{end+1} = 'fig9a-crosspred_masker';

    % c.Label.String = 'Deviance / trial benefit (adim)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure;
    figure('Position',[100 100 560 250]); % 560 420]);
    % il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');

    PC_all = (PC_pool(1:3,:)+PC_pool(4:6,:)+PC_pool(7:9,:))/3;
    PC_all = (PC_all(:,1:3)+PC_all(:,4:6)+PC_all(:,7:9)+PC_all(:,10:12))/4;

    x_var = 1:size(PC_all,2);
    y_var = 1:size(PC_all,1);
    imagesc(x_var,y_var,PC_all); hold on
    colormap('gray')
    caxis([0 10])
    set(gca,'FontSize',10);

    c = colorbar;
    c.Label.String = 'Percentage accuracy benefit (%)';

    for ii = 1:3
        for jj = 1:3
            text(ii-0.1,jj,sprintf('%.1f',PC_all(ii,jj)),'FontSize',10,'Color','b');
        end
    end
    % for ii = 5:8
    %     text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
    % end
    % for ii = 9:12
    %     text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
    % end
    % %%% Drawing a separator between participants:
    % plot(3.5*[1 1],[0 9.5],'Color','b');
    % plot(3.5*[1 1],[0 9.5],'Color','b');
    % 
    % plot(6.5*[1 1],[0 9.5],'Color','b');
    % plot(6.5*[1 1],[0 9.5],'Color','b');
    % 
    % plot(9.5*[1 1],[0 9.5],'Color','b');
    % plot(9.5*[1 1],[0 9.5],'Color','b');
    % 
    % plot([0 12.5],3.5*[1 1],'Color','b');
    % plot([0 12.5],6.5*[1 1],'Color','b');
    % plot([0 12.5],9.5*[1 1],'Color','b');
    % %%%
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'B.','Units','Normalized','FontSize',14,'FontWeight','Bold');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel('ACI from condition');
    ylabel('Data from condition');

    % text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

    h(end+1) = gcf;
    hname{end+1} = 'fig9b-crosspred_masker';

end
     
% hcrossPC = figure; hcrossDev = figure;
% for i = 1:N_maskers%length(results{1}.crosspred)
%     for k = 1:N_subjects
%         %plot cross PC
%         xvar = results{k}.crosspred(i).lambdas;
%         PC_test = results{k}.crosspred(i).PC_test;
%         Dev_test = results{k}.crosspred(i).Dev_test;
%         figure(hcrossPC)
%         yvar = mean(PC_test-PC_test(end,:),2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
%         evar = std(PC_test-PC_test(end,:),[],2);%/size(PC_test,2);%SEM
%         subplot(1,N_subjects,k);
%         errorbar(xvar,yvar,evar,'-','Color',Colour{i}); hold on;
%         set(gca,'XScale','log');
%         title([cfg_ACI{k}.Condition ', subject ' Subjects{k}]);
%         ylabel(['(cross)prediction percent accuracy benefit ' name_crosspred]);
%         xlabel('\lambda');
%         ylim([0 0.1])
%         %plot cross Dev
%         figure(hcrossDev)
%         yvar = mean(Dev_test-Dev_test(end,:),2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
%         evar = std(Dev_test-Dev_test(end,:),[],2);%/size(PC_test,2);%SEM
%         subplot(1,N_subjects,k);
%         errorbar(xvar,yvar,evar,'-','Color',Colour{i}); hold on;
%         set(gca,'XScale','log');
%         title([cfg_ACI{k}.Condition ', subject ' Subjects{k}]);
%         ylabel(['(cross)prediction deviance benefit ' name_crosspred]);
%         xlabel('\lambda');
%     end
%     legend(Masker_crosspred,'interpreter','none');
%     
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figs 10 and 11
if flags.do_fig10      || flags.do_fig8b      || flags.do_fig11_old || flags.do_fig12 || ...
   flags.do_fig3_suppl || flags.do_fig3b_suppl
    
    % All simulation results were stored in one output folder
    if isunix
        dir_results = '/home/alejandro/Documents/Databases/data/fastACI_data_z_tmp/20220715_sim_Q1_osses2022a/';
    else
        error('Leo put your folder here...')
        dir_results = [];
    end
    Subjects_ID = {'Results-S01-v1','Results-S02-v1','Results-S03-v1','Results-S04-v1', ...
                   'Results-S05-v1','Results-S06-v1','Results-S07-v1','Results-S08-v1', ...
                   'Results-S09-v1','Results-S10-v1','Results-S11-v1','Results-S12-v1'};
    model_str = 'osses2022a';
end

if flags.do_fig10
    
    N_noises = length(noise_str);
    
    Pos4 = [400 1200];
    N_subjs2plot = [3 9];
    
    for i_repeat = 1:2
        
        N_subj2plot = N_subjs2plot(i_repeat);
        figure('Position',[100 100 600 Pos4(i_repeat)]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
        il_tiledlayout(N_subj2plot,N_noises,'TileSpacing','tight');
        for i_subj = 1:N_subj2plot    
            dir_subj = [dir_results Subjects_ID{i_subj} filesep];
            for i_noise = 1:N_noises
                ACI_fname = [dir_subj 'ACI-' model_str '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4.mat'];

                ACI = [];
                cfg_ACI = [];
                results = [];
                info_toolbox = [];
                load(ACI_fname);

                if i_noise == 3
                    bColourbar = 'yes';
                else
                    bColourbar = 'no';
                end
                flags_opts = {'NfrequencyTicks',5,'colorbar',bColourbar}; 

                il_nexttile(N_noises*(i_subj-1)+i_noise);
                affichage_tf(ACI,'CInorm', cfg_ACI.t, cfg_ACI.f,flags_opts{:});
                if i_subj ~= N_subj2plot
                    xlabel('')
                    set(gca,'XTickLabel','');
                end
                if i_noise ~= 1
                    ylabel('');
                    set(gca,'YTickLabel','');
                end
                if i_subj == 1
                    title(noise_types_label{i_noise});
                    switch i_noise
                        case 1
                            text2show = 'A.';
                        case 2
                            text2show = 'B.';
                        case 3
                            text2show = 'C.';
                    end
                    text(0,1.15,text2show,'FontWeight','Bold','Units','Normalized','FontSize',14);
                end
                if i_noise == 3
                    Subj_ID = strsplit(Subjects_ID{i_subj},'-');
                    text2show = Subj_ID{2};
                    text(.7,.90,text2show,'FontWeight','Bold','Units','Normalized','FontSize',12);
                end
                
                disp('')
            end
        end
        Subjects_ID(1:N_subj2plot) = []; % remove the first six labels, so that it's easier to repeat this code
        
        h(end+1) = gcf;
        hname{end+1} = ['fig10-ACI-sim-set' num2str(i_repeat)];
    end
    disp('')
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig11_old || flags.do_fig8b || flags.do_fig12 || flags.do_fig3_suppl || ...
   flags.do_fig3b_suppl
    N_noises = length(noise_str);
    N_subjects = length(Subjects_ID);
    
    if flags.do_fig11_old
        % N_subj2plot = 6;
        figure('Position',[100 100 600 1200]);
        il_tiledlayout(N_noises,1,'TileSpacing','tight');
    end
    for i_noise = 1:N_noises
        PA_mean_chance = []; % emptied after having plotted every noise
        for i_subj = 1:N_subjects
            if i_noise == 1
                Subj_ID = strsplit(Subjects_ID{i_subj},'-');
                labels2use{i_subj} = Subj_ID{2};
            end
        
            dir_subj = [dir_results Subjects_ID{i_subj} filesep];
        
            ACI_fname_folder = [dir_subj 'ACI-' model_str '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4' filesep];
            Crossfile = [ACI_fname_folder 'Crosspred.mat'];
            
            cfg_ACI_here = load([ACI_fname_folder(1:end-1) '.mat']);
            res_here = cfg_ACI_here.results;
            %%% Save some data for group-level analysis
            idxlambda = res_here.idxlambda;
            
            % [cross,outs_cross] = Read_crosspred(Crossfile, idxlambda);
            % PA_mean_chance(i_subj,:) = outs_cross.PA_mean_re_chance;
            [cross,outs_cross] = Read_crosspred(Crossfile, idxlambda);
            PC_matrix(i_subj,:,i_noise)  = outs_cross.PA_mean_re_chance;
            Dev_matrix(i_subj,:,i_noise) = outs_cross.Dev_mean;
            Dev_matrix_plus_SEM(i_subj,:,i_noise) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
            Dev_matrix_min_SEM(i_subj,:,i_noise)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
                
            % num1 = res_here.FitInfo.Dev_test./res_here.FitInfo.CV.TestSize;
            % num2 = num1(end,:); % referential null ACI 
            % G_Devtrial_opt(:,i_subject,i_masker) = num1(idxlambda,:) - num2;
             
            num1 = res_here.FitInfo.PC_test;
            num2 = num1(end,:); % referential null ACI 
            G_PA_sim(:,i_subj,i_noise) = 100*(num1(idxlambda,:) - num2);
            
            %%%
            res = [];
            res.idxlambda = idxlambda;
            res.crosspred = cross;
            results{i_subj,i_noise} = res;
            %%%
            
            if flags.do_fig12
                % i_subj: ACI from
                % k: Data from
                for k = 1:N_subjects
                    PC_test  = cross(k).PC_test;
                    yvar = PC_test(idxlambda,:)-PC_test(end,:);
                    PC_matrix(i_subj,k,i_noise) = 100*mean(yvar);
                end
            end
        end
        
        if flags.do_fig11_old
            il_nexttile(i_noise)
            opts = [];
            opts.bColourbar = 1;
            opts.cmin = 0;
            opts.cmax = 25;
            opts.Labels = labels2use;
            my_image_plot(PC_matrix(:,:,i_noise),opts);

            if i_noise == 3
                xlabel('ACI_s_i_m from')
            else
                set(gca,'XTickLabel',[]);
            end
            ylabel('Waveforms from')
            % Pos = get(gcf,'Position');
            % Pos(3:4) = [700 550];
            % set(gcf,'Position',Pos);

            switch i_noise
                case 1
                    data.PA_mean_chance_white = PC_matrix(:,:,i_noise);
                    text2use = 'A. white';
                case 2
                    data.PA_mean_chance_bump = PC_matrix(:,:,i_noise);
                    text2use = 'B. bump';
                case 3
                    data.PA_mean_chance_MPS = PC_matrix(:,:,i_noise);
                    text2use = 'C. MPS';
            end
            title(text2use);
        end
    end
    
    if flags.do_fig11_old
        h(end+1) = gcf;
        hname{end+1} = 'fig11-crosspred';
        disp('')
    end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b || flags.do_fig3b_suppl
    %%% 4. Plotting the results:
	% 3 matrices of cross-prediction accuracy
    % xvar = results{1,1}.crosspred(1).lambdas;
    XTL = [];
    XT = 1:N_subjects;
    for ii = 1:N_subjects
        txt2use = 'S';
        if ii < 10 
            txt2use = [txt2use '0'];
        end
        txt2use = [txt2use num2str(ii)];
        YTL{ii} = txt2use;
        if mod(ii,2)==0
            XTL{ii} = txt2use;
        else
            XTL{ii} = '';
        end
    end

    figure('Position',[100 100 720 415]);
    il_tiledlayout(1,3,'TileSpacing','tight');

    for i_column = 1:3
        switch mod(i_column,3)
            case 1 % white
                Colour_here = masker_colours{1};
            case 2 % bump
                Colour_here = masker_colours{2};
            case 0 % MPS
                Colour_here = masker_colours{3};
        end
        Dev_matrix_here = Dev_matrix_plus_SEM(:,:,i_column);
        Dev_diag(:,i_column) = diag(Dev_matrix_here);
%             Dev_matrix_nodiag = Dev_matrix_here-repmat(Dev_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        il_nexttile(i_column)
        if flags.do_fig3b_suppl
            imagesc(1:N_subjects,1:N_subjects,Dev_matrix_here); hold on
%             imagesc(1:N_subjects,1:N_subjects,Dev_matrix(:,:,i_column)); hold on
            colormap('gray')
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            idx_total = 0;
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus_SEM(:,ii,i_column)<0);
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                        end
                        idx_total = idx_total + length(idx);
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
            fprintf('CVD - Noise cond=%s: idx=%.0f of %.0f with significant predictions\n',noise_types_label{i_noise},idx_total,N_subjects*N_subjects);
        end
        
        PC_matrix_here = PC_matrix(:,:,i_column);
        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8b
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([14 24.5])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            idx_total = 0;
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus_SEM(:,ii,i_column)<0);
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                        end
                        idx_total = idx_total + length(idx);
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
        end

        fprintf('PA - Noise cond=%s: idx=%.0f of %.0f with significant predictions\n',noise_types_label{i_noise},idx_total,N_subjects*N_subjects);
    end

    FS_here = 10;

    il_nexttile(1);
    set(gca,'XTickLabel',XTL);
    ylabel('Simulation data from');
    set(gca,'YTickLabel',YTL);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'D.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    il_nexttile(2);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',[]);
    text(0,1.05,'E.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    il_nexttile(3);
    text(0,1.05,'F.',opts_txt{:});
    % c = colorbar;
    % c.Label.String = 'Deviance / trial benefit (adim)';
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',[]);
    % set(c,'Limits',[-10.5 0])
    set(gca,'FontSize',FS_here);

    c = colorbar;

    for i_masker = 1:N_maskers
        nexttile(i_masker);
        ht = title(noise_types_label{i_masker});
        set(ht,'FontSize',FS_here);
        xlabel('ACI_s_i_m from');
    end
    h(end+1) = gcf;
    if flags.do_fig8b
        hname{end+1} = 'fig8b-crosspred_model_participant'; 
    end
    if flags.do_fig3b_suppl
        hname{end+1} = 'suppl-fig3b-crosspred_model_participant'; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig12
    between_noise = 1;
    if between_noise == 1
        % Build the exact matrix to plot
        for j = 1:N_subjects
            for k = 1:N_maskers
                PC_pool_here{j}(k,:) = transpose(squeeze(PC_matrix(j,k,:)));
                % Dev_pool_here{j}(k,:) = transpose(squeeze(Dev_matrix(j,k,:)));
            end
        end
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 1:4
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool = PC;
        PC_pool_diag_removed = PC_rem_tmp;
        % Dev_pool = Dev;
         
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 5:8
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool  = [PC_pool; PC];
        PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
        % Dev_pool = [Dev_pool; Dev];
        
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 9:12
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool = [PC_pool; PC];
        PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
        % Dev_pool = [Dev_pool; Dev];
    end
      
    if between_noise == 0
        XTL = [];
        XT = 1:N_subjects;
        for ii = 1:N_subjects
            txt2use = 'S';
            if ii < 10 
                txt2use = [txt2use '0'];
            end
            txt2use = [txt2use num2str(ii)];
            YTL{ii} = txt2use;
            if mod(ii,2)==0
                XTL{ii} = txt2use;
            else
                XTL{ii} = '';
            end
        end
   
        figure('Position',[100 100 800 600]);
        il_tiledlayout(2,3,'TileSpacing','tight');
         
        for i_column = 1:3
            switch mod(i_column,3)
                case 1 % white
                    Colour_here = masker_colours{1};
                case 2 % bump
                    Colour_here = masker_colours{2};
                case 0 % MPS
                    Colour_here = masker_colours{3};
            end
%             Dev_matrix_here = Dev_matrix(:,:,i_column);
%             Dev_diag(:,i_column) = diag(Dev_matrix_here);
%             Dev_matrix_nodiag = Dev_matrix_here-repmat(Dev_diag(:,i_column)',N_subjects,1); % matrix with no diagonal
             
            il_nexttile(i_column)
%             imagesc(1:N_subjects,1:N_subjects,Dev_matrix(:,:,i_column)); hold on
%             caxis([-10 0])
%             colormap('gray')
%             set(gca,'XTick',XT);
%             set(gca,'YTick',XT);
%             for ii = 1:N_subjects
%                 
%                 offxy = 0.5;
%                 offxy = offxy*.9;
%                 il_plot_square(ii,ii,Colour_here,offxy);
%                 % if Dev_diag(ii,i_column) ~= 0
%                 %     idx = find(Dev_matrix_nodiag(:,ii)<0);
%                 %     if ~isempty(idx)
%                 %         for jj = 1:length(idx)
%                 %             il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
%                 %         end
%                 %     end
%                 % end
%             end
             
            PC_matrix_here = PC_matrix(:,:,i_column);
%             PC_diag(:,i_column) = diag(PC_matrix_here);
%             PC_matrix_nodiag = PC_matrix_here-repmat(PC_diag(:,i_column)',N_subjects,1); % matrix with no diagonal
%             
            il_nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([14 24.5])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            
            for ii = 1:N_subjects
                
                offxy = 0.5;
                offxy = offxy*.9;
                il_plot_square(ii,ii,Colour_here,offxy);
                % if Dev_diag(ii,i_column) ~= 0
                %     idx = find(Dev_matrix_nodiag(:,ii)<0);
                %     if ~isempty(idx)
                %         for jj = 1:length(idx)
                %             il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                %         end
                %     end
                % end
            end
        end
        
        FS_here = 10;
        
        il_nexttile(1);
        set(gca,'XTickLabel',[]);
        ylabel('Data from');
        set(gca,'YTickLabel',YTL);
        opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
        text(0,1.05,'A.',opts_txt{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(2);
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        text(0,1.05,'B.',opts_txt{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(3);
        text(0,1.05,'C.',opts_txt{:});
        c = colorbar;
        c.Label.String = 'Deviance / trial benefit (adim)';
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        % set(c,'Limits',[-10.5 0])
        set(gca,'FontSize',FS_here);
        
        nexttile(4);
        ylabel('Simulation data from');
        set(gca,'YTickLabel',YTL);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        il_nexttile(5);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        il_nexttile(6);
        c = colorbar;
        c.Label.String = 'Percentage accuracy benefit (%)';
        % set(c,'Limits',[0 10.5])
        
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        for i_masker = 1:N_maskers
            nexttile(i_masker);
            ht = title(noise_types_label{i_masker});
            set(ht,'FontSize',FS_here);
            
            nexttile(i_masker+3);
            xlabel('ACI_s_i_m from');
        end
        h(end+1) = gcf;
        hname{end+1} = 'fig11-crosspred_model_participant';
        
        %%%
        PC_rem = [];
        for idx_noise = 1:3
            PC_rem{idx_noise} = PC_matrix(:,:,idx_noise)-repmat(diag(PC_matrix(:,:,idx_noise))',N_subjects,1);
            idx = find(PC_rem{idx_noise}>0); 
            length(idx)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Replicating the figure, but showing averages of cross predictions (between subjects)
        figure('Position',[100 100 700*2-200 550]);
        il_tiledlayout(2,2,'TileSpacing','tight');

        XL = [0.5 3.5]; 

        il_nexttile(1);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(2);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(3);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(4);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        for i_masker = 1:N_maskers
            for i_subject = 1:N_subjects

                offx = 0.05*(i_subject - ((N_subjects+2)/2-.5)); %  0.05*i_subject;

                idx_here = 1:N_subjects;
                idx_here(i_subject) = []; % removing the diagonal
                y2plot  = G_PA_sim(:,i_subject,i_masker); % 10 folds x 12 subjects x 
                y2plot2 = PC_rem{i_masker}(idx_here,i_subject);

                MCol = 'w';
                il_nexttile(1);
                errorbar(i_masker+offx, mean(y2plot), 1.64*sem(y2plot),'d', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                il_nexttile(3);
                Me = median(y2plot2);
                errU = max(y2plot2)-Me;
                errL = Me-min(y2plot2);
                errorbar(i_masker+offx, Me, errL, errU, 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                Benefit_re_diag_MI(i_subject,i_masker) = min(y2plot2);
                Benefit_re_diag_MA(i_subject,i_masker) = max(y2plot2);
                
                %%% Incorrect trials:
                y2plot  = nan([N_subjects 1]); 
                y2plot2 = nan([N_subjects 1]); 

                il_nexttile(2);
                Me = mean(y2plot);
                plot(i_masker+offx, Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                il_nexttile(4);
                Me = mean(y2plot2);
                plot(i_masker+offx,Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on            
            end
            y2plot  = mean(G_PA_sim(:,:,i_masker));
            il_nexttile(1);
            offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
            errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        end

        XLabel = 'Masker';
        % subplot(2,1,1)
        il_nexttile(1);
        title('All trials')
        grid on
        xlim(XL);
        % xlabel(XLabel); 
        ylim([4 26])
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', []); 
        ylabel('PA benefit (%)')
        xcoor = 0.01;
        ycoor = 0.94;
        txt_opts = {'Units','Normalized','FontWeight','Bold','FontSize',14};
        text(xcoor,ycoor,'D',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(2);
        title('Incorrect trials only')
        grid on
        xlim(XL);
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', []); 
        set(gca, 'YTickLabels',[]);
        text(xcoor,ycoor,'NaN',txt_opts{:});
        set(gca,'FontSize',FS_here);

        il_nexttile(3);
        grid on
        xlim(XL)
        ylim([-10 5])
        xlabel(XLabel); 
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', noise_types_label); 
        ylabel('PA benefit re. auto prediction (%)')
        text(xcoor,ycoor,'E',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(4);

        xlim(XL)
        ylim([-10 5])
        xlabel(XLabel); 
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', noise_types_label); 
        set(gca, 'YTickLabels',[]);
        % ylabel('Percent accuracy benefit (%)')
        text(xcoor,ycoor,'NaN',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        h(end+1) = gcf;
        hname{end+1} = ['fig11b-metrics-benefit-model-cross'];
        
        % fprintf('Delta PAsim values ranged between %.2f and %.2f\n', ...
        %     min(min(min(G_PA_sim))),max(max(max(G_PA_sim))));
        fprintf('Delta PAsim values ranged between %.1f and %.1f\n', ...
            min(min(PA_mean_chance)),max(max(PA_mean_chance)));
        fprintf('Benefit re. main diagonal between %.1f and %.1f (median=%.1f)\n', ...
            min(min(Benefit_re_diag_MI)), ...
            max(max(Benefit_re_diag_MA)), ...
            prctile([Benefit_re_diag_MI(:) Benefit_re_diag_MA(:)],50));
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure;
        % figure('Position',[100 100 1000 420]); 
        figure('Position',[100 100 560 420]); 
        il_tiledlayout(1,1,'TileSpacing','compact'); % 'tight');
         
        % il_nexttile(1);
        % x_var = 1:size(PC_pool,2);
        % y_var = 1:size(PC_pool,1);
        % imagesc(x_var,y_var,Dev_pool); hold on
        % colormap('gray')
        % caxis([-10 0])
          
        % c = colorbar;
        % c.Label.String = 'Deviance / trial benefit (adim)';
%         
%         for ii = 1:4
%             text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
%         end
%         for ii = 5:8
%             text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
%         end
%         for ii = 9:12
%             text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
%         end
%         %%% Drawing a separator between participants:
%         plot(3.5*[1 1],[0 9.5],'Color','b');
%         plot(3.5*[1 1],[0 9.5],'Color','b');
%         
%         plot(6.5*[1 1],[0 9.5],'Color','b');
%         plot(6.5*[1 1],[0 9.5],'Color','b');
%         
%         plot(9.5*[1 1],[0 9.5],'Color','b');
%         plot(9.5*[1 1],[0 9.5],'Color','b');
%         
%         plot([0 12.5],3.5*[1 1],'Color','b');
%         plot([0 12.5],6.5*[1 1],'Color','b');
%         plot([0 12.5],9.5*[1 1],'Color','b');
%         set(gca,'FontSize',10);
%         %%%
%         
%         disp('')
%         set(gca,'XTick',x_var)
%         set(gca,'YTick',y_var)
%         
%         XTL = repmat({'W','BP','MPS'},1,4);
%         YTL = repmat({'W','BP','MPS'},1,3);
%         set(gca,'XTickLabel',XTL);
%         set(gca,'YTickLabel',YTL);
%         
%         xlabel('ACI from condition');
%         ylabel('Data from condition');
%         
%         text(0,1.05,'D.','Units','Normalized','FontSize',14,'FontWeight','Bold');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % figure;
        il_nexttile(1);
        x_var = 1:size(PC_pool,2);
        y_var = 1:size(PC_pool,1);
        imagesc(x_var,y_var,PC_pool); hold on;
    	colormap('gray')
        % caxis([0 10])
        
        %%% Drawing a separator between participants:
        plot(3.5*[1 1],[0 9.5],'Color','b');
        plot(3.5*[1 1],[0 9.5],'Color','b');
        
        plot(6.5*[1 1],[0 9.5],'Color','b');
        plot(6.5*[1 1],[0 9.5],'Color','b');
        
        plot(9.5*[1 1],[0 9.5],'Color','b');
        plot(9.5*[1 1],[0 9.5],'Color','b');
        
        plot([0 12.5],3.5*[1 1],'Color','b');
        plot([0 12.5],6.5*[1 1],'Color','b');
        plot([0 12.5],9.5*[1 1],'Color','b');
        set(gca,'FontSize',10);
        %%%
        
        PC_pool_diag_removed = transpose(PC_pool_diag_removed);
        offxy = 0.5;
        offxy = offxy*.9;
        for jj = 1:size(PC_pool_diag_removed,1)
            idx_here = find(PC_pool_diag_removed(jj,:)>0);
            for ii = 1:length(idx_here)
                il_plot_square(jj,idx_here(ii),'m',offxy,'--');
            end
        end
        
        c = colorbar;
        c.Label.String = 'Percentage accuracy benefit (%)';
        
        for ii = 1:4
            text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 5:8
            text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 9:12
            text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        
        disp('')
        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
        
        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
        
        text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');
        
        h(end+1) = gcf;
        hname{end+1} = 'fig12a-crosspred_masker-sim';
        
        % c.Label.String = 'Deviance / trial benefit (adim)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure;
        figure('Position',[100 100 560 250]); % 560 420]);
        % il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');
         
        PC_all = (PC_pool(1:3,:)+PC_pool(4:6,:)+PC_pool(7:9,:))/3;
        PC_all = (PC_all(:,1:3)+PC_all(:,4:6)+PC_all(:,7:9)+PC_all(:,10:12))/4;
        
        x_var = 1:size(PC_all,2);
        y_var = 1:size(PC_all,1);
        imagesc(x_var,y_var,PC_all); hold on
        colormap('gray')
        caxis([14.5 23.5])
        set(gca,'FontSize',10);
        
        c = colorbar;
        c.Label.String = 'PA benefit (%)';
         
        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.1f',PC_all(ii,jj)),'FontSize',10,'Color','b');
            end
        end
        % %%%
        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
         
        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
         
        text(0,1.05,'B.','Units','Normalized','FontSize',14,'FontWeight','Bold');
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
         
        % text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

        h(end+1) = gcf;
        hname{end+1} = 'fig12b-crosspred_masker-sim';
        
        fprintf('Between-noise: Delta PAsim values ranged between %.1f and %.1f\n', ...
            min(min(PC_pool)),max(max(PC_pool)));
        % fprintf('Benefit re. main diagonal between %.1f and %.1f (median=%.1f)\n', ...
        %     min(min(Benefit_re_diag_MI)), ...
        %     max(max(Benefit_re_diag_MA)), ...
        %     prctile([Benefit_re_diag_MI(:) Benefit_re_diag_MA(:)],50));
        
        idx_here = find(PC_pool_diag_removed(:,:)>0); 
        fprintf('Number of cases where PA ''between'' was higher than the autoprediction: %.0f (out of %.0f)\n', ...
            length(idx_here),length(PC_pool_diag_removed(:))-3*N_subjects);
        
        idx_here_above2   = find(PC_pool_diag_removed(:,:)>=2); 
        idx_here_above1p5 = find(PC_pool_diag_removed(:,:)>=1.5); 
        
        idx_here_below1   = find(PC_pool_diag_removed(idx_here)<1);
        
        prctile(PC_pool_diag_removed(idx_here),90)
        prctile(PC_pool_diag_removed(idx_here),75)
        prctile(PC_pool_diag_removed(idx_here),50)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supplementary materials

% Inspired by ../Publications/Osses2022_JARO/ge20220401_Osses2022_JARO.m

if flags.do_fig1_suppl || flags.do_fig2_suppl
    % Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'};
    Gender     = {'M'  ,'M'  ,'F'  , 'M'  ,'F'  ,'M'  ,'M'  ,'F'  ,'M'  ,'F'  ,'M'  ,'M'};
    Language   = {'FR' ,'ES' ,'FR' , 'FR' ,'IT' ,'FR' ,'FR' , 'TR','FR' ,'FR' ,'ES' ,'FR'};
end

if flags.do_fig1_suppl 
    
    FontSize = 12;
    row4nan = [];
    
    f = [];
    
    Age        = [33    36     31    38    24    43    23    27   25     22    36    22]; mean(Age)
    
    Colours = {rgb('Brown'),rgb('Bisque'),rgb('NavajoWhite'),rgb('BurlyWood'),rgb('Peru'),rgb('Goldenrod'), ...
               rgb('RosyBrown'),rgb('Chocolate'),rgb('Sienna'),rgb('DarkGoldenrod'),rgb('Pink'),rgb('LightCoral')};
    
    for i_subj = 1:length(Subjects)
        dir_res = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep 'data_' ...
            Subjects{i_subj} filesep '1-experimental_results' filesep 'audiometry' filesep]; 
        
        file = Get_filenames(dir_res,'absthreshold_*.dat');
        if ~isempty(file)
            data_now = datread([dir_res file{1}]); % left ear: 'absthreshold_'  '_l.dat
            if isempty(f)
                f = transpose(data_now(:,1));
            end
            aud_l(i_subj,:) = transpose(data_now(:,2));
            
            data_now = datread([dir_res file{2}]); % right ear
            aud_r(i_subj,:) = transpose(data_now(:,2));
            
            disp('')
        else
            row4nan(end+1) = i_subj;
        end
    end

    % Average thresholds from freqs: 250, 500, 1000, 2000 and 4000 (BSA_RP_PTA_Final-2018-August.pdf)
    %   i.e., all frequencies tested, but excluding 8 kHz.
    f_idxs = 1:length(f);
    aud_avg(:,1) = mean(aud_l(:,f_idxs),2);
    aud_avg(:,2) = mean(aud_r(:,f_idxs),2);
    
    [~,aud_best_ear] = min(aud_avg,[],2);
    idxs = find(aud_best_ear==1);
    aud_all(idxs,:) = aud_l(idxs,:);
    
    idxs = find(aud_best_ear==2);
    aud_all(idxs,:) = aud_r(idxs,:);
    
    aud_Me_l = mean(aud_l);
    aud_Me_r = mean(aud_r);
    
    for i = 1:length(Subjects)
        aud_Me_across_f(i,1) = aud_avg(i,aud_best_ear(i));
    end
    % aud_Me_across_f_L = aud_avg(:,1);
    % aud_Me_across_f_R = aud_avg(:,2);
    [mi,mi_idx] = min(aud_Me_across_f);
    [ma,ma_idx] = max(aud_Me_across_f);
    fprintf('\tThe participants had average thresholds between %.1f (S%.0f) and %.1f dB~HL (S%.0f) or their best ear\n',mi,mi_idx,ma,ma_idx);
    
    %%% figure format:
    prefix = 'suppl-fig1-';
    Xlab = 'Frequency (Hz)';
    Ylab = 'Audiometric thresholds (dB HL)';
    YL = [-75 15]; % dB HL (inverted scale)
    XT = [125 250 500 1000 2000 4000 8000 16000];
    for i = 1:length(XT)
        if XT(i) ~= XT(end)
            XTL{i} = num2str(XT(i));
        else
            XTL{i} = 'avg';
        end
        
    end
    %%% End figure format

    figure; % ('Position',[100 100 1000 250*N_targets*correction]); 
    il_tiledlayout(2,1,'TileSpacing','tight'); % 'tight'); % none'); % 'Compact');
    
    %--- Fig. 1A ----------------------------------------------------------
    il_nexttile(1)
    idx2label = [];
    for i = 1:length(Subjects)
        if aud_best_ear(i)==1
            LW = 3;
            FaceColour = Colours{i};
            idx2label(end+1) = i;
            Style = '-';
        else
            LW = 1.5;
            FaceColour = 'w';
            Style = '--';
        end
        hpl(i) = semilogx(f,-aud_l(i,:),'-','Color',Colours{i},'LineStyle',Style,'LineWidth',LW); grid on, hold on;
        plot(XT(end),-aud_avg(i,1),'o','Color',Colours{i},'MarkerFaceColor',FaceColour,'LineWidth',LW,'MarkerSize',7);
    end
    plot(f,-aud_Me_l,'ko-','MarkerFaceColor','k','LineWidth',2);
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',[]);
    text(0.03,0.9,'A. Left ear','Units','Normalized','FontWeight','bold','FontSize',14);
    % set(gca,'XTickLabel',XTL);
    Ylabel(Ylab,FontSize);
    % Xlabel(Xlab,FontSize);

    ylim([-42 12]);
    xlim([200 24000])
    
    YT = -40:5:10;
    YT_label = -YT;
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YT_label);
 
    legend(hpl(idx2label),Subjects(idx2label),'Location','SouthEast');
    
    il_nexttile(2)
    
    idx2label = [];
    for i = 1:length(Subjects)
        if aud_best_ear(i)==2
            LW = 3;
            FaceColour = Colours{i};
            idx2label(end+1) = i;
            Style = '-';
        else
            LW = 1.5;
            FaceColour = 'w';
            Style = '--';
        end
        hpl(i) = semilogx(f,-aud_r(i,:),'-','Color',Colours{i},'LineStyle',Style,'LineWidth',LW); grid on; hold on;
        plot(XT(end),-aud_avg(i,2),'o','Color',Colours{i},'MarkerFaceColor',FaceColour,'LineWidth',LW,'MarkerSize',7);
    end
    plot(f,-aud_Me_r,'ko-','MarkerFaceColor','k','LineWidth',2);
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',[]);
    text(0.03,0.9,'B. Right ear','Units','Normalized','FontWeight','bold','FontSize',14);
    set(gca,'XTickLabel',XTL);
    Ylabel(Ylab,FontSize);
    Xlabel(Xlab,FontSize);

    ylim([-42 12]);
    xlim([200 24000])
    
    YT = -40:5:10;
    YT_label = -YT;
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YT_label);
    
    legend(hpl(idx2label),Subjects(idx2label),'Location','SouthEast');
    
    % ha_here(end+1) = gca;
    h(end+1) = gcf; % current handle
	hname{end+1} = [prefix 'PTA'];  
    
    Pos = get(gcf,'Position');
    Pos(4) = 600;
    set(gcf,'Position',Pos); % '570         754
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2_suppl
    
    thres_steady = nan([length(Subjects) 4]);
    list_steady  = nan([length(Subjects) 4]);
    thres_mod8   = nan([length(Subjects) 4]);
    list_mod8    = nan([length(Subjects) 4]);
    
    for i_subj = 1:length(Subjects)
        dir_res = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep 'data_' ...
            Subjects{i_subj} filesep '1-experimental_results' filesep 'intellitest' filesep]; 
        file1 = [dir_res 'Intellitest_' Subjects{i_subj} '_varspeech16.dat'];
        file2 = [dir_res 'Intellitest_' Subjects{i_subj} '_varspeech16_mod8.dat'];
        
        var = datread(file1);
        N = size(var,1);
        
        thres_steady(i_subj,1:N) = transpose(var(:,2));
        list_steady(i_subj,1:N)  = transpose(var(:,1));
        
        disp('')
        var = datread(file2);
        
        thres_mod8(i_subj,1:N) = transpose(var(:,2));
        list_mod8(i_subj,1:N)  = transpose(var(:,1));
        
    end
    
    thres_all     = nanmean(thres_steady,2);
    thres_mod_all = nanmean(thres_mod8  ,2);
    
    idx_FR    = find( strcmp(Language,'FR'));
    idx_Other = find(~strcmp(Language,'FR'));
    
    offx = 0.1;
    varx = 1; % ones(size(thres_all));
    vary = [thres_all thres_mod_all];
    figure('Position',[100 100 570 280]);
    hpl1 = plot([varx-offx varx+offx], vary(idx_FR,:)   , 'bo-','MarkerFaceColor','w'); hold on, grid on
    Me   = mean(vary(idx_FR,:));
    errL = sem(vary(idx_FR,:));
    errU = sem(vary(idx_FR,:));
    errorbar([varx-1.2*offx varx+1.2*offx], Me,errL,errU, 'bo','MarkerFaceColor','b','LineWidth',2);
    
    fprintf('Me     FR for steady=%.1f dB, mod8=%.1f\n',Me);
    
    hpl2 = plot(.5+[varx-offx varx+offx], vary(idx_Other,:), 'rs-','MarkerFaceColor','w');
    % plot(ones(size(thres_all))-offx, thres_all, 'rs-');
    Me   = mean(vary(idx_Other,:));
    errL = sem(vary(idx_Other,:));
    errU = sem(vary(idx_Other,:));
    errorbar([varx+.5-1.2*offx varx+.5+1.2*offx], Me,errL,errU, 'rs','MarkerFaceColor','r','LineWidth',2);
    
    fprintf('Me non-FR for steady=%.1f dB, mod8=%.1f\n',Me);
    
    xlim([0.5+.3 1.8]);
    set(gca,'XTick',[1-offx 1+offx 1.5-offx 1.5+offx])
    set(gca,'XTickLabel',{'steady','8-Hz','steady','8-Hz'})
    legend([hpl1(1) hpl2(1)],'French','Non-French');
    
    ylabel('SRT (dB)');
    xlabel('Condition');
    
    h(end+1) = gcf;
    hname{end+1} = 'suppl-fig2-Intellitest';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5_suppl
    prefix = 'suppl-fig5-';
    
    experiment_full = 'speechACI_Logatome-abda-S43M';
    experiment = strsplit(experiment_full,'-'); experiment = experiment{1};
    noise_here = noise_str{3}; % MPS noise
    subj_ID = 'S01';
    dir_ACI = [fastACI_paths('dir_data') experiment_full filesep subj_ID 'Results' filesep 'Results_ACI' filesep];
    
    if ~exist(dir_ACI,'dir')
        % Alejandro's local computer, all subjects are there:
        dir_ACI = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc_Leo/';
    end
    if ~exist(dir_ACI,'dir')
        error('None of the two ACI folder locations seem to contain the ACI of the participant...');
    end
    fname = ['ACI-' subj_ID '-' experiment '-' noise_here '-nob-gt-l1glm+pyrga-rev4.mat'];
    var = load([dir_ACI fname]);
    ACIs = var.results.ACIs;
    idxs = [1 8 14 20];
    
    figure('Position',[100 100 500 280]);
    il_tiledlayout(1,length(idxs),'TileSpacing','tight');
    
    NfrequencyTicks = 8;
    flags_extra = {'NfrequencyTicks',NfrequencyTicks,'colorbar','no'};
    
    t_idx = 10:70;
    for i = 1:length(idxs)
        il_nexttile(i);
        ACI_here = squeeze(ACIs(idxs(i),:,:));
        affichage_tf(ACI_here(:,t_idx),'CI',var.cfg_ACI.t(t_idx),var.cfg_ACI.f,flags_extra{:});
        if i~= 1
            ylabel('')
            set(gca,'YTickLabel',[]);
        end
        % set(gca,'FontSize',11);
        text(0.45,0.92,['\lambda_i_d_x=' num2str(idxs(i))],'Units','Normalized','FontWeight','Bold','FontSize',8);
    end
    h(end+1) = gcf; % current handle
	hname{end+1} = [prefix 'ACI-lambda'];
    
    disp('')
  
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.h = h;
data.hname = hname;
 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_add_formants(fname,cfg_ACI)

LW = 1;

%%% I still need to spot the exact values:
par_formants.timestep = 0.01; % positive timestep 0.01
par_formants.nformants = 5; % positive nformants 5

%%% Unsure:
% Formants
par_formants.maxformant = 5500; % positive maxformant 5500
par_formants.windowlength = 0.025;% 0.025 % positive windowlength 0.025
par_formants.dynamicrange = 30; % positive dynamic range 20

% F0
par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
%par_formants.pitchfloor = 100; % previous parameter value (14/10/2022)
par_formants.pitchfloor = 50; % positive pitch floor 100 (for f0)
par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

% Before 4/11/2021, I_min set to 40 dB:
par_formants.I_min = 59;%75; %, arbitrary value

affichage_tf_add_Praat_metrics_one_sound(fname,cfg_ACI,par_formants, '-', 'w',LW);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_To_dB(insig)

outsig = log(abs(insig));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [PC_tbt, Dev_tbt, PC_fold, Dev_fold] = il_tbtpred(cfg_ACI,results)
%   Version by Leo 
%
% N_lambdas = length(results.lambdas);
% N_folds = cfg_ACI.N_folds;
% 
% PC_all    = nan(N_lambdas,cfg_ACI.N_folds,cfg_ACI.N_trials);
% Dev_all   = nan(N_lambdas,cfg_ACI.N_folds,cfg_ACI.N_trials);
% 
% N = length(results.FitInfo.PC_test_t);
% PC_all = nan([N_lambdas, N_folds, N]);
% PC_all2 = nan(size(PC_all));
% Dev_all2 = nan(size(PC_all));
% 
% for i_lambda = 1:N_lambdas
%     
%     
%     for i_fold = 1:N_folds
%         PC_test_t  = squeeze(results.FitInfo.PC_test_t(i_lambda,i_fold,:));
%         Dev_test_t = squeeze(results.FitInfo.Dev_test_t(i_lambda,i_fold,:));
%         CVtest = find(results.FitInfo.CV.test(i_fold));
%         N_here = length(CVtest);
%         
%         if length(CVtest)~=length(Dev_test_t)
%             %fprintf(['unequal CVtest and MSEtest_t at fold # ' num2str(i_fold) '\n'])
%             if (length(CVtest)==length(Dev_test_t)-1) && (Dev_test_t(end)==0)
%                 PC_test_t = PC_test_t(1:end-1);
%                 Dev_test_t = Dev_test_t(1:end-1);
%             else
%                 error('Problem with the length of vectors Dev_test_t and PC_test_t')
%             end
%         end
%         PC_all(i_lambda,i_fold,CVtest)  = PC_test_t;
%         Dev_all(i_lambda,i_fold,CVtest) = Dev_test_t;
%         
%         PC_all2(i_lambda,i_fold,1:N_here)  = PC_test_t;
%         Dev_all2(i_lambda,i_fold,1:N_here) = Dev_test_t;
%     end
% end
% PC_tbt = squeeze(nanmean(PC_all,2));
% Dev_tbt = squeeze(nanmean(Dev_all,2));
% 
% PC_fold = nanmean(PC_all2,3); % across third dimension
% Dev_fold = nanmean(Dev_all2,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PC_fold, Dev_fold] = il_tbtpred_inc(cfg_ACI,results,data_passation)
% Adapted version by Alejandro

N_lambdas = length(results.lambdas);
N_folds = cfg_ACI.N_folds;

N = length(results.FitInfo.PC_test_t);
PC_all = nan([N_lambdas, N_folds, N]);
Dev_all = nan(size(PC_all));

if nargin >= 3
    idxs_inc  = find(data_passation.is_correct(cfg_ACI.idx_analysis)==0);
    idxs_corr = find(data_passation.is_correct(cfg_ACI.idx_analysis)==1); % correct indexes
else
    fprintf('All trials will be returned...')
end
                
for i_lambda = 1:N_lambdas
    for i_fold = 1:N_folds
        PC_test_t  = squeeze(results.FitInfo.PC_test_t(i_lambda,i_fold,:));
        Dev_test_t = squeeze(results.FitInfo.Dev_test_t(i_lambda,i_fold,:));
        CVtest = find(results.FitInfo.CV.test(i_fold));
        
        if nargin >= 3
            [idx_corr_this_fold,idx_corr_sort] = intersect(CVtest,idxs_corr);
        end
        % idx_inc_this_fold  = intersect(CVtest,idxs_inc);
        N_here = length(CVtest);
        
        if length(CVtest)~=length(Dev_test_t)
            %fprintf(['unequal CVtest and MSEtest_t at fold # ' num2str(i_fold) '\n'])
            if (length(CVtest)==length(Dev_test_t)-1) && (Dev_test_t(end)==0)
                PC_test_t = PC_test_t(1:end-1);
                Dev_test_t = Dev_test_t(1:end-1);
            else
                error('Problem with the length of vectors Dev_test_t and PC_test_t')
            end
        end
        
        PC_all(i_lambda,i_fold,1:N_here)  = PC_test_t;
        Dev_all(i_lambda,i_fold,1:N_here) = Dev_test_t;
        
        if nargin >= 3
            % setting back to NaN those trials within the fold that were correct:
            PC_all(i_lambda,i_fold,idx_corr_sort)  = nan; 
            Dev_all(i_lambda,i_fold,idx_corr_sort) = nan;
        end
    end
end

PC_fold = nanmean(PC_all,3); % across third dimension
Dev_fold = nanmean(Dev_all,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_plot_square(xcoor,ycoor,Colour_here,offxy,LineStyle)

if nargin < 5
    LineStyle = '-';
end
plot_opts = {'Color',Colour_here,'LineStyle',LineStyle};
% Plots a box:
plot([xcoor-offxy xcoor+offxy],(ycoor-offxy)*[1 1]      ,plot_opts{:}); % bottom side
plot((xcoor-offxy)*[1 1]      ,[ycoor-offxy ycoor+offxy],plot_opts{:}); % left side
plot((xcoor+offxy)*[1 1]      ,[ycoor-offxy ycoor+offxy],plot_opts{:}); % right side
plot([xcoor-offxy xcoor+offxy],(ycoor+offxy)*[1 1]      ,plot_opts{:}); % top side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_plot_square_XY(xcoor,ycoor,Colour_here,offx,offy,LineStyle,LW)

if nargin < 5
    LineStyle = '-';
end
plot_opts = {'Color',Colour_here,'LineStyle',LineStyle,'LineWidth',LW};
% Plots a box:
plot([xcoor-offx xcoor+offx],(ycoor-offy)*[1 1]      ,plot_opts{:}); % bottom side
plot((xcoor-offx)*[1 1]      ,[ycoor-offy ycoor+offy],plot_opts{:}); % left side
plot((xcoor+offx)*[1 1]      ,[ycoor-offy ycoor+offy],plot_opts{:}); % right side
plot([xcoor-offx xcoor+offx],(ycoor+offy)*[1 1]      ,plot_opts{:}); % top side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Me_diag,Me_nodiag,opts] = il_avg_diag_nodiag(insig)

Me_diag = mean(diag(insig));
opts.Me_diag_sem = 1.64*sem(diag(insig));

for i = 1:size(insig,1)
    insig(i,i) = NaN;
end
Me_nodiag = nanmean(insig(:));
insig = insig( find(~isnan(insig)) );
opts.Me_nodiag_sem = 1.64*sem(insig(:));
