function data = publ_osses2022b_JASA_figs(varargin)
% function data = publ_osses2022b_JASA_figs(varargin)
%
% 1. Description: Generates the figures
%
% % To display Fig. 1 of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig1'); % T-F representations /aba/, /ada/
%
% % To display Fig. 2a of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig2a'); % T-F representation of one realisation of each noise type
%
% % To display Fig. 2b of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig2b'); % Acoustic analysis of each noise type
%
% % To display Fig. 3 of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig3'); % Illustration of the Gaussian basis elements
%
% % To display Fig. 4 of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig4'); % SNR thresholds averaged across participants
%
% % To display Fig. 5 of Osses and Varnet, (2022, BioRxiv) use :::
%     publ_osses2022b_JASA_figs('fig5'); % Correct responses, dprime, criterion averaged across participants
%
% 'fig1' requires the speech samples from 'S01'
% 'fig2a','fig2b' require the noise samples from 'S01'
% 'fig3' does not require any additional data
% 'fig4','fig5' requires the savegame files from all ('S01'-'S12') participants. 
%        these data are taken directly from the publ_osses2022b data in fastACI.
%
% Author: Alejandro Osses
% See also: g20220716_get_paper_figures_JASA_ACI.m (from fastACI_sim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all, clc

if nargin == 0
    help publ_osses2022b_JASA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig1', ...
    'fig2a', ... % spectrogram white, bump, MPS noises
    'fig2b', ... % Other metrics to characterise white, bump, and MPS noises
    'fig3', ...  % Illustration of the Gaussian basis elements
    'fig4', ...
    'fig5_simu', ... % old name: fig4_simu
    'fig5', 'fig5_stats', ... % dprime and statistics
    'fig6', ...  % ACIs participants, L1109
    'fig7', ...
    'fig8', ...  % Cross predictions, L1608
    'fig8b', ... % Cross predictions 'between participant' for simulations
    'fig9', ...  % Cross predictions between noise
    'fig9b', ... % Cross predictions between noise for the simulations
    'fig10', ...
    'fig1_suppl', ... % Audiograms
    'fig2_suppl', ... % Intellitest, L3209 (on 29/09/2022)
    'fig3_suppl', ...
    'fig3b_suppl', ...
    'fig4_suppl', ...
    'fig4b_suppl', ...
    'fig5_suppl', ... % correlation plots
    'fig5b_suppl', ... % correlation plots, simulations    
    'fig6_suppl', ...
    'fig7_suppl'}; % different lambdas
definput.flags.plot={'plot','no_plot'};
definput.flags.local={'local','zenodo'};

definput.keyvals.dir_zenodo=[];
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

type_corr = 'Pearson'; % In the end, Pearson seems to be good enough
% type_corr = 'Spearman';

do_fig4_stats_N10 = 0; % by default no analysis of N=10 data
do_fig5_stats_N10 = 0; % by default no analysis of N=10 data
% End common variables
do_plot = flags.do_plot; % flags.do_plot into variable for easy debugging

bZenodo = flags.do_zenodo;
if bZenodo == 1
    fprintf('%s: The data from Zenodo will be loaded... \n',mfilename);
    if ~isempty(keyvals.dir_zenodo)
        fprintf('\t The data from Zenodo will be loaded from %s \n',keyvals.dir_zenodo);
        if exist(keyvals.dir_zenodo,'dir')
            bIs_dir = 1;
        else
            bIs_dir = 0;
        end
    else
        bIs_dir = 0;
    end
    if bIs_dir == 0
        msg_here = sprintf('\t You indicated that you want to proccess the data obtained from Zenodo. Please indicate the directory where it is...\n');
        fprintf(msg_here)
        keyvals.dir_zenodo = uigetdir([pwd filesep],msg_here);
        keyvals.dir_zenodo = [keyvals.dir_zenodo filesep];
    end
end

% colourbar_map = 'default';

colourbar_map = 'DG_jet';
% colourbar_map = 'DG_red_blue'; % Double gradient
% colourbar_map = 'SQ_red'; % Sequential red
flags_tf = {'colourbar_map',colourbar_map};

if bZenodo == 1
    dir_zenodo = keyvals.dir_zenodo;
    dir_ACI_exp = [dir_zenodo '03-Post-proc-data' filesep 'ACI_exp' filesep];
    dir_ACI_sim = [dir_zenodo '03-Post-proc-data' filesep 'ACI_sim' filesep]; % Copied from /20220715_sim_Q1_osses2022a/';
    
    dir_data = [dir_zenodo '01-Stimuli' filesep 'fastACI_data' filesep];
    dir_fastACI      = [dir_zenodo '02-Raw-data' filesep 'fastACI' filesep];
    dir_savegame_sim = [dir_zenodo '02-Raw-data' filesep 'ACI_sim' filesep];

else
    dir_ACI_exp = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc_Leo/'; % My run: dir_where = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc/';
        % Run-S01/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
        % Run-S02/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
        % ...
        % Run-S12/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
    dir_ACI_sim = '/home/alejandro/Documents/Databases/data/fastACI_data_z_tmp/20220715_sim_Q1_osses2022a/';
    dir_data = fastACI_paths('dir_data');
    dir_fastACI = fastACI_paths('dir_fastACI'); % fastACI_basepath
    dir_savegame_sim = [fastACI_dir_data 'speechACI_Logatome-abda-S43M' filesep 'osses2022a_debug_0715' filesep]; % '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/osses2022a_debug_0715/';

end
dir_savegame_exp = [dir_fastACI 'Publications' filesep 'publ_osses2022b' filesep]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
dir_subj_exp = [dir_data experiment filesep 'S01' filesep];

if flags.do_fig1 || flags.do_fig2a || flags.do_fig2b    
    %%% Checking the waveforms:
    if ~exist(dir_subj_exp,'dir')
        fprintf('The directory %s is not found on disk\n',dir_subj_exp);
        bRetrieve = input('    Do you want to retrieve this directory (this might take several minutes)? (1=yes; 0=no): ');
        
        if bRetrieve
            publ_osses2022b_utils('S01',[],'Check_sounds_for_a_participant');
        end
    end
    %%%
    
    dir_noise = [dir_subj_exp 'NoiseStim-white' filesep];
    
    basef = 8000;
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5, ...
        'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};
    
    CLim_base_e = [-8 -4.5];
    CLim = [-69.5 -39]; % obtained as 20*log10(exp(CLim_base_e))
    
    flags_extra = {'NfrequencyTicks',8,'colorbar','no'};
end

% if flags.do_fig1_old % Speech-in-noise representation
%     error('Removed from this script on 13/12/2022')
%     % Originally taken from: l20220602_Figures_ModulationGroup.m (by Leo)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1
    % Taken from: l20220602_Figures_ModulationGroup.m (by Leo)
    fname_aba = [dir_subj_exp 'speech-samples' filesep 'S43M_ab_ba.wav'];
    fname_ada = [dir_subj_exp 'speech-samples' filesep 'S43M_ad_da.wav'];
    [aba, fs] = audioread(fname_aba);
    [ada, fs] = audioread(fname_ada);
    
    %%% Plot targets + 1 example of noise
    [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
    [G_ada, fc, t, outs] = Gammatone_proc(ada, fs, flags_gamma{:});

    figure('Position',[100 100 500 200]); 
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile(1); 
    G_dB = To_dB(G_aba');
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
    
    %%% Parameters for Praat:
    par_formants = [];
    par_formants.timestep = 0.01; % positive timestep 0.01
    par_formants.nformants = 5; % positive nformants 5

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
    
    Add_formants(fname_aba,cfg_in,par_formants);
    %%%
    
    % flags_extra{end} = 'on';
    nexttile(2); 
    G_dB = To_dB(G_ada');
    affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); hold on;
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',XTL);
    
    caxis(CLim); 
    title('/ada/ target');
    ylabel('');
    set(gca,'YTickLabels',[]); %xlabel('')

    % fname1 = [cfg_ACI.dir_target files{1}];
	% fname2 = [cfg_ACI.dir_target files{2}];
           
    Add_formants(fname_ada,cfg_in,par_formants);
    
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
    fname_no1 = [dir_subj_exp 'NoiseStim-white'         filesep fname];
    fname_no2 = [dir_subj_exp 'NoiseStim-bumpv1p2_10dB' filesep fname];
    fname_no3 = [dir_subj_exp 'NoiseStim-sMPSv1p3'      filesep fname];
    
    [insig(:,1),fs] = audioread(fname_no1);
    insig(:,2) = audioread(fname_no2);
    insig(:,3) = audioread(fname_no3);
    
    figure('Position',[100 100 500 220]); 
    tiledlayout(1,3,'TileSpacing','tight');
    
    for i = 1:3
        %%% Plot targets + 1 example of noise
        [G, fc, t, outs] = Gammatone_proc(insig(:,i), fs, flags_gamma{:});
     
        dBFS = 100;
        floor2use = 0;
        G_dB(:,:,i) = max(20*log10(G)'+dBFS,floor2use); % To_dB(G');
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
        
        nexttile(i); 
        affichage_tf(G_dB(idx_f,idx_t,i), 'pow', t(idx_t), fc(idx_f), flags_extra{:}); hold on;
        caxis([0 dB_max])
        
        if i == 1
            text(-0.15,1.07,'A.','Units','Normalized','FontSize',12,'FontWeight','Bold');
        end
        if i ~= 1
            set(gca,'YTickLabel',[]);
            ylabel('');
        end
        if i == 2
            xlabel('Time (s)')
        else
            xlabel('')
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
	% N_sounds = 50; warning('temporal value'); % arbitrary choice for a quick plotting
    
    % Band levels are always computed
    do_kohlrausch2021_env = 1; % flags.do_fig2c; % Modulation spectrum
    % type_env = 'kohlrausch2021_env_noDC'; % as in the preprint
    type_env = 'kohlrausch2021_env'; % as in the new manuscript
    switch type_env
        case 'kohlrausch2021_env_noDC'
            % For preprint
            suff_DC = ''; 
            fmod_ylim = [10 50];
            suff_dB = '';
            
        case 'kohlrausch2021_env'
            % For JASA paper
            % With DC but also referenced to 0 dB
            
            suff_DC = '-with-DC'; 
            % fmod_ylim = [10 70];
            yoff = 70;
            % fmod_ylim = [10 70]-yoff; % showing the DC
            fmod_ylim = [20 45]-yoff; % not showing the DC
            suff_dB = ' re max.';
            
    end
    
	fmod_xlim = [0 60];
    fc_ylim = [35 65];
    percL = 25;
    percU = 75;
    
    % %%% Preprint:
    % fc_ylim = [30 70];
    % percL =  5;
    % percU = 95;
    %%%
    
    figure('Position',[100 100 650 450]); % before: Pos(4) = 650
    tiledlayout(2,3,'TileSpacing','tight');
    
    % i replaced by i_noise
	for i_noise = 1:length(noise_str)
        dir_where = [dir_subj_exp 'NoiseStim-' noise_str{i_noise} filesep];
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
                    type_env);
            end
 
            [outsig1,fc] = auditoryfilterbank(insig,fs);
            if j == 1
                t = (1:size(outsig1,1))/fs;
            end
            lvls(j,:) = rmsdb(outsig1) + dBFS; % These are the band levels
 
            if mod(j,50) == 1
                fprintf('\tProcessing sound %.0f of %.0f\n',j,N_sounds);
            end
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Band levels:
        L1me = prctile(lvls,50);
        L1perL = prctile(lvls,percL);
        L1perU = prctile(lvls,percU);
        
        nexttile(i_noise)
        semilogx(fc,L1perU,'Color',rgb('Gray')); hold on;
        plot(fc,L1perL,'Color',rgb('Gray'));
        plot(fc,L1me,'o-','Color',rgb('Maroon'),'MarkerFaceColor',rgb('Maroon')); grid on
               
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

        if do_kohlrausch2021_env
            idx_env = find(env_extra.f_env<=fmod_xlim(2));
            
            f_env  = env_extra.f_env(idx_env);
            % extra.fs_env = env_extra.fs_env;
            env_dB   = prctile(env_dB_full(idx_env,:),50,2);
            [DC_empirical(i_noise),idx_DC_empirical(i_noise)] = max(env_dB);
            if i_noise == 1
                DC = DC_empirical(i_noise);
            end
            switch type_env
                case 'kohlrausch2021_env'
                    env_dB = env_dB - DC;
                    env_dB_full = env_dB_full-DC;
            end
            env_dB_U = prctile(env_dB_full(idx_env,:),percU,2);
            env_dB_L = prctile(env_dB_full(idx_env,:),percL,2);
            
            nexttile(3+i_noise) % when there were 3 subpanels: (6+i_noise)
            plot(f_env,env_dB_U,'-','Color',rgb('Gray'),'LineWidth',2); hold on;
            plot(f_env,env_dB_L,'-','Color',rgb('Gray'),'LineWidth',2);
            plot(f_env,env_dB  ,'-','Color',rgb('Maroon'),'LineWidth',2); grid on
            
            if i_noise == 1
                data.f_env = f_env;
            end
            data.env_dB_U(:,i_noise) = env_dB_U(:);
            data.env_dB_L(:,i_noise) = env_dB_L(:);
            data.env_dB(:,i_noise)   = env_dB(:);
            data.DC_empirical(i_noise) = DC_empirical(i_noise);
            
            if i_noise == 1
                ylabel(sprintf('Envelope spectrum (dB %s)',suff_dB))
            else
                set(gca,'YTickLabel',[]);
            end
            if i_noise==2
                xlabel('Modulation frequency (Hz)')
            end

            if i_noise == 1
                text(0,1.05,'C.','Units','Normalize','FontWeight','Bold','FontSize',13);
            end
            % title(' ')
            % title(noise_types_label{i_noise});
            
            deltaf = 10;
            XT = deltaf:deltaf:fmod_xlim(2)-deltaf;
            set(gca,'XTick',XT);
            %%% Preprint:
            % set(gca,'YTick',fmod_ylim(1)+5:5:fmod_ylim(2)-5);
            set(gca,'YTick',fmod_ylim(1)+4:4:fmod_ylim(2)-4);
            xlim(fmod_xlim);
            ylim(fmod_ylim);

            % h(end+1) = gcf;
            % hname{end+1} = sprintf('Env-N-%.0f%s',N_sounds,suff);
        end
 
        % data.V_overall(:,i_noise) = V_overall(:);
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
        affichage_tf( squeeze(sumReWeight(2,:,:)), 'CI', cfg_ACI.t, cfg_ACI.f, flags_tf{:});
    else
        affichage_tf( sumReWeight, 'CI', cfg_ACI.t, cfg_ACI.f, flags_tf{:});
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

%  g20220720_behavstats_from_l20220325
%%% Analysis for 12 participants:
Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'}; 
N_subjects = length(Subjects);

if flags.do_fig4 || flags.do_fig5_simu || flags.do_fig5 || flags.do_fig5_stats
    %%% If analysis for the first 10 participants: 
    % % Subjects = {'S01','S02','S03','S04','S05','S07','S08','S09','S10','S11'};
    
    % % dir_experiment = [fastACI_dir_data experiment filesep]; % this is the default directory
    % dir_experiment = ['/home/alejandro/Desktop/fastACI_today/fastACI_data/' experiment filesep]; % /S01/Results/Results_final

    Markers = {'o','s','^'};
    LW = [1 1 2];

    for i_subject = 1:N_subjects
        subject = Subjects{i_subject};
        
        if flags.do_fig4 || flags.do_fig5 || flags.do_fig5_stats % Human listeners:
            dir_subject = [dir_savegame_exp 'data_' Subjects{i_subject} filesep '1-experimental_results' filesep];
        end
        if flags.do_fig5_simu
            % model = 'osses2022a_debug_0715';
            dir_subject = [dir_savegame_sim 'Results-' Subjects{i_subject} '-v1' filesep];
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
            file_savegame = [dir_results fname_results];
            load(file_savegame);
            %subj = Get_subjectname_from_dirname(cfg_game.dir_noise);

            outs_SNR = il_get_SNR(file_savegame,cfg_game); % New assessment of thresholds in an inline function
            
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
            
            for i_session = 1:N_sessions
                %%% Leo's way: Analysis only using the reversals: -------------
                % sessions are redefined as blocks of 400 trials
                idxi = 1+400*(i_session-1);
                idxf = 400*(i_session);
                [rev, idx] = Get_mAFC_reversals(SNR(idxi:idxf));
                
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
                % idxs2use_all = [idxs2use_all idxs2use]; % collating the idxs of the kept trials
                % SNR_here = SNR(idxs2use);
                
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
            medianSNR(i_subject, i_masker, :) = outs_SNR.SNR_me(:);    % median(SNR_here);
            percSNR_L(i_subject, i_masker, :) = outs_SNR.SNR_percL(:); % prctile(SNR_here,5);
            percSNR_U(i_subject, i_masker, :) = outs_SNR.SNR_percU(:); % prctile(SNR_here,95);
            idxs2use_all = outs_SNR.idxs2use_all;
            
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

        data.medianSNR = medianSNR;
        data.file_savegame = file_savegame;
    end
    
    if flags.do_fig4 || flags.do_fig5_simu
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

if flags.do_fig4 || flags.do_fig5_simu
    
    meanSNR_part_noise = mean(medianSNR(:,:,:),3); % average across test block
    data.meanSNR = meanSNR_part_noise;
    data.meanSNR_description = 'Average across test block for each participant and noise';
    for i_subject = 1:N_subjects
        meanSNR_participant(i_subject,:) = mean(squeeze(medianSNR(i_subject,:,:)));
        meanSNR_overall(i_subject) = mean(meanSNR_participant(i_subject,:),2);
    end
    [ma,ma_idx] = min(meanSNR_overall); % minimum is best performance
    [mi,mi_idx] = max(meanSNR_overall); % maximum is worst performance
    
    if do_plot
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
    end
end % end do_fig4

if flags.do_fig5 
    % Fig. 2: Performance as a function of SNR
    P_H_resp1 = 1-P_FA; % 'Hits' when the participant response was '1'
    P_H_resp2 = P_H; % Script3 assumes that the target sound is option '2'
    
    figure('Position',[100 100 700 600]);
    tiledlayout(2,2,'TileSpacing','tight');
    
    XL = [-18.5 -8.5];
    for i_masker = 1:N_maskers
        
        nexttile(1)
        
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
        text(0.05,0.92,'A. /aba/ trials (CR)',txt_options{:});
        
        nexttile(2)
        
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
        text(0.05,0.92,'B. /ada/ trials (H)',txt_options{:});
    end
    
    nexttile(3)
        
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

    nexttile(4)
        
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
    hname{end+1} = 'fig05-Rates'; 
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

if flags.do_fig6 || flags.do_fig7 || flags.do_fig8 || flags.do_fig9 || ...
   flags.do_fig3_suppl || flags.do_fig4_suppl || flags.do_fig5_suppl
    glmfct = 'l1glm';
    
    flags_base = {'gammatone', ... % gt
                  glmfct, ... % l1glm
                  'no_bias', ... % nob
                  'expvar_after_reversal',4, ... %'idx_trialselect', 1:4000 ... % tr4000
                  'pyramid_script','imgaussfilt', ...
                  'pyramid_shape',-1 ...
                  }; % the extra fields

    dir_where = dir_ACI_exp;
    % Location of the savegames:
               % /home/alejandro/Desktop/fastACI_today/fastACI_data/speechACI_Logatome-abda-S43M/S01/Results
    % dir_res = {'/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Publications/publ_osses2022b/',['1-experimental_results' filesep]};
    % dir_wav = '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/';

end

if flags.do_fig6 || flags.do_fig5_suppl || flags.do_fig5b_suppl
    
    N_noises = length(noise_str);
    flags_text = {'FontWeight','Bold','Units','Normalized'};
    
    if flags.do_fig6
        bAdd_thres = 1; % Option added on 10/10/2022
        if bAdd_thres 
            try
                flags_here = varargin(2:end); % excluding the fig number
                meanSNR = publ_osses2022b_JASA_figs('fig4','no_plot',flags_here{:});
            catch
                meanSNR = publ_osses2022b_JASA_figs('fig4','no_plot');
            end
            file_savegame = meanSNR.file_savegame;
            meanSNR = meanSNR.meanSNR; % recycling the variable
        end
    end
     
    Pos4 = [800 800];
    N_subjs2plot = [6 6];
    
    Group_ACIs = [];
    count_curr_subject = 1;
    for i_repeat = 1:2
        N_subj2plot = N_subjs2plot(i_repeat);
        if flags.do_fig6
            switch i_repeat
                case 1
                    suff = 'ac';
                case 2
                    suff = 'df';
            end
            figure('Position',[100 100 600 Pos4(i_repeat)]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
            tiledlayout(N_subj2plot,N_noises,'TileSpacing','tight');
        end
        
        for i_subj = 1:N_subj2plot    
            for i_noise = 1:N_noises
                ACI_fname = [dir_ACI_exp 'ACI-' Subjects{i_subj} '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4.mat'];
 
                ACI = [];
                cfg_ACI = [];
                results = [];
                info_toolbox = [];
                load(ACI_fname);
 
                % Copying the ACI to the group values:
                Group_ACIs{count_curr_subject,i_noise} = ACI;

                if flags.do_fig6
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
                        bColourbar = 1;
                    else
                        bColourbar = 0;
                    end
                    flags_opts = {'NfrequencyTicks',5,'colorbar',bColourbar}; 
                    flags_opts = [flags_opts flags_tf];
                
                    nexttile(N_noises*(i_subj-1)+i_noise);
                    plot_outs = affichage_tf(ACI,'CInorm', cfg_ACI.t, cfg_ACI.f,flags_opts{:}); hold on
                
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
                        il_plot_square_XY(xcoor,ycoor,'m',offx,offy,'--',4);
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
                        text(0,1.15,text2show,flags_text{:},'FontSize',14);
                    end
                    if i_noise == 3
                        text2show = Subjects{i_subj};
                        text(.7,.90,text2show,flags_text{:},'FontSize',12);
                    end
                    text2show = sprintf('%.1f',meanSNR(i_subj,i_noise));
                    text(.7,.75,text2show,flags_text{:},'FontSize',10,'Color',rgb('Gray'));

                    if bColourbar
                        hcl = plot_outs.tcolorbar;
                        set(hcl,'Ticks',[]);
                    end
                    if i_subj == 1
                        % Put 'aba'
                        if bColourbar
                            [xx,colour_max,colour_min] = Get_colourmap_rgb(colourbar_map);
                            text(1.1,1.15,'aba','FontSize',12,flags_text{:},'Color',colour_max);
                        end
                    end

                    if i_subj == N_subj2plot
                        if bColourbar
                            text(1.1,-.15,'ada','FontSize',12,flags_text{:},'Color',colour_min);
                        end
                    end

                end % end: do_fig6
            end
            count_curr_subject = count_curr_subject + 1; % incrementing the counter for the corresponding participant
        end
        Subjects(1:N_subj2plot) = []; % remove the first six labels, so that it's easier to repeat this code
        if flags.do_fig6
            meanSNR(1:N_subj2plot,:) = [];
        end
        
        if flags.do_fig6
            h(end+1) = gcf;
            hname{end+1} = ['fig06' suff '-ACI-exp'];
        end
    end
    
    if flags.do_fig6
        figure('Position',[100 100 2*600 250]);
        tiledlayout(1,N_noises,'TileSpacing','tight');
    end
    for i_noise = 1:N_noises
        if flags.do_fig6
            nexttile(i_noise);
        end
        G_ACI = Group_ACIs{1,i_noise}; % ACI from the first participant
        N_ACIs = length(Group_ACIs);
        for i_ACI = 2:N_ACIs
            G_ACI = G_ACI + Group_ACIs{i_ACI,i_noise};
        end
        G_ACI = G_ACI/N_ACIs;
        
        if flags.do_fig6
            % Using the cfg_ACI of the last participant:
            affichage_tf(G_ACI,'CI', cfg_ACI.t, cfg_ACI.f,flags_opts{:}); hold on
        
            dir_targets = [dir_subj_exp 'speech-samples' filesep]; % this is the folder of participant S01, arbitrarily chosen

            LineStyles = {'-','-.'};
            Colours = {'r','b'}; % {[0.6,0,0],[0,0,0.6]};
            LineWidth = 1.5;
            outs_aff = affichage_tf_add_Praat_metrics(dir_targets,cfg_ACI,[],LineStyles,Colours,LineWidth);
            if i_noise ~= N_maskers
                set(outs_aff.hl,'Visible','off'); % Removes the legend
            end

            title(noise_types_label{i_noise});
        
            switch i_noise
                case 1
                    text2show = 'G.';
                case 2
                    text2show = 'H.';
                case 3
                    text(1.2,1,'aba','FontSize',12,flags_text{:},'Color',colour_max);
                    text(1.2,0,'ada','FontSize',12,flags_text{:},'Color',colour_min);

                    text2show = 'I.';
            end
            text(0,1.08,text2show,'FontWeight','Bold','Units','Normalized','FontSize',14); 
            suff = 'gi';
        end % do_flag6
    end
    if flags.do_fig6
        h(end+1) = gcf;
        hname{end+1} = ['fig06' suff '-ACI-exp'];
    end
end % end do_fig6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flags.do_fig5_suppl || flags.do_fig5b_suppl
    rval = zeros([N_subjects N_subjects N_noises]);
    YL_PA = [0.05 1.05]; % range for experimental PA if correction for chance
    
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
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5_suppl
    for i_no = 1:N_noises
        for i_subj_aci = 1:N_subjects
            X = Group_ACIs{i_subj_aci,i_no};
            X = X(:);
            for i_subj_data = 1:N_subjects
                Y = Group_ACIs{i_subj_data,i_no};
                Y = Y(:);
                rval(i_subj_aci,i_subj_data,i_no) = corr(X,Y,'type',type_corr);
            end
        end
    end
    
    bDo_with_fig5b_suppl = input('Will you tile fig5_suppl right after? (1=yes; 0=no): ');
     
    if bDo_with_fig5b_suppl
        N_rows = 2;
        Pos4_here = 650;
        suff = 'suppl-fig5';
    else
        N_rows = 1;
        Pos4_here = 415;
        suff = '';
    end
    figure('Position',[100 100 720 Pos4_here]);
    tiledlayout(N_rows,3,'TileSpacing','tight'); % 'tight');
     
    idx_sig = [0 0 0]; % idx of significant ACIs
    for i_column = 1:3
        rval_here = rval(:,:,i_column);
        nexttile(i_column)
        imagesc(1:N_subjects,1:N_subjects,rval_here); hold on
        colormap('gray')
        caxis(YL_PA)
        set(gca,'XTick',XT);
        set(gca,'YTick',XT);
        if i_column ~= 1
            set(gca,'YTickLabel',[]);
        end
        if i_column == 3
            c = colorbar;
        end
    end
     
    FS_here = 10;
     
    nexttile(1);
    ylabel('ACI from');
    set(gca,'YTickLabel',YTL);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'A.',opts_txt{:});
    set(gca,'FontSize',FS_here);
     
    nexttile(2);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    text(0,1.05,'B.',opts_txt{:});
    set(gca,'FontSize',FS_here);
    
    nexttile(3);
    text(0,1.05,'C.',opts_txt{:});
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    
    % c.Label.String = 'Percentage accuracy benefit (%)';
    set(c,'Limits',YL_PA);
    fig_name = [suff '-corr_participant-' type_corr];

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
    [mean_val,sem_1p64_val, N_val] = il_get_offdiag_fig5_suppl(rval);
    
    fprintf('Type of correlation: %s\n',type_corr);
    fprintf('\tMean off diagonal for the noises: %.2f %.1f, %.2f\n',mean_val);
    fprintf('\tCI   off diagonal for the noises: %.2f %.2f, %.2f\n',sem_1p64_val);
    fprintf('\tN    off diagonal for the noises: %.1f %.1f, %.1f\n',N_val);
    %%%
end

bAll_participants = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7 || flags.do_fig8 || flags.do_fig9 || flags.do_fig3_suppl || ...
   flags.do_fig4_suppl
    % Loading cross predictions
    
    dir_res = {[dir_fastACI 'Publications' filesep 'publ_osses2022b' filesep],['1-experimental_results' filesep]};
           
    trialtype_analysis = 'total';

    if bAll_participants == 0 % bAll_participants = 0;
        warning('Temporal solution...');
        Subjects = {'S01','S02','S03','S04','S05','S07','S08','S09','S10','S11'};
        N_subjects = length(Subjects);
    end
    
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
                        fname_crosspred{ii} = [dir_ACI_exp 'ACI-' subject '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                        if exist(fname_crosspred{ii},'file')
                            bInclude(ii) = 1;
                        end
                    end
                end
            end
            %%%
            if flags.do_fig8 || flags.do_fig9 || flags.do_fig3_suppl || ...
                    flags.do_fig3b_suppl || flags.do_fig4_suppl || flags.do_fig4b_suppl 
                % 2.1. Checking whether you previously stored already the cross-prediction data:
                if flags.do_fig8 || flags.do_fig3_suppl 
                    if i_subject == 1 
                        bInclude = zeros(1,N_subjects);
                        cross_suffix = '';
                        fname_crosspred = [];
                        for ii = 1:N_subjects
                            fname_crosspred{ii} = [dir_ACI_exp 'ACI-' Subjects{ii} '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
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
                        fname_crosspred{ii} = [dir_ACI_exp 'ACI-' subject '-speechACI_Logatome-' noise_str{ii} '-nob-gt-' glmfct '+pyrga-rev4.mat'];
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
                
                if bAll_participants == 1
                    % 2.1. Checking whether you previously stored already the cross-prediction data:
                    [crosspred, fname_cross] = Check_if_crosspred(cfg_ACI_here.fnameACI,fname_crosspred,cross_suffix);
                else
                    crosspred = [];
                end

                bForceComplete = 0;

                if isempty(crosspred) || bForceComplete
                    % 2.2. If there is no cross-prediction in your computer, then
                    %      it obtains it
                    flags_ex = {'Data_matrix',Data_matrix_here}; % the extra fields
                    if bAll_participants == 1
                        flags_ex{end+1} = 'ACI_crosspred';
                        flags_ex{end+1} = fname_crosspred;
                    end
                    [col1,cfg_ACI_here,res,col4] = fastACI_getACI(fname_results, flags_base{:},flags_ex{:},flags_for_input{:});
                else
                    % (Part of 2.1): places 'crosspred' as a field of res. (local results)
                    if bAll_participants == 1
                        res.crosspred = crosspred;
                    end
                end

                %%% 3. It places the cross-prediction reults in a cell array:
                ACI{i_subject,i_masker}     = col1;
                cfg_ACI{i_subject,i_masker} = cfg_ACI_here; 
                results{i_subject,i_masker} = res;
                % idxlambda(i_subject) = res.idxlambda;

                if bAll_participants == 1
                    [cross,outs_cross] = Read_crosspred(fname_cross,res.idxlambda);
                    PC_matrix(i_subject,:,i_masker)  = outs_cross.PA_mean_re_chance;
                    Dev_matrix(i_subject,:,i_masker) = outs_cross.Dev_mean;
                    Dev_matrix_plus_SEM(i_subject,:,i_masker) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
                    Dev_matrix_min_SEM(i_subject,:,i_masker)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
                    
                    if i_masker == N_maskers
                        for ii_x = 1:N_maskers 
                            X = ACI{i_subject,ii_x};
                            X = X(:);
                            for ii_y = 1:N_maskers
                                Y = ACI{i_subject,ii_y}; Y = Y(:);
                                rval(i_subject,ii_x,ii_y) = corr(X,Y,'type',type_corr);
                            end
                        end
                    end
                end
                
                if flags.do_fig7 
                    res_here = results{i_subject,i_masker};
                    % PC and DEV for incorrect trial
                    outs_perf     = Read_ACI_performance_metrics(res_here,'test');
                    outs_perf_inc = Read_ACI_performance_metrics_inc(res_here,data_passation,cfg_ACI_here);
                    
                    % [PC_per_fold_inc,Dev_per_fold_inc] = il_tbtpred_inc(cfg_ACI_here,res_here,data_passation);

                    %%% Save some data for group-level analysis
                    idxlambda = res_here.idxlambda;
                    G_Devtrial_opt(:,i_subject,i_masker) = outs_perf.Dev;
                    G_PAtrial_opt(:,i_subject,i_masker)  = outs_perf.PA_re_chance;
                    G_OptDEV(i_subject,i_masker) = outs_perf.Dev_mean;
                    G_OptPA(i_subject,i_masker)  = outs_perf.PA_mean_re_chance;
                    
                    %%% Incorrect:
                    G_Devtrial_opt_inc(:,i_subject,i_masker) = outs_perf_inc.Dev;
                    G_PAtrial_opt_inc(:,i_subject,i_masker)  = outs_perf_inc.PA_re_chance;

                    G_OptDEV_inc(i_subject,i_masker) = outs_perf_inc.Dev_mean;
                    G_OptPA_inc(i_subject,i_masker)  = outs_perf_inc.PA_mean_re_chance;
                end % if flags.do_fig7
                
            end % end if bInclude
        end % end if i_subject
    end % end i_masker
    
end

if flags.do_fig7
    % Taken from l20220325_crossprediction_newversion.m and 
    %   g20220203_crossprediction.m
    
    %% function l20220502_analysis_pipeline_ACIstats
    figure('Position',[100 100 700*2-200 530]);
    tiledlayout(2,2,'TileSpacing','tight');
    
    factor_guess = 1/(1-.5);
    crit_sig = factor_guess*[0 0 1.3 2.39]; % significance criterion
    
    XL = [0.5 3.5]; 
    nexttile(1);
    plot(XL,crit_sig(1)*[1 1],'k--'); hold on; grid on
    
    nexttile(2);
    plot(XL,crit_sig(2)*[1 1],'k--'); hold on; grid on
    
    nexttile(3);
    plot(XL,crit_sig(3)*[1 1],'k--'); hold on; grid on
    
    nexttile(4);
    plot(XL,crit_sig(4)*[1 1],'k--'); hold on; grid on
    
    for i_masker = 1:N_maskers
        for i_subject = 1:N_subjects

            if i_subject == 9 || i_subject == 10
                disp('')
            end
            offx = 0.05*(i_subject - ((N_subjects+2)/2-.5)); %  0.05*i_subject;
            
            y2plot  = G_Devtrial_opt(:,i_subject,i_masker);
            y2plot2 = G_PAtrial_opt(:,i_subject,i_masker);
            
            nexttile(1);
            Me = mean(y2plot);
            EB = 1.64*sem(y2plot);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(1) && min(Me_SEM) <= crit_sig(1) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx, Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
            
            nexttile(3);
            Me = mean(y2plot2);
            EB = 1.64*sem(y2plot2);
            MCol = 'w';
            errorbar(i_masker+offx, Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on
            
            %%% Incorrect trials:
            y2plot  = G_Devtrial_opt_inc(:,i_subject,i_masker);
            y2plot2 = G_PAtrial_opt_inc(:,i_subject,i_masker);

            nexttile(2);
            Me = mean(y2plot);
            EB = 1.64*sem(y2plot);
            MCol = 'w';
            Me_SEM = [Me+EB Me-EB];
            if max(Me_SEM) >= crit_sig(2) && min(Me_SEM) <= crit_sig(2) % there is a zero crossing, non significant
                MCol = rgb('Gray');
            end
            errorbar(i_masker+offx, Me , EB , 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

            nexttile(4);
            Me = mean(y2plot2);
            EB = 1.64*sem(y2plot2);
            MCol = 'w';
            errorbar(i_masker+offx,Me, EB, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on            
        end

        y2plot  = G_OptDEV(:,i_masker);
        % subplot(2,1,1)
        nexttile(1);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        y2plot  = G_OptPA(:,i_masker);
        % subplot(2,1,2)
        nexttile(3);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        %%% Incorrect:
        y2plot  = G_OptDEV_inc(:,i_masker);
        nexttile(2);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        
        y2plot  = G_OptPA_inc(:,i_masker);
        nexttile(4);
        offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
        errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
    end

    XLabel = 'Masker';
    % subplot(2,1,1)
    nexttile(1);
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
    
    nexttile(2);
    title('Incorrect trials only')
    grid on
    xlim(XL);
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', []); 
    set(gca, 'YTickLabels',[]);
    text(xcoor,ycoor,'B',txt_opts{:});
    
    nexttile(3);
    grid on
    xlim(XL)
    ylim([-1 46])
    xlabel(XLabel); 
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', noise_types_label); 
    ylabel('Percent accuracy benefit (%)')
    text(xcoor,ycoor,'C',txt_opts{:});
    
    nexttile(4);
    
    xlim(XL)
    ylim([-1 46])
    xlabel(XLabel); 
    set(gca, 'XTick', 1:N_maskers); 
    set(gca, 'XTickLabels', noise_types_label); 
    set(gca, 'YTickLabels',[]);
    % ylabel('Percent accuracy benefit (%)')
    text(xcoor,ycoor,'D',txt_opts{:});
    
    h(end+1) = gcf;
    hname{end+1} = 'fig07-metrics-benefit';
    
    if bAll_participants == 0
        hname{end} = [hname{end} '-10-participants'];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if flags.do_fig8 || flags.do_fig3_suppl 
    % Code by Leo, taken from l20220325_crossprediction_newversion.m
 
    YL_PA = [-0.5 20.5]; % range for experimental PA if correction for chance
    
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
    
    bDo_with_fig8b = input('Will you tile fig8b right after? (1=yes; 0=no): ');
    
    if bDo_with_fig8b
        N_rows = 2;
        Pos4_here = 650;
        suff = '08';
    else
        N_rows = 1;
        Pos4_here = 415;
        suff = '8';
    end
    figure('Position',[100 100 720 Pos4_here]);
    tiledlayout(N_rows,3,'TileSpacing','tight'); % 'tight');

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
        Dev_diag(:,i_column) = diag(Dev_matrix_plus);
        
        if flags.do_fig3_suppl
            nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,Dev_matrix_plus); hold on
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
                    if bDo_with_fig8b
                        Y_start = 2.5;
                    else
                        Y_start = 2;
                    end
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
        % PC_matrix_nodiag = PC_matrix_here-repmat(PC_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8
            nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis(YL_PA)
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
                    if bDo_with_fig8b
                        Y_start = 2.5;
                    else
                        Y_start = 2;
                    end
                    
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
        end
    end
    PA_matrix_nan = PC_matrix;
    for ii = 1:3
        idxs = find(Dev_diag(:,ii)>=0);
        PA_matrix_nan(:,idxs) = nan;
    end
    fprintf('Delta PAexp values ranged between %.1f and %.1f\n', ...
        min(min(min(PA_matrix_nan))),max(max(max(PA_matrix_nan))));
    fprintf('ACIs leading to significant predictions using the CVD criterion: %.1f, %.1f, %.1f (white, bump, MPS)\n', ...
        idx_sig);
    
    fprintf('On-diagonal averages of: %.1f (white), %.1f (bump), %.1f (MPS) percent\n',Me_diag);
    fprintf('Off-diagonal averages of: %.1f (white), %.1f (bump), %.1f (MPS) percent\n',Me_nodiag);
    
    FS_here = 10;

    nexttile(1);
    ylabel('Data from');
    set(gca,'YTickLabel',YTL);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'A.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    nexttile(2);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    text(0,1.05,'B.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    nexttile(3);
    text(0,1.05,'C.',opts_txt{:});
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    if flags.do_fig8
        % c.Label.String = 'Percentage accuracy benefit (%)';
        set(c,'Limits',YL_PA);
        fig_name = ['fig' suff '-crosspred_participant'];
    end
    if flags.do_fig3_suppl
        % c.Label.String = 'Deviance / trial benefit (adim)';
        fig_name = 'suppl-fig3-crosspred_participant';
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig9 || flags.do_fig4_suppl
    
    YL_PA = [-0.5 20.5]; 
    
    if flags.do_fig9
        bDo_with_fig9b = input('Will you tile fig9b right after? (1=yes; 0=no): ');

        if bDo_with_fig9b
            N_rows = 1;
            N_cols = 2;
            Pos3_here = 580;
            Pos4_here = 250;
            suff = '09';
        else
            N_rows = 1;
            N_cols = 1;
            Pos4_here = 350;
            Pos4_here = 250;
            suff = '9b';
        end
    end
        
    for i_subj = 1:N_subjects % Participant number
        for i = 1:N_maskers % On-noise masker
            % res_here = results{i_subject,i_masker};
            % % PC and DEV for incorrect trial
            % outs_perf     = Read_ACI_performance_metrics(res_here,'test');

            idxlambda = results{i_subj,i}.idxlambda;
            res_here = results{i_subj,i};
            outs_perf = [];
            % outs_perf = Read_ACI_performance_metrics(res_here,'test');

            [xx,outs] = Read_crosspred(res_here, idxlambda);
            PA_matrix(i_subj,i,:) = outs.PA_mean_re_chance;
            % PA_matrix(i_subj,i,:) = outs.PA_mean_re_chance_no_corr;
            Dev_matrix(i_subj,i,:) = outs.Dev_mean;
            Dev_matrix_plus_SEM(i_subj,i,:) = outs.Dev_mean + outs.Dev_SEM;
        end
    end
        
    i_subj = 1:4;
    PA  = il_get_3x3(PA_matrix ,i_subj);
    Dev = il_get_3x3(Dev_matrix,i_subj);
    Dev_SEM = il_get_3x3(Dev_matrix_plus_SEM,i_subj);
    
    i_subj = 5:8;
    PA  = [PA ; il_get_3x3(PA_matrix ,i_subj)];
    Dev = [Dev; il_get_3x3(Dev_matrix,i_subj)];
    Dev_SEM = [Dev_SEM; il_get_3x3(Dev_matrix_plus_SEM,i_subj)];
    
    i_subj = 9:N_subjects;
    PA  = [PA ; il_get_3x3(PA_matrix ,i_subj)];
    Dev = [Dev; il_get_3x3(Dev_matrix,i_subj)];
    Dev_SEM = [Dev_SEM; il_get_3x3(Dev_matrix_plus_SEM,i_subj)];
    %%%
    
    if flags.do_fig4_suppl
        figure('Position',[100 100 500 800]); % 560 420]);
        tiledlayout(2,1,'TileSpacing','compact'); % 'tight');

        nexttile(1);
        x_var = 1:size(PA,2);
        y_var = 1:size(PA,1);
        imagesc(x_var,y_var,PA); hold on;
        colormap('gray')
        caxis(YL_PA)

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

        % PC_pool_diag_removed = transpose(PC_pool_diag_removed);
        offxy = 0.5;
        offxy = offxy*.9;
        for jj = 1:size(Dev_SEM,1)
            idx_here = find(Dev_SEM(jj,:)<0);
            for ii = 1:length(idx_here)
                il_plot_square(idx_here(ii),jj,'m',offxy,'--');
            end
        end

        c = colorbar;
        c.Label.String = 'PA benefit (%)';

        %%% Adding the Subject_ID to the panels
        for ii = 1:4
            text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 5:8
            text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 9:12
            text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        %%%
        
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
        hname{end+1} = 'suppl-fig4-crosspred_masker';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig9
        figure('Position',[100 100 Pos3_here Pos4_here]); 
        if bDo_with_fig9b
            tiledlayout(N_rows,N_cols,'TileSpacing','compact');
            nexttile(1);
        end

        bAll_participants = 1; %  bAll_participants = 0; % PA_sum
        if bAll_participants == 0
            idxs = [1 2 3 4 5 7 8 9 10 11];
            N_subjects = length(idxs);

            PC_matrix = PC_matrix(idxs,:,:);
            Dev_matrix = Dev_matrix(idxs,:,:);
            Dev_matrix_plus_SEM = Dev_matrix_plus_SEM(idxs,:,:);
        end

        PA_sum = zeros([3 3]);
        row_i = 1;
        row_f = 3;
        for i_subj = 1:4
            idxi = N_maskers*(i_subj-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end

        row_i = 4;
        row_f = 6;
        for i_subj = 5:8
            subj_off = 4;
            idxi = N_maskers*(i_subj-subj_off-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end

        row_i = 7;
        row_f = 9;
        for i_subj = 9:N_subjects
            subj_off = 8;
            idxi = N_maskers*(i_subj-subj_off-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end
        PA_sum = PA_sum/N_subjects;

        x_var = 1:size(PA_sum,2);
        y_var = 1:size(PA_sum,1);
        imagesc(x_var,y_var,PA_sum); hold on
        colormap('gray')
        caxis(YL_PA)
        set(gca,'FontSize',10);

        c = colorbar;
        if bDo_with_fig9b == 0
            c.Label.String = 'PA benefit (%)';
        end

        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.1f',PA_sum(jj,ii)),'FontSize',10,'Color','m');
            end
        end

        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)

        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);

        xlabel('ACI from condition');
        ylabel('Data from condition');

        if bAll_participants == 1
            text4label = 'A.';
        else
            text4label = 'G.';
        end
        text(-0.14,0.94,text4label,'Units','Normalized','FontSize',14,'FontWeight','Bold');
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xlabel('ACI from condition');
        ylabel('Data from condition');

        h(end+1) = gcf;
        hname{end+1} = ['fig' suff '-crosspred_masker'];
        
        if bAll_participants == 0
            hname{end} = [hname{end} '-10-participants'];
        end

        %%% Text
        fprintf('[...] the main diagonals gave overall Delta_PA values of %.1f, %.1f, and %.1f percent\n', ...
            PA_sum(1,1),PA_sum(2,2),PA_sum(3,3));
        fprintf('for the white, bump, and MPS noises (the same group values as in Fig. 7C)\n');
        
        diffe = PA_sum - repmat( transpose(diag(PA_sum)),3,1 );
        PA_sum_nodiag = PA_sum; 
        for ii = 1:3
            diffe(ii,ii) = nan; % setting the diagonal to NaN
            PA_sum_nodiag(ii,ii) = nan;
        end
        fprintf('that decreased (at least) by %.1f, %.1f, %.1f percent\n',-max(diffe));
        fprintf('that decreased (at most) by %.1f, %.1f, %.1f percent\n' ,-min(diffe));
        
        fprintf('Vertical direction: that decreased (at most) to %.1f, %.1f, %.1f percent\n',min(PA_sum_nodiag));
        
        %%% Finally I kept this one in the text:
        fprintf('Horizontal direction: that decreased (at most) to %.1f, %.1f, %.1f percent\n',min(PA_sum_nodiag,[],2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rval_me = squeeze(nanmean(rval));
        
        figure('Position',[100 100 Pos3_here Pos4_here]); 
        if bDo_with_fig9b
            tiledlayout(N_rows,N_cols,'TileSpacing','compact');
            nexttile(1);
        end
        imagesc(x_var,y_var,rval_me); hold on
        colormap('gray')
        caxis([-0.05 1.05])
        set(gca,'FontSize',10);
 
        c = colorbar;
 
        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.2f',rval_me(jj,ii)),'FontSize',10,'Color','m');
            end
        end

        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
 
        XTL = {'W','BP','MPS'};
        YTL = {'W','BP','MPS'};
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
 
        xlabel('ACI from condition');
        ylabel('ACI from condition');

        text4label = 'G.';
        text(-0.14,0.94,text4label,'Units','Normalized','FontSize',14,'FontWeight','Bold');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        h(end+1) = gcf;
        hname{end+1} = 'suppl-fig5b-cross_masker';
        
        if bAll_participants == 0
            hname{end} = [hname{end} '-10-participants'];
        end

    end
end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figs 10 and 11
if flags.do_fig10       || flags.do_fig8b || flags.do_fig9b || ...
   flags.do_fig3_suppl  || flags.do_fig3b_suppl|| flags.do_fig4b_suppl || ...
   flags.do_fig5b_suppl || flags.do_fig6_suppl
    
    % All simulation results were stored in one output folder
    if isunix
        dir_results = dir_ACI_sim;
    else
        error('Leo put your folder here...')
        dir_results = [];
    end
    Subjects_ID = {'Results-S01-v1','Results-S02-v1','Results-S03-v1','Results-S04-v1', ...
                   'Results-S05-v1','Results-S06-v1','Results-S07-v1','Results-S08-v1', ...
                   'Results-S09-v1','Results-S10-v1','Results-S11-v1','Results-S12-v1'};
    model_str = 'osses2022a';
end

if flags.do_fig9b || flags.do_fig10 || flags.do_fig5b_suppl || flags.do_fig6_suppl
    
    flags_text = {'FontWeight','Bold','Units','Normalized'};
    N_noises = length(noise_str);
    
    N_subjs2plot = [3 9];
    if flags.do_fig10 || flags.do_fig6_suppl
        Pos4 = [400 1200];
        bAdd_thres = 1; % Option added on 12/10/2022
    end
    
    if flags.do_fig10
        pref = 'fig10';
        idxs_repeat = 1;
    end
    if flags.do_fig6_suppl
        pref = 'suppl-fig6';
        idxs_repeat = [1 2]; % does all, but only plots '2'
    end
    if flags.do_fig5b_suppl || flags.do_fig9b
        bAdd_thres = 0;
        idxs_repeat = [1 2];
    end
    
    i_curr = 1;
    All_ACI_sim = [];
    for i_repeat = idxs_repeat
        N_subj2plot = N_subjs2plot(i_repeat);
        
        if (i_repeat == 1 && flags.do_fig10) || (i_repeat == 2 && flags.do_fig6_suppl)
            figure('Position',[100 100 600 Pos4(i_repeat)]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
            tiledlayout(N_subj2plot,N_noises,'TileSpacing','tight');
            bPlot = 1;
        else
            bPlot = 0;
        end
        
        for i_subj = 1:N_subj2plot    
            dir_subj = [dir_results Subjects_ID{i_subj} filesep];
            for i_noise = 1:N_noises
                
                if bAdd_thres 
                    dir_subj_save = [dir_savegame_sim Subjects_ID{i_subj} filesep];
                    file_savegame = Get_filenames(dir_subj_save,['savegame*' noise_str{i_noise} '*.mat']);
                    if length(file_savegame)==1
                        % Only one file should be retrieved
                        outs_SNR = il_get_SNR([dir_subj_save file_savegame{1}]);
                    else
                        error('Only 1 file_savegame should have been found')
                    end
                    meanSNR_sim(i_subj,i_noise) = mean(outs_SNR.SNR_me);
                    trial_nr(i_subj,i_noise) = length(outs_SNR.idxs2use_all);
                    % meanSNR = publ_osses2022b_JASA_figs('fig4','no_plot');
                    % meanSNR = meanSNR.meanSNR; % recycling the variable
                end
                
                ACI_fname = [dir_subj 'ACI-' model_str '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4.mat'];

                ACI = [];
                cfg_ACI = [];
                results = [];
                info_toolbox = [];
                load(ACI_fname);

                All_ACI_sim{i_curr,i_noise} = ACI; 
                
                if bPlot
                    if i_noise == 3
                        bColourbar = 'yes';
                    else
                        bColourbar = 'no';
                    end
                    flags_opts = {'NfrequencyTicks',5,'colorbar',bColourbar}; 
                    flags_opts = [flags_opts flags_tf];

                    nexttile(N_noises*(i_subj-1)+i_noise);
                    plot_outs = affichage_tf(ACI,'CInorm', cfg_ACI.t, cfg_ACI.f,flags_opts{:});
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
                        text(0,1.15,text2show,'FontSize',14, flags_text{:});
                    end
                    if i_noise == 3
                        Subj_ID = strsplit(Subjects_ID{i_subj},'-');
                        text2show = Subj_ID{2};
                        text(.7,.90,text2show,'FontSize',12, flags_text{:});
                        
                        if bColourbar
                            hcl = plot_outs.tcolorbar;
                            set(hcl,'Ticks',[]);
                        end
                        
                        if i_subj == 1
                            % Put 'aba'
                            if bColourbar
                                [xx,colour_max,colour_min] = Get_colourmap_rgb(colourbar_map);
                                text(1,1.15,'aba','FontSize',12,flags_text{:},'Color',colour_max);
                            end
                        end
                        
                        if i_subj == N_subj2plot
                            disp('')
                            if bColourbar
                                text(1,-.15,'ada','FontSize',12,flags_text{:},'Color',colour_min);
                            end
                        end
                        
                    end
                    text2show = sprintf('%.1f',meanSNR_sim(i_subj,i_noise));
                    text(.7,.75,text2show,flags_text{:},'FontSize',10,'Color',rgb('Gray'));
                end % end bPlot
                disp('')
            end
            i_curr = i_curr+1;
        end
        
        Subjects_ID(1:N_subj2plot) = []; % remove the first six labels, so that it's easier to repeat this code
        if flags.do_fig10 || flags.do_fig6_suppl
            if i_repeat == 1
                meanSNR_sim_all = [];
            end
            meanSNR_sim_all = [meanSNR_sim_all; meanSNR_sim];
            meanSNR_sim(1:N_subj2plot,:) = [];
        end
        
        if bPlot
            h(end+1) = gcf;
            hname{end+1} = [pref '-ACI-sim'];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5b_suppl % correlation simulations
    for i_no = 1:N_noises
        for i_subj_aci = 1:N_subjects
            X = All_ACI_sim{i_subj_aci,i_no};
            X = X(:);
            for i_subj_data = 1:N_subjects
                Y = All_ACI_sim{i_subj_data,i_no};
                Y = Y(:);
                rval(i_subj_aci,i_subj_data,i_no) = corr(X,Y,'type',type_corr);
            end
        end
    end
    
    for i_column = 1:3
        rval_here = rval(:,:,i_column);
        nexttile(i_column+3)
        imagesc(1:N_subjects,1:N_subjects,rval_here); hold on
        colormap('gray')
        caxis(YL_PA)
        set(gca,'XTick',XT);
        set(gca,'YTick',XT);
        if i_column ~= 1
            set(gca,'YTickLabel',[]);
        end
        if i_column == 3
            c = colorbar;
        end
    end
    
    FS_here = 10;
     
    nexttile(4);
    ylabel('ACISI from');
    set(gca,'YTickLabel',YTL);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'D.',opts_txt{:});
    set(gca,'FontSize',FS_here);
     
    nexttile(5);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
    text(0,1.05,'E.',opts_txt{:});
    set(gca,'FontSize',FS_here);
    
    nexttile(6);
    text(0,1.05,'F.',opts_txt{:});
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    
    % c.Label.String = 'Percentage accuracy benefit (%)';
    set(c,'Limits',YL_PA);
    
    set(gca,'FontSize',FS_here);
     
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',XTL);
    set(gca,'FontSize',FS_here);
     
    for i_masker = 1:N_maskers
        nexttile(i_masker+3);
        ht = title(noise_types_label{i_masker});
        set(ht,'FontSize',FS_here);
        xlabel('ACISI from');
    end
    
    %%%
    [mean_val,sem_1p64_val, N_val] = il_get_offdiag_fig5_suppl(rval);
    
    fprintf('Type of correlation: %s\n',type_corr);
    fprintf('\tMean off diagonal for the noises: %.2f %.1f, %.2f\n',mean_val);
    fprintf('\tCI   off diagonal for the noises: %.2f %.2f, %.2f\n',sem_1p64_val);
    fprintf('\tN    off diagonal for the noises: %.1f %.1f, %.1f\n',N_val);
    %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b || flags.do_fig9b || ... flags.do_fig3_suppl|| 
   flags.do_fig3b_suppl || flags.do_fig4b_suppl

    if flags.do_fig9b
        Subjects_ID = {'Results-S01-v1','Results-S02-v1','Results-S03-v1','Results-S04-v1', ...
                       'Results-S05-v1','Results-S06-v1','Results-S07-v1','Results-S08-v1', ...
                       'Results-S09-v1','Results-S10-v1','Results-S11-v1','Results-S12-v1'};
        warning('Temporal')
    end
    N_noises = length(noise_str);
    N_subjects = length(Subjects_ID);
    
    for i_noise = 1:N_noises
        PA_mean_chance = []; % emptied after having plotted every noise
        for i_subj = 1:N_subjects
            if i_noise == 1
                Subj_ID = strsplit(Subjects_ID{i_subj},'-');
                labels2use{i_subj} = Subj_ID{2};
            end
        
            dir_subj = [dir_results Subjects_ID{i_subj} filesep];
        
            ACI_fname_folder = [dir_subj 'ACI-' model_str '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4' filesep];
            if flags.do_fig9b == 1 || flags.do_fig4b_suppl == 1
                Crossfile = [ACI_fname_folder 'Crosspred-noise.mat'];
            else
                Crossfile = [ACI_fname_folder 'Crosspred.mat'];
            end
            
            cfg_ACI_here = load([ACI_fname_folder(1:end-1) '.mat']);
            res_here = cfg_ACI_here.results;
            %%% Save some data for group-level analysis
            idxlambda = res_here.idxlambda;
            
            [cross,outs_cross] = Read_crosspred(Crossfile, idxlambda);
            PC_matrix(i_subj,i_noise,:)  = outs_cross.PA_mean_re_chance;
            Dev_matrix(i_subj,i_noise,:) = outs_cross.Dev_mean;
            Dev_matrix_plus_SEM(i_subj,i_noise,:) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
            % Dev_matrix_min_SEM(i_subj,i_noise,:)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
            
            if flags.do_fig9b
                for ii_x = 1:N_maskers 
                    X = All_ACI_sim{i_subj,ii_x}; 
                    X = X(:);
                    for ii_y = 1:N_maskers
                        Y = All_ACI_sim{i_subj,ii_y}; Y = Y(:);
                        rval(i_subj,ii_x,ii_y) = corr(X,Y,'type',type_corr);
                    end
                end
            end
            
            %%%
            res = [];
            res.idxlambda = idxlambda;
            res.crosspred = cross;
            if flags.do_fig9b
                %%% Nothing to do
            else
                results{i_subj,i_noise} = res;
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b || flags.do_fig3b_suppl
    %%% 4. Plotting the results:
	% 3 matrices of cross-prediction accuracy
    % xvar = results{1,1}.crosspred(1).lambdas;
    
    YL_PA = [29 51]; % [14 24.5];
    
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

    bDo_with_fig8a = input('Do you want to tile this figure to fig8a? (1=yes; 0=no): ');
    
    if bDo_with_fig8a
        % No figure is created, the current figure is assumed to be fig8a
        offset_column = 3;
    else
        N_rows = 1;
        Pos4_here = 415;
        figure('Position',[100 100 720 Pos4_here]);
        tiledlayout(N_rows,3,'TileSpacing','tight');
        offset_column = 0;
    end

    for i_column = 1:3
        switch mod(i_column,3)
            case 1 % white
                Colour_here = masker_colours{1};
            case 2 % bump
                Colour_here = masker_colours{2};
            case 0 % MPS
                Colour_here = masker_colours{3};
        end
        Dev_matrix_here = squeeze(Dev_matrix_plus_SEM(:,i_column,:));
        Dev_diag(:,i_column) = diag(Dev_matrix_here);

        nexttile(i_column+offset_column)
        if flags.do_fig3b_suppl
            imagesc(1:N_subjects,1:N_subjects,Dev_matrix_here); hold on
            colormap('gray')
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            idx_total = 0;
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus_SEM(:,i_column,ii)<0);
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
        
        PC_matrix_here = squeeze(PC_matrix(:,i_column,:));
        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8b
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis(YL_PA)
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            idx_total = 0;
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus_SEM(:,i_column,ii)<0);
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

    nexttile(1+offset_column);
    if bDo_with_fig8a
        title(' ');
    end
    set(gca,'XTickLabel',XTL);
    ylabel('Simulation data from');
    set(gca,'YTickLabel',YTL);
    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    text(0,1.05,'D.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    nexttile(2+offset_column);
    if bDo_with_fig8a
        title(' ');
    end
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',[]);
    text(0,1.05,'E.',opts_txt{:});
    set(gca,'FontSize',FS_here);

    nexttile(3+offset_column);
    if bDo_with_fig8a
        title(' ');
    end
    text(0,1.05,'F.',opts_txt{:});
    % c = colorbar;
    % c.Label.String = 'Deviance / trial benefit (adim)';
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',[]);
    % set(c,'Limits',[-10.5 0])
    set(gca,'FontSize',FS_here);

    c = colorbar;

    for i_masker = 1:N_maskers
        nexttile(i_masker+offset_column);
        ht = title(noise_types_label{i_masker});
        set(ht,'FontSize',FS_here);
        xlabel('ACISI from');
    end
    if ~bDo_with_fig8a
        h(end+1) = gcf;
        if flags.do_fig8b
            hname{end+1} = 'fig8b-crosspred_model_participant'; 
        end
        if flags.do_fig3b_suppl
            hname{end+1} = 'suppl-fig3b-crosspred_model_participant'; 
        end
    end
    
    fprintf('White noise: PA between %.2f and %.1f\n',min(min(PC_matrix(:,:,1))), ...
        max(max(PC_matrix(:,:,1))));
    fprintf('Bump noise: PA between %.2f and %.1f\n',min(min(PC_matrix(:,:,2))), ...
        max(max(PC_matrix(:,:,2))));
    fprintf('MPS noise: PA between %.2f and %.1f\n',min(min(PC_matrix(:,:,3))), ...
        max(max(PC_matrix(:,:,3))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig9b || flags.do_fig4b_suppl % between noise
    YL_PA = [29 51]; % [14 24.5];

    bAll_participants = 1; %  bAll_participants = 0; % PA_sum
    if bAll_participants == 0
        idxs = [1 2 3 4 5 7 8 9 10 11];
        N_subjects = length(idxs);
        
        PC_matrix = PC_matrix(idxs,:,:);
        Dev_matrix = Dev_matrix(idxs,:,:);
        Dev_matrix_plus_SEM = Dev_matrix_plus_SEM(idxs,:,:);
    end
    
    PA = [];
    i_subj = 1:4;
    PA  = il_get_3x3(PC_matrix ,i_subj);
    Dev = il_get_3x3(Dev_matrix,i_subj);
    Dev_SEM = il_get_3x3(Dev_matrix_plus_SEM,i_subj);
    
    i_subj = 5:8;
    PA  = [PA ; il_get_3x3(PC_matrix ,i_subj)];
    Dev = [Dev; il_get_3x3(Dev_matrix,i_subj)];
    Dev_SEM = [Dev_SEM; il_get_3x3(Dev_matrix_plus_SEM,i_subj)];
    
    i_subj = 9:N_subjects;
    PA  = [PA ; il_get_3x3(PC_matrix ,i_subj)];
    Dev = [Dev; il_get_3x3(Dev_matrix,i_subj)];
    Dev_SEM = [Dev_SEM; il_get_3x3(Dev_matrix_plus_SEM,i_subj)];

    if flags.do_fig9b
        bDo_with_fig9a = input('Are you adding fig9b to fig9a? (1=yes; 0=no): ');

        if bDo_with_fig9a
            figure(1); warning('figure(1), temporal')
            nexttile(2);
        else
            figure('Position',[100 100 560 250]); % 560 420]);
        end

        PA_sum = zeros([3 3]);
        row_i = 1;
        row_f = 3;
        for i_subj = 1:4
            idxi = N_maskers*(i_subj-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end
    
        row_i = 4;
        row_f = 6;
        for i_subj = 5:8
            subj_off = 4;
            idxi = N_maskers*(i_subj-subj_off-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end

        row_i = 7;
        row_f = 9;
        for i_subj = 9:N_subjects
            subj_off = 8;
            idxi = N_maskers*(i_subj-subj_off-1)+1;
            idxf = idxi+N_maskers-1;
            PA_sum = PA_sum + PA(row_i:row_f,idxi:idxf);
        end
        PA_sum = PA_sum/N_subjects;
    
        x_var = 1:size(PA_sum,2);
        y_var = 1:size(PA_sum,1);
        imagesc(x_var,y_var,PA_sum); hold on
        colormap('gray')
        caxis(YL_PA)
        set(gca,'FontSize',10);

        xlabel('ACISI from condition');
        % ylabel(''); % ylabel('Data from condition');
        
        caxis(YL_PA)
        set(gca,'FontSize',10);

        c = colorbar;
        c.Label.String = 'PA benefit (%)';
        
        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.1f',PA_sum(jj,ii)),'FontSize',10,'Color','m');
            end
        end
    
        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)

        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);

        if bAll_participants == 1
            text4label = 'B.';
        else
            text4label = 'H.';
        end
        text(-0.14,0.94,text4label,'Units','Normalized','FontSize',14,'FontWeight','Bold');
        
        %%% Some text:
        diffe = PA_sum - repmat( transpose(diag(PA_sum)),3,1 );
        PA_sum_nodiag = PA_sum; 
        for ii = 1:3
            diffe(ii,ii) = nan; % setting the diagonal to NaN
            PA_sum_nodiag(ii,ii) = nan;
        end
        fprintf('The auto predictions were by %.1f, %.1f, %.1f percent\n',diag(PA_sum));
        fprintf('that decreased (at least) by %.1f, %.1f, %.1f percent\n',-max(diffe));
        fprintf('that decreased (at most) by %.1f, %.1f, %.1f percent\n' ,-min(diffe));
        fprintf('Vertical direction: that decreased (at most) to %.1f, %.1f, %.1f percent\n',min(PA_sum_nodiag));
        %%% Finally I kept this one in the text:
        fprintf('Horizontal direction: that decreased (at most) to %.1f, %.1f, %.1f percent\n',min(PA_sum_nodiag,[],2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rval_me = squeeze(nanmean(rval));
        
        figure(2); warning('Temporal figure(2)') 
        nexttile(2);

        imagesc(x_var,y_var,rval_me); hold on
        colormap('gray')
        caxis([-0.05 1.05])
        set(gca,'FontSize',10);
  
        c = colorbar;
  
        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.2f',rval_me(jj,ii)),'FontSize',10,'Color','m');
            end
        end

        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
 
        XTL = {'W','BP','MPS'};
        YTL = {'W','BP','MPS'};
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
 
        xlabel('ACISI from condition');
        ylabel('ACISI from condition');

        text4label = 'H.';
        text(-0.14,0.94,text4label,'Units','Normalized','FontSize',14,'FontWeight','Bold');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        h(end+1) = gcf;
        hname{end+1} = 'suppl-fig5b-cross_masker';
        
        if bAll_participants == 0
            hname{end} = [hname{end} '-10-participants'];
        end

    end % end flags.do_fig9b
    
    if flags.do_fig4b_suppl

        bDo_with_fig = input('Are you adding fig4b_suppl to fig4a_suppl? (1=yes; 0=no): ');

        if bDo_with_fig
            nexttile(2);
        else
            figure('Position',[100 100 560 250]); % 560 420]);
        end            

        x_var = 1:size(PA,2);
        y_var = 1:size(PA,1);
        imagesc(x_var,y_var,PA); hold on
        colormap('gray')
        caxis(YL_PA)
                
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
        
        %%% Adding the Subject_ID to the panels
        for ii = 1:4
            text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 5:8
            text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 9:12
            text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        %%%

        % PC_pool_diag_removed = transpose(PC_pool_diag_removed);
        offxy = 0.5;
        offxy = offxy*.9;
        for jj = 1:size(Dev_SEM,1)
            idx_here = find(Dev_SEM(jj,:)<0);
            for ii = 1:length(idx_here)
                il_plot_square(idx_here(ii),jj,'m',offxy,'--');
            end
        end        
        c = colorbar;
        c.Label.String = 'PA benefit (%)';

        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)

        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);

        xlabel('ACISI from condition');
        ylabel('Data from condition');

        text(0,1.05,'B.','Units','Normalized','FontSize',14,'FontWeight','Bold');
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
        dir_res = [dir_fastACI 'Publications' filesep 'publ_osses2022b' filesep 'data_' ...
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
    tiledlayout(2,1,'TileSpacing','tight'); % 'tight'); % none'); % 'Compact');
    
    %--- Fig. 1A ----------------------------------------------------------
    nexttile(1)
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
    
    nexttile(2)
    
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
        dir_res = [dir_fastACI 'Publications' filesep 'publ_osses2022b' filesep 'data_' ...
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
if flags.do_fig7_suppl
    prefix = 'suppl-fig7-';
    
    experiment_full = 'speechACI_Logatome-abda-S43M';
    experiment = strsplit(experiment_full,'-'); experiment = experiment{1};
    noise_here = noise_str{3}; % MPS noise
    subj_ID = 'S01';
    % dir_ACI = [fastACI_paths('dir_data') experiment_full filesep subj_ID filesep 'Results' filesep 'Results_ACI' filesep];
    
    if ~exist(dir_ACI_exp,'dir')
        % Alejandro's local computer, all subjects are there:
        dir_ACI_exp = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc_Leo/';
    end
    if ~exist(dir_ACI_exp,'dir')
        error('None of the two ACI folder locations seem to contain the ACI of the participant...');
    end
    fname = ['ACI-' subj_ID '-' experiment '-' noise_here '-nob-gt-l1glm+pyrga-rev4.mat'];
    var = load([dir_ACI_exp fname]);
    ACIs = var.results.ACIs;
    idxs = [1 8 14 20];
    
    figure('Position',[100 100 500 280]);
    tiledlayout(1,length(idxs),'TileSpacing','tight');
    
    NfrequencyTicks = 8;
    flags_extra = {'NfrequencyTicks',NfrequencyTicks,'colorbar','no'};
    flags_extra = [flags_extra flags_tf];
    t_idx = 10:70;
    for i = 1:length(idxs)
        nexttile(i);
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
% function tiledlayout(N,M,TileSpacing,TileSpacing_option)
% 
% if nargin < 4
%     TileSpacing_option = 'Compact';
% end
% bExist = exist('tiledlayout','file'); % tiledlayout.p
% if bExist
%     tiledlayout(N,M,TileSpacing,TileSpacing_option);
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function nexttile(N)
% 
% bExist = exist('nexttile','file'); % nexttile
% if bExist
%     if nargin == 0
%         nexttile;
%     else
%         nexttile(N);
%     end
% else
%     if N == 1
%         close;
%     end
%     figure;
% end

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [PC_fold, Dev_fold] = il_tbtpred_inc(cfg_ACI,results,data_passation)
% % Adapted version by Alejandro
% 
% N_lambdas = length(results.lambdas);
% N_folds = cfg_ACI.N_folds;
% 
% N = length(results.FitInfo.PC_test_t);
% PC_all = nan([N_lambdas, N_folds, N]);
% Dev_all = nan(size(PC_all));
% 
% if nargin >= 3
%     idxs_inc  = find(data_passation.is_correct(cfg_ACI.idx_analysis)==0);
%     idxs_corr = find(data_passation.is_correct(cfg_ACI.idx_analysis)==1); % correct indexes
% else
%     fprintf('All trials will be returned...')
% end
%                 
% for i_lambda = 1:N_lambdas
%     for i_fold = 1:N_folds
%         PC_test_t  = squeeze(results.FitInfo.PC_test_t(i_lambda,i_fold,:));
%         Dev_test_t = squeeze(results.FitInfo.Dev_test_t(i_lambda,i_fold,:));
%         CVtest = find(results.FitInfo.CV.test(i_fold));
%         
%         if nargin >= 3
%             [idx_corr_this_fold,idx_corr_sort] = intersect(CVtest,idxs_corr);
%         end
%         % idx_inc_this_fold  = intersect(CVtest,idxs_inc);
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
%         
%         PC_all(i_lambda,i_fold,1:N_here)  = PC_test_t;
%         Dev_all(i_lambda,i_fold,1:N_here) = Dev_test_t;
%         
%         if nargin >= 3
%             % setting back to NaN those trials within the fold that were correct:
%             PC_all(i_lambda,i_fold,idx_corr_sort)  = nan; 
%             Dev_all(i_lambda,i_fold,idx_corr_sort) = nan;
%         end
%     end
% end
% 
% PC_fold = nanmean(PC_all,3); % across third dimension
% Dev_fold = nanmean(Dev_all,3);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = il_get_SNR(file_savegame,cfg_game)

% data_passation = [];
% cfg_game = [];
load(file_savegame,'data_passation');

n_responses = data_passation.n_responses;
N_trials    = length(n_responses); % completed number of trials
% n_signal    = data_passation.n_targets(1:N_trials);
SNR         = data_passation.expvar;

resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
N_sessions  = length(resume_trial)-1;

% Fig. 2: Reversal analysis ---------------------------------------

if N_sessions > 10
    if isfield(cfg_game,'Subject_ID')
        subject = cfg_game.Subject_ID;
    else
        subject = '';
    end
    
    % In this case, the participants might have taken one or more pauses...
    idx_odd = find(mod(resume_trial(1:end-1)-1,400)~=0);
    fprintf('Participant %s has %.0f shorter sessions\n',subject,length(idx_odd));
    fprintf('\tI am going to merge that session with the previous one...\n',subject,length(idx_odd));
    for i_odd = 1:length(idx_odd)
        fprintf('\t(session %.0f started at trial %.0f)\n',idx_odd(i_odd),resume_trial(idx_odd));
    end
    resume_trial(idx_odd) = [];
    N_sessions = length(resume_trial)-1;
end
            
idxs2use_all = [];
for i_session = 1:N_sessions

    %%% Leo's way: Analysis only using the reversals: -------------
    % sessions are redefined as blocks of 400 trials
    idxi = 1+400*(i_session-1);
    idxf = 400*(i_session);
    [rev, idx] = Get_mAFC_reversals(SNR(idxi:idxf));
                
    % % Then uses median:
    % medianSNR_old(i_subject, i_masker, i_session) = median(rev);
    % percSNR_L_old(i_subject, i_masker, i_session) = prctile(rev,5);
    % percSNR_U_old(i_subject, i_masker, i_session) = prctile(rev,95);
    iscorr = data_passation.is_correct(idx);
    % perc_corr_old(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);

    %%% Alejandro's way:
    idxi = resume_trial(i_session);
    idxf = resume_trial(i_session+1)-1;

    idxi_rev = idx(4); % reversal number 4 (including it) and later
    idxi100 = idxf-100; % reversal number 4 (including it) and later

    idxs2use = idxi+idxi_rev:idxf;
    idxs2use_all = [idxs2use_all idxs2use]; % collating the idxs of the kept trials
    SNR_here = SNR(idxs2use);
    SNR_me(i_session)    = median(SNR_here);
    SNR_percL(i_session) = prctile(SNR_here,5);
    SNR_percU(i_session) = prctile(SNR_here,95);

    % iscorr = data_passation.is_correct(idxs2use);
    % perc_corr(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi:idxf);
    % perc_corr_all(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi_rev:idxf);
    % perc_corr_rev(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi100:idxf);
    % perc_corr_100(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % disp('')
end
outs = [];
outs.SNR_me = SNR_me;
outs.SNR_percL = SNR_percL;
outs.SNR_percU = SNR_percU;

outs.idxs2use_all = idxs2use_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PA_this_subj = il_get_3x3(PA_matrix,idxs)
% Hard coded... each row is one noise, each column is the cross validated noise

PA_this_subj = [];
for i_subj = idxs
    PA_tmp = squeeze(PA_matrix(i_subj,:,:));

    PA_this_subj(:,end+1) = transpose(PA_tmp(1,:));
    PA_this_subj(:,end+1) = transpose(PA_tmp(2,:));
    PA_this_subj(:,end+1) = transpose(PA_tmp(3,:));
end

if length(idxs) < 4
    PA_this_subj = [PA_this_subj nan([size(PA_this_subj,1) [4-length(idxs)]*3])];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mean_val,sem_val,N] = il_get_offdiag_fig5_suppl(rval)

N_noises = size(rval,3);
N = size(rval,1);
for i_no = 1:N_noises
    rval_here = squeeze(rval(:,:,i_no));
    for ii = 1:N
        rval_here(ii,ii) = NaN;
    end
    % rval_here = rval_here(:);
    rval_here = tril(rval_here);
    rval_here( find(rval_here == 0) ) = nan; % setting upper triangle to NaN
    
    rval_here = rval_here(:);
    rval_here( find(isnan(rval_here)) ) = [];
    mean_val(i_no) = nanmean(rval_here);
    sem_val(i_no) = 1.64*sem(rval_here);
    N(i_no) = length(rval_here);
end