function data = publ_osses2023a_JASA_EL_figs(varargin)
% function data = publ_osses2023a_JASA_EL_figs(varargin)
%
% 1. Description: Generates the figures
%
% % To display Fig. 1 of Osses et al. (2023, kerAMI) use :::
%     publ_osses2023a_JASA_EL_figs('fig1'); % Proportion of responses
%     publ_osses2023a_JASA_EL_figs('fig1','zenodo','dir_zenodo',dir_zenodo);
%
% % To display Fig. 2 of Osses et al. (2023, kerAMI) use :::
%     publ_osses2023a_JASA_EL_figs('fig2'); % Spectrograms and kernels for LAMI and LAPEL
%     publ_osses2023a_JASA_EL_figs('fig2','zenodo','dir_zenodo',dir_zenodo);
%
% % To display suppl. Fig. 1 of Osses et al. (2023, kerAMI) use :::
%     publ_osses2023a_JASA_EL_figs('fig1_suppl'); % Spectrograms and kernels for extra conditions
%     publ_osses2023a_JASA_EL_figs('fig1_suppl','zenodo','dir_zenodo',dir_zenodo);
%
% % To display suppl. Fig. 2 of Osses et al. (2023, kerAMI) use :::
%     publ_osses2023a_JASA_EL_figs('fig2_suppl'); % Target specific kernels
%     publ_osses2023a_JASA_EL_figs('fig2_suppl','zenodo','dir_zenodo',dir_zenodo);
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
data = [];

% close all, clc

if nargin == 0
    help publ_osses2023a_JASA_EL_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig1', ...
    'fig2', ...
    'fig1_suppl', ...
    'fig2_suppl'}; 
definput.flags.plot={'plot','no_plot'};
definput.flags.local={'local','zenodo'};

definput.keyvals.dir_zenodo=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';
 
experiment = 'segmentation';
% dir_exp = [dir_fastACI_results experiment filesep];

bZenodo = flags.do_zenodo;
bLocal = ~bZenodo;

if flags.do_local
    dir_data = fastACI_paths('dir_data');
else
    dir_zenodo = keyvals.dir_zenodo;
    dir_data = [dir_zenodo '02-Raw-data' filesep 'fastACI' filesep 'Publications' filesep 'publ_osses2023a' filesep];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1 || flags.do_fig2 || flags.do_fig1_suppl || flags.do_fig2_suppl
    %%%
    basef = 8000;
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5, ...
        'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};
    x_offset = 0.03;
    
    Condition_list = {  'LAMI','LAPEL', ... main conditions
                        'LACROCH','LALARM','LAMI_SHIFTED'};
    N_participants_expected = [16 18 5 5 3];
    
    if bZenodo
        dir_where_exp  = dir_data;
        dir_where_stim = [dir_zenodo '01-Stimuli' filesep];
        dir_where_post = [dir_zenodo '03-Post-proc-data' filesep];
        if ~exist(dir_where_post,'dir')
            mkdir(dir_where_post);
        end
    end
    if bLocal
        dir_where_exp  = [dir_data experiment filesep]; 
        % dir_where_stim = dir_where_exp;
        dir_where_post = [fastACI_dir_datapost experiment filesep];
    end
end

if flags.do_fig2 || flags.do_fig1_suppl || flags.do_fig2_suppl
    colourbar_map = 'hot_reversed'; % single gradient
    % colourbar_map = 'DG_jet'; % double gradient
    % colourbar_map = 'DG_red_blue'; % Double gradient
    % colourbar_map = 'SQ_red'; % Sequential red
    flags_extra = {'NfrequencyTicks',8,'colourbar_map',colourbar_map, ...
            'colorbar','no'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1
    Show_cell(Condition_list);

    idx_conds = 1:2;
    N_max = max(N_participants_expected(idx_conds));
    
    Perc = nan(N_max,length(idx_conds));
    Rf0 = nan([N_max,9,length(idx_conds)]); % memory allocation
    
    bin_step = 50;
    hist_bin_edges  = -175:bin_step:175;
    hist_bin_centre = -200:bin_step:200;
    
    for i_cond = idx_conds
        bInput = i_cond; % Always LAMI
        Cond = Condition_list{bInput};
        N_part_expected_cond = N_participants_expected(bInput);

        file_subj = il_get_filenames_cond(Cond,dir_where_exp);

        N_participants = length(file_subj);
        if N_participants ~= N_part_expected_cond
            warning('Not the same number of participants as expected!')
            pause(10);
        end
        for i_subject = 1:N_participants
            % dir_res = il_get_data_path('dir_res', bZenodo,file_subj{i_subject});
            
            dir_subj = [dir_where_exp file_subj{i_subject} filesep]; 
            if bZenodo == 0
                dir_res = [dir_subj 'Results' filesep];
            else
                dir_res = [dir_subj '1-experimental_results' filesep];
            end
            
            %%%
            if bLocal
                dir_res_post = [dir_where_post file_subj{i_subject} filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
                dir_res_post = [dir_res_post 'Results' filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            end
            if bZenodo
                dir_res_post = dir_where_post;
            end
            %%%
            %%% Load data
            file_save = Get_filenames(dir_res,'*savegame*.mat');
            if length(file_save)~=1
                error('More than one savegame found...');
            end
            savegame_path = [dir_res file_save{end}];

            cfg_game = [];
            data_passation = [];
            load(savegame_path);
            if bLocal
                dir_where_stim = dir_subj;
            end
            cfg_game.dir_target = [dir_where_stim 'speech-samples' filesep];
            %%%
            % Plot targets with segments (replaced by 'fig1a')
            if i_subject == 1
                N_segment = 8; % Prior knowledge
                t_edge = (0:N_segment)/10;
            end

            % Scores
            score(i_subject) = mean(data_passation.is_correct);
            bias(i_subject)  = mean(data_passation.n_responses);

            flags_in = {'no_plot','dir_out',dir_res_post}; % ,'dir_noise',dir_subj};

            if bLocal
                dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise = [dir_subj 'Stim-processed' filesep];
            end
            if bZenodo
                subj_id_here = strsplit(file_subj{i_subject},'_');
                subj_id_here = subj_id_here{end};
                dir_target = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise  = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'Stim-processed' filesep];
            end
            if exist(dir_noise,'dir')
                flags_in(end+1:end+2) = {'dir_noise',dir_noise};
            end
            
            if exist(dir_target,'dir')
                flags_in(end+1:end+2) = {'dir_target',dir_target};
            end

            % GLM-based revcorr
            flags_in_here = flags_in;
 
            bForce_dataload = 1;
            if bForce_dataload
                flags_in_here{end+1} = 'force_dataload';
            end

            [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savegame_path,'glm',flags_in_here{:});

            dim_f0 = 1; % dimension where the f0 is...
            switch i_cond
                case 1 % LAMI, point in f0 kernel with the highest weighting
                    segment_nr = 6; % segment number
                case 2
                    segment_nr = 5; % segment number
            end

            Data_f0_subset = Data_matrix(:,segment_nr,dim_f0);
            
            % %%% Correction by Alejandro: '<' replaced by '<='
            % Rate of responses for high- mid- and low-f0s:
            %   The responses could be ''l'aX'' (1) or ''la X'' (2). The average
            %   minus one indicate the rate for which 'la X' is chosen, because
            %   the closer to 0, the more l'aX was chosen and the closer to 1
            %   the more 'la X'.
            
            Perc(i_subject,i_cond) = 100*sum(data_passation.is_correct)/length(data_passation.is_correct);
            
            %%% Generating the histograms:
            n_responses_1_or_2 = data_passation.n_responses; % 1 for l'amie, 2 for la mie
            
            %%% Old plot showing la X proportions:
            % prefix_title = 'la X';
            % n_responses = n_responses_1_or_2 - 1; % 1 for l'amie, 0 for la mie
            %%%
            prefix_title = 'l''aX';
            n_responses = abs(n_responses_1_or_2 - 2); % 1 for l'amie, 0 for la mie
            %%%
            id_matrix = 1;
            idx_f0 = Data_f0_subset<= hist_bin_edges(id_matrix);
            Rf0(i_subject,id_matrix,i_cond) = 100*mean(n_responses(idx_f0));

            for id_matrix = 2:length(hist_bin_edges)
                idx_f0 = Data_f0_subset<= hist_bin_edges(id_matrix) & Data_f0_subset> hist_bin_edges(id_matrix-1);
                Rf0(i_subject,id_matrix,i_cond) = 100*mean(n_responses(idx_f0));
            end
            
            id_matrix = length(hist_bin_edges)+1;
            idx_f0 = Data_f0_subset>  hist_bin_edges(id_matrix-1);
            Rf0(i_subject,id_matrix,i_cond) = 100*mean( n_responses(idx_f0) );             
        end
    end
    for i_cond = idx_conds
        Cond = Condition_list{i_cond};
        Score_here = Perc(1:N_participants_expected(idx_conds),i_cond);
        Scores_mean(i_cond) = mean(Score_here);
        Scores_std(i_cond) = std(Score_here);        
        fprintf('%s: Mean score=%.1f (SD=%.1f)\n',Cond,Scores_mean(i_cond),Scores_std(i_cond));
        Me_seg = nanmean(Rf0(:,:,i_cond));
        fprintf('\tLeft- and right-most score in the histogram: %.1f and %.1f percent\n',Me_seg([1 end]));
    end

    figure('Position',[100 100 1000 300]);
    tiledlayout(1,2,'TileSpacing','tight');
    for i_cond = idx_conds
        nexttile(i_cond); 
        
        N_part_expected_cond = N_participants_expected(i_cond);
        idx_subj = 1:N_part_expected_cond;

        yvar = mean(Rf0(idx_subj,:,i_cond),1);
        errLU = 1.96*sem(Rf0(idx_subj,:,i_cond),[],1)/2;

        errorbar(hist_bin_centre, yvar, errLU); grid on
        xlim([-220 220])
        hold on
        
        gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
        plot([-250 250],[50 50],'k--');
        
        xvar = -220:220;
        mu   = median(xvar); % centred at mid point of xvar (0)
        sigma = 100;
        A    = 90; % arbitrary amplitude
        vo   = 0; % vertical offset
        area(xvar, gaus(xvar,mu,sigma,A,vo),'FaceColor', [0,0,1] , 'FaceAlpha', 0.1,'EdgeColor','none');

        ylim([0 100])
        xlabel('f0 shift (cents)')
        YT = 0:10:100;
        set(gca,'YTick',YT);
        if i_cond == 1
            ylabel(['Proportion of ' prefix_title ' responses (%)']);
            title('A. Condition LAMI')
        else
            set(gca,'YTickLabel','');
            title('B. Condition LAPEL')
        end
        % xlim([-250 250])
    end
    
    h(end+1) = gcf;
    hname{end+1} = 'fig1-hist';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2
    %%%
    for i_cond = 1:2
        switch i_cond
            case 1 % 'LAMI'
                Subject_ID = 'S001'; % first participant who conducted l'amie / la mie
                fname1 = 'l_amie.wav';
                fname2 = 'la_mie.wav';
                label1 = 'A. C''est l''amie';
                label2 = '   C''est la mie';
                t_seg1 = [.11  .305  .4507 0.6425]; sentence1 = {'C''est','l''a -','mie'}; % manually from Praat
                t_seg2 = [.104 .2935 .4507 0.71];   sentence2 = {'C''est','la','mie'}; % manually from Praat 
                
            case 2 %'LAPEL'
                Subject_ID = 'S021';
                fname1 = 'l_appel.wav';
                fname2 = 'la_pelle.wav';
                label1 = 'B. C''est l''appel';
                label2 = '   C''est la pelle';
                t_seg1 = [.105  .2818  .4567 0.77]; sentence1 = {'C''est','l''a -','ppel'}; % manually from Praat
                t_seg2 = [.11  .264 .463 0.77];     sentence2 = {'C''est','la','pelle'}; % manually from Praat 
                
        end
        
        if bLocal
            dir_subj = [dir_data experiment filesep Subject_ID filesep];
            dir_stim = [dir_subj 'speech-samples' filesep];
        end
        if bZenodo
            dir_stim = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep Subject_ID filesep 'speech-samples' filesep];
        end
        
        fname1 = [dir_stim fname1]; 
        fname2 = [dir_stim fname2]; 
            
        [insig1, fs] = audioread(fname1);
        [insig2, fs] = audioread(fname2);

        %%% Plot targets + 1 example of noise
        [G1, fc, t, outs] = Gammatone_proc(insig1, fs, flags_gamma{:});
        [G2, fc, t, outs] = Gammatone_proc(insig2, fs, flags_gamma{:});

        if i_cond == 1
            figure('Position',[100 100 1000 700]); % 600 500
            tiledlayout(4,2,'TileSpacing','tight');
        end
        nexttile(1+(i_cond-1));
        
        G_dB = To_dB(G1');
        if i_cond == 2
            flags_extra_here = flags_extra;
            flags_extra_here{end} = 'yes';
        else
            flags_extra_here = flags_extra;
        end
        affichage_tf(G_dB, 'pow', t, fc, flags_extra_here{:}); hold on;
        XTL = 0:.1:.8;
        XT  = XTL; XT(1) = min(t)/2;
        set(gca,'XTick',XT);
        set(gca,'XTickLabel','');
        if i_cond ~= 1
            set(gca,'YTickLabel','');
            ylabel('');
        end
        %%%
        max_dB = max(max(G_dB));  
        CLim = [max_dB-51 max_dB]; 
        %%%
        caxis(CLim); 
        
        xlabel('')
        if i_cond == 2
            hc = get(gca,'Colorbar');
            CB_Tick = get(hc,'Ticks'); % [-60 -50 -40 -30 -20]
            CB_Tick = [CB_Tick(1)-10 CB_Tick(1):10:CB_Tick(end)]; % 10 dB below
            
            CB_TickLabel = [];
            CB_TickLabel{1} = '<-50 dB';
            for i_tick = 2:length(CB_Tick)
                CB_TickLabel{i_tick} = sprintf('%.0f',CB_Tick(i_tick) - max(CB_Tick));
            end
            
            
            set(hc,'Ticks',CB_Tick);
            set(hc,'TickLabels',CB_TickLabel);
        end
        % title(label1); % title(['c''est l''' cfg_game.response_names{1}]);

        time_dly_compensation = 25e-3; % 30e-3; % 2*5.2e-3; % from Osses2019, pp. 1026
        YL_here = get(gca,'YLim');
        DR = YL_here(2)-YL_here(1); 
        t_seg = t_seg1 + time_dly_compensation; % because time delays are not compensated
        
        for i_sent = 1:length(sentence1)
            idxi = i_sent;
            idxf = i_sent+1;
            % Horizontal line:
            plot(t_seg(idxi:idxf), YL_here(2)*[1 1],'-','LineWidth',2,'Color','m');
            % Vertical lines:
            plot(t_seg(idxi)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
            plot(t_seg(idxf)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
        
            x_text = .25*(t_seg(idxf)-t_seg(idxi))+t_seg(idxi);
            y_text = YL_here(2) + .1*DR;
            text(x_text,y_text,sentence1{i_sent},'FontWeight','Bold','FontSize',11,'Color','m');
        end
        
        if i_cond == 1
            x_panel_label = -.2;
            y_panel_label = 1.15;
            text(x_panel_label,y_panel_label,'A.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end
        
        cfg_in = [];
        cfg_in.f = fc;
        cfg_in.t = t;
        bPlot_trajectories = 0;
        outs_from_Praat1 = Add_formants(fname1,cfg_in); % with default parameters

        %%% Adding vertical lines
        YL = get(gca,'YLim');
        for i = 1:length(XT)
            plot(XT(i)*[1 1],YL,'k--');
        end
        xlim([-x_offset max(XT)+x_offset])
        %%%

        % flags_extra{end} = 'on';
        nexttile(2+i_cond); 
        G_dB = To_dB(G2');
        affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); hold on;
        set(gca,'XTick',XT);
        set(gca,'XTickLabel',XTL);
        if i_cond ~= 1
            set(gca,'YTickLabel','');
            ylabel('');
        end
        caxis(CLim); 
        % title(label2);
        
        %%%
        YL_here = get(gca,'YLim');
        DR = YL_here(2)-YL_here(1); 
        t_seg = t_seg2 + time_dly_compensation; % because time delays are not compensated
        
        for i_sent = 1:length(sentence2)
            idxi = i_sent;
            idxf = i_sent+1;
            % Horizontal line:
            plot(t_seg(idxi:idxf), YL_here(2)*[1 1],'-','LineWidth',2,'Color','m');
            % Vertical lines:
            plot(t_seg(idxi)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
            plot(t_seg(idxf)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
        
            x_text = .25*(t_seg(idxf)-t_seg(idxi))+t_seg(idxi);
            y_text = YL_here(2) + .1*DR;
            text(x_text,y_text,sentence2{i_sent},'FontWeight','Bold','FontSize',11,'Color','m');
        end

        if i_cond == 1
            text(x_panel_label,y_panel_label,'B.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end
        
        outs_from_Praat2 = Add_formants(fname2,cfg_in); % with default parameters

        %%%
        YL = get(gca,'YLim');
        for i = 1:length(XT)
            plot(XT(i)*[1 1],YL,'k--');
        end
        xlim([-x_offset max(XT)+x_offset])
        %%%
        fprintf('Praat estimations:\n');
        f0_1 = outs_from_Praat1.f0{1};
        fprintf('\t Sound 1: %.1f Hz (%.1f,%.1f), mean=%.1f Hz\n',prctile(f0_1,50),prctile(f0_1,25),prctile(f0_1,75),nanmean(f0_1));
        f0_2 = outs_from_Praat2.f0{1};
        fprintf('\t Sound 2: %.1f Hz (%.1f,%.1f), mean=%.1f Hz\n',prctile(f0_2,50),prctile(f0_2,25),prctile(f0_2,75),nanmean(f0_2));
        % Sound 1: 206.3 Hz (195.5,239.0); % Praat 'Pitch listing': 216.1 Hz
        % Sound 2: 193.3 Hz (187.2,221.9); % Praat 'Pitch listing': 205.8 Hz

        Mean_of_mean = [nanmean(f0_1)+nanmean(f0_2)]/2;
        fprintf('\t The mean of both means is %.1f Hz\n',Mean_of_mean);

        % figure; plot(f0_1); hold on; plot(f0_2,'r--');
        %%%
    end
    
    Show_cell(Condition_list);
    
    for i_cond = 1:2
        bInput = i_cond; % Always LAMI
        Cond = Condition_list{bInput};
        N_part_expected_cond = N_participants_expected(bInput);

        file_subj = il_get_filenames_cond(Cond,dir_where_exp);
        
        N_participants = length(file_subj);
        if N_participants ~= N_part_expected_cond
            warning('Not the same number of participants as expected!')
            pause(10);
        end
        for i_subject = 1:N_participants
            dir_subj = [dir_where_exp file_subj{i_subject} filesep]; 
            if bZenodo == 0
                dir_res = [dir_subj 'Results' filesep];
            else
                dir_res = [dir_subj '1-experimental_results' filesep];
            end
            %%%
            if bLocal
                dir_res_post = [dir_where_post file_subj{i_subject} filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
                dir_res_post = [dir_res_post 'Results' filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            end
            if bZenodo
                dir_res_post = dir_where_post;
            end
            %%%
            %%% Load data
            file_save = Get_filenames(dir_res,'*savegame*.mat');
            if length(file_save)~=1
                error('More than one savegame found...');
            end
            savegame_path = [dir_res file_save{end}];

            cfg_game = [];
            data_passation = [];
            load(savegame_path)
            
            if bLocal
                dir_where_stim = dir_subj;
            end
            cfg_game.dir_target = [dir_where_stim 'speech-samples' filesep];
            % cfg_game.dir_target = [dir_where_exp filesep 'speech-samples' filesep];
            %%%
            % Plot targets with segments (replaced by 'fig1a')
            if i_subject == 1
                 N_segment = 8; % Prior knowledge
                t_edge = (0:N_segment)/10;
            end

            % Scores
            score(i_subject) = mean(data_passation.is_correct);
            bias(i_subject)  = mean(data_passation.n_responses);

            flags_in = {'no_plot','dir_out',dir_res_post};

            if bLocal
                dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise = [dir_subj 'Stim-processed' filesep];
            end
            if bZenodo
                subj_id_here = strsplit(file_subj{i_subject},'_');
                subj_id_here = subj_id_here{end};
                dir_target = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise  = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'Stim-processed' filesep];
            end
            
            % dir_noise = [dir_subj 'Stim-processed' filesep];
            if exist(dir_noise,'dir')
                flags_in(end+1:end+2) = {'dir_noise',dir_noise};
            end
            % dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
            if exist(dir_target,'dir')
                flags_in(end+1:end+2) = {'dir_target',dir_target};
            end

            % GLM-based revcorr
            flags_in_here = flags_in;

            bForce_dataload = 0;
            if bForce_dataload
                flags_in_here{end+1} = 'force_dataload';
            end

            [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savegame_path,'glm',flags_in_here{:});
            % [ACI_t1,cfg_ACI_t1,results_t1]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t1',flags_in{:});
            % [ACI_t2,cfg_ACI_t2,results_t2]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t2',flags_in{:});

            f0_kernel(:,i_subject) = ACI(:,1);
            time_kernel(:,i_subject) = ACI(:,2);
        end
        
        % Stats at the group level
        
        alpha = 0.005; % arbitrary p-value
        mYL = .35; % ylim for the plot
        
        [h_f0_kernel,p_f0_kernel,~,t_f0_kernel] = ttest(f0_kernel',0,'alpha', alpha);
        [h_time_kernel,p_time_kernel,~,t_time_kernel] = ttest(time_kernel',0,'alpha', alpha);
        
        %%%
        fprintf('In the %s condition, the f0 kernel reveals %.0f critical regions\n',Cond,sum(h_f0_kernel));
        dof = unique(t_f0_kernel.df); % normally, the number of participants minus 1
        for i_t = 1:length(h_f0_kernel)
            if h_f0_kernel(i_t) == 1
                tstat = t_f0_kernel.tstat(i_t);
                pval  = p_f0_kernel(i_t); 
                % Significant cases:
                fprintf('\t Time segment=%.0f: t(%.0f)=%.2f, p=%.5f\n',i_t,dof,tstat,pval);
            end
        end
        fprintf('Similarly for the time kernel, %.0f regions were critical\n',nansum(h_time_kernel)); 
        dof = unique(t_time_kernel.df); % normally, the number of participants minus 1
        for i_t = 1:length(h_time_kernel)
            if h_time_kernel(i_t) == 1
                tstat = t_time_kernel.tstat(i_t);
                pval  = p_time_kernel(i_t); 
                % Significant cases:
                fprintf('\t Time segment=%.0f: t(%.0f)=%.2f, p=%.5f\n',i_t,dof,tstat,pval);
            end
        end
        fprintf('\n\n');
        % the weights associated to the segment edges at 0.3 and 0.6 s \hlNew{were} found to be significant (segment at 0.3~s: $t(15)=3.80, p<0.005$; segment at 0.6~s: $t(15)=-3.58, p<0.005$).

        tile_start = 4;
        nexttile(1+tile_start+(i_cond-1)); 
        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on

        Me = mean(f0_kernel,2);
        errL = 1.96*sem(f0_kernel,[],2)/2;
        errU = 1.96*sem(f0_kernel,[],2)/2;
        errorbar(t_edge,Me,errL,errU,'k','LineWidth',2); hold on

        title(['Fundamental-frequency kernel (' cfg_game.response_names{1} '-' cfg_game.response_names{2} ')']); 
        ylabel('f_0 shift weight'); 
        % xlabel('Segment edge (s)')
        set(gca,'XTick',t_edge)
        set(gca,'XTickLabel','')
        if i_cond ~= 1
            set(gca,'YTickLabel','');
            ylabel('');
        end
        xlim([0-x_offset t_edge(end)+x_offset])
        grid on
        
        if i_cond == 1
            text(x_panel_label,y_panel_label,'C.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end
        
        %%% Asterisk when significantly different from zero:
        idxs = find(h_f0_kernel~=0 & ~isnan(h_f0_kernel));
        
        plot(t_edge(idxs),.9*mYL*h_f0_kernel(idxs),'k*');
        plot(t_edge(idxs),Me(idxs),'ko','MarkerFaceColor','r');
        %%%
        
        text(.85,.9,['N = ' num2str(i_subject)],'Units','normalized')
        % YL  = ylim; 
        
        ylim(mYL*[-1 1]);
        YT = -.4:.1:.4;
        set(gca,'YTick',YT);

        nexttile(2+tile_start+(i_cond));
        % plot(t_edge,zeros(size(t_edge)),'k--');
        % errorbar(t_edge,mean(time_kernel_t1,2),1.96*sem(time_kernel_t1,[],2)/2,'r','LineWidth',1); hold on
        % errorbar(t_edge,mean(time_kernel_t2,2),1.96*sem(time_kernel_t2,[],2)/2,'b','LineWidth',1); hold on
        
        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on
        
        Me = mean(time_kernel,2);
        errL = 1.96*sem(time_kernel,[],2)/2;
        errU = 1.96*sem(time_kernel,[],2)/2;
        errorbar(t_edge,Me,errL,errU,'k','LineWidth',2); hold on
        title(['Time kernel (' cfg_game.response_names{1} '-' cfg_game.response_names{2} ')']); 
        ylabel('Time-shift weight'); 
        xlabel('Segment edge (s)')
        set(gca,'XTick',t_edge)
        xlim([0-x_offset t_edge(end)+x_offset])
        grid on
        
        if i_cond == 1
            text(x_panel_label,y_panel_label,'D.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end
        
        %%% Asterisk when significantly different from zero:
        idxs = find(h_time_kernel~=0 & ~isnan(h_time_kernel));
        
        plot(t_edge(idxs),.9*mYL*h_time_kernel(idxs),'k*');
        plot(t_edge(idxs),Me(idxs),'ko','MarkerFaceColor','r');
        %%%
        
        ylim(mYL*[-1 1])
        set(gca,'YTick',YT);
        
        if i_cond ~= 1
            set(gca,'YTickLabel','');
            ylabel('');
        end
    end
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2-LAMI+LAPEL';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1_suppl
    idx_conds = [3 4 5]; % it skips the main condition LAMI
    XT_max = [.9 .9 .8];
    Label_panel = {'A.','B.','C.'};
    for i_cond = 1:length(idx_conds)
        bInput = idx_conds(i_cond);
        Cond = Condition_list{bInput};
        
        file_subj = il_get_filenames_cond(Cond,dir_where_exp);
        
        offx = 0;
        switch Cond
            case 'LACROCH' 
                wave1 = 'l_accroche.wav';
                wave2 = 'la_croche.wav';
                t_seg1 = [.1259  .2759 .4249 0.6756 0.81];  sentence1 = {'C''est','l''a ','c-cro - ','che'}; % manually from Praat
                t_seg2 = [.0959 .2738  .4099 0.7067 0.835]; sentence2 = {'C''est','la'   ,'cro - ','che'}; % manually from Praat 
                
            case 'LALARM'
                wave1 = 'l_alarme.wav';
                wave2 = 'la_larme.wav';
                t_seg1 = [0.0938  .2946  .4636 0.7269 0.85]; sentence1 = {'C''est','l''a -','lar - ','me'}; % manually from Praat
                t_seg2 = [0.1081 .3122 .4526 0.666 0.831];   sentence2 = {'C''est','la','lar - ','me'}; % manually from Praat 
                
            case 'LAMI_SHIFTED'
                offx = 0.02;
                wave1 = 'l_amie.wav';
                wave2 = 'la_mie.wav';     
                t_seg1 = [.11  .305  .4507 0.6425]; sentence1 = {'C''est','l''a -','mie   (LAMI-shifted)'}; % manually from Praat
                t_seg2 = [.104 .2935 .4507 0.71];   sentence2 = {'C''est','la','mie'}; % manually from Praat 
                
                t_for_kernels  = [0       0.1     0.2     0.3    0.4     0.5     0.6    0.7      0.8];
                kernel_f0_orig = [0.0266 -0.0217 -0.0158 -0.0763 0.00703 0.1996 -0.1362 0.01906 -0.02153];
                kernel_t_orig  = [0 0.0005 0.01557 0.0794 0.0807 0.0059 -0.1166 0.0076 0]; 
        end       
        file_subj = file_subj(:);

        N_participants = length(file_subj);
        if N_participants ~= N_participants_expected
            warning('Not the same number of participants as expected!')
            pause(10);
        end
        
        for i_subject = 1:N_participants
            Subject_ID = file_subj{i_subject};
            if bZenodo
                Subject_ID = strsplit(Subject_ID,'_');
                Subject_ID = Subject_ID{end};
            end
            dir_subj = [dir_where_exp file_subj{i_subject} filesep]; 
            if bZenodo == 0
                dir_res = [dir_subj 'Results' filesep];
            else
                dir_res = [dir_subj '1-experimental_results' filesep];
            end
            %%%
            if bLocal
                dir_res_post = [dir_where_post file_subj{i_subject} filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
                dir_res_post = [dir_res_post 'Results' filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            end
            if bZenodo
                dir_res_post = dir_where_post;
            end
            %%% Load data
            file_save = Get_filenames(dir_res,'*savegame*.mat');
            if length(file_save)~=1
                error('More than one savegame found...');
            end
                        
            % dir_subj = [dir_data experiment filesep Subject_ID filesep];
            % dir_res = [dir_subj 'Results' filesep];
            % dir_res_post = [dir_where_post Subject_ID filesep];
            % if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            % dir_res_post = [dir_res_post 'Results' filesep];
            % if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            % %%% Load data
            % file_save = Get_filenames(dir_res,'*savegame*.mat');
            % if length(file_save)~=1
            %     error('More than one savegame found...');
            % end
            savegame_path = [dir_res file_save{end}];
        
            cfg_game = [];
            data_passation = [];
            load(savegame_path);
            
            if bLocal
                dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise = [dir_subj 'Stim-processed' filesep];
            end
            if bZenodo
                subj_id_here = strsplit(file_subj{i_subject},'_');
                subj_id_here = subj_id_here{end};
                dir_target = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise  = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'Stim-processed' filesep];
            end
            % cfg_game.dir_target = dir_target;
            %%% End load data

            if i_subject == 1
                % N_segment = XT_max/.1; % Prior knowledge
                % t_edge = (0:N_segment)/10;

                % fname1 = [dir_subj 'speech-samples' filesep wave1];
                % fname2 = [dir_subj 'speech-samples' filesep wave2];
 
                if bLocal
                    dir_subj = [dir_data experiment filesep Subject_ID filesep];
                    dir_stim = [dir_subj 'speech-samples' filesep];
                end
                if bZenodo
                    dir_stim = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep Subject_ID filesep 'speech-samples' filesep];
                end

                fname1 = [dir_stim wave1]; 
                fname2 = [dir_stim wave2]; 
 
                [insig1, fs] = audioread(fname1);
                [insig2, fs] = audioread(fname2);
 
                %%% Plot targets + 1 example of noise
                [G1, fc, t, outs] = Gammatone_proc(insig1, fs, flags_gamma{:});
                [G2, fc, t, outs] = Gammatone_proc(insig2, fs, flags_gamma{:});

                if i_cond == 1
                    figure('Position',[100 100 2.2*600 750]); 
                    tiledlayout(4,3,'TileSpacing','tight');
                end
                curr_tile = i_cond;
                nexttile(curr_tile);
                
                if i_cond == 3
                    flags_extra_here = flags_extra;
                    flags_extra_here{end}='yes'; % setting colourbar to 'yes'
                else
                    flags_extra_here = flags_extra;
                end
                G_dB = To_dB(G1');
                affichage_tf(G_dB, 'pow', t, fc, flags_extra_here{:}); hold on;
                XTL = 0:.1:XT_max(i_cond);
                XT  = XTL+offx; 
                set(gca,'XTick',XT);
                set(gca,'XTickLabel','');
                
                if i_cond ~= 1
                    ylabel('');
                    set(gca,'YTickLabel',[]);
                end
                %%%
                max_dB = max(max(G_dB));  
                CLim = [max_dB-51 max_dB]; 
                xlabel('')
                caxis(CLim); 
                
                if i_cond == 3
                    hc = get(gca,'Colorbar');
                    % CB_Tick = get(hc,'Ticks'); % [-60 -50 -40 -30 -20]
                    CB_Tick = [CLim(1)+1:10:CLim(2)]; % 10 dB below

                    CB_TickLabel = [];
                    CB_TickLabel{1} = '<-50 dB';
                    for i_tick = 2:length(CB_Tick)
                        CB_TickLabel{i_tick} = sprintf('%.0f',CB_Tick(i_tick) - max(CB_Tick));
                    end

                    set(hc,'Ticks',CB_Tick);
                    set(hc,'TickLabels',CB_TickLabel);
                end
         
                time_dly_compensation = 15e-3; % 30e-3; % 2*5.2e-3; % from Osses2019, pp. 1026
                YL_here = get(gca,'YLim');
                DR = YL_here(2)-YL_here(1); 
                t_seg = t_seg1 + time_dly_compensation; % because time delays are not compensated
                for i_sent = 1:length(sentence1)
                    idxi = i_sent;
                    idxf = i_sent+1;
                    % Horizontal line:
                    plot(t_seg(idxi:idxf), YL_here(2)*[1 1],'-','LineWidth',2,'Color','m');
                    % Vertical lines:
                    plot(t_seg(idxi)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
                    plot(t_seg(idxf)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');

                    x_text = .25*(t_seg(idxf)-t_seg(idxi))+t_seg(idxi);
                    y_text = YL_here(2) + .1*DR;
                    text(x_text,y_text,sentence1{i_sent},'FontWeight','Bold','FontSize',11,'Color','m');
                end

                if i_cond == 1
                    x_panel_label = -.2;
                    y_panel_label = 1.15;
                    text(x_panel_label,y_panel_label,'A.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
                end
                
                %%% Adding vertical lines
                YL = get(gca,'YLim');
                for i = 1:length(XT)
                    plot(XT(i)*[1 1],YL,'k--');
                end
                xlim([-x_offset max(XT)+x_offset])
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curr_tile = i_cond+3;
                nexttile(curr_tile);
                
                G_dB = To_dB(G2');
                affichage_tf(G_dB, 'pow', t, fc, flags_extra{:}); hold on;
                % XT  = XTL; 
                set(gca,'XTick',XT);
                % set(gca,'XTickLabel','');
                if i_cond ~= 1
                    ylabel('');
                    set(gca,'YTickLabel',[]);
                end
                %%%
                max_dB = max(max(G_dB));  
                CLim = [max_dB-50 max_dB]; 
                %%%
                % xlabel('')
                caxis(CLim); 
                % title(sprintf('   C''est la %s',cfg_game.response_names{2})); 
                
                YL_here = get(gca,'YLim');
                DR = YL_here(2)-YL_here(1); 
                t_seg = t_seg2 + time_dly_compensation; % because time delays are not compensated

                for i_sent = 1:length(sentence2)
                    idxi = i_sent;
                    idxf = i_sent+1;
                    % Horizontal line:
                    plot(t_seg(idxi:idxf), YL_here(2)*[1 1],'-','LineWidth',2,'Color','m');
                    % Vertical lines:
                    plot(t_seg(idxi)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');
                    plot(t_seg(idxf)*[1 1], YL_here(2)*[1 .93],'-','LineWidth',2,'Color','m');

                    x_text = .25*(t_seg(idxf)-t_seg(idxi))+t_seg(idxi);
                    y_text = YL_here(2) + .1*DR;
                    text(x_text,y_text,sentence2{i_sent},'FontWeight','Bold','FontSize',11,'Color','m');
                end

                if i_cond == 1
                    text(x_panel_label,y_panel_label,'B.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
                end
                
                %%% Adding vertical lines
                YL = get(gca,'YLim');
                for i = 1:length(XT)
                    plot(XT(i)*[1 1],YL,'k--');
                end
                xlim([-x_offset max(XT)+x_offset])
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Getting Praat (not plotted because Add_formats has one output argument):
                cfg_in = [];
                cfg_in.f = fc;
                cfg_in.t = t;
                outs_from_Praat1 = Add_formants(fname1,cfg_in); % with default parameters
                outs_from_Praat2 = Add_formants(fname2,cfg_in); % with default parameters
                
                f0_all(1,i_cond) = nanmean(outs_from_Praat1.f0{1});
                f0_all(2,i_cond) = nanmean(outs_from_Praat1.f0{1});
                Mean_of_mean = [f0_all(1,:)+f0_all(2,:)]/2;
            end % if i_subj == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Scores
            score{i_cond}(i_subject) = mean(data_passation.is_correct);
            bias{i_cond}(i_subject)  = mean(data_passation.n_responses);

            flags_in = {'no_plot','dir_out',dir_res_post};
 
            % dir_noise = [dir_subj 'Stim-processed' filesep];
            if exist(dir_noise,'dir')
                flags_in(end+1:end+2) = {'dir_noise',dir_noise};
            end
            if exist(dir_target,'dir')
                flags_in(end+1:end+2) = {'dir_target',dir_target};
            end
 
            % GLM-based revcorr
            flags_in_here = flags_in;

            bForce_dataload = 0;
            if bForce_dataload
                flags_in_here{end+1} = 'force_dataload';
            end

            [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(savegame_path,'glm',flags_in_here{:});
            [ACI_t1,cfg_ACI_t1,results_t1]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t1',flags_in{:});
            [ACI_t2,cfg_ACI_t2,results_t2]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t2',flags_in{:});

            f0_kernel{i_cond}(:,i_subject) = ACI(:,1);
            time_kernel{i_cond}(:,i_subject) = ACI(:,2);
        end % if i_subject
        %%%
        f0_kernel_here = f0_kernel{i_cond};
        Me = mean(f0_kernel_here,2);
        errL = 1.96*sem(f0_kernel_here,[],2)/2;
        errU = 1.96*sem(f0_kernel_here,[],2)/2;
        
        %%%
        t_edge = (0:length(Me)-1)/10+offx;
        %%%
        
        nexttile(i_cond+2*3);

        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on
        switch Cond
            case 'LAMI_SHIFTED'
                plot(t_for_kernels,kernel_f0_orig,'o--','Color',0.5*[1 1 1],'LineWidth',2,'MarkerFaceColor',0.5*[1 1 1]);
        end
        
        errorbar(t_edge,Me,errL,errU,'k','LineWidth',2); hold on
         
        if i_cond == 2
            prefix = 'Fundamental-frequency kernel ';
        else
            prefix = '';
        end
        title(['   ' prefix '(' cfg_game.response_names{1} '-' cfg_game.response_names{2} ')']); 
        
        YT = -.4:.1:.4;
        set(gca,'YTick',YT);
        if i_cond == 1
            ylabel('f_0 shift weight'); 
        else
            set(gca,'YTickLabel','');
        end
        set(gca,'XTick',t_edge)
        set(gca,'XTickLabel','')
        xlim([0-x_offset t_edge(end)+x_offset])
        grid on
        % box on
        text(.85,.9,['N = ' num2str(N_participants)],'Units','normalized')
        % YL  = ylim; 
        mYL = 0.45; % 1.4*max(abs(YL));
        ylim(mYL*[-1 1]);
        % set(gca,'YTick',YT);
        if i_cond == 1
            text(x_panel_label,y_panel_label,'C.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end

        nexttile(i_cond+3*3);
        time_kernel_here = time_kernel{i_cond};
        
        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on
        switch Cond
            case 'LAMI_SHIFTED'
                plot(t_for_kernels,kernel_t_orig,'o--','Color',0.5*[1 1 1],'LineWidth',2,'MarkerFaceColor',0.5*[1 1 1]);
        end

        Me = mean(time_kernel_here,2);
        errL = 1.96*sem(time_kernel_here,[],2)/2;
        errU = 1.96*sem(time_kernel_here,[],2)/2;
        errorbar(t_edge,Me,errL,errU,'k','LineWidth',2); hold on
        if i_cond == 2
            prefix = 'Time kernel ';
        else
            prefix = '';
        end
        title(['   ' prefix '(' cfg_game.response_names{1} '-' cfg_game.response_names{2} ')']); 
        set(gca,'YTick',YT);
        if i_cond == 1
            ylabel('Time-shift weight'); 
        else
            set(gca,'YTickLabel','');
        end
        xlabel('Segment edge (s)')
        set(gca,'XTick',t_edge)
        xlim([0-x_offset t_edge(end)+x_offset])
        grid on
        % box on
        % text(.85,.8,['N = ' num2str(i_subject)],'Units','normalized')
        % YL = ylim; 
        % mYL=max(abs(YL));
        % ylim([-1.4 1.4]*mYL);
        ylim(mYL*[-1 1]);
        
        if i_cond == 1
            text(x_panel_label,y_panel_label,'D.','FontWeight','Bold','FontSize',14,'Unit','Normalized');
        end
    end % End cond
    h(end+1) = gcf;
    hname{end+1} = ['suppl-fig1-kernels-extra'];

    fprintf('\t The mean of both means is %.1f, %.1f, %.1f Hz\n',Mean_of_mean);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2_suppl
    Show_cell(Condition_list);
    
    YLs = [-0.3 0.4; -.3 .3; -.3 .5; -.5 .7];
    stepY = [.1 .1 .1 .2];
    lab = {'A.','B.','C.','D.'};
    
    for i_cond = 1:4
        bInput = i_cond; % Always LAMI
        Cond = Condition_list{bInput};
        N_part_expected_cond = N_participants_expected(bInput);

        file_subj = il_get_filenames_cond(Cond,dir_where_exp);
        
        N_participants = length(file_subj);
        if N_participants ~= N_part_expected_cond
            warning('Not the same number of participants as expected!')
            pause(10);
        end
        for i_subject = 1:N_participants
            Subject_ID = file_subj{i_subject};
            if bZenodo
                Subject_ID = strsplit(Subject_ID,'_');
                Subject_ID = Subject_ID{end};
            end
            dir_subj = [dir_where_exp file_subj{i_subject} filesep]; 
            if bZenodo == 0
                dir_res = [dir_subj 'Results' filesep];
            else
                dir_res = [dir_subj '1-experimental_results' filesep];
            end
            %%%
            if bLocal
                dir_res_post = [dir_where_post file_subj{i_subject} filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
                dir_res_post = [dir_res_post 'Results' filesep];
                if ~exist(dir_res_post,'dir'); mkdir(dir_res_post); end
            end
            if bZenodo
                dir_res_post = dir_where_post;
            end
            %%% Load data
            file_save = Get_filenames(dir_res,'*savegame*.mat');
            if length(file_save)~=1
                error('More than one savegame found...');
            end
            savegame_path = [dir_res file_save{end}];

            cfg_game = [];
            data_passation = [];
            load(savegame_path)
            % cfg_game.dir_target = [dir_where_exp filesep 'speech-samples' filesep];
            
            % Scores
            score(i_subject) = mean(data_passation.is_correct);
            bias(i_subject)  = mean(data_passation.n_responses);

            if bLocal
                dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise = [dir_subj 'Stim-processed' filesep];
            end
            if bZenodo
                subj_id_here = strsplit(file_subj{i_subject},'_');
                subj_id_here = subj_id_here{end};
                dir_target = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
                dir_noise  = [dir_where_stim 'fastACI_data' filesep 'segmentation' filesep subj_id_here filesep 'Stim-processed' filesep];
            end
            
            flags_in = {'no_plot','dir_out',dir_res_post};

            % dir_noise = [dir_subj 'Stim-processed' filesep];
            if exist(dir_noise,'dir')
                flags_in(end+1:end+2) = {'dir_noise',dir_noise};
            end
            % dir_target = [dir_subj 'speech-samples' filesep]; % not really used, but it removes one of the warnings...
            if exist(dir_target,'dir')
                flags_in(end+1:end+2) = {'dir_target',dir_target};
            end

            % GLM-based revcorr
            % % flags_in_here = flags_in;
            % % 
            % % bForce_dataload = 0;
            % % if bForce_dataload
            % %     flags_in_here{end+1} = 'force_dataload';
            % % end

            [ACI_t1,cfg_ACI_t1,results_t1]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t1',flags_in{:});
            [ACI_t2,cfg_ACI_t2,results_t2]    = fastACI_getACI(savegame_path,'glm','trialtype_analysis','t2',flags_in{:});

            f0_kernel_t1{i_cond}(:,i_subject) = ACI_t1(:,1);
            time_kernel_t1{i_cond}(:,i_subject) = ACI_t1(:,2);
            f0_kernel_t2{i_cond}(:,i_subject) = ACI_t2(:,1);
            time_kernel_t2{i_cond}(:,i_subject) = ACI_t2(:,2);
        end
        
        %%%
        % Plot targets with segments (replaced by 'fig1a')
        N_segment = size(f0_kernel_t1{i_cond},1)-1; % Prior knowledge
        t_edge = (0:N_segment)/10;
        XL = [0-x_offset 0.9+x_offset];
        
        if i_cond == 1
            figure('Position',[100 100 1000 750]); 
            tiledlayout(4,2,'TileSpacing','tight');
            
            x_panel_label = -.15;
            y_panel_label =  .92; % 1.15;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % f0 target1, target2:
        nexttile(1+2*(i_cond-1));
        
        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on
        
        f0_kernel_t1_here = f0_kernel_t1{i_cond};
        f0_kernel_t2_here = f0_kernel_t2{i_cond};
        Me1 = mean(f0_kernel_t1_here,2);
        errLU1 = 1.96*sem(f0_kernel_t1_here,[],2)/2;
        errorbar(t_edge,Me1,errLU1,'b','LineWidth',1); 
        errorbar(t_edge,mean(f0_kernel_t2_here,2),1.96*sem(f0_kernel_t2_here,[],2)/2,'r','LineWidth',1); 
        
        % Me = mean(f0_kernel,2);
        % errL = 1.96*sem(f0_kernel,[],2)/2;
        % errU = 1.96*sem(f0_kernel,[],2)/2;
        % errorbar(t_edge,Me,errL,errU,'k','LineWidth',2); hold on

        if i_cond == 1
            pref = 'Target-specific f_0 kernels';
            title(pref);
        else
            % pref = '';
            % title([pref 'f0 kernel (' cfg_game.response_names{1} '-' cfg_game.response_names{2} ')']); 
        end
        ylabel('Weight'); 
        
        if i_cond == 4
            xlabel('Segment edge (s)')
        else
            set(gca,'XTickLabel','');
        end
        
        set(gca,'XTick',t_edge)
        
        xlim(XL);
        grid on
        
        % if i_cond == 1
        text(x_panel_label,y_panel_label,sprintf('%s          %s - %s',lab{i_cond},cfg_game.response_names{1},cfg_game.response_names{2}),'FontWeight','Bold','FontSize',14,'Unit','Normalized');
        % end
        
        % text(.85,.9,['N = ' num2str(i_subject)],'Units','normalized')
        ylim(YLs(i_cond,:));
        YT = YLs(i_cond,1)+stepY(i_cond):stepY(i_cond):YLs(i_cond,2)-stepY(i_cond);
        set(gca,'YTick',YT);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time target1, target2:
        nexttile(2*(i_cond-1)+2);
        
        % Horizontal line:
        plot([min(t_edge) max(t_edge)],[0 0],'k--'); hold on
        
        time_kernel_t1_here = time_kernel_t1{i_cond};
        time_kernel_t2_here = time_kernel_t2{i_cond};
        errorbar(t_edge,mean(time_kernel_t1_here,2),1.96*sem(time_kernel_t1_here,[],2)/2,'b','LineWidth',1); 
        errorbar(t_edge,mean(time_kernel_t2_here,2),1.96*sem(time_kernel_t2_here,[],2)/2,'r','LineWidth',1); 
        
        if i_cond == 1
            pref = 'Target-specific time kernels';
            title(pref);
        end
        % ylabel('Weight'); 
        
        if i_cond == 4
            xlabel('Segment edge (s)')
        else
            set(gca,'XTickLabel','');
        end
        set(gca,'XTick',t_edge)
                
        set(gca,'YTickLabel','');
        ylabel('');
        
        xlim(XL);
        grid on
        
        text(.85,.9,['N = ' num2str(i_subject)],'Units','normalized')
        
        ylim(YLs(i_cond,:));
        set(gca,'YTick',YT);
    end
        
    h(end+1) = gcf;
    hname{end+1} = 'suppl-fig2-target-kernel';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.h = h;
data.hname = hname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function file_subj = il_get_filenames_cond(Cond, dir_where)

switch Cond
    case 'LAMI'
        fi1 = Get_filenames(dir_where,'*S00*');
        fi2 = Get_filenames(dir_where,'*S01*');
        file_subj = [fi1 fi2];
        
    case 'LAPEL'
        % SubjectFolders=[dir([dir_where '*S02*']); dir([dir_where '*S03*']); dir([dir_where '*S08*'])]; %LAPEL
        fi1 = Get_filenames(dir_where,'*S02*');
        fi2 = Get_filenames(dir_where,'*S03*');
        fi3 = Get_filenames(dir_where,'*S08*');
        file_subj = [fi1 fi2 fi3];
        
    case 'LACROCH' 
        % SubjectFolders=[dir([dir_where '*S04*']); dir([dir_where '*S06*'])]; %LACROCH
        fi1 = Get_filenames(dir_where,'*S04*');
        fi2 = Get_filenames(dir_where,'*S06*');
        file_subj = [fi1 fi2];
        
    case 'LALARM'
        file_subj = Get_filenames(dir_where,'*S05*');

    case 'LAMI_SHIFTED'
        file_subj = Get_filenames(dir_where,'*S07*'); % LAMI_shifted
        
end
file_subj = file_subj(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dir_path = il_get_data_path(path_type, bZenodo,dir_subj_name)

if bZenodo
    
else
    dir_data = fastACI_paths('dir_data');
end

switch path_type
    case 'dir_res'
        if bZenodo
            dir_path = '';
        else
            dir_path = '';
            
            dir_subj = [dir_where_exp dir_subj_name filesep];
            dir_res2 = [dir_subj 'Results' filesep];
        end
end
    
