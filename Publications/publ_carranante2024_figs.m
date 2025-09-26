function publ_carranante2024_figs(varargin)
% function publ_carranante2024_figs(varargin)
%
% % To display Fig. 1 of Carranante et al (2024) use :::
%     publ_carranante2024_figs('fig1');
%
% fig1: spectrogram of target sounds
% suppl_fig2: behavioral performances (and anovas)
% suppl_fig3: individual ACIs
% tab1: significant cues
% fig2: mean ACI for each condition
% fig3: crosspredictions
% suppl_fig1: audiograms
% suppl_fig4: difference of spectrograms
% suppl_fig5: simulated ACIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_carranante2024_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig1','suppl_fig2','suppl_fig3','fig2','fig3','tab1','suppl_fig1','suppl_fig2_suppl','suppl_fig3_suppl'};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

SubjectsABDA = {'S01','S08','S04','S06','S07','S12','S13'}; %abda
SubjectsABPA = {'S01','S04','S13','S20','S23','S33','S36'}; %abpa
SubjectsADGA = {'S01','S13','S16T','S18B','S19B','S21','S28'}; %adga
SubjectsADTA = {'S01','S13','S22','S24','S27','S39','S47'}; %adta
SubjectsAGKA = {'S01','S13','S32','S35','S38','S41','S42'}; %agka
SubjectsAKTA = {'S01','S13','S18B','S31','S34','S38B','S40B'}; %akta
SubjectsAPTA = {'S01B','S13','S16','S25B','S26','S28','S29'}; %apta
Subjects = unique([SubjectsABDA(:); SubjectsABPA(:); SubjectsADGA(:); SubjectsADTA(:); SubjectsAGKA(:); SubjectsAKTA(:); SubjectsAPTA(:)]);

Conditions = {'speechACI_Logatome-abda-S43M','speechACI_Logatome-adga-S43M','speechACI_Logatome-apta-S43M','speechACI_Logatome-akta-S43M','speechACI_Logatome-abpa-S43M','speechACI_Logatome-adta-S43M','speechACI_Logatome-agka-S43M'}; Condition_names = {'abda','adga','apta','akta','abpa','adta','agka'};
glmfct = 'l1glm';

N_subjects = length(Subjects);
max_N_subjects_condition = max([length(SubjectsABDA),length(SubjectsABPA),length(SubjectsADGA),length(SubjectsADTA),length(SubjectsAGKA),length(SubjectsAKTA),length(SubjectsAPTA)]);
masker = 'bumpv1p2_10dB';
N_conditions = length(Conditions);
conditioncolors = {'b','r','g','m','k','c','y'};
trialtype_analysis = 'total';
TF_type = 'gammatone';
expvar_after_reversal = 4;

resume_trial = [1,401,801,1201,1601,2001,2401,2801,3201,3601,4000];%resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
N_sessions  = length(resume_trial)-1;

flags_for_input = {TF_type,...
    'trialtype_analysis','total', ...'t1',...
    'N_folds', 10, ...
    'no_bias', ...
    'no_plot', ...
    'expvar_after_reversal',4, ...
    'add_signal',0, ...
    'apply_SNR',0, ...
    'skip_if_on_disk',1, ...
    'no_permutation', ...
    };

switch glmfct
    case 'l1glm'
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
        flags_for_input = [flags_for_input, ...
            'lambda',Lambdas, ...
            'pyramid_script','imgaussfilt', ...
            'pyramid_shape',-1, ...
            ];
    otherwise
        error('For this analysis, only the ''l1glm'' option was used.\n')
end

%dir_out = '';
%dir_out_figs = [pwd filesep];
fullpath = fastACI_paths('dir_data');


if flags.do_fig1 || flags.do_suppl_fig2_suppl
        
    path_aba = [fullpath 'speechACI_Logatome-abda-S43M\S01\speech-samples\S43M_ab_ba.wav'];
    [aba, fs] = audioread(path_aba);
    path_ada = [fullpath 'speechACI_Logatome-abda-S43M\S01\speech-samples\S43M_ad_da.wav'];
    [ada, fs] = audioread(path_ada);
    path_aga = [fullpath 'speechACI_Logatome-agka-S43M\S01\speech-samples\S43M_ag_ga.wav'];
    [aga, fs] = audioread(path_aga);
    path_apa = [fullpath 'speechACI_Logatome-abpa-S43M\S01\speech-samples\S43M_ap_pa.wav'];
    [apa, fs] = audioread(path_apa);
    path_ata = [fullpath 'speechACI_Logatome-adta-S43M\S01\speech-samples\S43M_at_ta.wav'];
    [ata, fs] = audioread(path_ata);
    path_aka = [fullpath 'speechACI_Logatome-agka-S43M\S01\speech-samples\S43M_ak_ka.wav'];
    [aka, fs] = audioread(path_aka);
    
    basef = 8000;
    flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5,'dboffset',100,'no_adt','binwidth',0.01, ...
        'no_outerear','no_middleear'};
    
    [G_aba, fc, t_spec, outs] = Gammatone_proc([aba], fs, flags_gamma{:});
    [G_ada, fc, t_spec, outs] = Gammatone_proc([ada], fs, flags_gamma{:});
    [G_aga, fc, t_spec, outs] = Gammatone_proc([aga], fs, flags_gamma{:});
    [G_apa, fc, t_spec, outs] = Gammatone_proc([apa], fs, flags_gamma{:});
    [G_ata, fc, t_spec, outs] = Gammatone_proc([ata], fs, flags_gamma{:});
    [G_aka, fc, t_spec, outs] = Gammatone_proc([aka], fs, flags_gamma{:});
    cfg_ACI.t = t_spec;
    cfg_ACI.f = fc;

    [~,color_aba,color_apa] = selectcolormap('speechACI_Logatome-abpa-S43M');
    [~,color_ada,color_ata] = selectcolormap('speechACI_Logatome-adta-S43M');
    [~,color_aga,color_aka] = selectcolormap('speechACI_Logatome-agka-S43M');
end

%% Fig 1. Stimuli

if flags.do_fig1
    h_fig1 = figure('Position', [100,100,700,450]); o=tiledlayout(2,3,'TileSpacing','compact')
    
    C_axis = [-9 -4.7];fontSize = 16;
    nexttile; affichage_tf(log(G_aba)', 'pow', t_spec, fc); caxis(C_axis); title('/aba/', 'Color', color_aba,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_aba, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); set(gca,'XTickLabels',{}); colorbar off
    nexttile; affichage_tf(log(G_ada)', 'pow', t_spec, fc); caxis(C_axis); title('/ada/', 'Color', color_ada,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_ada, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); set(gca,'XTickLabels',{}); ylabel(''); set(gca,'YTickLabels',{}); colorbar off
    nexttile; affichage_tf(log(G_aga)', 'pow', t_spec, fc); caxis(C_axis); title('/aga/', 'Color', color_aga,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_aga, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); set(gca,'XTickLabels',{}); ylabel(''); set(gca,'YTickLabels',{}); colorbar off
    nexttile; affichage_tf(log(G_apa)', 'pow', t_spec, fc); caxis(C_axis); title('/apa/', 'Color', color_apa,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_apa, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); colorbar off
    nexttile; affichage_tf(log(G_ata)', 'pow', t_spec, fc); caxis(C_axis); title('/ata/', 'Color', color_ata,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_ata, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); set(gca,'YTickLabels',{}); colorbar off
    nexttile; affichage_tf(log(G_aka)', 'pow', t_spec, fc); caxis(C_axis); title('/aka/', 'Color', color_aka,'FontSize', fontSize)
    hold on; affichage_tf_add_Praat_metrics_one_sound(path_aka, cfg_ACI, [], '-', 'w')
    xlabel(''); ylabel(''); set(gca,'YTickLabels',{}); colorbar off
    xlabel(o,'time (s)'); ylabel(o,'frequency (Hz)');
end
%% Fig 2. Behavioral performances and ANOVAs

if flags.do_suppl_fig2
    SNRlev = -20:1:-5;
    N_SNRlev = length(SNRlev);
    medianSNR = nan(N_conditions,max_N_subjects_condition,N_sessions);
    dprime = nan(N_conditions,max_N_subjects_condition,N_SNRlev);
    criterion = nan(N_conditions,max_N_subjects_condition,N_SNRlev);
    N_hist = nan(N_conditions,max_N_subjects_condition,N_SNRlev);
    
    hsuppl_fig2 = figure('Position', [100,100,800,400]); tiledlayout(2,2)
    
    for i_condition = 1:N_conditions
        
        experiment = Conditions{i_condition};
        colour_here = conditioncolors{i_condition};
        i_subject_in_condition = 1;
        
        for i_subject = 1:N_subjects
            
            subject = Subjects{i_subject};
            
            dir_subject = [fastACI_dir_data experiment filesep subject filesep];
            dir_results = [dir_subject 'Results' filesep];
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            
            if ~isempty(fname_results)
                fname_results = fname_results{end};
                fprintf(['Processing ' fname_results '\n'])
                load([dir_results fname_results]);
                cd(dir_results)
                
                n_responses = data_passation.n_responses;
                N_trials    = length(n_responses); % completed number of trials
                is_correct    = data_passation.is_correct;
                SNR         = data_passation.expvar;
                
                for i_session = 1:N_sessions
                    [rev,idxrev] = Get_mAFC_reversals(SNR(resume_trial(i_session):resume_trial(i_session+1)));
                    
                    medianSNR(i_condition,i_subject_in_condition,i_session) = median(rev(5:end));
                    PC(i_session) = mean(is_correct(idxrev(4):idxrev(end)));
                end
                meanPC(i_condition,i_subject_in_condition) = mean(PC);
                
                data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation,'histogram-group',0);
                %bin_centres = data_hist.bin_centres; N_bin = length(bin_centres);
                N_hist(i_condition,i_subject_in_condition,1:N_SNRlev) = data_hist.N_hist;
                P_H  = data_hist.H ./(data_hist.H +data_hist.M ); 
                P_FA = data_hist.FA./(data_hist.CR+data_hist.FA);

                dprime(i_condition,i_subject_in_condition,1:N_SNRlev) = norminv(P_H)-norminv(P_FA);           % Eq.  9 from Harvey2004
                criterion(i_condition,i_subject_in_condition,1:N_SNRlev) = -(norminv(P_H) + norminv(P_FA))/2; % Eq. 12 from Harvey2004
                
                i_subject_in_condition = i_subject_in_condition+1;
            end
        end

        % remove values if less than Nobs_min observations
        Nobs_min = 15;
        dprime(i_condition, :, sum(N_hist(i_condition,:,:)>=Nobs_min,2)<max_N_subjects_condition) = nan;
        criterion(i_condition, :, sum(N_hist(i_condition,:,:)>=Nobs_min,2)<max_N_subjects_condition) = nan;

        figure(hsuppl_fig2); nexttile(1,[1,2])
        lrev(i_condition) = errorbar(resume_trial(2:end)+30*i_condition,squeeze(nanmean(medianSNR(i_condition,:,:),2)),squeeze(std(medianSNR(i_condition,:,:),[],2,'omitnan')), [colour_here '-o']); hold on
        errorbar(resume_trial(end)+500+30*i_condition,squeeze(nanmean(median(medianSNR(i_condition,:,:),3),2)),squeeze(std(median(medianSNR(i_condition,:,:),3),[],2,'omitnan')), [colour_here 'o']); hold on
        
        nexttile([3])
        dprime(dprime==0)=nan;
        dprime(dprime==Inf)=nan;
        dprime(dprime==-Inf)=nan;
        errorbar(SNRlev+0.1*i_condition,squeeze(nanmean(dprime(i_condition,:,:),2))',squeeze(std(dprime(i_condition,:,:),[],2,'omitnan')),colour_here);
        hold on
        nexttile([4])
        criterion(criterion==0)=nan;
        criterion(criterion==Inf)=nan;
        criterion(criterion==-Inf)=nan;
        errorbar(SNRlev+0.1*i_condition,squeeze(nanmean(criterion(i_condition,:,:),2))',squeeze(std(criterion(i_condition,:,:),[],2,'omitnan')),colour_here);
        hold on
        %
        %         nexttile
        %         meanPC(meanPC==0)=nan;
        %         plot(i_condition,nanmean(meanPC(i_condition,:)),['o' colour_here]);
        %         errorbar(i_condition,nanmean(meanPC(i_condition,:)),std(meanPC(i_condition,:),[],2,'omitnan'),colour_here);
        %         hold on
    end
    
    figure(hsuppl_fig2); nexttile(1)
    xlim([0 5000]); xlabel('trial #'); ylabel('SNR [dB]'); grid on
    set(gca, 'XTick', [resume_trial 4500]);
    set(gca, 'XTickLabels', {'1' '401' '801' '1201' '1601' '2001' '2401' '2801' '3201' '3601' '4000' 'total'});
    legend(lrev,upper(Condition_names),'Location','NorthEastOutside')
    
    nexttile([3])
    xlim([-21 -6]); xlabel('SNR [dB]'); ylabel('d prime'); grid on
    
    nexttile([4])
    xlim([-21 -6]); xlabel('SNR [dB]'); ylabel('criterion'); grid on
    ylim([-0.8 0.8])
    %legend(Condition_names)
    
    %     figure(hPC)
    %     xlim([0.5 N_conditions+0.5]); set(gca, 'XTickLabel', Condition_names); xlabel('Condition')
    %     ylabel('% correct'); grid on

%     % remove values if less than Nobs_min observations
% Nobs_min = 10;
% dprime(i_condition, :, sum(N_hist(i_condition,:,:)>=Nobs_min,2)<max_N_subjects_condition) = nan;
% criterion(i_condition, :, sum(N_hist(i_condition,:,:)>=Nobs_min,2)<max_N_subjects_condition) = nan;
% dprime(:,:,SNRlev>=-13) = nan;
% criterion(:,:,SNRlev>=-13) = nan;

    % Simple ANOVA on SNR thresholds with factors Condition, Session and Condition*Session
    [idx_conditions, idx_subjects, idx_sessions] = ndgrid(1:N_conditions,1:max_N_subjects_condition,1:N_sessions);
    anovan(medianSNR(:), {idx_conditions(:), idx_sessions(:)},'model','interaction','varnames',{'Condition','Session'});
    
    % Simple ANOVA on dprime and criterion with factors Condition, SNR and Condition*SNR
        % remove values if less than Nobs_min observations
        %Nobs_min = 15;
    SNRlevANOVAidx = find(squeeze(sum(sum((N_hist>=Nobs_min),1),2))>=max_N_subjects_condition*N_conditions);
    %SNRlevANOVAidx = [min(SNRlevANOVAidx):max(SNRlevANOVAidx)];
    SNRlevANOVA = SNRlev(SNRlevANOVAidx); N_SNRlevANOVA = length(SNRlevANOVA);
    dprimeANOVA = dprime(:,:,SNRlevANOVAidx);
    criterionANOVA = criterion(:,:,SNRlevANOVAidx);
    [idx_conditions, idx_subjects, idx_SNRlev] = ndgrid(1:N_conditions,1:max_N_subjects_condition,1:N_SNRlevANOVA);
    anovan(dprimeANOVA(:), {idx_conditions(:), idx_SNRlev(:)},'model','interaction','varnames',{'Condition','SNR'});
    anovan(criterionANOVA(:), {idx_conditions(:), idx_SNRlev(:)},'model','interaction','varnames',{'Condition','SNR'});
    
end

%% Load ACI data

if flags.do_suppl_fig3 || flags.do_fig2 || flags.do_fig3 || flags.do_tab1
    
    for i_condition = 1:N_conditions
        
        experiment = Conditions{i_condition};
        colour_here = conditioncolors{i_condition};
        i_subject_in_condition = 1;
        
        for i_subject = 1:N_subjects
            
            subject = Subjects{i_subject};
            
            dir_subject = [fastACI_dir_data experiment filesep subject filesep];
            dir_results = [dir_subject 'Results' filesep];
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            
            if ~isempty(fname_results)
                fname_results = fname_results{end};
                fprintf(['Processing ' fname_results '\n'])
                load([dir_results fname_results]);
                cd(dir_results)
                
                cfg_games{i_condition,i_subject_in_condition} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
                
                Data_matrix = [];
                
                flags_here = {'dir_noise',cfg_games{i_condition,i_subject_in_condition}.dir_noise, ...
                    'dir_target',cfg_games{i_condition,i_subject_in_condition}.dir_target, ...
                    'skip_if_on_disk',1};
                
                switch glmfct
                    case {'lassoslow','lassoglmslow','l1glm'}
                        % calculate ACI lassoslow
                        [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:}, 'Data_matrix', Data_matrix);
                        
                        lambda = results{i_condition,i_subject_in_condition}.FitInfo.Lambda;
                        Dev_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.Dev_test./results{i_condition,i_subject_in_condition}.FitInfo.CV.TestSize;
                        PC_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.PC_test;
                        idxlambda(i_condition,i_subject_in_condition) = results{i_condition,i_subject_in_condition}.idxlambda;
                        Nfolds = results{i_condition,i_subject_in_condition}.FitInfo.CV.NumTestSets;
                        
                        if flags.do_suppl_fig3 || flags.do_tab1
                        % statistical thresholding 
                        % The thresholds are dependent on idxlambda. They
                        % were obtained by randomization (based on 1000
                        % permutations). 
                            B = results{i_condition,i_subject_in_condition}.B(:,results{i_condition,i_subject_in_condition}.idxlambda);
                            
                            switch idxlambda(i_condition,i_subject_in_condition)
                                case 12
                                    B(B<0.0240 & B>-0.0243) = 0;
                                case 13
                                    B(B<0.0142 & B>-0.0138) = 0;
                                case 14
                                    B(B<0.0043 & B>-0.0041) = 0;
                            end
                            [ACI] = Convert_lasso_B2ACI(B, cfg_ACI{i_condition,i_subject_in_condition}, 1, cfg_ACI{i_condition,i_subject_in_condition}.keyvals);
                        end
                    otherwise
                        % calculate ACI classic revcorr
                        [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
                end
                ACI_matrix(:,:,i_condition,i_subject_in_condition) = ACI;
                Subject_matrix{i_condition,i_subject_in_condition} = subject;
                i_subject_in_condition = i_subject_in_condition+1;
            end
        end
    end
end

%% Goodness-of-fit metrics

if flags.do_suppl_fig3 || flags.do_fig2 || flags.do_fig3 || flags.do_tab1
    for i_condition = 1:N_conditions
        DEVtemp = Dev_test(i_condition,:,:,:); DEVtemp(DEVtemp==0) = nan;
        PAtemp = PC_test(i_condition,:,:,:); PAtemp(PAtemp==0) = nan;%-PC_test(i_condition,:,end,:);
        for i_subject = 1:max_N_subjects_condition
            if idxlambda(i_condition,i_subject)>0
                DEVopt_temp = DEVtemp(1,i_subject,idxlambda(i_condition,i_subject),:);
                DEVnull_temp = DEVtemp(1,i_subject,end,:);
                DEVben(i_condition, i_subject, : ) = squeeze(DEVopt_temp-DEVnull_temp);
                
                PAopt_temp = PAtemp(1,i_subject,idxlambda(i_condition,i_subject),:);
                PAnull_temp = PAtemp(1,i_subject,end,:);
                PAben(i_condition, i_subject, : ) = squeeze(PAopt_temp-PAnull_temp);
            else
                DEVben(i_condition, i_subject, 1:10) = nan;
                PAben(i_condition, i_subject, 1:10) = nan;
            end
        end
    end
end

%% Fig 3. Individual ACIs

if flags.do_suppl_fig3
    
    for i_condition = 1:N_conditions
        
        hsuppl_fig3{i_condition} = figure('Position', [10,10,165,710]);
        tiledlayout(max_N_subjects_condition,1,'TileSpacing', 'compact');
        
        experiment = Conditions{i_condition};
        [cmap, color1, color2] = selectcolormap(experiment);
        
        for i_subject = 1:max_N_subjects_condition
            
            % plot ACI lassoslow, classic revcorr, lasso
            %idx_tile = i_condition+(i_subject-1)*N_conditions;
            figure(hsuppl_fig3{i_condition}); nexttile;%(idx_tile); %(i_subject,i_masker)
            affichage_tf(ACI_matrix(:,:,i_condition,i_subject),'CI', 'cfg', cfg_ACI{i_condition}, 'colorbar_map', cmap, 'NfrequencyTicks', 4); hold on
            %outs_aff = affichage_tf_add_Praat_metrics(cfg_games{i_condition}.dir_target,cfg_ACI{i_condition},[], {'-','-'},{color1,color2},1.5);
            %title(['ACI ' experiment(20:23) ', ' subject], 'interpreter','none');title('');
            
            text(0.95,0.95,Subject_matrix(i_condition,i_subject),'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')
            colorbar off
            %ylabel('');set(gca,'YTickLabels',[]);
            if ~(i_subject == max_N_subjects_condition)
                xlabel('');
                set(gca,'XTickLabels',[]);
            end
            if i_subject == 1
                title(upper(experiment(20:23)), 'interpreter','none');
            end
            
            % test if prediction accuracy is not significant
            if mean(DEVben(i_condition, i_subject, : ))+1.64*std(DEVben(i_condition, i_subject, : ))/sqrt(Nfolds) < 0
                ax = gca;
                ax.LineWidth = 1.5;
            end
        end
        
    end
    
end

% Number of significant weights in ROIs
if flags.do_tab1

    threshold_weight = 0.001;
    % define ROIs
    t_IVI = [0.2 0.25]; idx_IVI = find(cfg_ACI{1, 1}.t>=t_IVI(1) & cfg_ACI{1, 1}.t<=t_IVI(2));
    t_onset = [0.25 0.3]; idx_onset = find(cfg_ACI{1, 1}.t>=t_onset(1) & cfg_ACI{1, 1}.t<=t_onset(2));
    t_offset = [0.12 0.15]; idx_offset = find(cfg_ACI{1, 1}.t>=t_offset(1) & cfg_ACI{1, 1}.t<=t_offset(2));
    f_f0 = [100 170]; idx_f0 = find(cfg_ACI{1, 1}.f>=f_f0(1) & cfg_ACI{1, 1}.f<=f_f0(2));
    f_F1 = [350 550]; idx_F1 = find(cfg_ACI{1, 1}.f>=f_F1(1) & cfg_ACI{1, 1}.f<=f_F1(2));
    f_F2 = [1500 1900]; idx_F2 = find(cfg_ACI{1, 1}.f>=f_F2(1) & cfg_ACI{1, 1}.f<=f_F2(2));
    f_burst = [5000 7000]; idx_burst = find(cfg_ACI{1, 1}.f>=f_burst(1) & cfg_ACI{1, 1}.f<=f_burst(2));

    Ncue_F2onset = squeeze(any(any(abs(ACI_matrix(idx_F2,idx_onset,:,:))>threshold_weight,1),2));
    Ncue_F1onset = squeeze(any(any(abs(ACI_matrix(idx_F1,idx_onset,:,:))>threshold_weight,1),2));
    Ncue_f0IVI = squeeze(any(any(abs(ACI_matrix(idx_f0,idx_IVI,:,:))>threshold_weight,1),2));
    Ncue_F2offset = squeeze(any(any(abs(ACI_matrix(idx_F2,idx_offset,:,:))>threshold_weight,1),2));
    Ncue_F1offset = squeeze(any(any(abs(ACI_matrix(idx_F1,idx_offset,:,:))>threshold_weight,1),2));
    Ncue_burstIVI = squeeze(any(any(abs(ACI_matrix(idx_burst,idx_IVI,:,:))>threshold_weight,1),2));


fprintf(['\n*Table 1: Number of ACIs for which a given cue was found significant, in each condition*\n'])
fprintf(['(this table considers all ACIs, including non-significant ones)\n\n'])

fprintf(['\t\t\t' 'BD DG PT KT BP DT GK' '\n'])
fprintf(['F2 onset\t' num2str(sum(Ncue_F2onset,2)') '\n'])
fprintf(['F1 onset\t' num2str(sum(Ncue_F1onset,2)') '\n'])
fprintf(['burst IVI\t' num2str(sum(Ncue_burstIVI,2)') '\n'])
fprintf(['f0 IVI\t\t' num2str(sum(Ncue_f0IVI,2)') '\n'])
fprintf(['F2 offset\t' num2str(sum(Ncue_F2offset,2)') '\n'])
fprintf(['F1 offset\t' num2str(sum(Ncue_F1offset,2)') '\n'])
end

%% Fig 4. Mean ACIs

if flags.do_fig2
    
    for i_condition = 1:N_conditions
        hfig2{i_condition} = figure('Position', [100,100,400,300]);
        experiment = Conditions{i_condition};
        [cmap, color1, color2] = selectcolormap(experiment);
        
        affichage_tf(mean(ACI_matrix(:,:,i_condition,:),4),'CI', 'cfg', cfg_ACI{i_condition}, 'colorbar_map', cmap); hold on
        
        outs_aff = affichage_tf_add_Praat_metrics(cfg_games{i_condition}.dir_target,cfg_ACI{i_condition},[], {'-','-'},{color1,color2},1.5);
        
        title(upper(experiment(20:23)), 'interpreter','none');
        %cfg_ACI{i_condition}.keyvals.trialtype_analysis
        
        L = legend;
        L.String{1} = ['a' L.String{1}];
        L.String{2} = ['a' L.String{2}];
    end
end

%% Supplementary Fig 2. Diference of spectrograms

if flags.do_suppl_fig2_suppl

    for i_condition = 1:N_conditions
        switch i_condition
            case 1
                cmap = 'DG_red_blue';
                spec_diff = G_aba'-G_ada';
            case 2
                cmap = 'DG_blue_yellow';
                spec_diff = G_ada'-G_aga';
            case 3
                cmap = 'DG_red_green';
                spec_diff = G_apa'-G_ata';
            case 4
                cmap = 'DG_green_purple';
                spec_diff = G_aka'-G_ata';
            case 5
                cmap = 'DG_blue_purple';
                spec_diff = G_aba'-G_apa';
            case 6
                cmap = 'DG_cyan_purple';
                spec_diff = G_ada'-G_ata';
            case 7
                cmap = 'DG_yellow_cyan';
                spec_diff = G_aga'-G_aka';
        end

        hsuppl_fig2_suppl{i_condition} = figure('Position', [100,100,400,300]);
        experiment = Conditions{i_condition};
        [cmap, color1, color2] = selectcolormap(experiment);

        affichage_tf(spec_diff,'CI', 'cfg', cfg_ACI, 'colorbar_map', cmap); hold on

        %outs_aff = affichage_tf_add_Praat_metrics(cfg_games{i_condition}.dir_target,cfg_ACI);%,[], {'-','-'},{color1,color2},1.5);

        title(upper(experiment(20:23)), 'interpreter','none');
        %cfg_ACI{i_condition}.keyvals.trialtype_analysis

    end
end

%% Fig. 5 Auto- and cross-predictions
% only if glmfct = 'l1glm';

if flags.do_fig3
    Dev_test_cross = [];
    for i_condition = 1:N_conditions
        
        experiment = Conditions{i_condition};
        i_subject_in_condition = 1;
        
        for i_subject = 1:max_N_subjects_condition %N_subjects
            
            subject = Subject_matrix{i_condition,i_subject};
            
            dir_subject = [fastACI_dir_data experiment filesep subject filesep];
            dir_results = [dir_subject 'Results' filesep];
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            %fname_crosspred = [dir_subject 'Results' filesep 'Results_ACI' filesep 'ACI-' subject '-speechACI_Logatome-' masker '-nob-gt-l1glm-rev4.mat'];
            
            if ~isempty(fname_results)
                fname_results = fname_results{end};
                fprintf(['Processing ' fname_results '\n'])
                load([dir_results fname_results]);
                cd(dir_results)
                
                %             cfg_games{i_condition,i_subject_in_condition} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
                %
                %             Data_matrix = [];
                %
                %             flags_for_input = {TF_type, ...
                %                 'trialtype_analysis', trialtype_analysis,...
                %                 'N_folds', 10, ...
                %                 'dir_noise',cfg_games{i_condition,i_subject_in_condition}.dir_noise, ...
                %                 'dir_target',cfg_games{i_condition,i_subject_in_condition}.dir_target, ...
                %                 'add_signal',0, ...
                %                 'apply_SNR',0, ...
                %                 'skip_if_on_disk',1, ...
                %                 'no_permutation', ...
                %                 'no_bias', ...
                %                 'no_plot', ... 'idx_trialselect',idx_trialselect,
                %                 'expvar_after_reversal',expvar_after_reversal, ...
                %                 'pyramid_shape',-1, ...
                %                 'lambda',Lambdas ...
                %                 };
                
                flags_here = {'dir_noise',cfg_games{i_condition,i_subject_in_condition}.dir_noise, ...
                    'dir_target',cfg_games{i_condition,i_subject_in_condition}.dir_target, ...
                    'skip_if_on_disk',1};
                
                %%% 1. Obtaining the ACI of the current subject:
                [col1,cfg_ACI_here,res,Data_matrix_here] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:});
                crosspred = [];
                if isempty(Data_matrix_here)
                    % It means that Data_matrix_here was not loaded, maybe crosspred exists:
                    [fname_dir,fname] = fileparts(cfg_ACI_here.fnameACI);
                    
                    % 2.1. Checking whether you previously stored already the cross-prediction data:
                    clear fname_crosspred
                    suffix = '';
                    for ii_subject = 1:max_N_subjects_condition %N_subjects%
                        % TOFIX: problem with subject selection?
                        dir_subject_crosspred = [fastACI_dir_data experiment filesep Subject_matrix{i_condition,ii_subject} filesep];
                        fname_crosspred{ii_subject} = [dir_subject_crosspred 'Results' filesep 'Results_ACI' filesep 'ACI-' Subject_matrix{i_condition,ii_subject} '-speechACI_Logatome-' masker '-nob-gt-l1glm+pyrga-rev4.mat'];
                    end
                    crosspred = Check_if_crosspred(cfg_ACI_here.fnameACI,fname_crosspred,suffix);
                end
                if isempty(crosspred)
                    % 2.2. If there is no cross-prediction in your computer, then
                    %      it obtains it
                    flags = {'ACI_crosspred',fname_crosspred,'Data_matrix',Data_matrix_here}; % the extra fields
                    [col1,cfg_ACI_here,res,col4] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:}, flags{:});
                    crosspred = Check_if_crosspred(cfg_ACI_here.fnameACI,fname_crosspred,suffix);
                end
                
                % fill the crossprediction matrix
                for ii_subject = 1:max_N_subjects_condition%N_subjects%
                    % data from i_subject, ACI from ii_subject
                    % dir_subject_crosspred = [fastACI_dir_data experiment filesep Subject_matrix{i_condition,ii_subject} filesep];
                    Dev_test_cross(i_condition,i_subject,ii_subject) = mean(mean(crosspred(ii_subject).Dev_test_t(idxlambda(i_condition,ii_subject),:,:)-crosspred(ii_subject).Dev_test_t(end,:,:),2),3);
                    DEVben_cross(i_condition,i_subject,ii_subject,:) = mean(crosspred(ii_subject).Dev_test_t(idxlambda(i_condition,ii_subject),:,:)-crosspred(ii_subject).Dev_test_t(end,:,:),3);
                    PC_test_cross(i_condition,i_subject,ii_subject) = mean(mean(crosspred(ii_subject).PC_test_t(idxlambda(i_condition,ii_subject),:,:)-crosspred(ii_subject).PC_test_t(end,:,:),2),3);
                end
                i_subject_in_condition=i_subject_in_condition+1;
            end
        end
        %
        %     figure(hcrossPC{i_condition});
        %     imagesc(squeeze(2*PC_test_cross)); colorbar; caxis([0 0.2])
        %     title(['crosspred PA, ' Condition_names{i_condition}])
        %     set(gca, 'XTick', 1:max_N_subjects_condition); set(gca, 'XTickLabels', Subjects)
        %     set(gca, 'YTick', 1:max_N_subjects_condition); set(gca, 'YTickLabels', Subjects)
        %     xlabel('data from');     ylabel('ACI from')
        
%         figure%(hcrossPC{i_condition});
%         imagesc(squeeze(Dev_test_cross(i_condition,:,:))); colorbar; caxis([min(min(min(Dev_test_cross))) 0])
%         title(['crosspred deviance, ' Condition_names{i_condition}])
%         set(gca, 'XTick', 1:max_N_subjects_condition); set(gca, 'XTickLabels', Subject_matrix(i_condition,:))
%         set(gca, 'YTick', 1:max_N_subjects_condition); set(gca, 'YTickLabels', Subject_matrix(i_condition,:))
%         xlabel('data from');     ylabel('ACI from')
%         
    end

figure('Position', [100,100,1000,300]);
o = tiledlayout(1,N_conditions,'TileSpacing', 'compact');
for i_condition = 1:N_conditions
    nexttile
    Dev_test_cross_diag = diag(squeeze(Dev_test_cross(i_condition,:,:)));
    errorbar(1:7, Dev_test_cross_diag, std(DEVben(i_condition,:,:),[],3)/sqrt(Nfolds),'o'); hold on
    plot([0,8],[0,0],'k--');
    Dev_test_cross_nodiag = squeeze(Dev_test_cross(i_condition,:,:)) - diag(Dev_test_cross_diag);
    plot(1:7, sum(Dev_test_cross_nodiag)/length(Dev_test_cross_nodiag),'.'); hold on
    xlim([0 max_N_subjects_condition+1]); set(gca, 'XTick', 1:max_N_subjects_condition); set(gca, 'XTickLabels', Subject_matrix(i_condition,:))
    title(upper(Condition_names{i_condition}))
    ylim([-0.1, 0.04])
    if i_condition==1
        ylabel('cross-validated deviance');
    else
        set(gca, 'YTickLabels', {})
    end        
end
xlabel(o, 'Subject #')
end

if flags.do_suppl_fig1 
 
    %%% figure format:
    FontSize = 12;
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
    tiledlayout(1,2,'TileSpacing','tight'); % 'tight'); % none'); % 'Compact');
    
    f = [250 500 1000 2000 4000 8000];
    aud = [-5,0,0,10,15,30,5,10,30,30,40,40;10,10,5,5,0,0,15,15,10,10,5,0;5,10,0,0,10,5,10,10,5,0,15,0;10,10,5,5,10,-10,5,10,0,0,10,0;15,15,0,5,5,5,10,10,0,5,10,0;5,5,0,0,0,-5,10,5,0,0,0,0;15,10,0,0,-5,10,10,10,5,-5,-5,10;20,20,15,20,15,15,20,25,15,15,5,10;10,15,5,5,0,5,10,10,5,5,-10,15;10,5,5,0,5,15,15,10,0,5,5,5;5,10,5,-10,-5,5,10,10,10,-5,5,-5;10,10,10,5,5,15,10,5,10,10,0,15;15,10,15,10,10,0,5,10,15,15,0,-5;5,5,0,0,0,-5,10,5,0,0,10,0;5,5,5,0,0,-5,5,10,0,5,5,5;10,15,0,5,0,10,15,10,0,-5,0,-5;10,5,5,0,0,0,10,5,0,-5,-5,0;10,5,0,5,5,-5,10,10,5,0,0,0;10,5,5,-5,0,5,5,10,5,0,5,5;5,5,0,-5,5,5,5,-5,0,-10,5,-10;15,10,5,10,10,0,10,5,5,0,10,0;15,15,0,0,0,0,15,10,0,0,15,15;0,5,0,10,5,10,5,5,0,5,10,5;20,20,20,15,15,10,20,20,15,15,20,5;5,5,0,10,15,15,5,10,5,10,20,20;15,10,5,5,20,10,5,10,10,-5,10,15;5,5,5,5,15,10,5,10,10,10,15,15;10,0,-5,0,-5,5,10,5,0,5,0,15;10,10,5,0,5,15,10,5,10,5,10,5;10,10,5,15,5,5,10,10,5,10,5,10;10,10,5,5,5,10,10,10,5,5,0,10;10,10,5,5,15,5,15,10,-5,15,25,5;0,0,0,5,0,5,5,5,0,0,0,5;5,0,-10,-10,-10,-10,5,-10,-10,-10,-10,-10;0,0,0,0,0,0,0,0,0,0,0,0];
    Subjects = {'S01','S03','S04','S05','S06','S07','S12','S13','S16','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S31','S32','S33','S34','S35','S36','S37','S38','S39','S40','S41','S42','S43','S47'};    
    aud_l = aud(:,1:6);
    aud_r = aud(:,7:12);


    Colours = distinguishable_colors(N_subjects);

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

    nexttile(1)
    idx2label = [];
    for i = 1:length(Subjects)
        if aud_best_ear(i)==1
            LW = 1;%3;
            FaceColour = Colours(i,:);%Colours{i};
            idx2label(end+1) = i;
            Style = '-';
        else
            LW = 0.5;%1.5;
            FaceColour = 'w';
            Style = '--';
        end
        hpl(i) = semilogx(f+i*0.01*f-0.01*f*length(Subjects)/2,-aud_l(i,:),'-','Color',Colours(i,:),'LineStyle',Style,'LineWidth',LW); grid on, hold on;
        plot(XT(end),-aud_avg(i,1),'o','Color',Colours(i,:),'MarkerFaceColor',FaceColour,'LineWidth',LW,'MarkerSize',7);
    end
    plot(f,-aud_Me_l,'ko-','MarkerFaceColor','k','LineWidth',2);
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',[]);
    text(0.03,0.9,'A. Left ear','Units','Normalized','FontWeight','bold','FontSize',14);
    set(gca,'XTickLabel',XTL);
    Ylabel(Ylab,FontSize);
    Xlabel(Xlab,FontSize);

    ylim([-47 12]);
    xlim([200 24000])
    
    YT = -60:5:10;
    YT_label = -YT;
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YT_label);
 
    legend(hpl(idx2label),Subjects(idx2label),'Location','SouthWest','NumColumns',3);
    
    nexttile(2)
    
    idx2label = [];
    for i = 1:length(Subjects)
        if aud_best_ear(i)==2
            LW = 1;%3;
            FaceColour = Colours(i,:);
            idx2label(end+1) = i;
            Style = '-';
        else
            LW = 0.5;%1.5;
            FaceColour = 'w';
            Style = '--';
        end
        hpl(i) = semilogx(f+i*0.01*f-0.01*f*length(Subjects)/2,-aud_r(i,:),'-','Color',Colours(i,:),'LineStyle',Style,'LineWidth',LW); grid on; hold on;
        plot(XT(end),-aud_avg(i,2),'o','Color',Colours(i,:),'MarkerFaceColor',FaceColour,'LineWidth',LW,'MarkerSize',7);
    end
    plot(f,-aud_Me_r,'ko-','MarkerFaceColor','k','LineWidth',2);
    set(gca,'XTick',XT);
    set(gca,'XTickLabel',[]);
    text(0.03,0.9,'B. Right ear','Units','Normalized','FontWeight','bold','FontSize',14);
    set(gca,'XTickLabel',XTL);
    %Ylabel(Ylab,FontSize);
    Xlabel(Xlab,FontSize);

    ylim([-47 12]);
    xlim([200 24000])
    
    YT = -60:5:10;
    YT_label = [];%-YT;
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YT_label);
    
    legend(hpl(idx2label),Subjects(idx2label),'Location','SouthWest','NumColumns',3);
    
    % ha_here(end+1) = gca;
    h(end+1) = gcf; % current handle
	hname{end+1} = [prefix 'PTA'];  
    
    Pos = get(gcf,'Position');
    Pos(3) = 1000;
    Pos(4) = 500;
    set(gcf,'Position',Pos); % '570         754
    


% %% Load ACI data
% 
% if flags.do_suppl_fig3 || flags.do_fig2 || flags.do_fig3 || flags.do_tab1
%     
%     for i_condition = 1:N_conditions
%         
%         experiment = Conditions{i_condition};
%         colour_here = conditioncolors{i_condition};
%         i_subject_in_condition = 1;
%         
%         for i_subject = 1:N_subjects
%             
%             subject = Subjects{i_subject};
%             
%             dir_subject = [fastACI_dir_data experiment filesep subject filesep];
%             dir_results = [dir_subject 'Results' filesep];
%             
%             fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
%             
%             if ~isempty(fname_results)
%                 fname_results = fname_results{end};
%                 fprintf(['Processing ' fname_results '\n'])
%                 load([dir_results fname_results]);
%                 cd(dir_results)
%                 
%                 cfg_games{i_condition,i_subject_in_condition} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
%                 
%                 Data_matrix = [];
%                 
%                 flags_here = {'dir_noise',cfg_games{i_condition,i_subject_in_condition}.dir_noise, ...
%                     'dir_target',cfg_games{i_condition,i_subject_in_condition}.dir_target, ...
%                     'skip_if_on_disk',1};
%                 
%                 switch glmfct
%                     case {'lassoslow','lassoglmslow','l1glm'}
%                         % calculate ACI lassoslow
%                         [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:}, 'Data_matrix', Data_matrix);
%                         
%                         lambda = results{i_condition,i_subject_in_condition}.FitInfo.Lambda;
%                         Dev_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.Dev_test./results{i_condition,i_subject_in_condition}.FitInfo.CV.TestSize;
%                         PC_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.PC_test;
%                         idxlambda(i_condition,i_subject_in_condition) = results{i_condition,i_subject_in_condition}.idxlambda;
%                         Nfolds = results{i_condition,i_subject_in_condition}.FitInfo.CV.NumTestSets;
%                         
%                         if flags.do_suppl_fig3
%                         % statistical thresholding 
%                         % The thresholds are dependent on idxlambda. They
%                         % were obtained by randomization (based on 1000
%                         % permutations). 
%                             B = results{i_condition,i_subject_in_condition}.B(:,results{i_condition,i_subject_in_condition}.idxlambda);
%                             
%                             switch idxlambda(i_condition,i_subject_in_condition)
%                                 case 12
%                                     B(B<0.0240 & B>-0.0243) = 0;
%                                 case 13
%                                     B(B<0.0142 & B>-0.0138) = 0;
%                                 case 14
%                                     B(B<0.0043 & B>-0.0041) = 0;
%                             end
%                             [ACI] = Convert_lasso_B2ACI(B, cfg_ACI{i_condition,i_subject_in_condition}, 1, cfg_ACI{i_condition,i_subject_in_condition}.keyvals);
%                         end
%                     otherwise
%                         % calculate ACI classic revcorr
%                         [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
%                 end
%                 ACI_matrix(:,:,i_condition,i_subject_in_condition) = ACI;
%                 Subject_matrix{i_condition,i_subject_in_condition} = subject;
%                 i_subject_in_condition = i_subject_in_condition+1;
%             end
%         end
%     end
% end

    
end


%% Load simulated ACI data

if flags.do_suppl_fig3_suppl

    for i_condition = 1:N_conditions
        
        experiment = Conditions{i_condition};
        colour_here = conditioncolors{i_condition};
        i_subject_in_condition = 1;
        
        %for i_subject = 1:N_subjects
            
            subject = 'osses2022a'; % Subjects{i_subject};
            
            dir_subject = [fastACI_dir_data experiment filesep subject filesep];
            dir_results = [dir_subject 'Results' filesep];
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            
            if ~isempty(fname_results)
                fname_results = fname_results{end};
                fprintf(['Processing ' fname_results '\n'])
                load([dir_results fname_results]);
                cd(dir_results)
                
                cfg_games{i_condition,i_subject_in_condition} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
                
                Data_matrix = [];
                
                flags_here = {'dir_noise',cfg_games{i_condition,i_subject_in_condition}.dir_noise, ...
                    'dir_target',cfg_games{i_condition,i_subject_in_condition}.dir_target, ...
                    'skip_if_on_disk',1};
                
                switch glmfct
                    case {'lassoslow','lassoglmslow','l1glm'}
                        % calculate ACI lassoslow
                        [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:}, 'Data_matrix', Data_matrix);
                        
                        lambda = results{i_condition,i_subject_in_condition}.FitInfo.Lambda;
                        Dev_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.Dev_test./results{i_condition,i_subject_in_condition}.FitInfo.CV.TestSize;
                        PC_test(i_condition,i_subject_in_condition,:,:) = results{i_condition,i_subject_in_condition}.FitInfo.PC_test;
                        idxlambda(i_condition,i_subject_in_condition) = results{i_condition,i_subject_in_condition}.idxlambda;
                        Nfolds = results{i_condition,i_subject_in_condition}.FitInfo.CV.NumTestSets;
                        
                        if flags.do_suppl_fig3
                        % statistical thresholding 
                        % The thresholds are dependent on idxlambda. They
                        % were obtained by randomization (based on 1000
                        % permutations). 
                            B = results{i_condition,i_subject_in_condition}.B(:,results{i_condition,i_subject_in_condition}.idxlambda);
                            
                            switch idxlambda(i_condition,i_subject_in_condition)
                                case 12
                                    B(B<0.0240 & B>-0.0243) = 0;
                                case 13
                                    B(B<0.0142 & B>-0.0138) = 0;
                                case 14
                                    B(B<0.0043 & B>-0.0041) = 0;
                            end
                            [ACI] = Convert_lasso_B2ACI(B, cfg_ACI{i_condition,i_subject_in_condition}, 1, cfg_ACI{i_condition,i_subject_in_condition}.keyvals);
                        end
                    otherwise
                        % calculate ACI classic revcorr
                        [ACI,cfg_ACI{i_condition,i_subject_in_condition},results{i_condition,i_subject_in_condition}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
                end
                ACI_matrix(:,:,i_condition,i_subject_in_condition) = ACI;
                Subject_matrix{i_condition,i_subject_in_condition} = subject;
               % i_subject_in_condition = i_subject_in_condition+1;
            end
%        end
    end

    for i_condition = 1:N_conditions
        hsuppl_fig3_suppl{i_condition} = figure('Position', [100,100,400,300]);
        experiment = Conditions{i_condition};
        [cmap, color1, color2] = selectcolormap(experiment);
        
        affichage_tf(mean(ACI_matrix(:,:,i_condition,:),4),'CI', 'cfg', cfg_ACI{i_condition}, 'colorbar_map', cmap); hold on
        
        outs_aff = affichage_tf_add_Praat_metrics(cfg_games{i_condition}.dir_target,cfg_ACI{i_condition},[], {'-','-'},{color1,color2},1.5);
        
        title(upper(experiment(20:23)), 'interpreter','none');
        %cfg_ACI{i_condition}.keyvals.trialtype_analysis
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handle_tl = il_tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin)

handle_tl = [];
if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout','file'); % tiledlayout.p
if bExist
    if isempty(varargin)
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option);
    else
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin{:});
    end
else
    warning('We programmed this script to use more recent graphic options from MATLAB. It might be that you won''t be able to visualise these results as we have foreseen...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax = il_nexttile(N)

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

if nargout ~= 0
    ax = gcf;
end

function [cmap,color1,color2] = selectcolormap(experiment)
switch experiment
    case 'speechACI_Logatome-abda-S43M'
        cmap = 'DG_red_blue';
        color1 = [0.6,0,0];
        color2 = [0,0,0.6];
    case 'speechACI_Logatome-adga-S43M'
        cmap = 'DG_blue_yellow';
        color1 = [0,0,0.6];
        color2 = [0.6,0.6,0];
    case 'speechACI_Logatome-abpa-S43M'
        cmap = 'DG_red_green';
        color1 = [0.6,0,0];
        color2 = [0,0.6,0];
    case 'speechACI_Logatome-apta-S43M'
        cmap = 'DG_green_purple';
        color1 = [0,0.6,0];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-adta-S43M'
        cmap = 'DG_blue_purple';
        color1 = [0,0,0.6];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-akta-S43M'
        cmap = 'DG_cyan_purple';
        color1 = [0,0.6,0.6];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-agka-S43M'
        cmap = 'DG_yellow_cyan';
        color1 = [0.6,0.6,0];
        color2 = [0,0.6,0.6];
end

function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);
% Copyright 2010-2011 by Timothy E. Holy
  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end
  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
%end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
%end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
%end
