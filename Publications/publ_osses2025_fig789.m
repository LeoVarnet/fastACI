function [] = publ_osses2025_fig789(FigureName)
% function publ_osses2025_fig789(FigureName)
%%%%% Figures 7, 8 and 9 from fastACI paper, section 6.4 and 6.5 %%%%%
%
% % To display Fig. 7 of Osses et al use
%     publ_osses2025_fig789('fig7');
% % To display Fig. 8 of Osses et al use
%     publ_osses2025_fig789('fig8');
% % To display Fig. 9 of Osses et al use
%     publ_osses2025_fig789('fig9');
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FigureName
    case 'fig7'
        fig = 3;
    case 'fig8'
        fig = 2;
    case 'fig9'
        fig = 1;
end

%% parameters
h = []; hname = [];

Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'};% {'S01','S03','S04','S06','S07','S12','S13'}; %abda
N_subjects = length(Subjects);

switch fig 
    case 1
        GLMfcts =  {'classic_revcorr','l1glm','glmfitqp'};%
        GLMfcts_titles = {'correlation','glm_L1_GB','glm_L2'};
    case 2
        GLMfcts =  {'classic_revcorr'};%
        GLMfcts_titles = {'correlation'};
    case 3
        GLMfcts =  {'l1glm'};%
        GLMfcts_titles = {'glm_L1_GB'};
end
N_glmfcts = length(GLMfcts);

switch fig 
    case {1,2,3}
        Maskers = {'white','bumpv1p2_10dB','sMPSv1p3'};
        Maskers_titles = {'WN','BN','MPSN'};
end
N_maskers = length(Maskers);

% Conditions = {'speechACI_Logatome-abda-S43M','speechACI_Logatome-adga-S43M','speechACI_Logatome-apta-S43M','speechACI_Logatome-abpa-S43M','speechACI_Logatome-adta-S43M','speechACI_Logatome-agka-S43M','speechACI_Logatome-akta-S43M'};Condition_names = {'abda','adga','apta','abpa','adta','agka','akta'};
 experiment = 'speechACI_Logatome-abda-S43M';Condition_names = {'abda'};
% Conditions = {'speechACI_Logatome-adga-S43M'};Condition_names = {'adga'};
% experiment = 'speechACI_Logatome-abpa-S43M';Condition_names = {'abpa'};
% Conditions = {'speechACI_Logatome-apta-S43M'};Condition_names = {'apta'};
% Conditions = {'speechACI_Logatome-akta-S43M'};Condition_names = {'akta'};
% Conditions = {'speechACI_Logatome-adta-S43M'};Condition_names = {'adta'};
%experiment = 'speechACI_Logatome-agka-S43M';Condition_names = {'agka'};

trialtype_analysis = 'total';
TF_type = 'gammatone';
expvar_after_reversal = 4;

%% ACIs (and load ACI data)

clear cfg_games cfg_ACI results Dev_test PC_test results ACI_matrix

switch fig 
    case {1,2}
    N_conditions = 10;
    v_condition = 1:N_conditions;
    case 3
    N_conditions = 1;
    v_condition = 10;
end

randorder_trial = randperm(4000);

for i_glmfct = 1:N_glmfcts
    glmfct = GLMfcts{i_glmfct};%'l1glm';%'classic_revcorr';%

    switch glmfct
        case {'l1lm', 'l1glm'}
            N_lambda = 30;
            Lambdas = logspace(-4, -1, N_lambda);
            idx = find(Lambdas >= 10^-3);
            Lambdas = Lambdas(idx);
            glmflags = {'pyramid_shape',-1,'lambda',Lambdas};
        case {'glmfitqp'}
            glmflags = {'lambda0',1e3,'stepsize',10};
        otherwise
            glmflags = {};
    end

    for i_masker = 1:N_maskers
        masker = Maskers{i_masker};%'bumpv1p2_10dB';%'white';%'sMPSv1p3';%

        % 
        % hACIind = figure; hold on;
        % h(end+1) = hACIind;
        % hname{end+1} = ['ACI ' masker ' ' glmfct];
        % t = tiledlayout(N_subjects,N_conditions,'TileSpacing', 'tight');
        % title(t, hname{end}, 'Interpreter', 'none')

        for i_subject = 1:N_subjects
            for i_condition = v_condition
                Data_matrix = [];
                Ntrials = 400*i_condition;
                idx_trialselect = randorder_trial(1:Ntrials);
                %idx_trialselect = zeros(1,4000); idx_trialselect(trialselect) = 1;


                subject = Subjects{i_subject};

                dir_subject = [fastACI_dir_data experiment filesep subject filesep];
                dir_results = [dir_subject 'Results' filesep];

                fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);

                if ~isempty(fname_results)
                    fname_results = fname_results{end};
                    fprintf(['Processing ' fname_results '\n'])
                    load([dir_results fname_results]);
                    cd(dir_results)

                    cfg_games{i_condition,i_subject,i_masker,i_glmfct} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders

                    flags_for_input = {TF_type,...
                        'trialtype_analysis', trialtype_analysis,...
                        'idx_trialselect',idx_trialselect,...
                        'N_folds', 10, ...
                        'dir_noise',cfg_games{i_condition,i_subject,i_masker,i_glmfct}.dir_noise, ...
                        'dir_target',cfg_games{i_condition,i_subject,i_masker,i_glmfct}.dir_target, ...
                        'add_signal',0, ...
                        'apply_SNR',0, ...
                        'skip_if_on_disk',1, ...
                        'no_permutation', ...
                        'no_bias', ...
                        'no_plot', ... 'idx_trialselect',idx_trialselect,
                        'expvar_after_reversal',expvar_after_reversal, ...
                        };
                    flags_for_input = [flags_for_input glmflags];

                    switch glmfct
                        case {'lassoslow','lassoglmslow','l1glm'}
                            % calculate ACI lassoslow
                            [ACI,cfg_ACI{i_condition,i_subject,i_masker,i_glmfct},results{i_condition,i_subject,i_masker,i_glmfct}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

                            lambda = results{i_condition,i_subject,i_masker,i_glmfct}.FitInfo.Lambda;
                            Dev_test(i_condition,i_subject,i_masker,i_glmfct,:,:) = results{i_condition,i_subject,i_masker,i_glmfct}.FitInfo.Dev_test;
                            PC_test(i_condition,i_subject,i_masker,i_glmfct,:,:) = results{i_condition,i_subject,i_masker,i_glmfct}.FitInfo.PC_test;
                            idxlambda(i_condition,i_subject,i_masker,i_glmfct) = results{i_condition,i_subject,i_masker,i_glmfct}.idxlambda;
                            Nfolds = results{i_condition,i_subject,i_masker,i_glmfct}.FitInfo.CV.NumTestSets;
                            %statistical threshold (should be dependent on idxlambda)
                            B = results{i_condition,i_subject,i_masker,i_glmfct}.B(:,results{i_condition,i_subject,i_masker,i_glmfct}.idxlambda);
                            %
                            %                     switch idxlambda(i_condition,i_subject,i_masker,i_glmfct)
                            %                         case 12
                            %                             B(B<0.0240 & B>-0.0243) = 0;
                            %                         case 13
                            %                             B(B<0.0142 & B>-0.0138) = 0;
                            %                         case 14
                            %                             B(B<0.0043 & B>-0.0041) = 0;
                            %                     end
                            %                     [ACI] = Convert_lasso_B2ACI(B, cfg_ACI{i_condition,i_subject,i_masker,i_glmfct}, 1, cfg_ACI{i_condition,i_subject,i_masker,i_glmfct}.keyvals);

                        case {'glmfitqp'}
                             % calculate ACI glmfitqp
                            [ACI,cfg_ACI{i_condition,i_subject,i_masker,i_glmfct},results{i_condition,i_subject,i_masker,i_glmfct}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
                            finallambda = results{i_condition,i_subject,i_masker,i_glmfct}.finallambda;
                            lambdas = results{i_condition,i_subject,i_masker,i_glmfct}.lambdas;
                            if finallambda>lambdas(end)
                                ACI = zeros(size(ACI));
                            end
                        otherwise
                            % calculate ACI classic revcorr 
                            [ACI,cfg_ACI{i_condition,i_subject,i_masker,i_glmfct},results{i_condition,i_subject,i_masker,i_glmfct}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

                    end

                    % % plot ACI lassoslow, classic revcorr, lasso
                    % idx_tile = i_condition+(i_subject-1)*N_conditions;
                    % figure(hACIind); nexttile(idx_tile); %(i_subject,i_masker)
                    % affichage_tf(ACI,'CI', 'cfg', cfg_ACI{i_condition}); hold on
                    % 
                    % if i_subject == 1
                    %     title(['ACI ' experiment(20:23)], 'interpreter','none');
                    % end
                    % %title(['ACI ' experiment(20:23) ', ' subject], 'interpreter','none');title('');
                    % xlabel(''); ylabel('');set(gca,'XTick',[]);set(gca,'YTick',[]);colorbar off
                    % text(0.95,0.95,subject,'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')
                    % set(hACIind,'Name',glmfct);

                    ACI_matrix(:,:,i_condition,i_subject,i_masker,i_glmfct) = ACI;

                end
            end
        end
    end
end

%figure(hACIind);

%% Correlation plot

if fig<3
    maskercolor = {'b','g','r'};
    glmline = {'-o','--*','-.s'};

clear rho
switch fig
    case 1
        figure('Position',[100 100 500 500]);
    case 2
        figure('Position',[100 100 500 300]);
end

for i_glmfct = 1:N_glmfcts
    for i_masker = 1:N_maskers
        for i_subject = 1:N_subjects
            rho(:,i_subject,i_masker,i_glmfct) = corr(reshape(ACI_matrix(:,:,:,i_subject,i_masker,i_glmfct),64*86,10),reshape(ACI_matrix(:,:,end,i_subject,i_masker,i_glmfct),64*86,1));
        end

        % consider only points where more than 2/3 of ACI have converged
        last_before_converge = max(find(sum(isnan(rho(:,:,i_masker,i_glmfct)),2)>33*N_subjects/100)); % locate the last i_condition before 50% have converged
        rho(1:last_before_converge,:,i_masker,i_glmfct) = nan;

        %errorbar((1:10)*400,mean(rho(:,:,i_masker,i_glmfct),2),std(rho(:,:,i_masker,i_glmfct),[],2),[maskercolor{i_masker} glmline{i_glmfct}]); hold on
        plot((1:10)*400,nanmean(rho(:,:,i_masker,i_glmfct),2),[maskercolor{i_masker} glmline{i_glmfct}], 'LineWidth',1.5); hold on
    end
end
xlim([0.5 10.5]*400)
ylim([0 1])
ylabel('Pearson''s correlation')
xlabel('Number of trial')
grid on
plegend = plot(0,0,maskercolor{1},0,0,maskercolor{2},0,0,maskercolor{3},0,0,['k' glmline{1}],0,0,['k' glmline{2}],0,0,['k' glmline{3}]);
set(gca,'FontSize',13)
%legend(plegend,[Maskers GLMfcts],'interpreter','none')

%just Classic-revcorr
figure('Position',[100 100 500 300]);
for i_glmfct = 1
    for i_masker = 1:N_maskers
        %errorbar((1:10)*400,mean(rho(:,:,i_masker,i_glmfct),2),std(rho(:,:,i_masker,i_glmfct),[],2),[maskercolor{i_masker} glmline{i_glmfct}]); hold on
        
        plot((1:10)*400,mean(rho(:,:,i_masker,i_glmfct),2),[maskercolor{i_masker} glmline{i_glmfct}], 'LineWidth',1.5); hold on

    end
    xlim([0.5 10.5]*400)
    ylim([0.2 1])
    ylabel('Pearson''s correlation')
    xlabel('Number of trial')
    grid on
    %legend(Maskers,'interpreter','none')
end
end

%% cue-to-noise-ratio plot

if fig<3
    % define ROIs
    t_onset = [0.25 0.3]; idx_onset = find(cfg_ACI{1, 1}.t>=t_onset(1) & cfg_ACI{1, 1}.t<=t_onset(2));
    t_dummy = [0.55 0.8]; idx_dummy = find(cfg_ACI{1, 1}.t>=t_dummy(1) & cfg_ACI{1, 1}.t<=t_dummy(2));
    f_F2 = [1000 2000]; idx_F2 = find(cfg_ACI{1, 1}.f>=f_F2(1) & cfg_ACI{1, 1}.f<=f_F2(2));
    clear CNR CNRdB
    switch fig
    case 1
        figure('Position',[100 100 500 500]);
    case 2
        figure('Position',[100 100 500 300]);
    end
    for i_glmfct = 1:N_glmfcts
        for i_masker = 1:N_maskers
            for i_subject = 1:N_subjects
                squaredW_F2onset = squeeze(nanmean(nanmean(abs(ACI_matrix(idx_F2,idx_onset,:,i_subject,i_masker,i_glmfct)).^2,1),2));
                squaredW_dummy = squeeze(nanmean(nanmean(abs(ACI_matrix(:,idx_dummy,:,i_subject,i_masker,i_glmfct)).^2,1),2));
                CNR(:,i_subject,i_masker,i_glmfct) = squaredW_F2onset./squaredW_dummy;

                %CNR(isnan(CNR(:,i_subject,i_masker,i_glmfct)),i_subject,i_masker,i_glmfct)=0.001;
                % replace infinite ratio by very high
                CNR(isinf(CNR(:,i_subject,i_masker,i_glmfct)),i_subject,i_masker,i_glmfct)=10000;

                CNRdB(:,i_subject,i_masker,i_glmfct) = 10*log10(CNR(:,i_subject,i_masker,i_glmfct));
            end

            % consider only points where more than 2/3 of ACI have converged
            last_before_converge = max(find(sum(isnan(CNR(:,:,i_masker,i_glmfct)),2)>33*N_subjects/100)); % locate the last i_condition before 50% have converged
            CNR(1:last_before_converge,:,i_masker,i_glmfct) = nan;

            %errorbar((1:10)*400,mean(CNRdB(:,:,i_masker,i_glmfct),2),std(CNR(:,:,i_masker,i_glmfct),[],2),[maskercolor{i_masker} glmline{i_glmfct}]); hold on
            plot((1:10)*400,nanmean(CNR(:,:,i_masker,i_glmfct),2),[maskercolor{i_masker} glmline{i_glmfct}], 'LineWidth',1.5); hold on
        end
    end


    set(gca, 'YScale', 'log')
    grid on
    plot([0 11]*400,[1 1],'k-')
    xlim([0.5 10.5]*400)
    switch fig
        case 1
            ylim([0.7 50000]);
        case 2
            ylim([0.7 5]);
    end
    ylabel('Cue-to-noise ratio')
    xlabel('Number of trial')
    plegend = plot(0,0,maskercolor{1},0,0,maskercolor{2},0,0,maskercolor{3},0,0,['k' glmline{1}],0,0,['k' glmline{2}],0,0,['k' glmline{3}],'LineWidth',1.5);
    legend(plegend,[Maskers_titles GLMfcts_titles],'interpreter','none','Location','northwest')
    set(gca,'FontSize',13)
    %%
    % just classic_revcorr
    figure('Position',[100 100 500 300]);
    for i_glmfct = 1:1
        for i_masker = 1:N_maskers
            %errorbar((1:10)*400,mean(CNRdB(:,:,i_masker,i_glmfct),2),std(CNR(:,:,i_masker,i_glmfct),[],2),[maskercolor{i_masker} glmline{i_glmfct}]); hold on
            plot((1:10)*400,nanmean(CNR(:,:,i_masker,i_glmfct),2),[maskercolor{i_masker} glmline{i_glmfct}], 'LineWidth',1.5); hold on
        end
        set(gca, 'YScale', 'log')
        grid on
        plot([0 11]*400,[1 1],'k-')
        xlim([0.5 10.5]*400)
        ylim([0.7 10])
        ylabel('Cue-to-noise ratio')
        xlabel('Number of trial')
        legend(Maskers,'interpreter','none')

    end
end

if fig == 3
    %% Plot average ACIs side by side

    hACImean = figure('Position',[100 100 750 250]); %hold on;
    t = tiledlayout(1,N_maskers,'TileSpacing', 'compact');
    %title(t, 'average ACIs', 'Interpreter', 'none')

    figure(hACImean);

    for i_masker=1:N_maskers
        nexttile(i_masker); %(i_subject,i_masker)
        affichage_tf(squeeze(mean(ACI_matrix(:,:,10,:,i_masker),4)),'CI', 'cfg', cfg_ACI{i_condition}); %hold on
        title([Maskers_titles(i_masker)], 'interpreter','none');
        if i_masker > 1
            ylabel('');set(gca,'YTickLabels',{});
        end
        if i_masker == 3
            c_lim = caxis;
            colorbar('Ticks',c_lim,'TickLabels',{'ada','aba'});% change colorbar
        end
        if i_masker < 3
            colorbar off
        end
    end
end
end