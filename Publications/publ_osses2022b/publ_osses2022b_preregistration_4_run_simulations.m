function publ_osses2022b_preregistration_4_run_simulations
% function publ_osses2022b_preregistration_4_run_simulations
%
% Leo's analysis has 13 figures:
%   Use Fig #
%       1
%       2  - SNR thresholds as a function of session for all three noises
%       3
%   -   4  - ACI classical revcorr
%   -   5  - ACI lasso
%   y   6  - ACI lasso slow
%
%   -   9  - ACI incorrect classical revcorr
%   -   10 - ACI incorrect lasso
%   y   11 - ACI incorrect lasso slow
%
% My comments: 
%       - 'lasso' with 'lassoslow' should be the same
%       
% Author: Leo Varnet, adapted by Alejandro
% Copied from l20211201_analysis_pipeline_versie_Alejandro.m on 13/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('This script is not versioned yet, and has not been tested lately')
clc, close all
h = []; hname = [];

% white:         4000 trials
% bumpv1p2_10dB: 4000 trials
% sMPSv1p2:      4000 trials
% sMPSv1p1:      3200 trials
Subject_ID = {'osses2022a'};
Subject_ref = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14'};
N_Subjs_ref = length(Subject_ref);

% Subjects = {'dau1997'}; N_subjects = length(Subjects);
% Subjects = {'relanoiborra2019'}; N_subjects = length(Subjects);
% Subjects = {'king2019'}; N_subjects = length(Subjects);
% Subjects = {'maxwell2020'}; N_subjects = length(Subjects);
Maskers = {'white','bumpv1p2_10dB', 'sMPSv1p3'}; N_maskers = length(Maskers);
experiment = 'speechACI_Logatome-abda-S43M';
maskercolors = {'b','r',rgb('Green')};
trialtype_analysis = 'total';  % 'incorrect'
TF_type = 'gammatone';
% idx_trialselect = 1:1000; warning('Temporal')
expvar_after_reversal = 4;

% N_trials = 4000; idx_trialselect = 1:N_trials;
glmfct = 'lassoglmslow';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bPlot = 0;
if bPlot
    hrev = figure; hold on;    
    h(end+1) = hrev;
    hname{end+1} = 'SNR-reversals';

    hdprime = figure; hold on; il_set_position(hdprime,3,1);
    h(end+1) = hdprime;
    hname{end+1} = 'Confusion-H+CR';

    hhist   = figure; hold on; il_set_position(hhist,1,3);
    h(end+1) = hhist;
    hname{end+1} = 'Confusion-Hist';

    hACI = figure; tiledlayout(1,N_maskers); % 1 x 3
    il_set_position(hACI,1,N_maskers)
    h(end+1) = hACI{i_trialtype_analysis};
    hname{end+1} = ['ACI-' glmfct '-' trialtype_analysis];

    hMSEcurve = figure; hold on; % 1 x 1
    h(end+1) = hMSEcurve;
    hname{end+1} = ['MSE-' trialtype_analysis];

    hPAcurve = figure; hold on; % 1 x 1
    h(end+1) = hPAcurve;
    hname{end+1} = ['PA-' trialtype_analysis];

    lMSE = [];
    lPA = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_subject = [fastACI_dir_data experiment filesep Subject_ID filesep];

switch Subject_ID
    case 'osses2021'
        YL_fig1_2 = [-30 -10];
    otherwise
        YL_fig1_2 = [-20 0];
end

for i_Subj_ref = 1:N_Subjs_ref
    for i_masker = 1:N_maskers
        masker      = Maskers{i_masker};
        colour_here = maskercolors{i_masker};

        % loading data
        switch Subject_ID
            case {'osses2021','osses2022a','dau1997','king2019','maxwell2020','relanoiborra2019'}
                % Case for simulations:
                dirs = Get_filenames(dir_subject,'Run*');
                Show_cell(dirs);
                % bInput = input('Enter the number of the run you want to analyse: ');
                disp('the last folder will be used...')
                disp('pausing for 5 seconds, cancel (ctrl+c) if you want to abort');
                pause(5);
                bInput = length(dirs);

                dir_results = [dir_subject dirs{bInput} filesep 'Results' filesep];

            otherwise
                dirs = Get_filenames(dir_subject,'Results*');
        end

        fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
        fname_results = fname_results{end};
        load([dir_results fname_results]);
        cd(dir_results)

        n_responses = data_passation.n_responses;
        SNR         = data_passation.expvar;
        resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
        N_sessions  = length(resume_trial)-1;

        if bPlot
            %%% Preparing for Fig. 1 and 2:
            bPlot_in_Script3 = 0;
            data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation,'histogram-leo',bPlot_in_Script3);

            % Fig. 2: Reversal analysis ---------------------------------------
            bUseMean = 0;
            for i_session = 1:N_sessions
                rev = Get_mAFC_reversals(SNR(resume_trial(i_session):resume_trial(i_session+1)));
                if bUseMean
                    meanSNR(i_session)=mean(rev);
                    stdSNR(i_session)=std(rev);
                else
                    % Then uses median:
                    medianSNR(i_session) = median(rev);
                    percSNR_L(i_session) = prctile(rev,5);
                    percSNR_U(i_session) = prctile(rev,95);
                end
            end
            if bUseMean
                Me = meanSNR;
                errL = stdSNR/2;
                errU = stdSNR/2;
            else
                Me = medianSNR;
                errU = percSNR_U-medianSNR;
                errL = medianSNR-percSNR_L;
            end

            % plot reversal analysis
            x_offset = 10*(i_masker-1); % for visualisation purposes
            x_var = resume_trial(2:end)+x_offset;

            bAdd_whole_session = 1;
            if bAdd_whole_session
                dx = x_var(end)-x_var(end-1)+1; x_var(end+1) = x_var(end)+dx;
                Me(end+1) = nan; errL(end+1) = nan; errU(end+1) = nan; % idle variables (empty line)

                % For the whole session:
                rev = Get_mAFC_reversals(SNR);
                Me_all = median(rev);
                errL_all = Me_all-prctile(rev,5);
                errU_all = prctile(rev,95)-Me_all;

                x_offset = 60*(i_masker-1); % for visualisation purposes
                x_var(end+1) = x_var(end)+dx+x_offset;
                Me(end+1) = Me_all; errL(end+1) = errL_all; errU(end+1) = errU_all; % idle variables
            end

            figure(hrev)
            errorbar(x_var,Me,errL,errU,'o-','Color',colour_here,'MarkerFaceColor',colour_here);

            if i_masker == 1
                set(gca,'XTick',[x_var(1:2:end-3) resume_trial(end)]); grid on

                if bUseMean
                    title('SNR thresholds (\pm 0.5 SD)'); 
                else
                    title('SNR thresholds (perc 5-95)'); 
                end
                xlabel('trial #'); 
                ylabel('SNR (dB)');
            end
            if i_masker == N_maskers    
                ylim(YL_fig1_2)
                legend(Maskers, 'interpreter', 'none')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fig. 2: Performance as a function of SNR
        bin_centres = data_hist.bin_centres; 
        P_H  = data_hist.H ./(data_hist.H +data_hist.M );
        P_FA = data_hist.FA./(data_hist.CR+data_hist.FA);
        dprime = norminv(P_H)-norminv(P_FA);           % Eq.  9 from Harvey2004
        criterion = -(norminv(P_H) + norminv(P_FA))/2; % Eq. 12 from Harvey2004

        if bPlot
            % plot performance as a function of SNR
            P_H_resp2 = P_H; % Script3 assumes that the target sound is option '2'
            resp2_label = [data_hist.H_label '-trials (H)'];

            P_H_resp1 = 1-P_FA; % 'Hits' when the participant response was '1'
            resp1_label = [data_hist.CR_label '-trials (CR)'];

            figure(hdprime);
            subplot(3,1,1); 
            plot(bin_centres, P_H_resp2 ,'-o','Color',colour_here); hold on, grid on
            plot(bin_centres, P_H_resp1 ,'--o','Color',colour_here);
            if i_masker == N_maskers
                title('PC'); 
                xlabel('SNR (dB)'); 
                ylabel('rate correct answers');
                ylim([0.45 1])
                legend({resp2_label,resp1_label},'Location','SouthEast')
            end

            subplot(3,1,2)
            plot(bin_centres, dprime,'-o','Color',colour_here); hold on, grid on;
            if i_masker == 1
                title('d prime'); 
                xlabel('SNR (dB)'); 
                ylabel('d''');
            end
            if i_masker == N_maskers    
                legend(Maskers, 'interpreter', 'none','Location','SouthEast')
            end

            subplot(3,1,3)
            plot(bin_centres, criterion ,'-o','Color',colour_here); hold on, grid on;
            if i_masker == 1
                title('criterion'); 
                xlabel('SNR (dB)'); 
                ylabel('c');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N_hist = data_hist.N_hist;
            perc_correct = 100*data_hist.prop_correct;
            perc_bias    = 100*data_hist.bias_if_1;

            figure(hhist);
            subplot(1,3,i_masker)
            yyaxis left
            area(bin_centres,N_hist,'FaceAlpha',0.2,'EdgeColor','none');

            str4label = '(SNR in dB)';
            xlabel(['expvar ' str4label]);

            ylabel('Nr. of trials')
            yyaxis right

            ha = gca;
            % ha.YAxis(1).Color = 'k';
            ha.YAxis(2).Color = colour_here;

            hprop = plot(bin_centres,perc_correct,'-o','Color',colour_here,'MarkerFaceColor',colour_here);
            ylabel('Percentage correct (%)'); 
            ylim([0 100]); hold on

            hbias = plot(bin_centres,perc_bias, '--o','Color',colour_here,'MarkerFaceColor',colour_here);
            ylim([0 100])
            grid on

            yyaxis left
            ylim([0 950]);
            yyaxis right
            ylim([-3 103]);

            if i_masker == N_maskers
                legend([hprop, hbias],{'correct','bias (resp. was 1)'},'Location','NorthWest');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare ACI analysis
        cfg_game = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders

        Data_matrix = [];
        flags_for_input = {TF_type, ...
            'trialtype_analysis', trialtype_analysis,...
            'N_folds', 10, ...
            'dir_noise',cfg_game.dir_noise, ...
            'dir_target',cfg_game.dir_target, ...
            'add_signal',0, ...
            'apply_SNR',0, ...
            'skip_if_on_disk',1, ...
            'no_permutation', ...
            'no_bias', ...
            'no_plot','expvar_after_reversal',expvar_after_reversal, ... % 'idx_trialselect',idx_trialselect ...
            };

        switch glmfct
            case {'lassoslow','lassoglmslow'}
                % calculate ACI lassoslow
                [ACI,cfg_ACI,results, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

                % % Same, but no Data_matrix is requested
                % [ACI,cfg_ACI,results] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

                handle_here = hACI;

            otherwise
                % calculate ACI classic revcorr
                [ACI,cfg_ACI,results] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
                if strcmp(glmfct,'classic_revcorr')
                    handle_here = hACI_classic;
                end
                if strcmp(glmfct,'lasso')
                    handle_here = hACI_lasso;
                end
        end

        if bPlot
            % plot ACI lassoslow, classic revcorr, lasso
            figure(handle_here); nexttile; %(i_subject,i_masker)
            affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
            outs_aff = affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
            if i_masker ~= N_maskers
                set(outs_aff.hl,'Visible','off'); % Removes the legend
            end
            title(['ACI ' cfg_ACI.Condition ', ' cfg_ACI.keyvals.trialtype_analysis], 'interpreter','none')
            set(handle_here,'Name',glmfct);

            %MSE curve
            x_var = results.FitInfo.Lambda;
            y_var = results.FitInfo.MSEtest;
            idx_here = results.idxlambda;

            Me = mean(y_var,2);
            Std = std(results.FitInfo.MSEtest,[],2);

            figure(hMSEcurve)
            h_here = errorbar(x_var, Me, Std, 'Color',colour_here,'LineWidth',1); hold on, grid on

            plot(x_var(idx_here), Me(idx_here),'*','Color',colour_here,'LineWidth',2);
            lMSE(end+1) = h_here;

            if i_masker == N_maskers
                set(gca,'XScale','log'); 
                ylim([0.2 0.35])
                title('MSE'); 
                xlabel('\lambda'); 
                ylabel('MSE');
                legend(lMSE, Maskers, 'interpreter', 'none')
            end

            % Percent accuracy curve
            % Same 'x_var', but new 'y_var':
            y_var_bias = results.FitInfo.PCtest(end,:);
            y_var      = results.FitInfo.PCtest-y_var_bias;
            Me     = 100*mean(y_var,2);
            y_bars = 100*sem(y_var,[],2);

            figure(hPAcurve)
            h_here = errorbar(x_var, Me, y_bars, 'Color', colour_here,'LineWidth',1); hold on, grid on
            plot(x_var(idx_here), Me(idx_here),'*','Color',colour_here,'LineWidth',2);

            lPA(end+1) = h_here;
            if i_masker == N_maskers
                set(gca,'XScale','log'); 
                ylim([-10 30])
                title('Percentage of responses explained by ACI'); 
                xlabel('\lambda'); 
                ylabel('% accuracy difference');
                legend(lPA, Maskers, 'interpreter', 'none')
            end
        end
    end
end

if bPlot
    bSave = 1;

    if bSave == 1
        dir_out = fastACI_paths('dir_output');

        for i = 1:length(h)
            opts = [];
            opts.format = 'epsc';
            fname = [dir_out 'Subj-' Subject_ID '-' hname{i}];
            Saveas(h(i),fname,opts);

            opts.format = 'fig';
            Saveas(h(i),fname,opts);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_set_position(handle_fig,N_rows,N_cols)
% % Regular position:
% Pos34 = [ 570 420]; % for 1 x 1
% Pos34 = [ 800 420]; % for 1 x 2
% Pos34 = [1100 420]; % for 1 x 3
% Pos34 = [ 570 700]; % for 3 x 1

Pos = get(handle_fig,'Position');
switch N_rows
    case 1
        Pos(4) = 420;
    case 2
        % No figure with this config
    case 3
        Pos(4) = 700;
end
switch N_cols
    case 1
        Pos(3) = 570;
    case 2
        Pos(3) = 800;
    case 3
        Pos(3) = 1100;
end
set(handle_fig,'Position',Pos);
