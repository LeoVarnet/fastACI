function [] = publ_osses2025_fig6()
% function publ_osses2025_fig6()
%%%%% Figure 6 from fastACI paper %%%%%
%
% % To display Fig. 6 of Osses et al use
%     publ_osses2025_fig6();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all, clear all

N_ACI = 3;

glmfct = 'l1glm';%'classic_revcorr';% 

trialtype_analysis = 'total';
TF_type = 'gammatone';
expvar_after_reversal = 4;

switch glmfct
    case {'l1lm', 'l1glm'}
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
    otherwise
        Lambdas = [];
end

%% ACIs (and load ACI data)

clear cfg_games cfg_ACI results Dev_test PC_test results ACI_matrix

hACIind = figure('Position', [100 100 600 500]); hold on;
tiledlayout(3,N_ACI,'TileSpacing', 'compact');

for i_ACI = 1:N_ACI

    switch i_ACI
        case 1
            subject = 'S01';
            masker = 'white';%'bumpv1p2_10dB';%'sMPSv1p3';%
            experiment = 'speechACI_Logatome-abda-S43M'; 
            Condition_name = 'ABDA22';
        case 2
            subject = 'SLeo';
            masker = 'bumpv1p2_5dB';%'bumpv1p2_10dB';%'white';%'sMPSv1p3';%
            experiment = 'speechACI_Logatome-abda-S41F'; 
            Condition_name = 'ABDA24';
        case 3
            subject = 'SLeo';
            masker = 'white';%'sMPSv1p3';%'bumpv1p2_10dB';%
            experiment = 'speechACI_varnet2013'; 
            Condition_name = 'ABDA13';
    end

    Data_matrix = [];
        dir_subject = [fastACI_dir_data experiment filesep subject filesep];
        dir_results = [dir_subject 'Results' filesep];

        fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);

        if ~isempty(fname_results)
            fname_results = fname_results{end};
            fprintf(['Processing ' fname_results '\n'])
            load([dir_results fname_results]);
            cd(dir_results)

            cfg_games{i_ACI} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders

            flags_for_input = {TF_type, ...
                'trialtype_analysis', trialtype_analysis,...
                'N_folds', 10, ...
                'dir_noise',cfg_games{i_ACI}.dir_noise, ...
                'dir_target',cfg_games{i_ACI}.dir_target, ...
                'add_signal',0, ...
                'apply_SNR',0, ...
                'skip_if_on_disk',1, ...
                'no_permutation', ...
                'no_bias', ...
                'no_plot', ...
                'expvar_after_reversal',expvar_after_reversal, ...
                'pyramid_shape',-1, ...
                'lambda',Lambdas ...
                };

            switch glmfct
                case {'lassoslow','lassoglmslow','l1glm'}
                    % calculate ACI lassoslow
                    [ACI,cfg_ACI{i_ACI},results{i_ACI}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

                    lambda = results{i_ACI}.FitInfo.Lambda;
                    Dev_test(i_ACI,:,:) = results{i_ACI}.FitInfo.Dev_test;
                    PC_test(i_ACI,:,:) = results{i_ACI}.FitInfo.PC_test;
                    idxlambda(i_ACI) = results{i_ACI}.idxlambda;
                    Nfolds = results{i_ACI}.FitInfo.CV.NumTestSets;
                    %statistical threshold (should be dependent on idxlambda)
                    B = results{i_ACI}.B(:,results{i_ACI}.idxlambda);
                   
                otherwise
                    % calculate ACI classic revcorr
                    [ACI,cfg_ACI{i_ACI},results{i_ACI}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

            end

            % plot ACI lassoslow, classic revcorr, lasso
            figure(hACIind); nexttile(3+(i_ACI-1)*N_ACI); %(i_subject,i_masker)
            affichage_tf(ACI,'CI', 'cfg', cfg_ACI{i_ACI}); hold on
            switch i_ACI
                case 2
                    c_axis = [-6.5 6.5]*1e-3;
                    clim(c_axis);
            end
            %title(['ACI ' experiment], 'interpreter','none');
            %xlabel(''); ylabel('');set(gca,'XTick',[]);set(gca,'YTick',[]);colorbar off
            %text(0.95,0.95,subject,'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')
            set(hACIind,'Name',glmfct);

            if i_ACI~=3
                xlabel('');
                %set(gca,'XTickLabels',[]);colorbar off
            end
            ylabel('');set(gca,'YTickLabels',[]);
            c_lim = clim;
            colorbar(gca, ...
               'XTickLabel',{'ada','aba'}, ...
               'XTick', c_lim)

            cfg_ACI{i_ACI} = Update_dir_main([fastACI_dir_data experiment filesep],cfg_ACI{i_ACI});

            files = Get_filenames(cfg_ACI{i_ACI}.dir_target,'*.wav');
            fname1 = [cfg_ACI{i_ACI}.dir_target files{1}];
            fname2 = [cfg_ACI{i_ACI}.dir_target files{2}];
            [T1,fs] = audioread(fname1);
            T2      = audioread(fname2);

            if fs==48000
                T1 = resample(T1,16000,fs);
                T2 = resample(T2,16000,fs);
                fs = 16000;
                % %downsample to fs = 16000
                % T1 = decimate(T1,3);
                % T2 = decimate(T2,3);
                % fs = 48000/3;
            end

            undersampling = 100;
            basef = 8000;
            flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5,'dboffset',100,'no_adt','binwidth',undersampling/fs, ...
                'no_outerear','no_middleear'};

            [G1]              = Gammatone_proc(T1, fs, flags_gamma{:});
            [G2, fc, t, outs] = Gammatone_proc(T2, fs, flags_gamma{:});

            G2 = log(G2'/(max(max(G2))));
            G1 = log(G1'/(max(max(G1))));

            bColour_bar = 'no'; % 'yes';
            switch i_ACI
                case 1
                    c_axis = [-5 0];%[0 0.7];
                case 2
                    c_axis = [-5 0];%[0 0.7];
                case 3
                    c_axis = [-5 0];%[0 0.7];
            end

            opts_colourbar = {'NfrequencyTicks', 8, 'colorbar', bColour_bar,'caxis',c_axis};
            figure(hACIind); nexttile(1+(i_ACI-1)*N_ACI);
            affichage_tf(G1, 'CI', t, fc, opts_colourbar{:}); hold on; % caxis([0 0.01]);
            affichage_tf_add_Praat_metrics_one_sound(fname1,cfg_ACI{i_ACI},[], '-', 'w');
            if i_ACI~=3
                xlabel('');
                %set(gca,'XTickLabels',[]);colorbar off
            end

            % show labels for formants and f0
            switch i_ACI
                case 1
                    text(0.6, 550, 'f_0','Color','w')
                    text(0.47, 2949, 'F1','Color','w')
                    text(0.47, 4086, 'F2','Color','w')
                    text(0.47, 5348, 'F3','Color','w')
                    text(0.47, 6320, 'F4','Color','w')
            end
    
            figure(hACIind); nexttile(2+(i_ACI-1)*N_ACI);
            affichage_tf(G2, 'CI', t, fc, opts_colourbar{:}); hold on; % caxis([0 0.01]);
            affichage_tf_add_Praat_metrics_one_sound(fname2,cfg_ACI{i_ACI},[], '-', 'w');
            ylabel('');set(gca,'YTick',[]);
            if i_ACI~=3
                xlabel('');
                %set(gca,'XTickLabels',[]);colorbar off
            end
        end

        % display formants on ACI 
        figure(hACIind); nexttile(3+(i_ACI-1)*N_ACI); %(i_subject,i_masker)
        affichage_tf_add_Praat_metrics_one_sound(fname1,cfg_ACI{i_ACI},[], '-', 'r');
        affichage_tf_add_Praat_metrics_one_sound(fname2,cfg_ACI{i_ACI},[], '-', 'b');
            

end

figure(hACIind);

end