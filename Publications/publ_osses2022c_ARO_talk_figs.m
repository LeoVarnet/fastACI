function publ_osses2022c_ARO_talk_figs(varargin)
% function [h,hname] = publ_osses2022c_ARO_talk_figs(varargin)
%
% See also the local file (fastACI_sim repo): g20210301_recreating_varnet2013.m
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% experiment = 'speechACI_varnet2013';
% model = 'osses2022a';
% g20211207_calibrating_the_model(experiment,model)

h = [];
hname = [];

definput.keyvals.dir_out = [];
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

dir_out = keyvals.dir_out;

dir_out_figs = fastACI_paths('dir_output'); % where figures and results will be stored

noise_types = {'white','SSN'};
Subjects = {'osses2022a','SLV','SAO'}; %{'SLV','SAO','osses2022a'};
Subjects_lab = {'osses2022a','S01','S02'}; % {'S01','S02','osses2022a'};

do_figACI = 1;
do_figKernels = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the ACIs
for i_subject = 1:length(Subjects)
    subj = Subjects{i_subject};
    for i_noise = 1:length(noise_types)
        noise_type = noise_types{i_noise};
        
        outs = publ_osses2022c_utils(subj,noise_type,'Get_flags');
        flags_for_input = outs.flags_for_input;
        
        switch subj
            case 'osses2022a'
                switch noise_type
                    case 'white'
                        dir_target = []; % outs.dir_target;
                        dir_noise  = []; % outs.dir_noise;
                        % fname_results = [fastACI_dir_data 'speechACI_varnet2013' filesep 'osses2022a' filesep 'Results' filesep 'savegame_2022_02_02_14_29_osses2022a_speechACI_varnet2013_white.mat'];
                        folder2look = {'Results-run-1-white'};
                    case 'SSN'
                        dir_target = []; % outs.dir_target;
                        dir_noise  = []; % outs.dir_noise;
                        % fname_results = [fastACI_dir_data 'speechACI_varnet2013' filesep 'osses2022a' filesep 'Results' filesep 'savegame_2022_02_02_19_35_osses2022a_speechACI_varnet2013_SSN.mat'];
                        folder2look = {'Results-run-1-SSN'};
                end
                dir_subj = [fastACI_dir_data 'speechACI_varnet2013' filesep subj filesep];
                dir_local = Get_filenames(dir_subj,[folder2look{1} '*']);
                if length(dir_local)>1
                    Show_cell(dir_local);
                    bInput = input('Please indicate the folder that should contain the results to process: ');
                    dir_local = dir_local{bInput};
                elseif ~isempty(dir_local)
                    dir_local = dir_local{1};
                else
                    % if not found at all
                    dir_local = [];
                end
                dir_local = [dir_subj dir_local filesep];
                fname_results = Get_filenames(dir_local,['savegame*' noise_type '.mat']); % [fastACI_dir_data 'speechACI_varnet2013' filesep subj filesep folders{1} filesep 'savegame_*_SSN.mat'];
                if ~isempty(fname_results)
                    fname_results = [dir_local fname_results{1}];
                end
                if ~exist(fname_results,'file')
                    % If no savegame file is found, the simulations will be run
                    publ_osses2022c_ARO_talk_1_sim(noise_type);
                
                    % Now we look again for a savegame file, and now should be found:
                    fname_results = Get_filenames(dir_local,['savegame*' noise_type '.mat']); % [fastACI_dir_data 'speechACI_varnet2013' filesep subj filesep folders{1} filesep 'savegame_*_SSN.mat'];
                    fname_results = [dir_local fname_results{1}];
                end
                
            otherwise
                outs = publ_osses2022c_utils(subj,noise_type,'Get_filenames');
                dir_target = outs.dir_target;
                dir_noise  = outs.dir_noise;
                fname_results = outs.fname_results;
        end
        
        %%% Extra flags if needed:
        if ~isempty(dir_noise)
            flags_for_input{end+1} = 'dir_noise';
            flags_for_input{end+1} = dir_noise;
        end
        if ~isempty(dir_target)
            flags_for_input{end+1} = 'dir_target';
            flags_for_input{end+1} = dir_target;
        end
        flags_for_input{end+1} = 'dir_out';
        flags_for_input{end+1} = dir_out;
        
        bCross_prediction = 0;
        if ~strcmp(subj,'osses2022a')
            % fcross = [fastACI_dir_data 'speechACI_varnet2013' filesep 'osses2022a' filesep 'Results' filesep 'Results_ACI' filesep 'ACI-osses2022a-speechACI_varnet2013-' noise_type '-nob-gt-l1glm-rev4.mat'];
            if ~isempty(dir_out)
                fcross = [dir_out 'ACI-osses2022a-speechACI_varnet2013-' noise_type '-nob-gt-l1glm-rev4.mat'];
            else
                fcross = [fastACI_dir_data 'speechACI_varnet2013' filesep 'osses2022a' ...
                    filesep 'Results' filesep 'Results_ACI' filesep 'ACI-osses2022a-speechACI_varnet2013-' noise_type '-nob-gt-l1glm-rev4.mat'];
            end
            
            if exist(fcross,'file')
                flags_for_input(end+1) = {'ACI_crosspred'};
                flags_for_input(end+1) = {fcross};
                bCross_prediction = 1;
            end
        end
        
        [ACI,cfg_ACI,results] = fastACI_getACI(fname_results,flags_for_input{:});
        if bCross_prediction
            crosspred = Check_if_crosspred(cfg_ACI.fnameACI,fcross);
        end
        if do_figACI
            h(end+1) = gcf;
            hname{end+1} = ['ACI-' Subjects_lab{i_subject} '-' noise_type];
            text4title = [Subjects_lab{i_subject} ': noise type=' noise_type];
            title(text4title);
            
            bAdd_formants = 0;
            if bAdd_formants
                affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
                outs_aff = affichage_tf_add_Praat_metrics(dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
            else
                if i_subject == 1 && i_noise == 1
                    warning('Temporal: Formants by-passed')
                end
            end
        end
        
        if do_figKernels
            YL_ker = [-10 90];
            XL_ker = [0 0.65]; % arbitrary
            XT_ker = [0.05:.05:0.6];
            YT_ker = -10:10:90;
            bApply_SNR = 1;
            if bApply_SNR == 1
                suff_SNR = '+SNR';
            else
                suff_SNR = '';
            end
            bMFB_with_LPF = 0;
            if bMFB_with_LPF
                suff_LPF = '';
            else
                suff_LPF = '-noLP';
            end
            model = 'osses2022a';
            fname_kern = sprintf('Kernels-%s-%s-model-%s%s%s.mat',Subjects_lab{i_subject},noise_type,model,suff_SNR,suff_LPF);
            fname_kern_full = sprintf('%s%s',dir_out_figs,fname_kern);
            
            bRun = ~exist(fname_kern_full,'file');
            
            dBFS = 100;
            
            if bRun
                kernel_H1 = [];  kernel_H2 = [];  kernel_H3 = [];
                kernel_M1 = [];  kernel_M2 = [];  kernel_M3 = [];
                kernel_CR1 = []; kernel_CR2 = []; kernel_CR3 = [];
                kernel_FA1 = []; kernel_FA2 = []; kernel_FA3 = [];
                
                kernel_H_avg1 = [];  kernel_H_avg2 = [];  kernel_H_avg3 = []; 
                kernel_M_avg1 = [];  kernel_M_avg2 = [];  kernel_M_avg3 = [];
                kernel_CR_avg1 = []; kernel_CR_avg2 = []; kernel_CR_avg3 = [];
                kernel_FA_avg1 = []; kernel_FA_avg2 = []; kernel_FA_avg3 = [];

                out_kernels = []; % struct to store the outputs
            
                idx = find(strcmp(flags_for_input,'idx_trialselect'),1,'first');
                if ~isempty(idx)
                    trials2load = max(flags_for_input{idx+1});
                else
                    trials2load = 4000; % max(trials_partial);
                end
                opts = [];
                opts.dir_noise = dir_noise;
                opts.dir_target = dir_target;
                [cfg_game, data_passation, ListStim] = Convert_ACI_data_type(fname_results,opts);
                cfg_ACI.dir_noise  = dir_noise; % make sure of the latest folder
                cfg_ACI.dir_target = dir_target; % make sure of the latest folder

                if isempty(dir_noise)
                    dir_noise = cfg_game.dir_noise;
                end
                if isempty(dir_target)
                    dir_target = cfg_game.dir_target;
                end
                files = Get_filenames(dir_noise,'*.wav');
                files = files(cfg_game.stim_order); % sorted as presented
                
                files = files(1:trials2load);
                n_responses = data_passation.n_responses(1:trials2load);
                is_correct  = data_passation.is_correct(1:trials2load);

                idxs_correct   = find(is_correct==1);
                idxs_incorrect = find(is_correct==0);
                idxs_target    = find(n_responses==1);
                idxs_nontarget = find(n_responses==2);

                idxs_target_H  = intersect(idxs_target,idxs_correct);
                idxs_target_M  = intersect(idxs_target,idxs_incorrect);
                idxs_target_CR = intersect(idxs_nontarget,idxs_correct);
                idxs_target_FA = intersect(idxs_nontarget,idxs_incorrect);

                idxs_labelled = zeros(size(n_responses));
                idxs_labelled(idxs_target_H) = 1;
                idxs_labelled(idxs_target_M) = 2;
                idxs_labelled(idxs_target_CR) = 3;
                idxs_labelled(idxs_target_FA) = 4;
            
                fi = Get_filenames(dir_target,'*.wav');
                in_aba = audioread([dir_target fi{1}]);
                in_ada = audioread([dir_target fi{2}]);
                
                for i_wav = 1:length(files)
                    fprintf('Processing sound %.1f\n',i_wav);
                    [insig,fs] = audioread([dir_noise files{i_wav}]);
                    if i_subject == 1 && strcmp(noise_type,'white')
                        if i_wav == 1
                            fprintf('TEMPORAL for subj=%s: Noises are scaled to the correct level...\n',Subjects{i_subject});
                        end
                        lvl_here = 65;
                        lvl_bef = rmsdb(insig);
                        insig = scaletodbspl(insig,lvl_here,dBFS);
                        lvl_aft = rmsdb(insig);
                        gain2use = 10^((lvl_aft-lvl_bef)/20); 
                    else
                        gain2use = 1;
                    end
                    if bApply_SNR
                        if data_passation.n_targets(i_wav) == 1
                            speech_here = gain2use*in_aba;
                        else
                            speech_here = gain2use*in_ada;
                        end
                        speech_here = 10^(data_passation.expvar(i_wav)/20)*speech_here;
                        
                        insig = insig+speech_here;
                    end

                    switch model
                        case 'osses2022a'
                            subfs = 16000;
                            flags_model = {'mfb','subfs',subfs,'mfb_osses2022a'};
                            [outsig,fc,mfc] = osses2022a(insig,fs,flags_model{:});

                            idx2use = [8 17 19]; % approx at 400, 1500, and 2000 Hz.
                            out_kernels.model = model;
                            out_kernels.flags_model = flags_model;
                            out_kernels.fc = fc(idx2use);
                            out_kernels.mfc = mfc;

                            outsig = il_sum_mfb(outsig,idx2use,bMFB_with_LPF);
                            idx2use = 1:length(idx2use);

                            out_kernels.idx2use = idx2use;
                    end

                    switch idxs_labelled(i_wav)
                        case 1
                            kernel_H1(:,end+1) = outsig(:,idx2use(1));
                            kernel_H2(:,end+1) = outsig(:,idx2use(2));
                            kernel_H3(:,end+1) = outsig(:,idx2use(3));
                        case 2
                            kernel_M1(:,end+1) = outsig(:,idx2use(1));
                            kernel_M2(:,end+1) = outsig(:,idx2use(2));
                            kernel_M3(:,end+1) = outsig(:,idx2use(3));
                        case 3
                            kernel_CR1(:,end+1) = outsig(:,idx2use(1));
                            kernel_CR2(:,end+1) = outsig(:,idx2use(2));
                            kernel_CR3(:,end+1) = outsig(:,idx2use(3));
                        case 4
                            kernel_FA1(:,end+1) = outsig(:,idx2use(1));
                            kernel_FA2(:,end+1) = outsig(:,idx2use(2));
                            kernel_FA3(:,end+1) = outsig(:,idx2use(3));
                        otherwise
                            warning('Skipping')
                    end
                end
                disp('')
            
                perc = [5 10 25 50 75 90 95];
                out_kernels.perc = perc;
                out_kernels.perc_description = 'Percentiles that were tested';

                for j = 1:length(perc)
                    kernel_H_avg1(:,j) = prctile(kernel_H1,perc(j),2);
                    kernel_H_avg2(:,j) = prctile(kernel_H2,perc(j),2);
                    kernel_H_avg3(:,j) = prctile(kernel_H3,perc(j),2);

                    kernel_M_avg1(:,j) = prctile(kernel_M1,perc(j),2);
                    kernel_M_avg2(:,j) = prctile(kernel_M2,perc(j),2);
                    kernel_M_avg3(:,j) = prctile(kernel_M3,perc(j),2);

                    kernel_CR_avg1(:,j) = prctile(kernel_CR1,perc(j),2);
                    kernel_CR_avg2(:,j) = prctile(kernel_CR2,perc(j),2);
                    kernel_CR_avg3(:,j) = prctile(kernel_CR3,perc(j),2);

                    kernel_FA_avg1(:,j) = prctile(kernel_FA1,perc(j),2);
                    kernel_FA_avg2(:,j) = prctile(kernel_FA2,perc(j),2);
                    kernel_FA_avg3(:,j) = prctile(kernel_FA3,perc(j),2);
                end

                t_ms = 1000*(1:size(kernel_H1,1))/subfs;
                out_kernels.t_ms = t_ms;

                out_kernels.kernel_H_avg1 = kernel_H_avg1;
                out_kernels.kernel_H_avg2 = kernel_H_avg2;
                out_kernels.kernel_H_avg3 = kernel_H_avg3;

                out_kernels.kernel_M_avg1 = kernel_M_avg1;
                out_kernels.kernel_M_avg2 = kernel_M_avg2;
                out_kernels.kernel_M_avg3 = kernel_M_avg3;

                out_kernels.kernel_CR_avg1 = kernel_CR_avg1;
                out_kernels.kernel_CR_avg2 = kernel_CR_avg2;
                out_kernels.kernel_CR_avg3 = kernel_CR_avg3;

                out_kernels.kernel_FA_avg1 = kernel_FA_avg1;
                out_kernels.kernel_FA_avg2 = kernel_FA_avg2;
                out_kernels.kernel_FA_avg3 = kernel_FA_avg3;

                out_kernels.idxs_labelled = idxs_labelled;
                out_kernels.idxs_labelled_description = '1=Hit, 2=Miss, 3=Correct rejection, 4=False alarm';
                
                save(fname_kern_full,'out_kernels');
            else
                out_kernels = []; % it will be loaded again
                load(fname_kern_full);
                perc = out_kernels.perc;
                
                kernel_H_avg1 = out_kernels.kernel_H_avg1;
                kernel_H_avg2 = out_kernels.kernel_H_avg2;
                kernel_H_avg3 = out_kernels.kernel_H_avg3;
                
                kernel_M_avg1 = out_kernels.kernel_M_avg1;
                kernel_M_avg2 = out_kernels.kernel_M_avg2;
                kernel_M_avg3 = out_kernels.kernel_M_avg3;
                
                kernel_CR_avg1 = out_kernels.kernel_CR_avg1;
                kernel_CR_avg2 = out_kernels.kernel_CR_avg2;
                kernel_CR_avg3 = out_kernels.kernel_CR_avg3;
                
                kernel_FA_avg1 = out_kernels.kernel_FA_avg1;
                kernel_FA_avg2 = out_kernels.kernel_FA_avg2;
                kernel_FA_avg3 = out_kernels.kernel_FA_avg3;
                
                t_ms = out_kernels.t_ms;
            end
            idx_me = find(perc==50,1,'first');
            t = t_ms/1000;
            
            %%% Only the kernels:
            offy = 20;
            
            figure;
            plot(t,kernel_H_avg1(:,idx_me),'r'); hold on;
            plot(t,kernel_M_avg1(:,idx_me),'k'); hold on;
            
            plot(t,offy+kernel_H_avg2(:,idx_me),'r'); 
            plot(t,offy+kernel_M_avg2(:,idx_me),'k'); 
            
            plot(t,2*offy+kernel_H_avg3(:,idx_me),'r'); 
            plot(t,2*offy+kernel_M_avg3(:,idx_me),'k'); 
                     
            title(sprintf('target present; subj=%s; noise=%s',Subjects_lab{i_subject},noise_type))
            
            h(end+1) = gcf;
            hname{end+1} = ['target-present-' Subjects_lab{i_subject} '-' noise_type];
            xlabel('Time (s)');
            
            ylim(YL_ker); grid on
            xlim(XL_ker);
            set(gca,'XTick',XT_ker);
            set(gca,'YTick',YT_ker);
            set(gca,'YTickLabel',[]);
            
            %%%
            if i_subject == 1
                if strcmp(noise_type,'SSN')
                    figure;
                    
                    idxs = 1:2:size(kernel_H_avg2,2);
                    plot(t,offy+kernel_H_avg2(:,idxs),'Color',.7*[1 1 1],'LineWidth',2);  hold on
                    plot(t,offy+kernel_H_avg2(:,idx_me),'r','LineWidth',2); 
                    % plot(t,offy+kernel_M_avg2(:,idx_me),'k','LineWidth',2); 

                    title(sprintf('target present + perc; subj=%s; noise=%s',Subjects_lab{i_subject},noise_type))

                    h(end+1) = gcf;
                    hname{end+1} = ['target-present+perc-HIT-' Subjects_lab{i_subject} '-' noise_type];
                    xlabel('Time (s)');

                    ylim(YL_ker); grid on
                    xlim(XL_ker);
                    set(gca,'XTick',XT_ker);
                    set(gca,'YTick',YT_ker);
                    set(gca,'YTickLabel',[]);
                    
                    %%%
                    figure;
                    idxs = 1:2:size(kernel_H_avg2,2);
                    plot(t,offy+kernel_M_avg2(:,idxs),'Color',.7*[1 1 1],'LineWidth',2);  hold on
                    plot(t,offy+kernel_M_avg2(:,idx_me),'k','LineWidth',2); 
                    
                    title(sprintf('target present + perc; subj=%s; noise=%s',Subjects_lab{i_subject},noise_type))

                    h(end+1) = gcf;
                    hname{end+1} = ['target-present+perc-MISS-' Subjects_lab{i_subject} '-' noise_type];
                    xlabel('Time (s)');

                    ylim(YL_ker); grid on
                    xlim(XL_ker);
                    set(gca,'XTick',XT_ker);
                    set(gca,'YTick',YT_ker);
                    set(gca,'YTickLabel',[]);
                end
            end
            %%%
            
            figure;
            plot(t,kernel_CR_avg1(:,idx_me),'b'); hold on;
            plot(t,kernel_FA_avg1(:,idx_me),'k'); hold on;
            
            plot(t,offy+kernel_CR_avg2(:,idx_me),'b'); 
            plot(t,offy+kernel_FA_avg2(:,idx_me),'k'); 
            
            plot(t,2*offy+kernel_CR_avg3(:,idx_me),'b'); 
            plot(t,2*offy+kernel_FA_avg3(:,idx_me),'k'); 
            xlabel('Time (s)');
            
            title(sprintf('target absent; subj=%s; noise=%s',Subjects_lab{i_subject},noise_type))
            h(end+1) = gcf;
            hname{end+1} = ['target-absent-' Subjects_lab{i_subject} '-' noise_type];
            
            ylim(YL_ker); grid on
            xlim(XL_ker);
            set(gca,'XTick',XT_ker);
            set(gca,'YTick',YT_ker);
            set(gca,'YTickLabel',[]);
           
            if i_subject == 1
                figure;
                fig_kernel(i_noise) = gcf;
                Style_here = 'r-';
                LW = 1;
            elseif i_subject == 2
                figure(fig_kernel(i_noise));
                Style_here = 'm--';
                LW = 2;
            else
                figure(fig_kernel(i_noise));
                Style_here = 'g:';
                LW = 3;
            end
                
            offy = 10;
            plot(t,kernel_H_avg1(:,idx_me)-kernel_M_avg1(:,idx_me),Style_here,'LineWidth',LW); hold on;
            plot(t,offy+(kernel_H_avg2(:,idx_me)-kernel_M_avg2(:,idx_me)),Style_here,'LineWidth',LW);
            plot(t,2*offy+(kernel_H_avg3(:,idx_me)-kernel_M_avg3(:,idx_me)),Style_here,'LineWidth',LW);
             
            title(sprintf('target-present kernel; noise=%s',noise_type))
            
            % ylim(YL_ker); 
            grid on
            xlim(XL_ker);
            set(gca,'XTick',XT_ker);
            set(gca,'YTick',[0:5:20]);
            set(gca,'YTickLabel',[]);

            if i_subject == 1
                h(end+1) = gcf;
                hname{end+1} = ['target-present-kernel-' noise_type];
                xlabel('Time (s)');
            end
            %%%%%
            
            if i_subject == 1
                figure;
                fig_kernel_abs(i_noise) = gcf;
                Style_here = 'b-';
                LW = 1;
            elseif i_subject == 2
                figure(fig_kernel_abs(i_noise));
                Style_here = 'm--';
                LW = 2;
            else
                figure(fig_kernel_abs(i_noise));
                Style_here = 'g:';
                LW = 3;
            end    
            plot(t,kernel_CR_avg1(:,idx_me)-kernel_FA_avg1(:,idx_me),Style_here,'LineWidth',LW); hold on;
            plot(t,offy+(kernel_CR_avg2(:,idx_me)-kernel_FA_avg2(:,idx_me)),Style_here,'LineWidth',LW);
            plot(t,2*offy+(kernel_CR_avg3(:,idx_me)-kernel_FA_avg3(:,idx_me)),Style_here,'LineWidth',LW);
             
            title(sprintf('target-absent kernel; noise=%s',noise_type))
            
            if i_subject == 1
                h(end+1) = gcf;
                hname{end+1} = ['target-absent-kernel-' noise_type];
                xlabel('Time (s)');
            end
            
            grid on
            xlim(XL_ker);
            set(gca,'XTick',XT_ker);
            set(gca,'YTick',[0:5:20]);
            set(gca,'YTickLabel',[]);

            disp('')
        end
        
        %% for i = 1:length(crosspred)
        Colour_here = 'b';

        xvar = results.lambdas;
        PC_test = results.FitInfo.PC_test;
        yvar = 100*mean(PC_test,2);
        
        figure;
        semilogx(xvar,yvar,'-','Color',Colour_here); hold on; grid on;
        idx = results.idxlambda;
        plot(xvar(idx),yvar(idx),'o','LineWidth',2,'MarkerFaceColor','b');
        
        
        if ~strcmp(subj,'osses2022a')
            Colour_here = 'g';
            xvar = crosspred.lambdas;
            PC_test = crosspred.PC_test;
            yvar = 100*mean(PC_test,2); % avoiding to write down 'results{1}.crosspred.PC_test' twice

            semilogx(xvar,yvar,'--','Color',Colour_here,'LineWidth',2); hold on;
            idx = crosspred.idxlambda;
            plot(xvar(idx),yvar(idx),'s','LineWidth',2,'MarkerFaceColor',Colour_here);
        end
        title(sprintf('%s: %s',subj,noise_type));
        % ylabel(['(cross)prediction accuracy benefit ' name_crosspred]);
        xlabel('\lambda');
        ylim([44 64]);
        ylabel('Prediction accuracy');
        
        plot([min(xvar) max(xvar)],50*[1 1],'k:','LineWidth',2);
        
        set(gca,'YTick',46:2:62);
        h(end+1) = gcf;
        hname{end+1} = ['Crosspred-' subj '-' noise_type];
        
        disp('')

        % legend(Masker_crosspred,'interpreter','none');
        % end 
    end
    
end

bSave = input('1=to save, 0 to cancel: ');
if bSave
    for i = 1:length(h)
        opts = [];
        opts.format = 'epsc';
        Saveas(h(i),[dir_out_figs hname{i}],opts);

        opts.format = 'png';
        Saveas(h(i),[dir_out_figs hname{i}],opts);

        opts.format = 'fig';
        Saveas(h(i),[dir_out_figs hname{i}],opts);
    end
end
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_sum_mfb(insig,idx2use,bWith_LPF)

insig = insig(idx2use);
outsig = [];

if bWith_LPF
    idxi_MFB = 1;
else
    idxi_MFB = 2;
end
for i = 1:length(insig)
    N_mod = size(insig{i},2);
    if bWith_LPF == 0
        N_mod = N_mod-1; % excluding the LPF
    end
    N_samples = size(insig{i},1);
    out_tmp = zeros([N_samples 1]);
    for j = idxi_MFB:N_mod
       out_tmp = out_tmp + insig{i}(:,j);
    end
    outsig(:,i) = out_tmp/N_mod;
end
