function [h, hname] = publ_osses2022c_ARO_poster_figs(varargin)
% function [h,hname] = publ_osses2022c_ARO_poster_figs(varargin)
%
% It generates the figures from Osses and Varnet (2022, ARO poster).
% You need to specify the figure number.
%
% % Characterisation of the background noises: 
% %   You can plot all panels of Fig. 2 by using::
%   publ_osses2022c_ARO_poster_figs('fig2');
%
% % Auditory Classification Images (ACIs): 
% %   You can plot all panels of Fig. 3 by using::
%   publ_osses2022c_ARO_poster_figs('fig3');
%
% %   or you can plot the ACIs of participant S1 or S2 only by using:
%   publ_osses2022c_ARO_poster_figs('fig3a');
%   publ_osses2022c_ARO_poster_figs('fig3b');
%
% % Correlation between partial ACIs and full ACIs:
% %   You can plot the correlations for both participants by using::
%   publ_osses2022c_ARO_poster_figs('fig4');
%
% %   or you can plot the correlations for participant S1 or S2 only by using::
%   publ_osses2022c_ARO_poster_figs('fig4a');
%   publ_osses2022c_ARO_poster_figs('fig4b');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.type={'fig2a','fig2b','fig2c','fig2d','fig3','fig3a','fig3b', ...
    'fig4','fig4a','fig4b'};
definput.keyvals.dir_out = [];
[flags,keyvals]  = ltfatarghelper({},definput,varargin);
dir_out = keyvals.dir_out; % relevant for 'fig3' and 'fig4'

if nargin == 0
    close all
    clc
end

h = [];
ha = [];
hname = [];

experiment = 'speechACI_Logatome-abda-S43M';
noise_types = {'white','sMPSv1p3','bumpv1p2_10dB'};
noise_types_label = {'WN','MPS','BP'};

Subjects = {'SLV','SAO'};
Subjects_id = [1 2];
    
do_fig2 = flags.do_fig2a || flags.do_fig2b || flags.do_fig2c || flags.do_fig2d;

if do_fig2 || flags.do_fig3a || flags.do_fig4a
    idx = 1;
    Subjects = Subjects(idx);
    Subjects_id = Subjects_id(idx);
end
if flags.do_fig3b || flags.do_fig4b
    idx = 2;
    Subjects = Subjects(idx);
    Subjects_id = Subjects_id(idx);
end
%%% First, we check whether the participants' data are stored on disk:
for i = 1:length(Subjects)
    for j = length(noise_types):-1:1 % reversed order (in case we need to remove one noise condition)
        outs = publ_osses2022c_ARO_poster_utils(Subjects{i},noise_types{j},'Get_filenames');
        
        if outs.bGenerate_stimuli == 1
            
            warning('Re-generate waveforms first...');
            % if exist(outs.dir_noise,'dir')
            Subject_id_here = ['S' num2str(Subjects_id(i))];
            dir_crea = [fastACI_basepath 'Publications' filesep 'publ_osses2022c' filesep 'data_' Subject_id_here filesep '0-init' filesep];
            file_crea = Get_filenames(dir_crea,['cfgcrea*' noise_types{j} '.mat']);
            
            if length(file_crea) == 1
                [~,bReproducible] = fastACI_experiment_init_from_cfg_crea([dir_crea file_crea{1}]);    
                
                if bReproducible == 0
                    noise = noise_types{j};
                    dir_subj = [fastACI_dir_data  experiment filesep Subjects{i} filesep];
                    dir_crea_dst = [dir_subj 'Results' filesep];
                    if exist(dir_crea_dst,'dir')
                        bContinue = input(['Do you have the sounds for participant ' Subjects{i} ', condition ' noise '? (1=yes, 0=no): ']);
                        if bContinue
                            disp(['Please place the sounds under ' dir_subj 'NoiseStim-' noise filesep ' and press any button to continue (press ctrl+c to abort)']);
                            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                            pause;
                            
                            copyfile([dir_crea file_crea{1}],dir_crea_dst);
                            dir_save = [fileparts(dir_crea(1:end-1)) filesep '1-experimental_results' filesep];
                            file_save = Get_filenames(dir_save,['savegame*' noise_types{j} '.mat']);

                            % This file is from the repository: it is always there
                            copyfile([dir_save file_save{1}],dir_crea_dst);
                        end
                    end
                    
                    warning('The sounds for Subject %s, cond=%s, don''t seem reproducible. Skippping this condition.',Subjects{i},noise_types{j});
                    noise_types(j) = [];
                    noise_types_label(j) = [];
                end
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig2
    % Migrated from g20220201_characterising_noises_AROposter.m
    
    Colours = {rgb('Gray'),'k',rgb('Maroon')};
    
    if flags.do_fig2a 
        fname_sound = 'Noise_00001.wav';
        % cfg_affi = [];
        for i = 1:length(noise_types)
            outs = publ_osses2022c_ARO_poster_utils('SLV',noise_types{i},'Get_filenames');
            sound2read = [outs.dir_noise fname_sound];
            [insig,fs] = audioread(sound2read);

            L = 512;
            L_overlap = round(0.9*L);
            L_f = L;

            win = hamming(L);
            [~,f_spec,t_spec,P] = spectrogram(insig,win,L_overlap,L_f,fs);
            P_dB = 10*log10(abs(P));
            max_dB = max(max(P_dB));
            P_dB = P_dB-max_dB;
            P_dB(1)   =-65; % 'min value': Trick to fix the ylimits
            P_dB(end) =  0; % 'max value': Trick to fix the ylimits

            figure;
            plot_stft(t_spec,f_spec,P_dB);
            title4plot = sprintf('%s, Noise0001.wav: Short-time Fourier transform',noise_types_label{i});
            title(title4plot);

            xlim([t_spec(1) .5])
            set(gca,'XTick',.1:.1:.4);
            % ACI: 64x34 (64 freqs)
            il_Post_figure;

            h(end+1) = gcf;
            hname{end+1} = ['STFT-' noise_types_label{i}];
        end
    end
    
    if flags.do_fig2b || flags.do_fig2c || flags.do_fig2d
        N_sounds = 5000;
        opts = [];    
        
        % Band levels are always computed
        do_kohlrausch2021_env = flags.do_fig2c; % Modulation spectrum
        do_V                  = flags.do_fig2d; % V metric

        fmod_xlim = [0 200];
        fmod_ylim = [10 50];
        fc_ylim = [30 70];
        V_lim = [-10 0];
        %%%

        for i = 1:length(noise_types)
            outs = publ_osses2022c_ARO_poster_utils('SLV',noise_types{i},'Get_filenames');
            dir_where = outs.dir_noise;
            suff = noise_types_label{i};
            opts.Colour = Colours{i};

            dBFS = 100;

            [~,dir2check1] = fileparts(dir_where(1:end-1));
            % dir_where = [dir_where filesep];

            files1 = Get_filenames(outs.dir_noise,'*.wav');
            files1 = files1(1:N_sounds);

            lvls = [];
            for j = 1:N_sounds
                file = files1{j};
                [insig,fs] = audioread([outs.dir_noise file]);

                if do_kohlrausch2021_env
                    [env_dB(:,j),xx,env_extra] = Get_envelope_metric(insig,fs,'kohlrausch2021_env_noDC');
                end

                % if do_varnet2017_env
                %     % Removed, have a look at g20220201...
                % end

                [outsig1,fc] = auditoryfilterbank(insig,fs);
                t = (1:size(outsig1,1))/fs;

                lvls(j,:) = rmsdb(outsig1) + dBFS;

                for i_fc = 1:length(fc)
                    % if do_W
                    %     % Removed, have a look at g20220201...
                    % end
                    if do_V
                        [V1(j,i_fc),description,yenv1] = Get_envelope_metric(outsig1(:,i_fc),fs,'V');
                    end

                end

                if mod(j,50) == 1
                    fprintf('\tProcessing sound %.0f of %.0f\n',j,N_sounds);
                end
            end

            % if do_varnet2017_env
            %     % Removed, have a look at g20220201...
            % end

            if do_kohlrausch2021_env
                extra.f_env  = env_extra.f_env;
                extra.fs_env = env_extra.fs_env;
                extra.env_dB_U = prctile(env_dB,95,2);
                extra.env_dB   = prctile(env_dB,50,2);
                extra.env_dB_L = prctile(env_dB, 5,2);

                figure;
                plot(extra.f_env,extra.env_dB,'-','Color',Colours{i},'LineWidth',2); hold on, grid on
                plot(extra.f_env,extra.env_dB_U,'-','Color',rgb('Gray'));
                plot(extra.f_env,extra.env_dB_L,'-','Color',rgb('SlateGray'));
                ylabel('Envelope spectrum (dB)')
                xlabel('Modulation frequency (Hz)')

                il_Post_figure;

                deltaf = 25;
                XT = deltaf:deltaf:extra.fs_env/2-deltaf;
                for i_xt = 1:length(XT)
                    if mod(i_xt,2) == 1
                        XTL{i_xt} = '';
                    else
                        XTL{i_xt} = num2str(XT(i_xt));
                    end
                end
                set(gca,'XTick',XT);
                set(gca,'XTickLabel',XTL);

                % Modulation spectrum
                xlim(fmod_xlim);
                ylim(fmod_ylim);

                title(sprintf('Envelope for %.0f sounds',N_sounds))

                h(end+1) = gcf;
                hname{end+1} = sprintf('Env-%s-N-%.0f%s',dir2check1,N_sounds,suff);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if flags.do_fig2b
                L1me = prctile(lvls,50);
                errLL1 = L1me - prctile(lvls,5);
                errUL1 = prctile(lvls,95)-L1me;

                figure;
                semilogx(fc,L1me,'o-','Color',Colours{i},'MarkerFaceColor',Colours{i}); hold on; grid on
                errorbar(fc,L1me,errLL1,errUL1,'-','Color',Colours{i});
                ylabel('Band level (dB)')
                xlabel('Frequency (Hz)')

                XT = [125 250 500 1000 2000 4000 8000];
                set(gca,'XTick',XT);

                il_Post_figure;

                ylim(fc_ylim);

                title(sprintf('rmsdb for %.0f sounds',N_sounds))
                h(end+1) = gcf;
                hname{end+1} = sprintf('BL-%s-N-%.0f%s',dir2check1,N_sounds,suff);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if do_V
                %%% Metric V
                V1me = prctile(V1,50);
                errL1 = V1me - prctile(V1,5);
                errU1 = prctile(V1,95)-V1me;

                figure;
                semilogx(fc,V1me,'o-','Color',Colours{i},'MarkerFaceColor',Colours{i}); hold on; grid on
                errorbar(fc,V1me,errL1,errU1,'-','Color',Colours{i});

                ylabel('V (dB)')
                xlabel('Frequency (Hz)')

                XT = [125 250 500 1000 2000 4000 8000];
                set(gca,'XTick',XT);

                ylim(V_lim);

                title(sprintf('Metric V for %.0f sounds',N_sounds))
                h(end+1) = gcf;
                hname{end+1} = sprintf('V-metric-%s-N-%.0f%s',dir2check1,N_sounds,suff);
            end

            if flags.do_fig2b
                lvls_broadband = sum_dB_power(lvls);
                [prctile(lvls_broadband,5) prctile(lvls_broadband,50) prctile(lvls_broadband,95)]
            end
        end    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any of the remaining figures:
if flags.do_fig3 || flags.do_fig3a || flags.do_fig3b || flags.do_fig4 || flags.do_fig4a || flags.do_fig4b
    
    bPlot_ACIs      = ~(flags.do_fig3 + flags.do_fig3a + flags.do_fig3b==0);
    bDo_partialACIs = flags.do_fig4 || flags.do_fig4a || flags.do_fig4b;    

    %%% Specific ACI configuration:
    TF_type = 'gammatone';
    glmfct  = 'l1glm'; 

    N_lambda = 30;
    Lambdas = logspace(-4, -1, N_lambda);
    idx = find(Lambdas >= 10^-3);
    Lambdas = Lambdas(idx);
    %%%

    flags_for_input = {TF_type, ...
                       glmfct, ... % 'dir_noise',dir_noise, 'dir_target',dir_target, ...
                       'trialtype_analysis', 'total', ...
                       'add_signal',0, ...
                       'apply_SNR',0, ...
                       'skip_if_on_disk',1, ...
                       'expvar_after_reversal', 4, ...
                       'lambda', Lambdas, ...
                       'no_permutation', 'no_bias', ...
                       'pyramid_script','imresize', ...
                       'pyramid_shape',0}; 
    if bPlot_ACIs == 0
        flags_for_input{end+1} = 'no_plot';
    end
    flags_for_input{end+1} = 'dir_out';
    flags_for_input{end+1} = dir_out;

    do_fig_formatting = 1;

    if do_fig_formatting
        Pos3 = 500;
        Pos4 = 600;
        XL = [0 0.5];
        XT = 0.1:.1:.4;
        FS = 18;
    end

    if bDo_partialACIs
        LW = 2;
        Markers = {'o-','s--','d:'};
        Colours = {rgb('Gray'),rgb('Maroon'),'k'};
    end

    for i_subject = 1:length(Subjects)
        subj = Subjects{i_subject};
        subj_here = ['S0' num2str(Subjects_id(i_subject))];

        for i_noise = 1:length(noise_types)
            noise_type = noise_types{i_noise};

            Data_matrix = []; % init
            dir_noise = []; % init
            dir_target = []; % init
            
            dir_subj = [fastACI_dir_data experiment filesep subj filesep];
            dir_res = [dir_subj 'Results' filesep];
            fname_results = Get_filenames(dir_res,['savegame*' noise_type '.mat']);
            if length(fname_results) ~= 1
                error('More than one file was found...')
            end
            fname_results = [dir_res fname_results{1}];

            [ACI,cfg_ACI,results, Data_matrix] = fastACI_getACI(fname_results,flags_for_input{:});

            if bPlot_ACIs
                h(end+1) = gcf;
                hname{end+1} = [subj_here '-ACI-' noise_type];
                ha(end+1) = gca;
                il_Post_figure;

                if do_fig_formatting
                    title('');

                    text4title = sprintf('%s: %s',subj_here,noise_types_label{i_noise});
                    text(.05,.95,text4title,'Units','Normalized','FontSize',FS,'FontWeight','Bold')

                    Pos = get(h(end),'Position');
                    Pos(3) = Pos3; 
                    Pos(4) = Pos4; 
                    set(h(end),'Position',Pos);

                    xlim(XL);
                    set(gca,'XTick',XT);
                    set(gca,'FontSize',FS);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if bDo_partialACIs
                %%% Loading Data_matrix upfront: It will only be reloaded when
                trials_partial = 400:400:3600;
                lambda_opt = results.lambdas(results.idxlambda);

                flags_here = flags_for_input;
                idx = find(strcmp(flags_here,'lambda'));
                if ~isempty(idx)
                    flags_here{idx+1} = lambda_opt; % Replacing the lambda for the optimal lambda
                end

                ACI_all = [];
                for i = 1:length(trials_partial)

                    idx_trialselect = 1:trials_partial(i);

                    idx = find(strcmp(flags_here,'idx_trialselect'));
                    if ~isempty(idx)
                        % Then the field exists and needs to be replaced:
                        flags_here{idx+1} = idx_trialselect; % Replacing the lambda for the optimal lambda
                    else
                        % Then the keyval does not exist and can be appended to the list of parametrs:
                        flags_here{end+1} = 'idx_trialselect';
                        flags_here{end+1} = idx_trialselect;
                    end
                    fnameACI = fastACI_getACI_fname(fname_results, flags_here{:});

                    Data_partial = [];
                    if ~exist(fnameACI,'file')
                        % Then the ACI still needs to be generated
                        if isempty(Data_matrix)
                            % Loading the data matrix:
                            trials2load = max(trials_partial);
                            [cfg_game, data_passation, ListStim] = Convert_ACI_data_type(fname_results);
                            cfg_ACI_here = cfg_ACI;
                            cfg_ACI_here.N = trials2load;
                            cfg_ACI_here.N_trialselect = trials2load;
                            Data_matrix = fastACI_getACI_dataload(cfg_ACI_here, ListStim,cfg_game);
                        else
                            % Nothing to do, because then data matrix exists already
                        end
                        Data_partial = Data_matrix(idx_trialselect,:,:);
                        % calculate ACI lassoslow
                        % ACI_partial(:,:,i) = fastACI_getACI(fname_results, flags_here{:}, ...
                        %     'Data_matrix', Data_partial); 
                    end
                    ACI_partial = fastACI_getACI(fname_results, flags_here{:}, 'no_plot',...
                        'Data_matrix', Data_partial); 

                    corr_curve_pearson(i_noise,i,i_subject) = corr(ACI_partial(:),ACI(:),'type','Pearson');
                    corr_curve_spearman(i_noise,i,i_subject) = corr(ACI_partial(:),ACI(:),'type','Spearman');
                end
            end        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        if bDo_partialACIs
            figure;
            % LW = 2;
            % Markers = {'o-','s--','d:'};
            % Colours = {rgb('Gray'),rgb('Maroon'),'k'};
            for i_noise = 1:length(noise_types)
                plot(trials_partial,corr_curve_pearson(i_noise,:,i_subject), ...
                    Markers{i_noise}, 'Color', Colours{i_noise}, ...
                    'MarkerFaceColor',Colours{i_noise},'LineWidth',LW); hold on, grid on
            end
            text4title = sprintf('%s: Pearson correlation',subj_here);
            text(.05,.95,text4title,'Units','Normalized','FontSize',FS,'FontWeight','Bold');

            h(end+1) = gcf;
            hname{end+1} = [subj_here '-correlation'];
            ha(end+1) = gca;
            il_Post_figure;

            legend(noise_types_label,'Location','SouthEast')

            Pos = get(h(end),'Position');
            Pos(3) = 1.4*Pos3; 
            Pos(4) = Pos4; 
            set(h(end),'Position',Pos);

            xlim([0 4000]);
            set(gca,'XTick',trials_partial);
            XT_here = [];
            for i_tick = 1:length(trials_partial)
                if mod(i_tick,2)==1
                    XT_here{i_tick} = num2str(trials_partial(i_tick));
                else
                    XT_here{i_tick} = '';
                end
            end
            set(gca,'XTickLabel',XT_here);
            set(gca,'YTick',0:.1:1)
            xlabel('Trials used to assess the partial ACIs')

            ylabel('Correlation value');
            set(gca,'FontSize',FS);

            ylim([0 1])

            disp('')
        end
    end
end

if nargout == 0
    dir_out_figs = fastACI_paths('dir_output');

    for i = 1:length(h)
        opts = [];
        opts.format = 'epsc';
        Saveas(h(i),[dir_out_figs hname{i}],opts);

        opts.format = 'fig';
        Saveas(h(i),[dir_out_figs hname{i}],opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_Post_figure

pause_len = 1;
fprintf('Pausing for %.1f seconds\n',pause_len);
pause(pause_len)