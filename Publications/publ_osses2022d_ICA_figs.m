function publ_osses2022d_ICA_figs(varargin)
% function l20220303_figs_for_CFA_version_Alejandro(varargin)
%
% % To display Fig. 2 of Osses, Lorenzi, and Varnet (2022, ICA) use :::
%     publ_osses2022d_ICA_figs('fig2');
%
% % To display Table. 2 of Osses, Lorenzi, and Varnet (2022, ICA) use :::
%     publ_osses2022d_ICA_figs('table2');
%
% Original name: l20220303_figs_for_CFA_version_Alejandro.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc

if nargin == 0
    help publ_osses2022d_ICA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2','table2'};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if ~isunix
    % Leo's local names
    experiments = { {'modulationACI', 'white', 'S4'}, ...
                    {'modulationACI', 'white', 'S1'}, ...
                    {'speechACI_varnet2013', 'white', 'SLeo'}, ...
                    {'speechACI_varnet2013', 'white', 'SAO'}, ...
                    {'speechACI_varnet2013', 'SSN', 'S01'}, ...
                    {'speechACI_varnet2013', 'SSN', 'S02'}, ...
                    {'speechACI_Logatome-abda-S43M', 'white', 'SLV'}, ...
                    {'speechACI_Logatome-abda-S43M', 'white', 'SAO'}};
else
    % Alejandro's local names
    experiments = { {'modulationACI', 'white', 'S4'}, ...
                    {'modulationACI', 'white', 'S1'}, ...
                    {'speechACI_varnet2013', 'white', 'SLeo'}, ...
                    {'speechACI_varnet2013', 'white', 'SAO-5000-trials'}, ... % 'SAO' for Leo
                    {'speechACI_varnet2013', 'SSN', 'osses2021c_S01'}, ... % 'S01' for Leo
                    {'speechACI_varnet2013', 'SSN', 'osses2021c_S02'}, ... % 'S02' for Leo
                    {'speechACI_Logatome-abda-S43M', 'white', 'SLV'}, ...
                    {'speechACI_Logatome-abda-S43M', 'white', 'SAO'}};
end

glmfct = 'l1glm'; %'classic_revcorr';%

switch glmfct
    case 'l1glm'
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
end

DimCI = 'gammatone';
flags_for_input = {'trialtype_analysis','total', ...'t1',...
    'N_folds', 10, ...
    'no_permutation', ...
    'no_bias', ...
    'no_plot', ... 'expvar_after_reversal',expvar_after_reversal, ...
    'lambda',Lambdas, ...
    'pyramid_script','imresize', ... % with the newer settings, we need to specify this explicitly
    'pyramid_shape',0 ...
    };

dir_out = '';
dir_out_figs = [pwd filesep];
if iswindows
    fullpath = 'C:\Users\LeoVarnet\ownCloud\Data\fastACI_data\';
else
    dir_out = '/home/alejandro/Documents/Databases/data/fastACI_data_TST/publ_varnet2022b_CFA/';
    fullpath = fastACI_paths('dir_data');
    flags_for_input{end+1} = 'dir_out';
    flags_for_input{end+1} = dir_out;
    dir_out_figs = [dir_out 'Figures-new' filesep]; mkdir(dir_out_figs);
end

if flags.do_table2
    
    lvltot_with_sil = [];
    for k = 1:2
        switch k
            case 1
                %%% ABDA13 + ABDA21:
                dir2look = [fastACI_paths('dir_data') 'speechACI_varnet2013' filesep];% '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_varnet2013/';
                dir_whereF_A = [dir2look 'SLeo' filesep 'speech-samples' filesep]; % Aba_F.txt'];
                idx_start = 5;
                idx1 = 15; % manually observed
                idx2 = 28; % manually observed
                idx_end = 42;
            case 2
                %%% ABDA22:
                dir2look = [fastACI_paths('dir_data') 'speechACI_Logatome-abda-S43M' filesep];
                dir_whereF_A = [dir2look 'SLV' filesep 'speech-samples' filesep]; 
                idx_start = 5;
                idx1 = 15; % manually observed
                idx2 = 28; % manually observed
                idx_end = 55;
                %%%
        end
        filesF  = Get_formants_from_dir(dir_whereF_A);
        filesI  = Get_intensity_from_dir(dir_whereF_A);
        filesf0 = Get_f0_from_dir(dir_whereF_A);
        filesW = Get_filenames(dir_whereF_A,'*.wav');
        
        Nsounds = length(filesF);
        % minIforF = 75; % Minimum intensity to be added to the plot (relative value)

        I_min = 59;
        pars_here = [];
        count = 1;
        
        for i = 1:Nsounds
            %%% Loading the Praat results
            [t_f0{i},f0{i}] = Get_f0_from_txt([dir_whereF_A filesf0{i}]);
            [t_I{i} ,I{i}]  = Get_intensity_from_txt([dir_whereF_A filesI{i}]);
            [t_F{i} ,F{i}]  = Get_formants_from_txt([dir_whereF_A filesF{i}]); 
            
            minIforF = max(max(I{i})-20, I_min); % In case the user requests
            idxsNaN = find(I{i}<minIforF | isnan(I{i}));
            F{i}(idxsNaN,:) = nan;
            f0{i}(idxsNaN) = nan;
            %%%
            
            idx_all = [1:idx1 idx2:size(f0{i}(idx2:end),1)];
            me = [nanmean(f0{i}(1:idx1)) nanmean(f0{i}(idx2:end)) nanmean(f0{i}(idx_all))];
            st = [nanstd(f0{i}(1:idx1)) nanstd(f0{i}(idx2:end)) nanstd(f0{i}(idx_all))];
            pars_here(count,1:2:6) = me;
            pars_here(count,2:2:6) = st;
            count = count+1;
            
            for j = 1:size(F{i},2)
                me = [nanmean(F{i}(1:idx1,j)) nanmean(F{i}(idx2:end,j)) nanmean(F{i}(idx_all,j))];
                st = [nanstd(F{i}(1:idx1,j))  nanstd(F{i}(idx2:end,j))  nanstd(F{i}(idx_all,j))];
                pars_here(count,1:2:6) = me;
                pars_here(count,2:2:6) = st;
                count = count+1;
            end

            [insig,fs] = audioread([dir_whereF_A filesW{i}]);
            
            dBFS = 100; % 100 dB FS (-3 to 'correct for the silence'), 93.61; % 100;
            idxi_sa = round(idx_start/100*fs);
            idx1_sa = round(idx1/100*fs);
            idxf_sa = round(idx_end  /100*fs);
            
            lvl1 = rmsdb(insig(idxi_sa:idx1_sa)+eps)+dBFS;
            lvl2 =  rmsdb(insig(idx1_sa+1:idxf_sa)+eps)+dBFS;
            lvltot =  rmsdb(insig(idxi_sa:idxf_sa)+eps)+dBFS;
            lvltot_with_sil(end+1) =  rmsdb(insig+eps)+dBFS;
            
            me = [lvl1 lvl2 lvltot];
            pars_here(count, 1:2:6) = me;
            count = count+1;
            
            fprintf('Sound %s\n',filesf0{i});
            if i == Nsounds
                var2latex(pars_here);
            end
            % t_here = t_f0{i};
            % figure; 
            % plot(t_here,f0{i}); hold on;
            % plot(t_here,F{i},'b');
        end
        disp('')
        
    end
end

if flags.do_fig2
    % figure('Position', [50 50 650 650])
    Pos = [50 50 750 650];
    figure('Position', Pos)
    tspacing_style = 'none'; % 'Compact'
    tlayout = tiledlayout(4,4,'TileSpacing',tspacing_style, 'Padding',tspacing_style');

    for i_experiment = 1:length(experiments)
        nexttile

        dir_expe = [fullpath experiments{i_experiment}{1} filesep experiments{i_experiment}{3} filesep];
        switch experiments{i_experiment}{1}
            case 'modulationACI'
                dir_results = [dir_expe ];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
            otherwise
                dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        end
        % %%% Old code from Leo:
        % % (Leo, note that your code fails to detect the correct savegame if you 
        % %  have multiple conditions in the same folder - as white+MPS+bump).
        % cd(dir_results)
        filtextra = '';
        switch experiments{i_experiment}{1}
            case 'modulationACI'
                filtextra = 'swap*';
                % filtextra = ''; warning('No swap...')
                bAdd_formants = 0;
            case {'speechACI_varnet2013','speechACI_Logatome-abda-S43M'}
                filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
                bAdd_formants = 1;
        end

        % D = dir([dir_results 'savegame_*' filtextra '.mat']);% fname_results = D(end).name;
        % fname_results = D(end).name;
        % %%% End old code Leo

        %%% New code Alejandro:
        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        if length(files) > 1 
            Show_cell(files);
            bInput = input(['Choose from the above list. The expected choice is related to ''' experiments{i_experiment}{2} ''' noise: ']);
            files = files(bInput);
        end
        if isempty(files)
            files = Get_filenames(dir_results,['savegame_*.mat']);
            if length(files)>1
                Show_cell(files);
                bInput = input(['Choose from the above list. The expected choice is related to ''' experiments{i_experiment}{2} ''' noise: ']);
                files = files(bInput);
            end
        end

        switch experiments{i_experiment}{1}
            case 'modulationACI'
                % dir_data = fastACI_dir_data; % dir_data of this computer
                flags_here = {'dir_out',dir_out,'dir_data',fastACI_dir_data};
                switch experiments{i_experiment}{3}
                    case 'S4'
                        Subject_here = 'SA';
                    case 'S1'
                        Subject_here = 'SB';
                end
                [xx,xx,data] = publ_varnet2022b_CFA('fig2_mod22',Subject_here,flags_here{:});
                close; % closes the last figure

                ACI = data.ACI;
                cfg_ACI = data.cfg_ACI;
                results = data.results;

            otherwise
                fname_results = [dir_results files{1}];
                %%% End new code Alejandro
                [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, DimCI, glmfct, flags_for_input{:});
        end
        thres(i_experiment) = prctile(results.data_passation.expvar,50); % median

        if mod(i_experiment,2)==1
            % switch experiments{i_experiment}{3}
            %     % Every time participant Leo is found:
            %     case {'SLV', 'S01', 'SLeo', 'S4'}
                % plot targets
                files = Get_filenames(cfg_ACI.dir_target,'*.wav');
                fname1 = [cfg_ACI.dir_target files{1}];
                fname2 = [cfg_ACI.dir_target files{2}];
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

                G2 = G2'/(max(max(G2)));
                G1 = G1'/(max(max(G1)));

                bColour_bar = 'no'; % 'yes';
                switch experiments{i_experiment}{1}
                    case 'modulationACI'  
                        c_axis = [0.7 1];
                    case 'speechACI_varnet2013'
                        c_axis = [0 0.5];
                    case 'speechACI_Logatome-abda-S43M'
                        c_axis = [0 1];
                end
                opts_colourbar = {'NfrequencyTicks', 8, 'colorbar', bColour_bar,'caxis',c_axis};
                switch experiments{i_experiment}{1}
                    case 'modulationACI'  
                        % Inverting: placing the AM target first:
                        affichage_tf(G2, 'pow', t, fc, opts_colourbar{:}); hold on; % caxis([0 0.01]);
                    otherwise
                        affichage_tf(G1, 'pow', t, fc, opts_colourbar{:}); hold on; % caxis([0 0.01]);
                end
                if bAdd_formants
                    LW = 1;

                    %%% I still need to spot the exact values:
                    par_formants.timestep = 0.01; % positive timestep 0.01
                    par_formants.nformants = 4; % positive nformants 5

                    %%% Unsure:
                    % Formants
                    par_formants.maxformant = 5500; % positive maxformant 5500
                    par_formants.windowlength = 0.025;% 0.025 % positive windowlength 0.025
                    par_formants.dynamicrange = 30; % positive dynamic range 20

                    % F0
                    par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
                    par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
                    par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

                    % Before 4/11/2021, I_min set to 40 dB:
                    par_formants.I_min = 59;%75; %, arbitrary value

                    affichage_tf_add_Praat_metrics_one_sound(fname1,cfg_ACI,par_formants, '-', 'w',LW);
                end
                nexttile
                switch experiments{i_experiment}{1}
                    case 'modulationACI' 
                        affichage_tf(G1, 'pow', t, fc, opts_colourbar{:}); hold on; %, caxis([0 0.01]);
                    otherwise
                        affichage_tf(G2, 'pow', t, fc, opts_colourbar{:}); hold on; %, caxis([0 0.01]);
                end
                set(gca,'YTickLabel',[]);
                ylabel('');
                if bAdd_formants
                    affichage_tf_add_Praat_metrics_one_sound(fname2,cfg_ACI,[], '-', 'w',LW);
                end
                nexttile
        end

        % display ACI
        if mod(i_experiment,2) == 1
            flags_extra = {'colorbar','no'};
        else
            flags_extra = {'colorbar','yes'};
        end
        set(gca,'YTickLabel',[]);
        switch experiments{i_experiment}{1}
            case 'modulationACI'
                ACI_to_plot = ACI; %-squeeze(results.ACI);
                warning('Temporal')
            otherwise
                ACI_to_plot = ACI; % squeeze(results.ACI);
        end
        outs = affichage_tf(ACI_to_plot,'CI', 'cfg', cfg_ACI, 'NfrequencyTicks', 8, flags_extra{:}); hold on % figure; imagesc(cfg_ACI.t,1:length(cfg_ACI.f),ACI_to_plot)
        if isfield(outs,'tcolorbar')
            set(outs.tcolorbar,'Ticks',[]); % removes the 'Ticks'
            %title('ACI','interpreter','none')
        end
        set(gca,'YTickLabel',[]);
        ylabel('');

    end

    thres = [thres(1:2:end); thres(2:2:end)]; % SA and SB now in different rows

    FS_thres = 8;

    S01_Lab = 'SA';
    S02_Lab = 'SB';
    YLab = 'Frequency (Hz)';
    posX = 0.95; % position, normalised units
    posY = 0.9;  % position, normalised units
    x_cb    = 1.1;
    y_cb(2) = 1.05; % position for colourbar, normalised units
    y_cb(1) = -.05; % position for colourbar, normalised units

    flags_common = {'Unit','Normalized','FontWeight','bold','HorizontalAlignment', 'right'};
    ax=nexttile(1);colormap(ax,'hot'); xlabel('');ylabel(YLab); text(posX,posY,'mod. tone'  ,'Color','white',flags_common{:});
    ax=nexttile(2);colormap(ax,'hot'); xlabel('');ylabel('');   text(posX,posY,'unmod. tone','Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=nexttile(3);
    xlabel('');ylabel('');
    i_subj = 1;
    i_exp = 1;
    text(posX,posY    , S01_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    ax=nexttile(4);
    xlabel('');ylabel('');
    i_subj = 2;
    text(posX,posY    , S02_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    text(x_cb,y_cb(2),'mod'  ,'Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'unmod','Color','blue','Units','Normalized');

    % set(gca,'YTickLabel',[]);
    ax=nexttile(5);colormap(ax,'hot'); 
    xlabel('');ylabel(YLab); 
    text(posX,posY,'aba'        ,'Color','white',flags_common{:});
    ax=nexttile(6);colormap(ax,'hot'); 
    xlabel('');ylabel('');   
    text(posX,posY,'ada'        ,'Color','white',flags_common{:});

    % set(gca,'YTickLabel',[]);
    ax=nexttile(7);
    xlabel('');ylabel('');   
    i_subj = 1;
    i_exp  = 2; 
    text(posX,posY     , S01_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    ax=nexttile(8);
    xlabel('');ylabel('');   
    i_subj = 2;
    text(posX,posY    , S02_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ada','Color','blue','Units','Normalized');

    % set(gca,'YTickLabel',[]);
    ax=nexttile(9);colormap(ax,'hot'); xlabel('');ylabel(YLab); text(posX,posY,'aba'        ,'Color','white',flags_common{:});
    ax=nexttile(10);colormap(ax,'hot');xlabel('');ylabel('');   text(posX,posY,'ada'        ,'Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=nexttile(11);                   
    xlabel('');ylabel('');   
    text(posX,posY, S01_Lab                     ,flags_common{:});
    i_subj = 1;
    i_exp  = 3; 
    text(posX,posY     , S01_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    ax=nexttile(12);
    xlabel('');ylabel('');
    text(posX,posY, S02_Lab                     ,flags_common{:});
    i_subj = 2;
    text(posX,posY     , S02_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ada','Color','blue','Units','Normalized');

    ax=nexttile(13);colormap(ax,'hot');           ylabel(YLab); text(posX,posY,'aba'        ,'Color','white',flags_common{:});
    ax=nexttile(14);colormap(ax,'hot');           ylabel('');   text(posX,posY,'ada'        ,'Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=nexttile(15);
    ylabel('');   
    text(posX,posY, S01_Lab                     ,flags_common{:});
    i_subj = 1;
    i_exp  = 4; 
    text(posX,posY     , S01_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    ax=nexttile(16);
    ylabel('');   
    text(posX,posY, S02_Lab                     ,flags_common{:});
    i_subj = 2;
    text(posX,posY     , S02_Lab, flags_common{:});
    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ada','Color','blue','Units','Normalized');

    nexttile(14);
    text(0.52,0.08,'f_0','Color','white','Units','Normalized');
    text(0.52,0.35,'F_1','Color','white','Units','Normalized');
    text(0.52,0.48,'F_2','Color','white','Units','Normalized');
    text(0.52,0.62,'F_3','Color','white','Units','Normalized');
    text(0.52,0.75,'F_4','Color','white','Units','Normalized');

    tAlignment = 'left';
    FS = 10;
    ax=nexttile(1);  title('MOD22' ,'FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(5);  title('ABDA13','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(9);  title('ABDA21','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(13); title('ABDA22','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);

    fout = [dir_out_figs 'Figure2'];
    % print(fout,'-dpdf')
    opts = [];
    opts.format = 'epsc'; % Pos34 = [650 800]; Pos = get(gcf,'Position'); Pos(3:4) = Pos34; set(gcf,'Position',Pos);
    Saveas(gcf,fout,opts);

    opts = [];
    opts.format = 'fig';
    Saveas(gcf,fout,opts);

    % %% Plot two examples of WN and SSN
    % f_limits = [0 4050]; % Hz, the exact bin is at 4048 Hz (see the paper, p. 4, top)
    % t_limits = [0 0.3425]; % s, first 0.34 s were accounted for (see the paper, p. 3)
    % 
    % %%% ABDA13 / ABDA21
    % dir_target = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Stimuli/varnet2013/';
    % files = {'Aba.wav','Ada.wav'};
    % %%
    % 
    % % %%% ABDA22
    % % dir_target = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Stimuli/Logatome/';
    % % files = {'S43M_ab_ba_eq.wav','S43M_ad_da_eq.wav'};
    % % %%%
    % 
    % opts = [];
    % opts.window = 'hamming';
    % 
    % [T_dB,f_spec,t_spec] = Time_frequency_converter(dir_target,files,length(files),opts);
    % %%%
    % tf_s = t_limits(end); 
    % ff   = f_limits(end); 
    % 
    % idx_t = find(round(100*t_spec)/100 <= tf_s);
    % idx_f = find(round(f_spec) <= ff);
    % f_spec = f_spec(idx_f);
    % t_spec = t_spec(idx_t);
    % T_dB   = T_dB(idx_f,idx_t,:);
    % %%%
    % 
    % max_dB = max(max(max(T_dB)));
    % 
    % figure;
    % subplot(1,3,1)
    % plot_stft(t_spec,f_spec,T_dB(:,:,1)-max_dB);
    % title('/aba/');
    % 
    % % figure;
    % subplot(1,3,2)
    % plot_stft(t_spec,f_spec,T_dB(:,:,2)-max_dB);
    % title('/ada/');
    % 
    % % figure;
    % subplot(1,3,3)
    % plot_stft(t_spec,f_spec,T_dB(:,:,1)-T_dB(:,:,2));
    % title('/aba/-/ada/');
    % 
    % Pos = get(gcf,'Position');
    % Pos(3) = 2.5*Pos(3); % widening the figure
    % set(gcf,'Position',Pos);
    % 
    % h(end+1) = gcf;
    % hname{end+1} = 'fig1-aba-ada-diff';
end