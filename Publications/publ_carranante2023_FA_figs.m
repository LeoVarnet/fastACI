function publ_carranante2023_FA_figs(varargin)
% function publ_carranante2023_FA_figs(varargin)
%
% This script generates the figures 2 to 4 from the publication by 
%   Carranante, G.; Giavazzi, M.; and Varnet, L. (2023). Auditory reverse 
%   correlation applied to the study of place and voicing: Four new phoneme-
%   discrimination tasks. To be presented at Forum Acusticum.
% This script reads the savegame (MAT) files for participant S01 and S13. 
%   the files are attempted to be retrieved from the local fastACI_dir_data
%   path (where the noises are expected to be in the folder 'NoiseStim-bumpv1p2_10dB',
%   these folders are strictly needed for Fig. 4). Otherwise, the savegame
%   files will be retrieved from the fastACI repository (Publications/publ_
%   carranante2023/). However, we (the fastACI team) have not enabled yet 
%   the automatic waveform re-generation. This means that in the latter case
%   only 'fig2' and 'fig3' will be recreated.
%   As a temporal solution, contact the fastACI team to obtain the noise 
%   waveforms to obtain the corresponding ACIs.
%
% % To display Fig. 2 of Carranante and Varnet (2023) use :::
%     publ_carranante2023_FA_figs('fig2');
%
% % To display Fig. 3 of Carranante and Varnet (2023) use :::
%     publ_carranante2023_FA_figs('fig3');
%
% % To display Fig. 4 of Carranante and Varnet (2023) use :::
%     publ_carranante2023_FA_figs('fig4');
% 
% % If the data are locally stored at a different location, specify the keyval
% %   'dir_data'. For instance, for figure 4:
%     dir_data = 'C:\Users\LeoVarnet\ownCloud\Data\fastACI_data\';
%     publ_carranante2023_FA_figs('fig4','dir_data',dir_data);
%
% Author: Leo Varnet
% See also: publ_osses2022d_ICA_figs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc

if nargin == 0
    help publ_carranante2023_FA_figs;
    return
end

% data = [];
% h = [];
% hname = [];

definput.flags.type={'missingflag','fig2','fig3','fig4'};
definput.keyvals.dir_data=[];
% definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% Leo's local names
experiments = { ...
    {'speechACI_Logatome-abda-S43M', 'bumpv1p2_10dB', 'S01'}, ...
    {'speechACI_Logatome-abda-S43M', 'bumpv1p2_10dB', 'S13'}, ...
    {'speechACI_Logatome-adga-S43M', 'bumpv1p2_10dB', 'S01'}, ...
    {'speechACI_Logatome-adga-S43M', 'bumpv1p2_10dB', 'S13'}, ...
    {'speechACI_Logatome-apta-S43M', 'bumpv1p2_10dB', 'S01'}, ...
    {'speechACI_Logatome-apta-S43M', 'bumpv1p2_10dB', 'S13'}, ...
    {'speechACI_Logatome-abpa-S43M', 'bumpv1p2_10dB', 'S01'}, ...
    {'speechACI_Logatome-abpa-S43M', 'bumpv1p2_10dB', 'S13'}, ...
    {'speechACI_Logatome-adta-S43M', 'bumpv1p2_10dB', 'S01'}, ...
    {'speechACI_Logatome-adta-S43M', 'bumpv1p2_10dB', 'S13'}};

glmfct = 'l1glm'; 
% glmfct = 'classic_revcorr'; % for quick prototyping (but not used in the publication)

switch glmfct
    case 'l1glm'
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
        flags_for_input = {'trialtype_analysis','total', ...'t1',...
            'N_folds', 10, ...
            'no_bias', ...
            'no_plot', ...
            'expvar_after_reversal',4, ...
            'lambda',Lambdas, ...
            'pyramid_script','imgaussfilt' ...
            };
    % case 'classic_revcorr'
    %     flags_for_input = {'trialtype_analysis','total', ...'t1',...
    %         'N_folds', 10, ...
    %         'no_bias', ...
    %         'no_plot', ...
    %         'expvar_after_reversal',4, ...
    %         };
end

conditioncolors = {'b','r','g','m','k'};
subjectsymbol = {'o','square'};
DimCI = 'gammatone';

% dir_out = '';
% dir_out_figs = [pwd filesep]; % current folder

if ~isempty(keyvals.dir_data)
    dir_data_fullpath = keyvals.dir_data;
else
    % if empty, the dir_data is set to the automatic data folder for fastACI:
    dir_data_fullpath = fastACI_paths('dir_data');
end

dir_publ = [fastACI_basepath 'Publications' filesep 'publ_carranante2023' filesep];

% if iswindows
%     % Leo's local folder
%     dir_data_fullpath = 'C:\Users\LeoVarnet\ownCloud\Data\fastACI_data\'; % Leo's local folder
% else
%     % Alejandro's local folder:
%     dir_out = '/home/alejandro/Documents/Databases/data/fastACI_data_TST/publ_varnet2022b_CFA/';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2
    Pos = [50 50 450 200];
    figure('Position', Pos)
    for i_experiment = 1:length(experiments)
        subj_id = experiments{i_experiment}{3};
        dir_expe = [dir_data_fullpath experiments{i_experiment}{1} filesep subj_id filesep];
        dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        
        filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        
        if isempty(files)
            % Then the savegames are not locally: we will load the stored results then...
            dir_results = [dir_publ 'data_' subj_id filesep '1-experimental_results' filesep];
            files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        end
        
        if ~isempty(files)
            % Now it should always work
            fname_results = [dir_results files{1}];
            load(fname_results,'data_passation');
            r = Get_mAFC_reversals(data_passation.expvar);
            SNRthres(i_experiment) = median(r(5:end));
            SNRsd(i_experiment) = std(r(5:end));
        end
    end
    Xpos = [0.9 1.1 1.9 2.1 2.9 3.1 3.9 4.1 4.9 5.1];
    for i_experiment = 1:length(experiments)
        errorbar(Xpos(i_experiment),SNRthres(i_experiment),SNRsd(i_experiment)/2,'LineStyle','none','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1});hold on
    end
    %title('SNR thresholds (\pm 0.5 SD)');
    set(gca,'XTick',1:5)
    set(gca,'XTickLabels',{'ABDA22','ADGA23','APTA23','ABPA23','ADTA23'})
    xlabel('experiment'); 
    ylabel('SRT (dB)');
    legend({'','','','','','','','','S01','S13'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig3
    Pos = [50 50 450 600];
    figure('Position', Pos);
    tlayout = tiledlayout(5,2,'TileSpacing','compact', 'Padding','compact');
    
    for i_experiment = 1:length(experiments)
        subj_id = experiments{i_experiment}{3};
        cond = experiments{i_experiment}{1};
        dir_expe = [dir_data_fullpath cond filesep subj_id filesep];
        dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        
        filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        
        if isempty(files)
            % Then the savegames are not locally: we will load the stored results then...
            dir_results = [dir_publ 'data_' subj_id filesep '1-experimental_results' filesep];
            files = Get_filenames(dir_results,['savegame_*' cond '*' filtextra '.mat']);
        end
        
        if ~isempty(files)
            fname_results = [dir_results files{1}];
            cfg_game = []; % loaded next
            data_passation = []; % loaded next
            load(fname_results);
            
            data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation,'histogram-leo',0);
            
            bin_centres = data_hist.bin_centres;
            P_H  = data_hist.H ./(data_hist.H +data_hist.M );
            P_FA = data_hist.FA./(data_hist.CR+data_hist.FA);
            PC   = (data_hist.H + data_hist.CR) ./(data_hist.H + data_hist.M + data_hist.CR + data_hist.FA);
            bias = (data_hist.M + data_hist.CR) ./(data_hist.H + data_hist.M + data_hist.CR + data_hist.FA);
            dprime = norminv(P_H)-norminv(P_FA);           % Eq.  9 from Harvey2004
            criterion = -(norminv(P_H) + norminv(P_FA))/2; % Eq. 12 from Harvey2004
            
            % plot performance as a function of SNR
            P_H_resp2 = P_H; % Script3 assumes that the target sound is option '2'
            resp2_label = [data_hist.H_label '-trials (H)'];
            
            P_H_resp1 = 1-P_FA; % 'Hits' when the participant response was '1'
            resp1_label = [data_hist.CR_label '-trials (CR)'];
            
            nexttile(1+floor((i_experiment-1)/2)*2);
            plot(bin_centres, dprime ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
            ylabel('d'''); ylim([0 3]); xlim([-20 -8]); 
            %plot([-20 -8], [0 0] ,'-','Color',[0.8 0.8 0.8]); grid on; hold on
            %set(gca, 'YTick', [0 0.5 1])
            
            nexttile(2+floor((i_experiment-1)/2)*2);
            plot(bin_centres, criterion ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
            ylabel('c'); ylim([-0.7 0.7]); xlim([-20 -8]); 
            %plot([-20 -8], [0 0] ,'-','Color',[0.8 0.8 0.8]); grid on; hold on
            %set(gca, 'YTick', [0 0.5 1])
        end
    end
    
    nexttile(9); 
    xlabel('SNR (dB)');
    
    nexttile(10); 
    xlabel('SNR (dB)');
    
    tAlignment = 'left';
    FontSize = 10;
    ax=nexttile(1); 
    title('ABDA22','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(3); 
    title('ADGA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(5); 
    title('APTA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(7); 
    title('ABPA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(9); 
    title('ADTA23','FontSize',FontSize,'FontWeight','bold');
    set(ax,'TitleHorizontalAlignment',tAlignment);
        
    legend({'S01','S13'})
end

if flags.do_fig4
    % figure('Position', [50 50 650 650])
    Pos = [50 50 750 950];
    figure('Position', Pos)
    tlayout = tiledlayout(5,4,'TileSpacing','compact', 'Padding','none');

    for i_experiment = 1:length(experiments)
        % curr_tile = curr_tile + 1;
        nexttile; % (curr_tile);
    
        dir_expe = [dir_data_fullpath experiments{i_experiment}{1} filesep experiments{i_experiment}{3} filesep];
        dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        
        filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
        bAdd_formants = 1;

        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        if length(files) > 1 
            Show_cell(files);
            bInput = input(['Choose from the above list. The expected choice is related to ''' experiments{i_experiment}{2} ''' noise: ']);
            files = files(bInput);
        end
        if isempty(files)
            files = Get_filenames(dir_results,'savegame_*.mat');
            if length(files)>1
                Show_cell(files);
                bInput = input(['Choose from the above list. The expected choice is related to ''' experiments{i_experiment}{2} ''' noise: ']);
                files = files(bInput);
            end
        end

        if ~isempty(files)
            fname_results = [dir_results files{1}];
            %%% End new code Alejandro
            [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, DimCI, glmfct, flags_for_input{:});

            thres(i_experiment) = prctile(results.data_passation.expvar,50); % median

            if mod(i_experiment,2)==1

                k = strfind(cfg_ACI.dir_target,'fastACI_data');

                files = Get_filenames([dir_data_fullpath cfg_ACI.dir_target(k+13:end)],'*.wav');
                fname1 = [dir_data_fullpath cfg_ACI.dir_target(k+13:end) '\' files{1}];
                fname2 = [dir_data_fullpath cfg_ACI.dir_target(k+13:end) '\'  files{2}];
                [T1,fs] = audioread(fname1);
                T2      = audioread(fname2);

                undersampling = 100;
                basef = 8000;
                flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5,'dboffset',100,'no_adt','binwidth',undersampling/fs, ...
                    'no_outerear','no_middleear'};

                [G1]              = Gammatone_proc(T1, fs, flags_gamma{:});
                [G2, fc, t, outs] = Gammatone_proc(T2, fs, flags_gamma{:});

                G2 = G2'/(max(max(G2)));
                G1 = G1'/(max(max(G1)));

                bColour_bar = 'no'; % 'yes';
                c_axis = [0 1];
                opts_colourbar = {'NfrequencyTicks', 8, 'colorbar', bColour_bar,'caxis',c_axis};
                affichage_tf(G1, 'pow', t, fc, opts_colourbar{:}); hold on; % caxis([0 0.01]);

                if bAdd_formants
                    LW = 1;
                    par_formants = {};
                    outs=affichage_tf_add_Praat_metrics_one_sound(fname1,cfg_ACI,par_formants, '-', 'w',LW);
                end
                nexttile
                affichage_tf(G2, 'pow', t, fc, opts_colourbar{:}); hold on; %, caxis([0 0.01]);
                set(gca,'YTickLabel',[]);
                ylabel('');
                if bAdd_formants
                    outs=affichage_tf_add_Praat_metrics_one_sound(fname2,cfg_ACI,[], '-', 'w',LW);
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

            ACI_to_plot = ACI; % squeeze(results.ACI);

            outs = affichage_tf(ACI_to_plot,'CI', 'cfg', cfg_ACI, 'NfrequencyTicks', 8, flags_extra{:}); hold on % figure; imagesc(cfg_ACI.t,1:length(cfg_ACI.f),ACI_to_plot)
            if isfield(outs,'tcolorbar')
                set(outs.tcolorbar,'Ticks',[]); % removes the 'Ticks'
                %title('ACI','interpreter','none')
            end
            set(gca,'YTickLabel',[]);
            ylabel('');
        end
    end

    FS_thres = 8;

    S01_Lab = 'S1';
    S02_Lab = 'S13';
    YLab = 'Frequency (Hz)';
    posX = 0.95; % position, normalised units
    posY = 0.9;  % position, normalised units
    x_cb    = 1.07;
    y_cb(2) = 1.07; % position for colourbar, normalised units
    y_cb(1) = -.05; % position for colourbar, normalised units

    flags_common = {'Unit','Normalized','FontWeight','bold','HorizontalAlignment', 'right'};
    ax=nexttile(1);
    colormap(ax,'hot'); 
    xlabel('');
    ylabel(YLab); 
    text(posX,posY,'aba','Color','white',flags_common{:});
    
    ax=nexttile(2);
    colormap(ax,'hot'); 
    xlabel('');
    ylabel('');   
    text(posX,posY,'ada','Color','white',flags_common{:});
    
    ax=nexttile(3);
    xlabel('');ylabel('');
    i_subj = 1;
    i_exp = 1;
    text(posX,posY, S01_Lab, flags_common{:});
    
    ax=nexttile(4);
    xlabel('');ylabel('');
    i_subj = 2;
    text(posX,posY, S02_Lab, flags_common{:});

    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ada','Color','blue','Units','Normalized');

    ax=nexttile(5);
    colormap(ax,'hot'); 
    xlabel('');
    ylabel(YLab); 
    text(posX,posY,'ada','Color','white',flags_common{:});
    
    ax=nexttile(6);
    colormap(ax,'hot'); 
    xlabel('');
    ylabel('');   
    text(posX,posY,'aga','Color','white',flags_common{:});

    ax=nexttile(7);
    xlabel('');ylabel('');   
    i_subj = 1;
    i_exp  = 2; 
    text(posX,posY, S01_Lab, flags_common{:});

    ax=nexttile(8);
    xlabel('');ylabel('');   
    i_subj = 2;
    text(posX,posY, S02_Lab, flags_common{:});

    text(x_cb,y_cb(2),'ada','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'aga','Color','blue','Units','Normalized');

    ax=nexttile(9);
    colormap(ax,'hot'); 
    xlabel(''); % removing the xlabel
    ylabel(YLab); 
    text(posX,posY,'apa','Color','white',flags_common{:});
    
    ax=nexttile(10);
    colormap(ax,'hot');
    xlabel(''); % removing the xlabel
    ylabel(''); % removing the ylabel
    text(posX,posY,'ata','Color','white',flags_common{:});
    
    ax=nexttile(11);                   
    xlabel('');ylabel('');   
    text(posX,posY, S01_Lab,flags_common{:});
    i_subj = 1;
    i_exp  = 3; 
    text(posX,posY, S01_Lab, flags_common{:});

    ax=nexttile(12);
    xlabel('');ylabel('');
    text(posX,posY, S02_Lab, flags_common{:});
    i_subj = 2;
    text(posX,posY, S02_Lab, flags_common{:});

    text(x_cb,y_cb(2),'apa','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ata','Color','blue','Units','Normalized');

    ax=nexttile(13);
    colormap(ax,'hot'); 
    xlabel('');
    ylabel(YLab); 
    text(posX,posY,'aba','Color','white',flags_common{:});
    
    ax=nexttile(14);
    colormap(ax,'hot');
    xlabel('');
    ylabel(''); 
    text(posX,posY,'apa','Color','white',flags_common{:});
    
    ax=nexttile(15);                   
    xlabel('');
    ylabel('');   
    text(posX,posY, S01_Lab, flags_common{:});
    i_subj = 1;
    i_exp  = 3; 
    text(posX,posY, S01_Lab, flags_common{:});

    ax=nexttile(16);
    xlabel('');
    ylabel('');
    text(posX,posY, S02_Lab, flags_common{:});
    i_subj = 2;
    text(posX,posY, S02_Lab, flags_common{:});

    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'apa','Color','blue','Units','Normalized');

    ax=nexttile(17);colormap(ax,'hot'); 
    ylabel(YLab); 
    text(posX,posY,'ada','Color','white',flags_common{:});
    
    ax=nexttile(18);colormap(ax,'hot'); 
    ylabel('');   
    text(posX,posY,'ata','Color','white',flags_common{:});
    
    ax=nexttile(19);
    ylabel('');   
    text(posX,posY, S01_Lab ,flags_common{:});
    i_subj = 1;
    i_exp  = 4; 
    text(posX,posY, S01_Lab, flags_common{:});

    ax=nexttile(20);
    ylabel('');   
    text(posX,posY, S02_Lab, flags_common{:});
    i_subj = 2;
    text(posX,posY, S02_Lab, flags_common{:});

    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'ada','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ata','Color','blue','Units','Normalized');

    nexttile(18);
    text(0.52,0.08,'f_0','Color','white','Units','Normalized');
    text(0.52,0.3 ,'F_1','Color','white','Units','Normalized');
    text(0.52,0.48,'F_2','Color','white','Units','Normalized');
    text(0.52,0.62,'F_3','Color','white','Units','Normalized');
    text(0.52,0.75,'F_4','Color','white','Units','Normalized');

    tAlignment = 'left';
    FontSize = 10;
    ax=nexttile(1);  
    title('ABDA22','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(5);  
    title('ADGA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(9);  
    title('APTA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(13); 
    title('ABPA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);
    
    ax=nexttile(17); 
    title('ADTA23','FontSize',FontSize,'FontWeight','bold'); 
    set(ax,'TitleHorizontalAlignment',tAlignment);

    % fout = [dir_out_figs 'Figure2'];
    % opts = [];
    % opts.format = 'epsc'; % Pos34 = [650 800]; Pos = get(gcf,'Position'); Pos(3:4) = Pos34; set(gcf,'Position',Pos);
    % opts = [];
    % opts.format = 'fig';
end