function publ_carranante2023_FA_figs(varargin)
% function publ_carranante2023_FA_figs(varargin)
%
% % To display Fig. 1 of Carranante and Varnet (2023) use :::
%     publ_osses2022d_ICA_figs('fig1');
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc

if nargin == 0
    help publ_carranante2023_FA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2','fig3','fig4',};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

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

glmfct = 'l1glm'; %'classic_revcorr';%

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
    case 'classic_revcorr'
        flags_for_input = {'trialtype_analysis','total', ...'t1',...
            'N_folds', 10, ...
            'no_bias', ...
            'no_plot', ...
            'expvar_after_reversal',4, ...
            };
end

conditioncolors = {'b','r','g','m','k'};
subjectsymbol = {'o','square'};
DimCI = 'gammatone';

dir_out = '';
dir_out_figs = [pwd filesep];
if iswindows
    fullpath = 'C:\Users\LeoVarnet\ownCloud\Data\fastACI_data\'; % Leo's local folder
else
    % Alejandro's local folder:
    dir_out = '/home/alejandro/Documents/Databases/data/fastACI_data_TST/publ_varnet2022b_CFA/';
    fullpath = fastACI_paths('dir_data');
end

if exist(dir_out,'dir')
    flags_for_input{end+1} = 'dir_out';
    flags_for_input{end+1} = dir_out;
    
    dir_out_figs = [dir_out 'Figures-new' filesep]; 
    if ~exist(dir_out_figs,'dir')
        mkdir(dir_out_figs);
    end
else
    dir_out = ''; % making sure that a non-existing dir_out is set to empty (i.e., defaults will be used)
end

if flags.do_fig2
    Pos = [50 50 450 200];
    figure('Position', Pos)
    for i_experiment = 1:length(experiments)
        dir_expe = [fullpath experiments{i_experiment}{1} filesep experiments{i_experiment}{3} filesep];
        dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        
        filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        
        if ~isempty(files)
        fname_results = [dir_results files{1}];
        load(fname_results)
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
    xlabel('experiment'); ylabel('SRT (dB)');
    legend({'','','','','','','','','S01','S13'})
end


if flags.do_fig3
    Pos = [50 50 450 600];
    figure('Position', Pos);
    tlayout = il_tiledlayout(5,2,'TileSpacing','compact', 'Padding','compact');
    
    for i_experiment = 1:length(experiments)
        dir_expe = [fullpath experiments{i_experiment}{1} filesep experiments{i_experiment}{3} filesep];
        dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
        
        filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
        files = Get_filenames(dir_results,['savegame_*' filtextra '.mat']);
        
         % (curr_tile);
        
        if ~isempty(files)
            fname_results = [dir_results files{1}];
            load(fname_results)
            
            
            data_hist = Script3_AnalysisComplex_functions(cfg_game,data_passation,'histogram-leo',0);
            
            bin_centres = data_hist.bin_centres;
            P_H  = data_hist.H ./(data_hist.H +data_hist.M );
            P_FA = data_hist.FA./(data_hist.CR+data_hist.FA);
            PC = (data_hist.H + data_hist.CR) ./(data_hist.H + data_hist.M + data_hist.CR + data_hist.FA);
            bias = (data_hist.M + data_hist.CR) ./(data_hist.H + data_hist.M + data_hist.CR + data_hist.FA);
            dprime = norminv(P_H)-norminv(P_FA);           % Eq.  9 from Harvey2004
            criterion = -(norminv(P_H) + norminv(P_FA))/2; % Eq. 12 from Harvey2004
            
            % plot performance as a function of SNR
            P_H_resp2 = P_H; % Script3 assumes that the target sound is option '2'
            resp2_label = [data_hist.H_label '-trials (H)'];
            
            P_H_resp1 = 1-P_FA; % 'Hits' when the participant response was '1'
            resp1_label = [data_hist.CR_label '-trials (CR)'];
            
%             il_nexttile(1+floor((i_experiment-1)/2)*2);
%             plot(bin_centres, PC ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
%             ylabel('percent correct'); ylim([0.5 1]); xlim([-20 -10]); 
%             
%             il_nexttile(2+floor((i_experiment-1)/2)*2);
%             plot(bin_centres, bias ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
%             ylabel('bias'); ylim([0 1]); plot([-20 -10], [0.5 0.5] ,'-','Color',[0.8 0.8 0.8]); grid on; hold on
%             set(gca, 'YTick', [0 0.5 1])
               
            il_nexttile(1+floor((i_experiment-1)/2)*2);
            plot(bin_centres, dprime ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
            ylabel('d'''); ylim([0 3]); xlim([-20 -8]); 
            %plot([-20 -8], [0 0] ,'-','Color',[0.8 0.8 0.8]); grid on; hold on
            %set(gca, 'YTick', [0 0.5 1])
            
            il_nexttile(2+floor((i_experiment-1)/2)*2);
            plot(bin_centres, criterion ,'-','Color',conditioncolors{ceil(i_experiment/2)},'Marker',subjectsymbol{mod(i_experiment+1,2)+1}); grid on; hold on
            ylabel('c'); ylim([-0.7 0.7]); xlim([-20 -8]); 
            %plot([-20 -8], [0 0] ,'-','Color',[0.8 0.8 0.8]); grid on; hold on
            %set(gca, 'YTick', [0 0.5 1])
        end
    end
    
    il_nexttile(9); xlabel('SNR (dB)');
    il_nexttile(10); xlabel('SNR (dB)');
    
    tAlignment = 'left';
    FS = 10;
    ax=nexttile(1);  title('ABDA22','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(3);  title('ADGA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(5);  title('APTA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(7); title('ABPA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(9); title('ADTA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
%     ax=nexttile(2); set(gca, 'YTickLabel', {'da', 0, 'ba'});
%     ax=nexttile(4); set(gca, 'YTickLabel', {'ga', 0, 'da'});
%     ax=nexttile(6); set(gca, 'YTickLabel', {'ta', 0, 'pa'});
%     ax=nexttile(8); set(gca, 'YTickLabel', {'pa', 0, 'ba'});
%     ax=nexttile(10); set(gca, 'YTickLabel', {'ta', 0, 'da'});
        
    legend({'S01','S13'})
end

if flags.do_fig4
    % figure('Position', [50 50 650 650])
    Pos = [50 50 750 950];
    figure('Position', Pos)
    tlayout = il_tiledlayout(5,4,'TileSpacing','compact', 'Padding','none');

    % curr_tile = 0;
    for i_experiment = 1:length(experiments)
        % curr_tile = curr_tile + 1;
        il_nexttile; % (curr_tile);
    
        dir_expe = [fullpath experiments{i_experiment}{1} filesep experiments{i_experiment}{3} filesep];
%         switch experiments{i_experiment}{1}
%             case 'modulationACI'
%                 dir_results = [dir_expe ];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
%             otherwise
                dir_results = [dir_expe 'Results' filesep];%['C:\Users\LSP005\ownCloud\Data\Projet ModulationACI\AM_4Hz_75ms_m\A_king2019' filesep participant];
%         end
        % %%% Old code from Leo:
        % % (Leo, note that your code fails to detect the correct savegame if you 
        % %  have multiple conditions in the same folder - as white+MPS+bump).
        % cd(dir_results)

                filtextra = [experiments{i_experiment}{2} '*']; % using the condition as filter
                bAdd_formants = 1;

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

        if ~isempty(files)
        fname_results = [dir_results files{1}];
        %%% End new code Alejandro
        [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, DimCI, glmfct, flags_for_input{:});

        thres(i_experiment) = prctile(results.data_passation.expvar,50); % median

        if mod(i_experiment,2)==1
            % switch experiments{i_experiment}{3}
            %     % Every time participant Leo is found:
            %     case {'SLV', 'S01', 'SLeo', 'S4'}
                % plot targets
                
                k = strfind(cfg_ACI.dir_target,'fastACI_data');
                
                files = Get_filenames([fullpath cfg_ACI.dir_target(k+13:end)],'*.wav');
                fname1 = [fullpath cfg_ACI.dir_target(k+13:end) '\' files{1}];
                fname2 = [fullpath cfg_ACI.dir_target(k+13:end) '\'  files{2}];
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
%                     %%% I still need to spot the exact values:
%                     par_formants.timestep = 0.01; % positive timestep 0.01
%                     par_formants.nformants = 4; % positive nformants 5
% 
%                     %%% Unsure:
%                     % Formants
%                     par_formants.maxformant = 5500; % positive maxformant 5500
%                     par_formants.windowlength = 0.025;% 0.025 % positive windowlength 0.025
%                     par_formants.dynamicrange = 30; % positive dynamic range 20
% 
%                     % F0
%                     par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
%                     par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
%                     par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)
% 
%                     % Before 4/11/2021, I_min set to 40 dB:
%                     par_formants.I_min = 59;%75; %, arbitrary value

                    outs=affichage_tf_add_Praat_metrics_one_sound(fname1,cfg_ACI,par_formants, '-', 'w',LW);
                end
                il_nexttile
                affichage_tf(G2, 'pow', t, fc, opts_colourbar{:}); hold on; %, caxis([0 0.01]);
                set(gca,'YTickLabel',[]);
                ylabel('');
                if bAdd_formants
                    outs=affichage_tf_add_Praat_metrics_one_sound(fname2,cfg_ACI,[], '-', 'w',LW);
                end
                il_nexttile
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

    
  %  thres = [thres(1:2:end); thres(2:2:end)]; % SA and SB now in different rows

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
    ax=il_nexttile(1);colormap(ax,'hot'); xlabel('');ylabel(YLab); text(posX,posY,'aba','Color','white',flags_common{:});
    ax=il_nexttile(2);colormap(ax,'hot'); xlabel('');ylabel('');   text(posX,posY,'ada','Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(3);
    xlabel('');ylabel('');
    i_subj = 1;
    i_exp = 1;
    text(posX,posY    , S01_Lab, flags_common{:});
   % text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    ax=il_nexttile(4);
    xlabel('');ylabel('');
    i_subj = 2;
    text(posX,posY    , S02_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ada','Color','blue','Units','Normalized');

    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(5);colormap(ax,'hot'); 
    xlabel('');ylabel(YLab); 
    text(posX,posY,'ada'        ,'Color','white',flags_common{:});
    ax=il_nexttile(6);colormap(ax,'hot'); 
    xlabel('');ylabel('');   
    text(posX,posY,'aga'        ,'Color','white',flags_common{:});

    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(7);
    xlabel('');ylabel('');   
    i_subj = 1;
    i_exp  = 2; 
    text(posX,posY     , S01_Lab, flags_common{:});
%   text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    ax=il_nexttile(8);
    xlabel('');ylabel('');   
    i_subj = 2;
    text(posX,posY    , S02_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres, flags_common{:});

    text(x_cb,y_cb(2),'ada','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'aga','Color','blue','Units','Normalized');

    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(9);colormap(ax,'hot'); xlabel('');ylabel(YLab); text(posX,posY,'apa'        ,'Color','white',flags_common{:});
    ax=il_nexttile(10);colormap(ax,'hot');xlabel('');ylabel('');   text(posX,posY,'ata'        ,'Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(11);                   
    xlabel('');ylabel('');   
    text(posX,posY, S01_Lab                     ,flags_common{:});
    i_subj = 1;
    i_exp  = 3; 
    text(posX,posY     , S01_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});
%
    ax=il_nexttile(12);
    xlabel('');ylabel('');
    text(posX,posY, S02_Lab                     ,flags_common{:});
    i_subj = 2;
    text(posX,posY     , S02_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});


    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'apa','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ata','Color','blue','Units','Normalized');
% set(gca,'YTickLabel',[]);
    ax=il_nexttile(13);colormap(ax,'hot'); xlabel('');ylabel(YLab); text(posX,posY,'aba'        ,'Color','white',flags_common{:});
    ax=il_nexttile(14);colormap(ax,'hot');xlabel('');ylabel('');   text(posX,posY,'apa'        ,'Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(15);                   
    xlabel('');ylabel('');   
    text(posX,posY, S01_Lab                     ,flags_common{:});
    i_subj = 1;
    i_exp  = 3; 
    text(posX,posY     , S01_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});
%
    ax=il_nexttile(16);
    xlabel('');ylabel('');
    text(posX,posY, S02_Lab                     ,flags_common{:});
    i_subj = 2;
    text(posX,posY     , S02_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'aba','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'apa','Color','blue','Units','Normalized');

    

    ax=il_nexttile(17);colormap(ax,'hot');           ylabel(YLab); text(posX,posY,'ada'        ,'Color','white',flags_common{:});
    ax=il_nexttile(18);colormap(ax,'hot');           ylabel('');   text(posX,posY,'ata'        ,'Color','white',flags_common{:});
    % set(gca,'YTickLabel',[]);
    ax=il_nexttile(19);
    ylabel('');   
    text(posX,posY, S01_Lab                     ,flags_common{:});
    i_subj = 1;
    i_exp  = 4; 
    text(posX,posY     , S01_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    ax=il_nexttile(20);
    ylabel('');   
    text(posX,posY, S02_Lab                     ,flags_common{:});
    i_subj = 2;
    text(posX,posY     , S02_Lab, flags_common{:});
%    text(posX,posY-0.15, sprintf('thres\n%.1f dB',thres(i_subj,i_exp)),'FontSize',FS_thres,flags_common{:});

    % set(gca,'YTickLabel',[]);
    text(x_cb,y_cb(2),'ada','Color','red' ,'Units','Normalized');
    text(x_cb,y_cb(1),'ata','Color','blue','Units','Normalized');

    nexttile(18);
    text(0.52,0.08,'f_0','Color','white','Units','Normalized');
    text(0.52,0.3,'F_1','Color','white','Units','Normalized');
    text(0.52,0.48,'F_2','Color','white','Units','Normalized');
    text(0.52,0.62,'F_3','Color','white','Units','Normalized');
    text(0.52,0.75,'F_4','Color','white','Units','Normalized');

    tAlignment = 'left';
    FS = 10;
    ax=nexttile(1);  title('ABDA22' ,'FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(5);  title('ADGA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(9);  title('APTA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(13); title('ABPA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);
    ax=nexttile(17); title('ADTA23','FontSize',FS,'FontWeight','bold'); set(ax,'TitleHorizontalAlignment',tAlignment);

    fout = [dir_out_figs 'Figure2'];
    % print(fout,'-dpdf')
    opts = [];
    opts.format = 'epsc'; % Pos34 = [650 800]; Pos = get(gcf,'Position'); Pos(3:4) = Pos34; set(gcf,'Position',Pos);
%    Saveas(gcf,fout,opts);

    opts = [];
    opts.format = 'fig';
%    Saveas(gcf,fout,opts);
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