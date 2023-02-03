function publ_osses2020a_JASA_EL_figs(varargin)
% function publ_osses2020a_JASA_EL_figs(varargin)
%
% 1. Description:
%       This function recreates the data processing and figures from the 
%       study by Osses, McLachlan, and Kohlrausch (2020, JASA-EL).
%       To run this script it is needed that you download the Zenodo dataset
%       located at: doi.org/10.5281/zenodo.7601160. You need to specify the
%       location of the downloaded directories into the variable dir_data
%       (see below).
%
% 2. Stand-alone example:
%       % Location in Alejandro's computer:
%       dir_data = '/home/alejandro/Documents/Databases/data/Osses-et-al-2020-JASA-EL/';
%       publ_osses2020a_JASA_EL_figs('table1','dir_data',dir_data); % instrument levels
%       publ_osses2020a_JASA_EL_figs('fig1'  ,'dir_data',dir_data); % Simulated reverberance
%       publ_osses2020a_JASA_EL_figs('fig2a' ,'dir_data',dir_data); % with no format
%       publ_osses2020a_JASA_EL_figs('fig2b' ,'dir_data',dir_data); % with no format
%       publ_osses2020a_JASA_EL_figs('table4','dir_data',dir_data); % one-way ANOVA
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2018
% Original name: r20180412_ch06_binaural_additional_plots.m, g2022_r20180412_proc_Glen.m
%
% Created on    : 12/04/2018
% Last update on: 03/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
end

h = [];
ha = [];
hname = [];

definput.import={'amt_cache'};
definput.flags.type={'fig1','fig2a','fig2b','table1','table4'};
definput.keyvals.dir_data = [];
[flags,keyvals]  = ltfatarghelper({},definput,varargin);
 
% clc, close all

%%%
if isempty(keyvals.dir_data)
    error('You need to specify the ''external folder'' where the raw data and waveforms are: keyvals.dir_data is a compulsory input');
end
dir_data = keyvals.dir_data;
dir_audio_sim = [dir_data '04-Stimuli-9s-for-simulations' filesep];

% dir_audio2use = dir_audio_exp;
dir_audio2use = dir_audio_sim; dBFS = 107.7;

%% Obtaining or loading the data:
Prev_model     = [];
Prev_model_max = [];
Prev_model_min = [];

fname = 'fig1_simulations_binaural_model';
c = amt_cache('get',fname,flags.cachemode);

if ~isempty(c)
    bRun = 0;

    Prev_model = c.Prev_model;
    Prev_model_max = c.Prev_model_max;
    Prev_model_min = c.Prev_model_min;
    Leq_A = c.Leq_A;
    Leq_A_max = c.Leq_A_max;
    Leq_Z = c.Leq_Z;
    Leq_Z_max = c.Leq_Z_max;

    instr_id  = c.instr_id;
    filter_id = c.filter_id;

    room_id = c.room_id;
    room_idx = c.room_idx;
else
    %%% Common parameters:
    instr_id  = {'b01_vio','b07_db','b08_fl','b09_pic','b16_contra_bsn','b17_horn','b20_tr','b22_timp'};
    filter_id = {'b01_vio','b07_db','b08_fl','b09_pic','b16_con'       ,'b17_hrn' ,'b20_tr','b22_tim'};
    N_instr = length(instr_id);
    %%%

    bRun = 1;
    subfs = 11025;
    % flags_model = {'subfs',subfs,'dboffset',dBFS,'publ_osses2017'}; % 'no_binaural','exclude'};
    flags_model = {'subfs',subfs,'dboffset',dBFS,'exclude'};

    %%% Using the closest configuration to that used by van Dorp:
    % flags_model = {'subfs',subfs,'dboffset',dBFS,'publ_vandorpschuitman2013'}; warning('temporal')
end  
                
if bRun
    for i = 1:N_instr
        dir_audio_cur = [dir_audio2use instr_id{i} filesep];
        files = Get_filenames(dir_audio_cur,[filter_id{i} '*.wav']);

        pRev = [];
        for j = 1:length(files)
            fprintf('\t Assessing Prev for %s\n',files{j});

            if i == 1
                idx = findstr(files{j},'_');
                label_room{j} = files{j}(idx(end)+1:end-4);
            end

            filename = [dir_audio_cur files{j}];

            [insig,fs] = audioread(filename);
            [output, info] = vandorpschuitman2013(insig,fs,flags_model{:});
            Prev_model(i,j) = median(output.prev_frame); % par.pRev
            Prev_model_min(i,j) = prctile(output.prev_frame,25);
            Prev_model_max(i,j) = prctile(output.prev_frame,75);

            % If you want to check... the results for the two first
            %   instruments should be:
            %  %     % A        Aabs      B          Babs     C         Cabs1     Cabs2     D          
            %  %   Prev_model = [ ...            
            %  %      11.1895    9.7612   10.7617    7.8050   15.9289   11.9120   12.9275   13.2082    % vio
            %  %      16.5423   16.4591   18.2064   16.4924   21.9031   17.6238   18.5414   20.1603    % db

            Leq_A(i,j)     = output.insig_level(1);
            Leq_A_max(i,j) = output.insig_level(2);
            Leq_Z(i,j)     = output.insig_level(3);
            Leq_Z_max(i,j) = output.insig_level(4);
            % numFrames(i,j) = length(outputs.par.FL_frame);
        end
    end

    instr_id  = {'b01_vio','b07_db','b08_fl','b09_pic','b16_contra_bsn','b17_horn','b20_tr','b22_timp'};
    filter_id = {'b01_vio','b07_db','b08_fl','b09_pic','b16_con'       ,'b17_hrn' ,'b20_tr','b22_tim'};

    i2sort = [2 4 1 6 3 7 8 5]; % maybe do a better job with the file naming
    room_id = label_room(i2sort);
    room_idx = 1:length(filter_id);

    % [Prev_model_stored,Prev_model_min_stored,Prev_model_max_stored] = il_get_PREV_9s_NEW;

    %%% Saving the cache:
    c = [];
    c.Prev_model = Prev_model(:,i2sort);
    c.Prev_model_max = Prev_model_max(:,i2sort);
    c.Prev_model_min = Prev_model_min(:,i2sort);
    c.Leq_A = Leq_A(:,i2sort);
    c.Leq_A_max = Leq_A_max(:,i2sort);
    c.Leq_Z = Leq_Z(:,i2sort);
    c.Leq_Z_max = Leq_Z_max(:,i2sort);

    c.instr_id  = instr_id;
    c.filter_id = filter_id;

    c.room_id = room_id;
    c.room_idx = room_idx;

    amt_cache('set',fname,c);
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Starting with the processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_table1
    LAFmax_90s = [77.1 73.8 88.7 74.2 64.7 78.5 84.4 78.2]'; % hard coded, from Osses et al. (2017)
    LZFmax_90s = [80.4 97.1 90.5 76.9 74.3 85.2 86.4 95.5]'; % hard coded, from Osses et al. (2017)

    var2show = [mean(Leq_A,2) mean(Leq_A_max,2) mean(Leq_Z,2) mean(Leq_Z_max,2) round(10*mean(Leq_Z_max,2))/10-LZFmax_90s round(10*mean(Leq_A_max,2))/10-LAFmax_90s];
    var2show(:,end+1) = [var2show(:,3)-var2show(:,1)];

    var2show
    % var2latex(round(var2show*10)/10)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1
    
    roomsLbls = room_id; % {'A_a_b_s','B_a_b_s','A','C_a_b_s_1','B','C_a_b_s_2','D','C'};
    LineStyle = {'kp'       ,'kv','r<','mo','b^','kd','b>','ks'};
    FaceColor = {rgb('Teal'),'y' ,'r' ,'w' ,'b' ,'k' ,rgb('CadetBlue') ,rgb('YellowGreen')};
    LineWidth = [ 1         , 1  , 2  , 2  , 2  , 2  , 1  , 1];
    % instrLbls = {'Vio','Viola','Cello','Cla','Bsn','FrHrn','Ti','Flute','Picc','Oboe','Trum','Tri','DBass','CBsn'}; 
    instrLbls = {'Vio','Viola','Cello','DBass','Flute','Picc','Oboe','Cla','Bsn','CBsn','FrHrn','Trum','Ti','Tri'}; 
    
    M = length(roomsLbls);
    
    FontSize = 14;
    
    %     T304corr(:,j) = transpose(mean( T30toti(instrGrps(1,j):instrGrps(2,j),:) ,1));
    %     EDT4corr(:,j) = transpose(mean( EDTtoti(instrGrps(1,j):instrGrps(2,j),:) ,1));
    
    Rev_ans     = Prev_model; % il_get_PREV_9s_NEW;
    % Rev_ans_10s = il_get_PREV_10s;
    Rev_ans_old = local_get_PREV_90s; % 8 x 14
    
    offsetx = 0.10;
    % x2plot = -2.5*offsetx:offsetx:2.5*offsetx;
    x2plot = -(M/2-0.5)*offsetx:offsetx:(M/2-0.5)*offsetx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx2plot = 1:8;
    idx2plot_old = [1 4 5 6 10 11 12 13]; %  Previous knowledge... % 1:14;
    % idx2plot = [1 2 3 9 11 13 5 7 12 4 10];
    % limits4trends = [6 9];
    
    %%%
    % Corr: between simulations:
    for j = 1:8
        Rev_old = Rev_ans_old(idx2plot_old(j),:)';
        
        [rp(j,1),rp(j,2)]   = corr(Rev_old,Rev_ans(idx2plot(j),:)','type','Pearson');
        [rs(j,1),rs(j,2)]   = corr(Rev_old,Rev_ans(idx2plot(j),:)','type','Spearman');
        
        % [rp_with10(j,1),rp_with10(j,2)]   = corr(Rev_ans_old(idx2plot_old(j),:)',Rev_ans_10s(idx2plot(j),:)','type','Pearson');
        % [rs_with10(j,1),rs_with10(j,2)]   = corr(Rev_ans_old(idx2plot_old(j),:)',Rev_ans_10s(idx2plot(j),:)','type','Spearman');
    end
    %%%
    
    N = length(idx2plot);
    Nold = length(idx2plot_old);
    ColorOld = rgb('LightGray');
    
    % maxEstimate = max( max(Rev_ans(:,idx2plot)) );
    
    step = 1;
    varX = step:step:N*step;
    
    figure;
    % Current data for 'Room 1':
    handle_room = plot(varX+x2plot(1),Rev_ans(idx2plot,1),LineStyle{1},'MarkerFaceColor',FaceColor{1},'LineWidth',LineWidth(1)); grid on; hold on
    %%%
    plot(varX+x2plot(1),Rev_ans_old(idx2plot_old,1),LineStyle{1}(2),'Color',ColorOld,'MarkerFaceColor',ColorOld,'LineWidth',LineWidth(1)); grid on
    for j = 2:M
        plot(varX+x2plot(j),Rev_ans_old(idx2plot_old,j),LineStyle{j}(2),'Color',ColorOld,'MarkerFaceColor',ColorOld,'LineWidth',LineWidth(j)); 
    end
    %%%
    
    for j = 2:M
        handle_room(end+1) = plot(varX+x2plot(j),Rev_ans(idx2plot,j),LineStyle{j},'MarkerFaceColor',FaceColor{j},'LineWidth',LineWidth(j)); 
    end
    
    % text(1*step+1.8*x2plot(1)-to_subtract(1), 2,sprintf('T30: \nEDT:'),'FontSize',FontSize);
    ylabel('Reverberance P_R_E_V [MU]','FontSize',FontSize);
    xlabel('Instrument','FontSize',FontSize);
    % Title(sprintf('Odeon orchestra for %.d-s long samples',dur));
    
    ha(end+1) = gca;
    h(end+1) = gcf;
    hname{end+1} = 'Odeon-Prev-sim4experiment';
    
    Pos = get(h(end),'Position');
    Pos(3) = 2.5*Pos(3);
    Pos(4) = 320;
    set(h(end),'Position',Pos);
    set(ha(end),'XTick',(1:N)*step);
    set(ha(end),'XTickLabel',instrLbls(idx2plot_old))
    set(ha(end),'YTick',2:2:26)
    
    linkaxes(ha,'xy');
    set(ha,'FontSize',FontSize)
    xval_min = min(varX)-step+0.2;
    xval_max = max(varX)+step-0.3; % 0.6
    xlim([xval_min xval_max]) 
    
    yval_min = 3;
    yval_max = 26.95;
    ylim([yval_min yval_max])
        
    set(gca,'XGrid','off')
    TickLength = get(gca,'TickLength');
    TickLength(1) = 0;
    set(gca,'TickLength',TickLength)
    
    for i = 1.5:N
        plot([i i]*step, [0 round(yval_max)*1.05],'Color','k')
    end
    
    leg2show = [];
    for i = 1:N
        leg2show{i} = ['Room ' num2str(room_idx(i))];
    end
    hleg = legend([handle_room],leg2show,'Location','SouthEast','FontSize',10);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig2a || flags.do_fig2b
    % Loading stored data:
    res = []; % fake initialisation to browse the variable more easily in MATLAB

    instr_id = {'b01_vio','b07_db','b08_fl','b09_pic','b16_contra_bsn','b17_horn','b20_tr','b22_timp'};
    room_id_not_sorted = {'A','Aabs','B','Babs','C','Cabs1','Cabs2','D'};
    idxEDT = [2 4 1 6 3 7 8 5];
    room_id  = room_id_not_sorted(idxEDT); %  {'Aabs','Babs','A','Cabs1','B','Cabs2','D','C'};

    file2load = [dir_data '03-Results-summary' filesep 'raw-data.mat'];
    load(file2load);
    %%%
                        
    h = [];
    hname = [];
    
    leg4plot = [];
    for j = 1:8
        leg4plot{end+1} = instr_id{j}(5:end);
    end
    
    idx = [1:8]; lab = 'Instr';
    
    if flags.do_fig2a
        bWithinRoom = 0;
        bWithinInstrument = 1; % between rooms
    end
    if flags.do_fig2b
        bWithinRoom = 1; % between instruments
        bWithinInstrument = 0;
    end
    if bWithinRoom == 1
        for j = 1:8
            figure;
            exp2eval=sprintf('h(end+1) = il_plot_res(res.sc_room%.0f(:,idx),res.sc_room%.0f_name{1});',j,j);
            eval(exp2eval);

            xlabel(lab)
            ylabel('Scores')
            ylim([-0.05 1.05])
            xlim([0.5 8.5])
            set(gca,'XTick',1:8);
            set(gca,'XTickLabel',leg4plot);
        end
    end
    
    if bWithinRoom == 0
        idx = [1:8]; lab = 'Room Nr.'; 
        room_lab = room_id; 

        figure; 
        [h(end+1),Me,errL,errU] = local_plot_res(res.sc_instr1(:,idx),res.sc_instr1_name{1});
        xlabel(lab)
        ylabel('Scores')
        ylim([-0.05 1.05])
        xlim([0.5 8.5])
        set(gca,'XTick',1:8);
        set(gca,'XTickLabel',room_lab);
        hname{end+1} = ['PREVexp-instr-' num2str(1)];

        for j = 2:8
            figure;
            exp2eval=sprintf('[h(end+1),Me(end+1,:),errL(end+1,:),errU(end+1,:)] = il_plot_res(res.sc_instr%.0f(:,idx),res.sc_instr%.0f_name{1});',j,j);
            eval(exp2eval);
            xlabel(lab)
            ylabel('Scores')
            ylim([-0.05 1.05])
            xlim([0.5 8.5])
            set(gca,'XTick',1:8);
            set(gca,'XTickLabel',room_lab);
            hname{end+1} = ['PREVexp-instr-' num2str(1)];
        end
    end
    
    disp('')
    % bCorr = 1;
    % if bCorr
    % 
    %     rp = zeros(8,1);
    %     pp = zeros(8,1);
    %     rs = zeros(8,1);
    %     ps = zeros(8,1);
    %     rp_max = zeros(8,1);
    %     pp_max = zeros(8,1);
    %     rs_max = zeros(8,1);
    %     ps_max = zeros(8,1);
    % 
    %     Prev_max_max = max(max(Prev_model));
    % 
    %     Prev_model     = Prev_model(:,idxEDT);
    %     Prev_model_max = Prev_model_max(:,idxEDT);
    %     Me = Me(:,idxEDT);
    %     room_lab = room_lab(idxEDT);
    % 
    %     bUseMax = input('0 = Use avg. 90 s; 1 = use max 90s: ');
    %     switch bUseMax
    %         case 0
    %             Prev_model_90 = il_get_PREV_90s;
    %         case 1
    %             [xx,Prev_model_90] = il_get_PREV_90s;  
    %     end
    %     idxs2use = [1 4 5 6 10 11 12 13];
    %     Prev_model_90 = transpose(Prev_model_90(idxs2use,:)); % now each row is one of the used instruments
    % 
    %     for j = 1:8
    %         xdata = Prev_model(j,:)';
    %         ydata = Me(j,:)';
    %         [rp(j),pp(j)] = corr(xdata,ydata,'Type','Pearson');
    %         [rs(j),ps(j)] = corr(xdata,ydata,'Type','Spearman');
    % 
    %         [rp_max(j),pp_max(j)] = corr(Prev_model_max(j,:)',ydata,'Type','Pearson');
    %         [rs_max(j),ps_max(j)] = corr(Prev_model_max(j,:)',ydata,'Type','Spearman');
    % 
    %         %%%
    %         % 90-s long data aseem not to be so representative
    %         xdata_90 = Prev_model_90(j,:)';
    % 
    %         [rp2(j),pp2(j)] = corr(xdata_90,ydata,'Type','Pearson');
    %         [rs2(j),ps2(j)] = corr(xdata_90,ydata,'Type','Spearman');
    % 
    %         % model with model
    %         [rpm(j),ppm(j)] = corr(xdata,xdata_90,'Type','Pearson');
    %         [rsm(j),psm(j)] = corr(xdata,xdata_90,'Type','Spearman');
    %         %%%
    % 
    %         % Prev_model_no = il_norm_model(xdata,min(ydata),max(ydata));
    %         ydata_no = il_norm_model(ydata,min(xdata),max(xdata));
    % 
    %         errL = xdata-Prev_model_min(j,:)';
    %         errU = Prev_model_max(j,:)'-xdata;
    % 
    %         errL = max(errL,0);
    % 
    %         figure;
    %         subplot(1,3,1)
    %         errorbar(1:8, xdata,errL,errU,'ro--','LineWidth',2,'MarkerFaceColor','r'); hold on
    %         plot(1:8, ydata_no ,'ks-');
    %         xlim([0.5 8.5]);
    %         % set(gca,'YTick',.1:.1:.9);
    %         set(gca,'XTick',1:8);
    %         set(gca,'XTickLabel',room_lab);
    %         set(gca,'FontSize',14);
    %         legend('P_R_E_V_,_s_i_m','P_R_E_V_,_m_e_d (scaled)','Location','NorthWest')
    %         Title(name2figname(instr_id{j}));
    %         y1 = floor(min(xdata-errL))-1;
    %         y2 = ceil(max(xdata+errU))+1;
    %         ylim([y1 y2])
    %         set(gca,'YTick',y1+1:y2-1);
    %         Ylabel('Reverberance P_R_E_V [MU]')
    %         Xlabel('Room')
    %         % text('')
    % 
    %         opts = [];
    %         opts.bLegend = 0;
    % 
    %         subplot(1,3,2)
    %         corr_plot(xdata,ydata,'Pearson',opts);
    %         Xlabel('P_R_E_V_,_s_i_m [MU]');
    %         Ylabel('P_R_E_V_,_e_x_p');
    %         ylim([-0.1 1.1])
    %         xlim([floor(min(xdata))-1 ceil(max(xdata))+1])
    %         [r,p] = corr(xdata,ydata,'type','Pearson');
    %         text2print = sprintf('r_p=%.2f\n p_v_a_l=%.3f',r,p);
    %         % xlim_this = get(gca,'XLim');
    %         text(min(xdata),.95,text2print,'FontSize',14)
    % 
    %         subplot(1,3,3)
    %         corr_plot(xdata,ydata,'Spearman', opts);
    %         Xlabel('P_R_E_V_,_s_i_m (ordinal)');
    %         Ylabel('P_R_E_V_,_e_x_p (ordinal)');
    %         set(gca,'XTick',1:8);
    %         xlim([0.5 8.5]);
    %         ylim([0.5 8.5]);
    % 
    %         [r,p] = corr(xdata,ydata,'type','Spearman');
    %         text2print = sprintf('r_s=%.2f\n p_v_a_l=%.3f',r,p);
    %         % xlim_this = get(gca,'XLim');
    %         text(1,7.5,text2print,'FontSize',14)
    % 
    %         set(gcf,'Position',[40 250 1300 420]);
    %         h(end+1) = gcf;
    %         hname{end+1} = ['PREV-within-instr-' name2figname(instr_id{j})];
    % 
    %     end
    % 
    %     % instr_id = {'b01_vio','b07_db','b08_fl','b09_pic','b16_contra_bsn','b17_horn','b20_tr','b22_timp'};
    %     idxTrend1 = [1 6 8]; % 1 = vio / 6 = FrHrn / Ti
    %     idxTrend2 = [3 4 7]; % 3 = fl / 4 = pic / 7 = tru
    %     idxTrend3 = [2 5]; % 2 db / 5 = Cbsn
    %     Prev1 = Prev_model(idxTrend1,:)'; Me1 = Me(idxTrend1,:)';
    %     Prev2 = Prev_model(idxTrend2,:)'; Me2 = Me(idxTrend2,:)';
    %     Prev3 = Prev_model(idxTrend3,:)'; Me3 = Me(idxTrend3,:)';
    %     [rp_grand,pp_grand] = corr(Prev_model(:),Me(:),'type','Pearson');
    %     [rs_grand,ps_grand] = corr(Prev_model(:),Me(:),'type','Spearman');
    %     [rp_trend(1),pp_trend(1)] = corr(Prev1(:),Me1(:),'type','Pearson');
    %     [rs_trend(1),ps_trend(1)] = corr(Prev1(:),Me1(:),'type','Spearman');
    %     [rp_trend(2),pp_trend(2)] = corr(Prev2(:),Me2(:),'type','Pearson');
    %     [rs_trend(2),ps_trend(2)] = corr(Prev2(:),Me2(:),'type','Spearman');
    %     [rp_trend(3),pp_trend(3)] = corr(Prev3(:),Me3(:),'type','Pearson');
    %     [rs_trend(3),ps_trend(3)] = corr(Prev3(:),Me3(:),'type','Spearman');
    %     % [rp_trend(4),pp_trend(4)] = corr(Prev3(:,2),Me3(:,2),'type','Pearson');
    %     % [rs_trend(4),ps_trend(4)] = corr(Prev3(:,2),Me3(:,2),'type','Spearman');
    % 
    %     figure;
    %     subplot(1,2,1)
    %     corr_plot(Prev_model(:),Me(:),'Pearson'); grid on
    %     xlabel('P_R_E_V (model)')
    %     ylabel('P_R_E_V (exp)')
    %     ylim([0.1 1.1])
    % 
    %     subplot(1,2,2)
    %     corr_plot(Prev_model(:),Me(:),'Spearman'); grid on
    %     xlabel('P_R_E_V ordinal (model)')
    %     ylabel('P_R_E_V ordinal (exp)')
    %     % ylim([0.1 1.1])
    % 
    % 
    %     disp('The most important result: Prev,exp with Prev,sim (within instrument, 10-s sounds): ')
    %     var2latex(round(100*[rp pp rs ps])/100)
    % 
    %     if bUseMax
    %         disp('Prev,exp with Prev,sim-MAX (but sim with 90-s sounds): ')
    %     else
    %         disp('Prev,exp with Prev,sim-AVG (but sim with 90-s sounds): ')
    %     end
    %     var2latex(round(100*[rp2' pp2' rs2' ps2'])/100)
    % 
    %     % var2latex(round(100*[rp_max pp_max rs_max ps_max])/100)
    %     disp('')
    % 
    % end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_table4
    var = load([dir_data '03-Results-summary' filesep 'raw-data.mat']);
    
    sc = median(var.res.sc_room1);
    
    instr_id = [];
    room_id = [];
    
    instr_id{end+1} = var.res.sc_instr1_name{1};
    room_id{end+1}  = var.res.sc_room1_name{1};
    for i = 2:8
        exp2eval = sprintf('sc = [sc; median(var.res.sc_room%.0f)];',i);
        eval(exp2eval);
        
        exp2eval = sprintf('room_id{end+1}=var.res.sc_room%.0f_name{1};',i);
        eval(exp2eval); 
        
        exp2eval = sprintf('instr_id{end+1}=var.res.sc_instr%.0f_name{1};',i);
        eval(exp2eval); 
    end
    
    % room_id_not_sorted = {'A','Aabs','B','Babs','C','Cabs1','Cabs2','D'};
    % i2sort = [2 4 1 6 3 7 8 5]; 
    % sc = sc(i2sort,:);
    % room_id = room_id(i2sort);
    
    %%%
    % p = anova1(sc);
    
    sc
    
    % figure;
    for i = 1:8
        % % Two-way ANOVA
        % exp2eval = sprintf('[p,out1,out2] = anova2(var.res.sc_room%.0f); title(room_id{%.0f});',i,i);
        % One-way ANOVA
        exp2eval = sprintf('[p,out1,out2] = anova1(var.res.sc_room%.0f); close; title(room_id{%.0f});',i,i);
        eval(exp2eval);
        
        disp('')
    end
    %%%
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rev_ans,Rev_max] = local_get_PREV_90s
% Run: m20180125_working4paper_JASAEL, stop at L456

% instrGrps =
%      1     5     6     7     8     9    10    12    14    16    17    20    22    23
%      4     5     6     7     8     9    11    13    15    16    19    21    22    23
%      Vio   Viola Cello DBass Flute Picc Oboe  Cla   Bsn   CBsn  FrHrn Trum  Ti    Tri

Rev_max = [... %Val_max =
   15.3123   15.5387   16.4999   17.7971   14.7098   10.3776   18.9353   10.8187   13.3986   14.1480   11.0260    9.4866   20.6316    7.3789
   14.7606   14.5307   16.7058   16.6621   11.0512    6.3250   11.4375    8.3918   12.8943   14.9747   10.3996    9.9741   21.7808    4.8116
   15.1953   15.8669   16.5716   17.3327   12.1565    9.9984   16.7652    9.9594   13.7856   12.4320   12.0357   10.4269   23.3926    5.9422
   15.0798   16.0993   16.8114   18.3957   11.2913    7.9011   14.3171   10.0024   14.0498   13.8579   10.4140   10.3218   22.7612    5.5451
   16.2496   16.3392   17.8086   19.3745   11.5448    8.3544   14.8733   11.8647   14.8429   13.5560   14.1621   10.8467   23.3549    6.2536
   16.2041   16.8635   16.9949   19.1675   12.2131    9.3332   15.7600   11.0598   15.8956   14.1200   13.2287   11.2837   24.2728    6.2586
   15.9700   17.4143   18.6044   19.7364   10.6692    7.6146   15.4062   11.9550   16.5248   16.4242   13.6768   10.8033   25.1084    2.5385
   19.5533   18.7086   20.4532   22.3169   16.4522   12.4841   19.4874   14.8749   18.8637   17.0142   15.9769   13.5794   26.9083    7.7374];
Rev_max = Rev_max';

Rev_ans = [...
   12.2633   13.3466   13.6059   15.0666   11.1275    8.6675   11.8549    6.9005   10.2818    6.7791    7.0471    6.4068   18.2858    7.0157
   11.2291   12.1128   13.2556   12.2480    8.5834    5.3881    7.8385    5.7077    9.5469   10.1649    5.9754    5.2814   19.6861    4.5113
   12.0586   13.4499   14.3517   14.4881    9.6291    7.7944   11.0116    7.1857   10.1041    5.1655    7.4187    5.8944   20.8324    5.6011
   11.9042   13.4929   14.9818   15.0623    7.8576    6.6162    9.6531    7.0827   10.8633    4.9754    6.1630    5.4206   19.7431    4.7448
   13.0040   14.2512   15.8884   14.6228    8.5442    7.0653    9.9658    7.9807   11.1994    5.9801    7.4101    5.0831   20.8737    5.3358
   13.1479   14.4507   15.9189   16.1827    8.9512    7.3701   11.1287    7.9857   11.9317    5.2550    6.9983    6.4100   20.8191    5.4219
   12.9816   15.6077   16.3784   14.3101    7.2913    5.9567   10.1072    7.8547   12.1946   12.4093    6.9099    5.0751   21.5170    1.2823
   16.4841   17.2998   18.5668   20.7709   11.8134   10.5564   13.3692   10.8360   14.6186    9.3039    9.8252    7.4184   24.8661    7.1700];
Rev_ans = Rev_ans';
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,Me,errL,errU] = local_plot_res(scores,label4plot);

Me = median(scores);
errL = Me-prctile(scores,25);
errU = prctile(scores,75)-Me;

N = size(scores,1);
M = size(scores,2);

errorbar(1:M,Me,errL,errU,'s-'); 
title(sprintf('%s (N=%.0f)',label4plot,N),'interpreter','none');

h = gcf;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function xdata = il_norm_model(xdata,Min_val,Max_val)
% 
% Min_x = min(xdata);
% Max_x = max(xdata);
% xdata = (xdata-Min_x)/(Max_x-Min_x); % now scaled between 0 and 1
% xdata = xdata*(Max_val-Min_val)+Min_val;