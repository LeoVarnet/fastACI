function [h,hname] = publ_varnet2022a_JASA_figs(varargin)
% function [h,hname] = publ_varnet2022a_JASA_figs(varargin)
%
% Generates the figures from the publication by Varnet and Lorenzi (2022,
%   JASA). This script requires some preprocessed MAT files.
%
% % To display Fig. 2 of Varnet and Lorenzi (2022, JASA) use :::
%     publ_varnet2022a_JASA_figs('fig2');
%
% % To display Fig. 3 of Varnet and Lorenzi (2022, JASA) use :::
%     publ_varnet2022a_JASA_figs('fig3');
%
% % To display Fig. 4 of Varnet and Lorenzi (2022, JASA) use :::
%     publ_varnet2022a_JASA_figs('fig4');
%
% % To display Fig. 1 of the supplementary materials of Varnet and Lorenzi 
% %     (2022, JASA) use :::
%     publ_varnet2022a_JASA_figs('fig1_suppl');
%
% % To display Fig. 2 of the supplementary materials of Varnet and Lorenzi 
% %     (2022, JASA) use :::
%     publ_varnet2022a_JASA_figs('fig2_suppl');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_varnet2022a_JASA_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2','fig3','fig3a','fig3c','fig3d','fig4','fig1_suppl','fig2_suppl'};
definput.keyvals.models=[];
definput.flags.publ = {'varnet2022a_JASA','do_varnet2022b_CFA'}; % the first one is the default
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%%% Checking whether the data are stored locally:
bAre_stored_locally = publ_varnet2022a_JASA_0_checkdata;
if ~isempty(find(bAre_stored_locally==0))
    fprintf('%s: Please make sure you donwload the experimental data from Zenodo\n',upper(mfilename));
    fprintf('\tas instructed in the script publ_varnet2022a_JASA_0_checkdata.m\n');
end
%%%

data = il_load_data(flags,keyvals);
N_subjects = data.N_subjects;

plot_models = 0;

% if N_subjects == 8
    colorcode = {   [0 0.447 0.741]; ...     % Colour for S1
                    [0.85 0.325 0.098]; ...  %            S2
                    [0.929 0.694 0.125]; ... %            S3
                    [0.494 0.184 0.556]; ... %            S4
                    [0.466 0.674 0.188]; ... %            S5
                    [0.301 0.745 0.933]; ... %            S6
                    [0.635 0.078 0.184]; ... %            S7
                    [1 0.6 1]};              %            S8
    legendnames = 'S' + string(1:8); % Make sure that the labels are S1 to S8
    
    CI_FaceAlpha = 0.15;
    Data_LineWidth = 2;%1.2;
    IT_LineWidth = 2; 
    IT_LineStyle = ':';
    IT_Color = 'b'; % Ideal template colour
% else
%     error('The number of participants should be 8, please check why MATLAB cannot find all the raw data...')
% end

%% do_fig2: H, CR rates ---------------------------------------------------
if flags.do_fig2
    % if plot_models == 1
    %     error('Continue here...')
    %     Pos = [100 100 700 400];
    % else
    Pos = [100 100 800 400];
    % end
    figure('Position', Pos);

    hist_H = data.hist_H;
    hist_M = data.hist_M;
    hist_CR = data.hist_CR;
    hist_FA = data.hist_FA;
    hist_N = data.hist_N;
    m_levels = data.m_levels;
    
    TPmatrix = hist_H./(hist_H+hist_M); 
    TPmatrix(hist_N<150)=nan;
    TAmatrix = hist_CR./(hist_CR+hist_FA); 
    TAmatrix(hist_N<150)=nan;
    Nmatrix = hist_N; %Nmatrix(hist_N<150)=nan;

    t=tiledlayout('flow','TileSpacing','compact');
    for i_subject = 1:N_subjects

        nexttile(i_subject); % subplot(2,3,i_subject)

        yyaxis left
        plot(m_levels*2,TPmatrix(i_subject,:),'-o','LineWidth',1.5,'Color',[colorcode{i_subject}]); hold on;
        plot(m_levels*2,TAmatrix(i_subject,:),':x','LineWidth',1.5,'Color',[colorcode{i_subject}]);
        ylim([0 1])
        if plot_models == 1
            xlim([-40 -1])
        else
            xlim([-14 -3]*2)
        end

        xlabel(t,'m (dB)')
        ylabel(t,'correct response rate')
        
        yyaxis right
        area(m_levels*2,Nmatrix(i_subject,:),'FaceColor',[colorcode{i_subject}], ...
            'FaceAlpha',0.3*(rgb2lightness(reshape(colorcode{i_subject},[1,1,3])))/100,'EdgeColor','none');
        ylim([0 1000])

        grid on
        if (plot_models==0 && mod(i_subject,4)==0) || (plot_models==1 && mod(i_subject,3)==0)
            ylabel('N of trials')
        else
            set(gca,'YTickLabels',[]);
        end
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = [0.5 0.5 0.5];

        title(legendnames(i_subject),'interpreter','none')

    end
    legend({'H','CR'},'Location','SouthEast')
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2';
end

%% do_fig3: plot temporal revcorr - group plot ----------------------------
if flags.do_fig3 
    
    tE = data.tE;
    ideal_template = data.ideal_template;
    CIt_h = data.CIt_h;
    CIt_l = data.CIt_l;
    CIt_all = data.CIt_all;
    CI1t_all = data.CI1t_all;
    CI1t_l = data.CI1t_l;
    CI1t_h = data.CI1t_h;
    CI2t_l = data.CI2t_l;
    CI2t_h = data.CI2t_h;
    CI2t_all = data.CI2t_all;
    
    %figure('Name', 'Temporal kernel','Position', [100 100 600 450]); 
    tl = tiledlayout(3,1,'TileSpacing','compact');%subplot(3,1,1)

    if plot_models == 1
        Fideal = 1.25;
    else
        Fideal = 1e-4;
    end

    nexttile(1)
    hold on
    if plot_models == 1
        plot(tE, 1.5*Fideal*ideal_template'/max(max(ideal_template)), IT_Color, 'linewidth', IT_LineWidth, 'LineStyle', IT_LineStyle)%, @plot, 1, size(undersmplE,1))
        % plot(tE, mean(CIt_real)/std(mean(CIt_real)),'Color','k', 'linewidth', Data_LineWidth);
        % plot(tE, (CIt_real)/std(mean(CIt_real)),'Color',[0 0 0 0.2],'LineWidth',0.5);    
    else
        plot(tE, 1.5*Fideal*ideal_template'/max(max(ideal_template)), IT_Color, 'linewidth', IT_LineWidth, 'LineStyle', IT_LineStyle)%, @plot, 1, size(undersmplE,1))
        fill([tE fliplr(tE)], [mean(CIt_h) fliplr(mean(CIt_l))]',':k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    end

    for i_subject = 1:N_subjects
        if plot_models == 1
           CIt_all(i_subject,:) = CIt_all(i_subject,:)/std(CIt_all(i_subject,:));
           plot(tE, CIt_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth)
        else
            plot(tE, CIt_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth)
        end
    end
    % plot(tE, mean(CIt_all)','k-', 'linewidth', 2)
    %fill([tE fliplr(tE)], [mean(CIt_all)+0.5*std(CIt_all) fliplr(mean(CIt_all)-0.5*std(CIt_all))]', 'k', 'FaceAlpha',0.1,'LineStyle','none');
    plot(tE, zeros(size(tE)), 'k', 'linewidth', 0.5)%, @plot, 1, size(undersmplE,1))

    if plot_models == 1
        ylabel('normalized weights');
    else
        ylabel('weights');
    end

    ylim([-3.5 3.5]*Fideal); xlim([tE(1) tE(end)])
    title('general kernel'); %title('general kernel');
    xlabel('time (s)');
    %grid on
    %box on
    %set(gca,'XTickLabels',[])
    %%
    nexttile(2)%subplot(3,1,2)
    hold on
    if plot_models == 1
        % plot(tE, mean(CI2t_real)/std(mean(CI2t_real)),'Color','k', 'linewidth', Data_LineWidth);
        % plot(tE, (CI2t_real)/std(mean(CI2t_real)),'Color',[0 0 0 0.2],'LineWidth',0.5);    
    else
        plot(tE, 1.5*Fideal*ideal_template'/max(max(ideal_template)), IT_Color, 'linewidth', IT_LineWidth)%, @plot, 1, size(undersmplE,1))
        fill([tE fliplr(tE)], [mean(CI2t_h) fliplr(mean(CI2t_l))]',':k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    end

    for i_subject = 1:N_subjects
        if plot_models == 1
            %CI2t_all(i_subject,:) = CI2t_all(i_subject,:)/std(CI2t_all(i_subject,:));
        end
        plot(tE, CI2t_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth)
    end
    %plot(tE, mean(CI2t_all)','k-', 'linewidth', 2)
    if plot_models == 1
        ylabel('normalized weights');
    else
        ylabel('weights');
    end
    title('target-present kernel');
    ylim([-3.5 3.5]*Fideal); xlim([tE(1) tE(end)])
    set(gca,'XTickLabels',[])
    grid on
    box on

    %xlabel('time (s)');
    h_here = [];
    nexttile(3)%subplot(3,1,3)%subplot(3,1,2)figure
    hold on
    if plot_models == 1
        % h(end+1) = plot(tE, mean(CI1t_real)/std(mean(CI1t_real)),'Color','k', 'linewidth', Data_LineWidth);
        % plot(tE, (CI1t_real)/std(mean(CI1t_real)),'Color',[0 0 0 0.2],'LineWidth',0.5);    
    else
        h_here(end+1) = plot(tE, 1.5*Fideal*ideal_template'/max(max(ideal_template)), IT_Color, 'linewidth', IT_LineWidth);%, @plot, 1, size(undersmplE,1))
        fill([tE fliplr(tE)], [mean(CI1t_h) fliplr(mean(CI1t_l))]',':k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    end

    for i_subject = 1:N_subjects
        if plot_models == 1
            %CI1t_all(i_subject,:) = CI1t_all(i_subject,:)/std(CI1t_all(i_subject,:));
        end
        h_here(end+1) = plot(tE, CI1t_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth);
    end
    %plot(tE, mean(CI1t_all)','k-', 'linewidth', 2)
    if plot_models == 1
        ylabel('normalized weights');
    else
        ylabel('weights');
    end
    title('target-absent kernel');
    ylim([-3.5 3.5]*Fideal);xlim([tE(1) tE(end)])
    xlabel('time (s)');
    grid on
    box on

    %nexttile([3 1]);
    if plot_models == 1
        lgd=legend(h_here, ['S average', legendnames], 'Location', 'eastoutside', 'Interpreter', 'none');
    else
        lgd=legend(h_here, ['IT',legendnames], 'Location', 'eastoutside', 'Interpreter', 'none');
    end
    lgd.Layout.Tile = 'east';
    lgd.Layout.TileSpan = [3 4];

    h(end+1) = gcf;
    hname{end+1} = 'fig3';
end

%% do_fig4: plot spectral revcorr -----------------------------------------
if flags.do_fig4
    
    fE = data.fE;
    CIf_all = data.CIf_all;
    CIf_l = data.CIf_l;
    CIf_h = data.CIf_h;
    
    CI2f_all = data.CI2f_all;
    CI2f_l = data.CI2f_l;
    CI2f_h = data.CI2f_h;
        
    CI1f_all = data.CI1f_all;
    CI1f_l   = data.CI1f_l;
    CI1f_h   = data.CI1f_all;
    figure('Name', 'Fourier kernel','Position', [100 100 600 450]); 
    tl = tiledlayout(3,1,'TileSpacing','compact'); % subplot(3,1,1)

    if plot_models == 1
        arrowheight = 0.01;
    else
        arrowheight = 0.002;
    end

    nexttile(1)
    hold on
    if plot_models == 1
        % plot(fE, mean(4*CIf_real),'Color','k')
        % fill([fE fliplr(fE)], [mean(4*CIf_real)+std(4*CIf_real) fliplr(mean(4*CIf_real)-std(4*CIf_real))]', 'k', 'FaceAlpha',0.1,'LineStyle','none');
    end
    for i_subject = 1:N_subjects
        plot(fE, CIf_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth)
    end
    %plot(fE, max(CIf_all(:))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))), 'r', 'linewidth', 2)%, @plot, 1, size(undersmplE,1))
    %plot(fE, mean(CIf_all)','k-', 'linewidth', 2)
    fill([fE fliplr(fE)], [mean(CIf_h) fliplr(mean(CIf_l))]', 'k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    plot([4 4],[0 arrowheight],'-^','Color',IT_Color, 'linewidth', IT_LineWidth, 'MarkerIndices', [2], 'MarkerSize', 4)
    set(gca, 'XTickLabels', [])
    ylabel('weight amplitude'); title('general kernel'); 
    if plot_models == 1
        ylim([0 0.07])
    else
        ylim([0 0.027])
    end
    grid on
    box on

    nexttile(2)
    hold on
    if plot_models == 1
        % plot(fE, mean(4*CI2f_real),'Color','k')
        % fill([fE fliplr(fE)], [mean(4*CI2f_real)+std(4*CI2f_real) fliplr(mean(4*CI2f_real)-std(4*CI2f_real))]', 'k', 'FaceAlpha',0.1,'LineStyle','none');
    end

    for i_subject = 1:N_subjects
        plot(fE, CI2f_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth)
    end
    %plot(fE, 0.02*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))), 'r', 'linewidth', 2)%, @plot, 1, size(undersmplE,1))
    %plot(fE, mean(CI2f_all)','k-', 'linewidth', 2)
    fill([fE fliplr(fE)], [mean(CI2f_h) fliplr(mean(CI2f_l))]', 'k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    plot([4 4],[0 arrowheight],'-^','Color',IT_Color, 'linewidth', IT_LineWidth, 'MarkerIndices', [2], 'MarkerSize', 4)
    ylabel('weight amplitude');
    title('target-present kernel'); 
    if plot_models == 1
        ylim([0 0.07])
    else
        ylim([0 0.027])
    end
    %xlabel('frequency (Hz)'); 
    set(gca, 'XTickLabels', [])
    %legend({'ideal observer' D.name}, 'Location', 'eastoutside', 'Interpreter', 'none');
    grid on
    box on

    nexttile(3)
    hold on
    h_here = [];
    if plot_models == 1
    %     h(end+1) = plot(fE, mean(4*CI1f_real),'Color','k')
    %     fill([fE fliplr(fE)], [mean(4*CI1f_real)+std(4*CI1f_real) fliplr(mean(4*CI1f_real)-std(4*CI1f_real))]', 'k', 'FaceAlpha',0.1,'LineStyle','none');
    end
    for i_subject = 1:N_subjects
        h_here(end+1) = plot(fE, CI1f_all(i_subject,:)','Color',[colorcode{i_subject}],'Linewidth',Data_LineWidth);
    end
    %plot(fE, 0.02*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))), 'r', 'linewidth', 2)%, @plot, 1, size(undersmplE,1))
    %plot(fE, mean(CI1f_all)','k-', 'linewidth', 2)
    fill([fE fliplr(fE)], [mean(CI1f_h) fliplr(mean(CI1f_l))]', 'k', 'FaceAlpha',CI_FaceAlpha,'LineStyle','none');
    plot([4 4],[0 arrowheight],'-^','Color',IT_Color, 'linewidth', IT_LineWidth, 'MarkerIndices', [2], 'MarkerSize', 4)
    ylabel('weight amplitude');
    title('target-absent kernel'); 
    if plot_models == 1
        ylim([0 0.07])
    else
        ylim([0 0.027])
    end
    xlabel('frequency (Hz)'); 
    grid on
    box on
    % 
    % if isyes(plot_models)
    %     legend(h,['S average', legendnames], 'Location', 'eastoutside', 'Interpreter', 'none');
    % else
    lgd = legend(h_here,[legendnames], 'Location', 'eastoutside', 'Interpreter', 'none');
    
    lgd.Layout.Tile = 'east';
    lgd.Layout.TileSpan = [3 4];
    
    h(end+1) = gcf;
    hname{end+1} = 'fig4';
end

%% fig1_suppl: Audiograms -------------------------------------------------
if flags.do_fig1_suppl
    %%% Loading the data:
    f = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
    Audiog = [...
         10,  5, 10,  5, 10, 15, 20, 20, 35, 35, 10,  5, 10, 10, 15, 20, 20, 20, 25, 35; ... % S1
         15, 10,  5,  5, 10,  5,  5, 20, 25, 10,  5, 10, 10, 10, 10, -5, 10, 10, 25, 15; ...
         10, 10, 10,  5, 10,  0,  0,  5, 25, 15, 10,  5, 10, 10,  5,  5,  5, 10, 20,  5; ...
         -5,  0,  0,  0,  5,  0, -5, 15, 30,  0,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan; ... % S4
          5,  0, -5,  0,  0,  0,  5, 20, 10, 0,  5,  0, -5,  0,  0,  0,  5, 20, 15,  5; ...
         15, 10, 10, 10, 10,  5,  5,  5, 25,  5, 10,  5,  5,  0,  0,  5, -5,  5, 15, -5; ...
         10, 10,  5,  5, 10, 15, 10, 15, 25, 10,  5, 10, 10,  5,  5, 10, 10, 15, 20, 15; ...
          5,  0,  5,  5,  5,  5, -5,-10,  5,  5,  5,  0,  0,  0,  5,  5,  5, -5, 15, 10; ...
          10, 10, 10, 10, 10,  5, 15,  5, 25, 15, 10,  5, 10, 10,  5, 10,  5,  0, 20, 35]; % Audiogram of the excluded participant
    Audiog = Audiog(1:N_subjects,:);
    
    figure('Position', [100 100 750 300]); 

    for i_subject = 1:N_subjects
        subplot(1,2,1);
        semilogx(f,Audiog(i_subject,1:10),'Color',[colorcode{i_subject} 0.7]);hold on;
        subplot(1,2,2); 
        semilogx(f,Audiog(i_subject,11:end),'Color',[colorcode{i_subject} 0.7]);hold on; 
    end
    subplot(1,2,1);
    semilogx(f,mean(Audiog(:,1:10)),'k','LineWidth',2);
    set(gca,'Ydir','reverse');
    title('Left ear');
    xlabel('frequency (Hz)');
    xlim([f(1) f(end)]);
    ylabel('dB HL');
    
    YL = [-15 65];
    %%% Drawing a grey rectangle:
    basef = 1000; % Hz
    coor_x = audtofreq(freqtoaud(basef)+[-.5 +.5]); % +/-.5 ERB_N from 1000 Hz
    coor_x = [coor_x coor_x([2 1])];
    coor_y = [YL(1) YL(1) YL(2) YL(2)];
    colour_here = [0.5 0.5 0.5]; % grey
    
    fill(coor_x,coor_y,colour_here,'EdgeColor','none','FaceAlpha',0.2)
    ha = gca; % handle axis
    
    subplot(1,2,2); 
    semilogx(f,nanmean(Audiog(:,11:end)),'k','LineWidth',2);
    set(gca,'Ydir','reverse');
    title('Right ear');
    xlabel('frequency (Hz)');
    xlim([f(1) f(end)]);
    ylabel('dB HL');
    fill(coor_x,coor_y,colour_here,'EdgeColor','none','FaceAlpha',0.2)
    
    ha(end+1) = gca; % handle axis
    
    linkaxes(ha,'xy');
    ylim(YL);
    XT = [500,1000,2000,5000];
    set(ha, 'XTick', XT);
    set(ha,'XTickLabels',XT);
    
    legend(legendnames, 'Location', 'southwest', 'Interpreter', 'none'); 
    
    h(end+1) = gcf;
    hname{end+1} = 'Audiograms.png';
end

%% do_fig2_suppl: Behavior ------------------------------------------------
if flags.do_fig2_suppl
    
    m = data.m;
    trialnum = data.trialnum;
    bias = data.bias;
    m_reject = data.m_reject;
    bias_reject = data.bias_reject;
    
    figure('Position', [100 100 550 350]); 

    tiledlayout('flow','TileSpacing','compact');
    for i_subject = 1:N_subjects
        nexttile(1) % subplot(2,1,1);
        plot(trialnum, m(i_subject,:)*2,'Color',[colorcode{i_subject}]); 
        hold on;
        nexttile(2) % subplot(2,1,2);
        plot(trialnum, bias(i_subject,:)-1,'Color',[colorcode{i_subject}]); 
        hold on;
    end

    nexttile(1) % subplot(2,1,1); 
    plot(trialnum, m_reject*2,'k:'); 
    ylabel('m (dB)'); xlim([trialnum(1) trialnum(end)]); 
    ylim([-30 0]) %xlabel(' trial #');
    set(gca,'XTickLabels',[])
    title('performance level')

    nexttile(2) % subplot(2,1,2); 
    plot(trialnum, bias_reject-1,'k:');
    xlabel(' trial #'); 
    ylabel('rate of ''AM'' answer'); 
    xlim([trialnum(1) trialnum(end)]); 
    ylim([0 1]); 
    title('response bias')
    
    h(end+1) = gcf;
    hname{end+1} = 'fig2suppl';
end
%%% End of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = il_load_data(flags,keyvals)

dir_data = [fastACI_paths('dir_data') 'modulationACI' filesep];
% files = Get_filenames(dir_data,'S*');

if flags.do_fig3a
    files = {'S1'}; % S_AO
elseif flags.do_fig3c
    files = {'S3'}; % LL
elseif flags.do_fig3d
    files = {'S4'}; % S_LV2
else
    % All participants:
    files = {'S1','S2','S3','S4','S5','S6','S7','S8','Srej'};
    % (1): S1
    % (2): S2
    % (3): S3
    % (4): S4
    % (5): S5
    % (6): S6
    % (7): S7
    % (8): S8
    % (9): Srej
end

N_subjects = length(files);
if N_subjects == 9
    if strcmp(files{end},'Srej')
        % The rejected participant
        files = files(1:end-1);
        N_subjects = length(files);
    end
end
data.N_subjects = N_subjects;

for i_subject = N_subjects:-1:1 % First participant read at last...
    dir_local = [dir_data files{i_subject} filesep];
    
    if flags.do_fig2
        %%% Loading the data from 'Behavior.mat'
        CR = []; % Correct rejection, variable from 'Behavior.mat'
        FA = []; % False alarm,       variable from 'Behavior.mat'
        H = [];  % Hit rate,          variable from 'Behavior.mat' 
        M = [];  % Miss rate,         variable from 'Behavior.mat'
        N_m = [];
        % PC_targetabsent = [];
        % PC_targetpresent = [];
        file2load = [dir_local 'Behavior.mat'];
        if exist(file2load,'file')
            load(file2load);
        else
            publ_varnet2022a_utils(files{i_subject},'Get_Behavior',flags);
            load(file2load);
        end
        
        data.hist_N(i_subject,:) = N_m;
        data.hist_H(i_subject,:) = H;
        data.hist_M(i_subject,:) = M;
        data.hist_FA(i_subject,:) = FA;
        data.hist_CR(i_subject,:) = CR;  
        data.m_levels = m_edge(1:end-1)+diff(m_edge)/2;
    end
    
    if flags.do_fig2_suppl
        bias_windowed = [];
        % m_edge = [];
        m_windowed = [];
        trialnum = [];
        load([dir_local 'Behavior'])
        %%% End of loading
        
        data.m(i_subject,:) = m_windowed;
        data.trialnum = trialnum;
        data.bias(i_subject,:) = bias_windowed;
        % PC(i_subject,:) = (PC_targetabsent+PC_targetpresent)/2;
    end
    
    if flags.do_fig3 || flags.do_fig3a || flags.do_fig3c || flags.do_fig3d
        %%% Loading the data from 'CIt'
        CI = [];
        CI1 = [];
        CI1rand_ci = [];
        CI2 = [];
        CI2rand_ci = [];
        CIrand_ci = [];
        ideal_template = [];
        tE = [];
        file2load = [dir_local 'CIt.mat'];
        if exist(file2load,'file')
            load(file2load);
        else
            publ_varnet2022a_utils(files{i_subject},'Get_CIt',flags);
            load(file2load);
        end
        %%% End loading
        data.tE = tE;
        data.ideal_template = ideal_template;
        data.CIt_l(i_subject,:) = CIrand_ci(1,:)';
        data.CIt_h(i_subject,:) = CIrand_ci(2,:)';
        data.CIt_all(i_subject,:) = CI;
        data.CI1t_all(i_subject,:) = CI1;
        data.CI2t_l(i_subject,:) = CI2rand_ci(1,:)';
        data.CI2t_h(i_subject,:) = CI2rand_ci(2,:)';
        data.CI2t_all(i_subject,:) = CI2;
        data.CI1t_l(i_subject,:) = CI1rand_ci(1,:)';
        data.CI1t_h(i_subject,:) = CI1rand_ci(2,:)';
    end
    
    if flags.do_fig4
        %%% Loading the data in 'CIf':
        CIfft = [];
        CIfftrand_ci = [];
        CI1fft = [];
        CI1fftrand_ci = [];
        CI2fft = [];
        CI2fftrand_ci = [];
        
        fE = [];
        % ideal_templatecfft = [];
        file2load = [dir_local 'CIf'];
        if exist(file2load,'file')
            load(file2load)
        else
            publ_varnet2022a_utils(files{i_subject},'Get_CIf',flags);
        end
        %%% End loading
         
        CIf_all(i_subject,:) = abs(CIfft);
        CIf_h(i_subject,:) = CIfftrand_ci(2,:)';
        CIf_l(i_subject,:) = CIfftrand_ci(1,:)';
        
        CI1f_all(i_subject,:) = abs(CI1fft);
        CI1f_h(i_subject,:) = CI1fftrand_ci(2,:)';
        CI1f_l(i_subject,:) = CI1fftrand_ci(1,:)';
        
        CI2f_all(i_subject,:) = abs(CI2fft);
        CI2f_l(i_subject,:) = CI2fftrand_ci(1,:)';
        CI2f_h(i_subject,:) = CI2fftrand_ci(2,:)';
        
        data.fE = fE;
        data.CIf_all = CIf_all;
        data.CIf_l = CIf_l;
        data.CIf_h = CIf_h;
        data.CI1f_all = CI1f_all;
        data.CI1f_h = CI1f_h;
        data.CI1f_l = CI1f_l;
        data.CI2f_all = CI2f_all;
        data.CI2f_l = CI2f_l;
        data.CI2f_h = CI2f_h;
    end
    % %%% Loading the data in 'metrics.mat'
    % Aftarget_CIfft = [];
    % Aftarget_CI1fft = [];
    % Aftarget_CI2fft = [];
    % % Aftarget_CI1fftboot: [1×200 double]
    % % Aftarget_CI1fftrand: [1×200 double]
    % 
    % % Aftarget_CI2fftboot: [1×200 double]
    % % Aftarget_CI2fftrand: [1×200 double]
    % 
    % % Aftarget_CIfftboot: [1×200 double]
    % % Aftarget_CIfftrand: [1×200 double]
    % freqmaxA_CI1fftboot = [];
    % % freqmaxA_CI1fftrand: [1×200 double]
    % freqmaxA_CI2fftboot = [];
    % % freqmaxA_CI2fftrand: [1×200 double]
    % freqmaxA_CIfftboot = [];
    % % freqmaxA_CIfftrand: [1×200 double]
    % maxA_CIfft = [];
    % maxA_CI1fft = [];
    % maxA_CI2fft = [];
    % freqmaxA_CIfft = [];
    % freqmaxA_CI1fft = [];
    % freqmaxA_CI2fft = [];
    % % maxA_CI1fftboot: [1×200 double]
    % % maxA_CI1fftrand: [1×200 double]
    % % maxA_CI2fftboot: [1×200 double]
    % % maxA_CI2fftrand: [1×200 double]
    % % maxA_CIfftboot: [1×200 double]
    % % maxA_CIfftrand: [1×200 double]
    % phaseftarget_CIfft = [];
    % phaseftarget_CI1fft = [];
    % phaseftarget_CI2fft = [];
    % phaseftarget_CI1fftboot = [];
    % % phaseftarget_CI1fftrand: [1×200 double]
    % phaseftarget_CI2fftboot = [];
    % % phaseftarget_CI2fftrand: [1×200 double]
    % phaseftarget_CIfftboot = [];
    % % phaseftarget_CIfftrand: [1×200 double]
    % tdecay_CI = [];
    % tdecay_CI1 = [];
    % tdecay_CI2 = [];
    % tdecay_CIboot = [];
    % tdecay_CI1boot = [];
    % tdecay_CI2boot = [];
    % % tdecay_CI1rand = [];
    % % tdecay_CI2rand: [1×200 double]
    % % tdecay_CIrand: [1×200 double]    
    % load([dir_local 'metrics']);

    
    % tdecayCI_all(i_subject) = tdecay_CI;
    % tdecayCI1_all(i_subject) = tdecay_CI1;
    % tdecayCI2_all(i_subject) = tdecay_CI2;
    % % AftCI_all(i_subject) = Aftarget_CIfft;
    % % AftCI1_all(i_subject) = Aftarget_CI1fft;
    % % AftCI2_all(i_subject) = Aftarget_CI2fft;
    % P=unwrap([-pi/2*ones(size(phaseftarget_CIfft)); phaseftarget_CIfft],[],1);
    % pftCI_all(i_subject) = P(2,:);
    % P=unwrap([-pi/2*ones(size(phaseftarget_CI1fft)); phaseftarget_CI1fft],[],1);
    % pftCI1_all(i_subject) = P(2,:);
    % P=unwrap([-pi/2*ones(size(phaseftarget_CI2fft)); phaseftarget_CI2fft],[],1);
    % pftCI2_all(i_subject) = P(2,:);
    % % maxACI_all(i_subject) = maxA_CIfft;
    % % maxACI1_all(i_subject) = maxA_CI1fft;
    % % maxACI2_all(i_subject) = maxA_CI2fft;
    % fmaxACI_all(i_subject) = freqmaxA_CIfft; % Used for do_fig5
    % fmaxACI1_all(i_subject) = freqmaxA_CI1fft;
    % fmaxACI2_all(i_subject) = freqmaxA_CI2fft;
    % fmaxACIboot_all(i_subject,:) = freqmaxA_CIfftboot;
    % 
    % tdecayCIboot_all(i_subject,:) = tdecay_CIboot;
    % tdecayCI1boot_all(i_subject,:) = tdecay_CI1boot;
    % tdecayCI2boot_all(i_subject,:) = tdecay_CI2boot;
    % % AftCIboot_all(i_subject,:) = Aftarget_CIfftboot;
    % % AftCI1boot_all(i_subject,:) = Aftarget_CI1fftboot;
    % % AftCI2boot_all(i_subject,:) = Aftarget_CI2fftboot;
    % P=unwrap([-pi/2*ones(size(phaseftarget_CIfftboot)); phaseftarget_CIfftboot],[],1);
    % pftCIboot_all(i_subject,:) = P(2,:);
    % % % pftCIboot_all(i_subject,:) = phaseftarget_CIfftboot;
    % P=unwrap([-pi/2*ones(size(phaseftarget_CI1fftboot)); phaseftarget_CI1fftboot],[],1);
    % pftCI1boot_all(i_subject,:) = P(2,:);
    % % % pftCI1boot_all(i_subject,:) = phaseftarget_CI1fftboot;
    % P=unwrap([-pi/2*ones(size(phaseftarget_CI2fftboot)); phaseftarget_CI2fftboot],[],1);
    % pftCI2boot_all(i_subject,:) = P(2,:);
    % % % pftCI2boot_all(i_subject,:) = phaseftarget_CI2fftboot;
    % % maxACIboot_all(i_subject,:) = maxA_CIfftboot;
    % % maxACI1boot_all(i_subject,:) = maxA_CI1fftboot;
    % % maxACI2boot_all(i_subject,:) = maxA_CI2fftboot;
    % fmaxACI1boot_all(i_subject,:) = freqmaxA_CI1fftboot;
    % fmaxACI2boot_all(i_subject,:) = freqmaxA_CI2fftboot;
    % % 
    % % tdecayCIrand_all(i_subject,:) = tdecay_CIrand;
    % % tdecayCI1rand_all(i_subject,:) = tdecay_CI1rand;
    % % tdecayCI2rand_all(i_subject,:) = tdecay_CI2rand;
    % % AftCIrand_all(i_subject,:) = Aftarget_CIfftrand;
    % % AftCI1rand_all(i_subject,:) = Aftarget_CI1fftrand;
    % % AftCI2rand_all(i_subject,:) = Aftarget_CI2fftrand;
    % % pftCIrand_all(i_subject,:) = phaseftarget_CIfftrand;
    % % pftCI1rand_all(i_subject,:) = phaseftarget_CI1fftrand;
    % % pftCI2rand_all(i_subject,:) = phaseftarget_CI2fftrand;
    % % maxACIrand_all(i_subject,:) = maxA_CIfftrand;
    % % maxACI1rand_all(i_subject,:) = maxA_CI1fftrand;
    % % maxACI2rand_all(i_subject,:) = maxA_CI2fftrand;
    % % fmaxACIrand_all(i_subject,:) = freqmaxA_CIfftrand;
    % % fmaxACI1rand_all(i_subject,:) = freqmaxA_CI1fftrand;
    % % fmaxACI2rand_all(i_subject,:) = freqmaxA_CI2fftrand;
    % 
    % if do_figX
    %     error('Continue here')
    %     file2load = [dir_local 'corr.mat'];
    %     if exist(file2load,'file') % from Script3_AnalysisComplex.m
    %     else
    %         % If the file is not on disk then it is created
    %         savefile = Get_filenames(dir_local, 'savegame*');
    %         savefile = [dir_local savefile{1}];
    %         il_Script3_AnalysisComplex(savefile);
    %     end
    %     cor1 = [];
    %     cor2 = [];
    %     cor1_boot = [];
    %     cor2_boot = [];
    %     cor1_rand = [];
    %     cor2_rand = [];
    %     load(file2load)
    % 
    %     cor1_all(i_subject) = cor1;
    %     cor2_all(i_subject) = cor2;
    %     cor1boot_all(i_subject,:) = cor1_boot;
    %     cor2boot_all(i_subject,:) = cor2_boot;
    %     cor1rand_all(i_subject,:) = cor1_rand;
    %     cor2rand_all(i_subject,:) = cor2_rand;
    % end
    % 
    % if do_figY
    %     CIp = prctile(phaseftarget_CIfftboot(:),50);
    %     CI1p = prctile(phaseftarget_CI1fftboot(:),50);
    %     CI2p = prctile(phaseftarget_CI2fftboot(:),50);
    % 
    %     CIp_all(i_subject,:) = CIp;
    %     CI1p_all(i_subject,:) = CI1p;
    %     CI2p_all(i_subject,:) = CI2p;
    % end
end

% %%% Info about the rejected participant:
% dir_local = [dir_data 'Srej' filesep];
% load([dir_local 'Behavior'],'m_windowed','bias_windowed');
data.m_reject = m_windowed;
data.bias_reject = bias_windowed;
