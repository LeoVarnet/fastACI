function pres_osses2022_02_AABBA(varargin)
% functiong pres_osses2022_02_AABBA(varargin)
%
% Generates some figures from the presentation titled 'Paris II - Simulating 
%   the perception of soundscapes, speech-, AM-, and FM-sounds' given in the
%   context of the AABBA meeting celebrated in Vienna on 3 and 4 February 2022.
%
% At this moment, this script requires that you already have the simulation
%    results stored locally in your computer.
%
% % To display the figure in page 11 of the final presentation (PDF) use :::
%     % ACI of participant S2 from Osses & Varnet (2021, DAGA):
%     pres_osses2022_02_AABBA('fig_page11');
%
% % To display the top figure in page 13 use :::
%     pres_osses2022_02_AABBA('fig_page13_top'); % ACIs for models and the human listener
%
% % To display the bottom figure in page 13 use :::
%     pres_osses2022_02_AABBA('fig_page13_bottom'); % thresholds and correlation plot
%
% Author: Alejandro Osses
% Original name: g20220202_analysis_pipeline_sim_AABBA.m (fastACI_sim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help pres_osses2022_02_AABBA;
    return
end

close all

h = [];
hname = [];

definput.flags.type={'missingflag','fig_page11','fig_page13_top','fig_page13_bottom'};
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%%% Common variables:
experiment_full = 'speechACI_Logatome-abda-S43M';
%%%    

if ~isempty(keyvals.dir_out)
    flags_here = {'dir_out',keyvals.dir_out};
else
    flags_here = {};
end
if flags.do_fig_page11
    %%% ACI for participant S2 from Osses2021c (DAGA paper):
    data = publ_osses2021c_DAGA_2_figs('fig1b',flags_here{:});
    
    for i = 1:length(data.h)
        if strfind(data.hname{i},'ACI-')
            h(end+1) = data.h(i);
            hname{end+1} = data.hname{i};
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment = 'speechACI_Logatome-abda-S43M';
 
subjects   = {'dau1997','king2019','relanoiborra2019','osses2021','osses2022a','SLV'};
noise_type = 'white';    
N_subjects = length(subjects);

if flags.do_fig_page13_top
    % g20220113_analysis_pipeline_versie_Alejandro
    figure; 
    il_tiledlayout(1,N_subjects);
    
    tile_current = 0;
    for i_subj = 1:N_subjects
        tile_current = tile_current + 1;
        
        subj_here = subjects{i_subj};
        [h1,hname1,outs1] = pres_osses2022_02_AABBA_utils(experiment_full, subj_here, noise_type,flags_here{:});
        
        il_nexttile(tile_current); 
        affichage_tf(outs1.ACI,'CI', 'cfg', outs1.cfg_ACI); hold on
        outs_aff = affichage_tf_add_Praat_metrics(outs1.cfg_game.dir_target, ...
            outs1.cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);

        if i_subj ~= N_subjects
            set(outs_aff.hl,'Visible','off'); % Removes the legend
        end
        title(subj_here);
        h(end+1) = gcf; % every time obtains the current figure handle
    end
    
    h = unique(h);
    
    Pos = get(gcf,'Position');
    if length(h) == 1 
        % then all figures are tiled
        Pos(3:4) = [1500 420];
    else
        Pos(3:4) = [250 420];
    end
    for i_h = 1:length(h)
        set(h(i_h),'Position',Pos);
        if i_h ~= 1
            suff = num2str(i_h);
        else
            suff = '';
        end
        hname{end+1} = ['fig13-top' suff];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig_page13_bottom
    lab4plot = []; % for do_fig_page13_bottom
    
    for i_subj = 1:N_subjects
        subj_here = subjects{i_subj};
        [h1,hname1,outs1] = pres_osses2022_02_AABBA_utils(experiment_full, subj_here, noise_type,flags_here{:});
        
        ACIs(:,i_subj) = outs1.ACI(:);
        lab4plot{end+1} = subj_here;
        
        expvar = outs1.data_passation.expvar;
        expvar_me(i_subj) = median(expvar);
        expvar_L(i_subj) = prctile(expvar,25);
        expvar_U(i_subj) = prctile(expvar,75);
    end
    
    errL = expvar_me - expvar_L;
    errU = expvar_U - expvar_me;
    idx = 1:6; % [2 3 6 4 5 1];
    x_var = 1:N_subjects;

    figure;
    errorbar(x_var,expvar_me(idx),errL(idx),errU(idx),'bo','LineWidth',2);

    ylabel('SNR (dB)')
    xlim([0 length(subjects)+1]);
    set(gca,'XTick',x_var);
    set(gca,'XTickLabel',subjects(idx));
    try
        % If you use an old MATLAB version you won't get an error
        set(gca,'XTickLabelRotation',35);
    end
    ylim([-40 0])    
    set(gca,'YTick',-35:5:0);
    grid on

    h(end+1) = gcf;
    hname{end+1} = 'thres';

    Pos = get(gcf,'Position');
    Pos(4) = .7*Pos(4);
    set(gcf,'Position',Pos);

    %%%
    type = 'Spearman';
    
    idx = [6 1:2 4:5 3]; % to put SLV first
    opts = [];
    % opts.corrmat = corr_ma;
    opts.title = type;
    opts.labels = lab4plot(idx);
    opts.label_FontSize = 10;
    opts.outline = 1;
    opts.type = type;
    [h_corrmat, h_colorbar, corrmat] = my_plot_corrmat(ACIs(:,idx),opts); 
    
    h(end+1) = h_corrmat;
    hname{end+1} = ['corr-' type];
end

bSave = input('Do you want to save the figures? (1=yes; 0=no): ');

if bSave == 1
    dir_out = fastACI_paths('dir_output');

    for i = 1:length(h)
        opts = [];
        opts.format = 'epsc';
        fname = [dir_out hname{i}];
        Saveas(h(i),fname,opts);

        opts.format = 'png';
        Saveas(h(i),fname,opts);
        
        opts.format = 'fig';
        Saveas(h(i),fname,opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_tiledlayout(N,M,TileSpacing,TileSpacing_option)

if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout','file'); % tiledlayout.p
if bExist
    tiledlayout(N,M,TileSpacing,TileSpacing_option);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_nexttile(N)

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