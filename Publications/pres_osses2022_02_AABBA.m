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
% definput.keyvals.models=[];
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
    
if flags.do_fig_page13_top % Old: bPlot_ACI
    % g20220113_analysis_pipeline_versie_Alejandro
    
    Pos3 = 350;
    
    noise_type = 'white';
    N_subjects = length(subjects);
    
    figure; 
    tiledlayout(1,N_subjects);
    for i_subj = 1:N_subjects
        
        subj_here = subjects{i_subj};
        [h1,hname1,outs1] = pres_osses2022_02_AABBA_utils(experiment_full, subj_here, noise_type,flags_here{:});
        
        % figure(handle_here);
        nexttile;
        affichage_tf(outs1.ACI,'CI', 'cfg', outs1.cfg_ACI); hold on
        outs_aff = affichage_tf_add_Praat_metrics(outs1.cfg_game.dir_target, ...
            outs1.cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
        
        if i_subj ~= N_subjects
            set(outs_aff.hl,'Visible','off'); % Removes the legend
        end
        title(subj_here)
    end
    Pos = get(gcf,'Position');
    Pos(3:4) = [1500 420];
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig13-top';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig_page13_bottom % (old 'bPlot_thres_corr' variable)
    
    experiment      = strsplit(experiment_full,'-'); experiment = experiment{1};
    lab4plot = [];

    noise_type = 'white'; % 'sMPSv1p3';

    for i = 1:length(subjects)
        subject    = subjects{i};
 
        dir_subject = [fastACI_dir_data experiment_full filesep subject filesep];
 
        fname_suff = ['ACI-' subject '-' experiment '-' noise_type '-nob-gt-l1glm.mat'];
        fname_save_filt = ['savegame_*' subject '_' experiment '*' noise_type '*.mat'];
 
        switch subject
            case 'SLV'
                % The only normal-hearing participant...
                dirs = Get_filenames(dir_subject,'Results*');
                bInput = length(dirs);

                dir_res = [dirs{bInput} filesep];
                fi = Get_filenames([dir_subject dir_res 'Results_ACI' filesep], [fname_suff(1:end-4) '*rev*.mat']);
                fname_suff = fi{1};
                disp('')

            otherwise
                % Then it is a model:
                dirs = Get_filenames(dir_subject,'Run*');
                Show_cell(dirs);
                disp('the last folder will be used...')
                bInput = length(dirs);

                dir_res = [dirs{bInput} filesep 'Results' filesep];
        end
        dir_results = [dir_subject dir_res 'Results_ACI' filesep];

        fname = [dir_results  fname_suff];
        var = load(fname);

        ACIs(:,i) = var.ACI(:);
        lab4plot{end+1} = subject;

        dir_results = [dir_subject dir_res];
        fi = Get_filenames(dir_results,fname_save_filt);
        fname = [dir_results  fi{1}];
        var = load(fname);

        expvar = var.data_passation.expvar;
        expvar_me(i) = median(expvar);
        expvar_L(i) = prctile(expvar,25);
        expvar_U(i) = prctile(expvar,75);

    end

    errL = expvar_me - expvar_L;
    errU = expvar_U - expvar_me;
    idx = [2 3 6 4 5 1];
    x_var = 1:length(subjects);

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
    
    opts = [];
    % opts.corrmat = corr_ma;
    opts.title = type;
    opts.labels = lab4plot;
    opts.label_FontSize = 10;
    opts.outline = 1;
    opts.type = type;
    [h_corrmat, h_colorbar, corrmat] = my_plot_corrmat(ACIs,opts); 
    
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