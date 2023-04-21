function data = pres_osses2023_04_WASdag(varargin)
% functiong data = pres_osses2023_04_WASdag(varargin)
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
% See also: publ_osses2022b_JASA_figs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % close all, clc
% 
% if nargin == 0
%     help publ_osses2022b_JASA_figs;
%     return
% end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig1', ... % ACIs using 3 types GLM fitting algorithms
    'fig2a', ... % ACIs using 2 types GLM fitting algorithms
    'fig2b', ...
    'fig3a', ...
    'fig3b'};
definput.flags.plot={'plot','no_plot'};
definput.flags.local={'local','zenodo'};
 
definput.keyvals.dir_zenodo=[];
definput.keyvals.dir_out=[];
 
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

dir_out = keyvals.dir_out;

% Common variables:
noise_str         = {'bumpv1p2_10dB'};

experiment        = 'speechACI_Logatome-abda-S43M';
% N_maskers = length(noise_str);

bZenodo = flags.do_zenodo;
if bZenodo == 1
    fprintf('%s: The data from Zenodo will be loaded... \n',mfilename);
    if ~isempty(keyvals.dir_zenodo)
        fprintf('\t The data from Zenodo will be loaded from %s \n',keyvals.dir_zenodo);
        if exist(keyvals.dir_zenodo,'dir')
            bIs_dir = 1;
        else
            bIs_dir = 0;
        end
    else
        bIs_dir = 0;
    end
    if bIs_dir == 0
        msg_here = sprintf('\t You indicated that you want to proccess the data obtained from Zenodo. Please indicate the directory where it is...\n');
        fprintf(msg_here)
        keyvals.dir_zenodo = uigetdir([pwd filesep],msg_here);
        keyvals.dir_zenodo = [keyvals.dir_zenodo filesep];
    end
end

% % colourbar_map = 'default';
 
colourbar_map = 'DG_jet';
flags_tf = {'colourbar_map',colourbar_map};
 
if bZenodo == 1
    dir_zenodo = keyvals.dir_zenodo;
    dir_ACI_exp = [dir_zenodo '03-Post-proc-data' filesep 'ACI_exp' filesep];
%     dir_ACI_sim = [dir_zenodo '03-Post-proc-data' filesep 'ACI_sim' filesep]; % Copied from /20220715_sim_Q1_osses2022a/';
%     
%     dir_data = [dir_zenodo '01-Stimuli' filesep 'fastACI_data' filesep];
%     dir_fastACI      = [dir_zenodo '02-Raw-data' filesep 'fastACI' filesep];
%     dir_savegame_sim = [dir_zenodo '02-Raw-data' filesep 'ACI_sim' filesep];
% 
else
    if ~isempty(dir_out)
        dir_ACI_exp = dir_out;
    else
        % directory in Alejandro's computer:
        dir_ACI_exp = '/home/alejandro/Documents/Databases/data/fastACI_datapost/speechACI_Logatome-abda-S43M/S04/';
    end
end
% dir_savegame_exp = [dir_fastACI 'Publications' filesep 'publ_osses2022b' filesep]; 

if flags.do_fig1
    Subjects = {'S04'};
    glmfct2test = {'classic_revcorr','glmfitqp','l1glm'}; 
end
if flags.do_fig2a
    Subjects = {'S04'};
    glmfct2test = {'classic_revcorr','l1glm'}; 
end
if flags.do_fig2b
    Subjects = {'S08'};
    glmfct2test = {'classic_revcorr','l1glm'}; 
end

if flags.do_fig3a
    Subjects = {'S04'};
    glmfct2test = {'classic_revcorr','l1glm'}; 
end
if flags.do_fig3b
    Subjects = {'S08'};
    glmfct2test = {'classic_revcorr','l1glm'}; 
end

% glmfct2test = {'classic_revcorr','l1glm'}; 

% Subjects = {'S08'}; 
N_subjects = length(Subjects);
 
N_ACIs = length(glmfct2test);
i_subj = 1;
i_noise = 1;

%%%
dir_results = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep ...
                    'data_' Subjects{i_subj} filesep '1-experimental_results' filesep];
        
fname_results = Get_filenames(dir_results,['savegame_*' noise_str{i_noise} '.mat']);
fname_results = fname_results{end};
%%%

if flags.do_fig1 || flags.do_fig2a || flags.do_fig2b || flags.do_fig3a || flags.do_fig3b
    for i = 1:N_ACIs
        glmfct = glmfct2test{i};
        switch glmfct
            case 'l1glm'
                fname_suffix{i} = ['ACI-' Subjects{i_subj} '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4.mat'];
            case 'classic_revcorr'
                fname_suffix{i} = ['ACI-' Subjects{i_subj} '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-classic_revcorr+pyrga-rev4.mat'];
            case 'glmfitqp'
                fname_suffix{i} = ['ACI-' Subjects{i_subj} '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-glmfitqp+pyrga-rev4.mat'];
        end
        
        bCalculate = ~exist([dir_ACI_exp fname_suffix{i}],'file');
        if bCalculate
            trialtype_analysis = 'total';   % prereg, Sec. 6.1, after Eq. (1)
            TF_type = 'gammatone';          % prereg, Sec. 6.2.1

            flags_in = {TF_type, ...
                glmfct, ...
                'no_bias', ... % nob
                'expvar_after_reversal',4,'dir_out',dir_out}; %'idx_trialselect', 1:4000 ... % tr4000};
            [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI([dir_results fname_results],flags_in{:});
        end
        
    end
            
end

% dir_where = dir_ACI_exp;
% Location of the savegames:
           % /home/alejandro/Desktop/fastACI_today/fastACI_data/speechACI_Logatome-abda-S43M/S01/Results
% dir_res = {'/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Publications/publ_osses2022b/',['1-experimental_results' filesep]};
% dir_wav = '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/';

bAdd_formants = 1;
if bAdd_formants
    var = load([dir_results fname_results]);
    
    dir_here = [fastACI_dir_data experiment filesep 'S01' filesep 'speech-samples' filesep];
    warning('temporal dir_here')
    fname_aba = [dir_here var.cfg_game.filename_target{1}];
    fname_ada = [dir_here var.cfg_game.filename_target{2}];
    disp('')
end

if flags.do_fig1 || flags.do_fig2a || flags.do_fig2b
    N_noises = length(noise_str);
    flags_text = {'FontWeight','Bold','Units','Normalized'};

    % bAdd_thres = 1; % Option added on 10/10/2022
    % if bAdd_thres 
    %     try
    %         flags_here = varargin(2:end); % excluding the fig number
    %         meanSNR = publ_osses2022b_JASA_figs('fig4','no_plot',flags_here{:});
    %     catch
    %         meanSNR = publ_osses2022b_JASA_figs('fig4','no_plot');
    %     end
    %     file_savegame = meanSNR.file_savegame;
    %     meanSNR = meanSNR.meanSNR; % recycling the variable
    % end

    figure('Position',[100 100 800 350]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
    tiledlayout(N_noises,N_ACIs,'TileSpacing','tight');
    
    for i_aci = 1:N_ACIs
        nexttile(i_aci);
        ACI_fname = [dir_ACI_exp fname_suffix{i_aci}];
  
        ACI = [];
        cfg_ACI = [];
        results = [];
        info_toolbox = [];
        load(ACI_fname);

        if i_aci == N_ACIs
            bColourbar = 1;
        else
            bColourbar = 0;
        end
        flags_opts = {'NfrequencyTicks',9,'colorbar',bColourbar}; 
        flags_opts = [flags_opts flags_tf];
                 
        plot_outs = affichage_tf(ACI,'CInorm', cfg_ACI.t, cfg_ACI.f,flags_opts{:}); hold on
        
        if bAdd_formants
            basef = 8000;
            flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',0.5, ...
                'dboffset',100,'no_adt','binwidth',0.01,'no_outerear','no_middleear'};

            [aba, fs] = audioread(fname_aba);
            [ada, fs] = audioread(fname_ada);

            %%% Plot targets + 1 example of noise
            [G_aba, fc, t, outs] = Gammatone_proc(aba, fs, flags_gamma{:});
            [G_ada, fc, t, outs] = Gammatone_proc(ada, fs, flags_gamma{:});

            cfg_in = [];
            cfg_in.f = fc;
            cfg_in.t = t;

            %%% Parameters for Praat:
            par_formants = [];
            par_formants.timestep = 0.01; % positive timestep 0.01
            par_formants.nformants = 5; % positive nformants 5

            % Formants
            par_formants.maxformant = 5500; % positive maxformant 5500
            par_formants.windowlength = 0.025;% 0.025 % positive windowlength 0.025
            par_formants.dynamicrange = 30; % positive dynamic range 20

            % F0
            par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
            par_formants.pitchfloor = 50; % positive pitch floor 100 (for f0)
            par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

            % Before 4/11/2021, I_min set to 40 dB:
            par_formants.I_min = 59;%75; %, arbitrary value

            Add_formants(fname_aba,cfg_in,par_formants,'r');
            Add_formants(fname_ada,cfg_in,par_formants,'b');
            %%%
        end
        
        if i_aci ~= 1
            set(gca,'YTickLabel','');
            ylabel('');
        end
        
        if bColourbar
            hcl = plot_outs.tcolorbar;
            set(hcl,'Ticks',[]);
        end
        % Put 'aba'
        if bColourbar
            [xx,colour_max,colour_min] = Get_colourmap_rgb(colourbar_map);
            text(1.05,1.05,'aba','FontSize',12,flags_text{:},'Color',colour_max);
        end

        if bColourbar
            text(1.05,-.05,'ada','FontSize',12,flags_text{:},'Color',colour_min);
        end
        title(glmfct2test{i_aci},'interpreter','none');

    end
    h(end+1) = gcf;
    if flags.do_fig1
        fig_pref = 'fig01';
    end
    if flags.do_fig2a
        fig_pref = 'fig02a';
    end
    if flags.do_fig2b
        fig_pref = 'fig02b';
    end
    
    hname{end+1} = [fig_pref '-ACI-exp'];

end % end do_fig1, do_fig2

if flags.do_fig3a || flags.do_fig3b
    N_noises = length(noise_str);
    flags_text = {'FontWeight','Bold','Units','Normalized'};

    % figure('Position',[100 100 1200 350]); % set(gcf,'Position',[100 100 600 Pos4(i_repeat)])
    % tiledlayout(N_noises,N_ACIs,'TileSpacing','tight');
    
    for i_aci = 1:N_ACIs
        % nexttile(i_aci);
        suff = '';
        ACI_fname = [dir_ACI_exp fname_suffix{i_aci}];
  
        ACI = [];
        cfg_ACI = [];
        results = [];
        info_toolbox = [];
        load(ACI_fname);

        switch glmfct2test{i_aci}
            case 'classic_revcorr'
                PA(i_aci,:) = results.FitInfo.PC;
            case 'l1glm'
                idxlambda = results.idxlambda;
                PA(i_aci,:) = results.FitInfo.PC_test(idxlambda,:); % for the 10 folds
        end
    end
    Me_PA = mean(100*PA,2);
    Std_PA = std(100*PA,[],2);        
    
    figure;
    errorbar([1 2],Me_PA,Std_PA,'bo','LineWidth',2,'MarkerFaceColor','b');
    xlim([0 3]); grid on;
    ylim([49.5 60.5]);
    
    Pos = get(gcf,'Position');
    Pos(3:4) = [170 420];
    set(gcf,'Position',Pos);
    
    ylabel('PA (%)');
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabels','');
    disp('')
    
    h(end+1) = gcf;
    if flags.do_fig3a
        hname{end+1} = 'fig03a-PA';
    end
    if flags.do_fig3b
        hname{end+1} = 'fig03b-PA';
    end
    

end % end do_fig7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.h = h;
data.hname = hname;