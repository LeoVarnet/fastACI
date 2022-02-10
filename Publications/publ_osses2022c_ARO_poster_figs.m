function [h, hname] = publ_osses2022c_ARO_poster_figs(varargin)
% function [h,hname] = publ_osses2022c_ARO_poster_figs(varargin)
%
% It generates the figures from Osses and Varnet (2022, ARO poster).
% You need to specify the figure number.
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
% %   You can plot the correlations for both participants by using:
%   publ_osses2022c_ARO_poster_figs('fig4');
%
% %   or you can plot the correlations for participant S1 or S2 only by using:
%   publ_osses2022c_ARO_poster_figs('fig4a');
%   publ_osses2022c_ARO_poster_figs('fig4b');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.type={'fig2','fig3','fig3a','fig3b','fig4','fig4a','fig4b'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

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

if flags.do_fig3 || flags.do_fig4
    Subjects = {'SLV','SAO'};
    Subjects_id = [1 2];
end
if flags.do_fig2 || flags.do_fig3a || flags.do_fig4a
    Subjects = {'SLV'};
    Subjects_id = 1;
end
if flags.do_fig3b || flags.do_fig4b
    Subjects = {'SAO'};
    Subjects_id = 2;
end
%%% First, we check whether the participants' data are stored on disk:
for i = 1:length(Subjects)
    for j = length(noise_types):-1:1 % reversed order (in case we need to remove one noise condition)
        outs = publ_osses2022c_ARO_poster_utils(Subjects{i},noise_types{j},'Get_filenames');
        if outs.bGenerate_stimuli == 1
            error('Re-generate waveforms first...')
        end
    end
end

%%% Specific ACI configuration:
TF_type = 'gammatone';
glmfct  = 'lassoglmslow'; 

N_lambda = 30;
Lambdas = logspace(-4, -1, N_lambda);
idx = find(Lambdas >= 10^-3);
Lambdas = Lambdas(idx);
%%%

bPlot_ACIs      = ~(flags.do_fig3 + flags.do_fig3a + flags.do_fig3b==0);
bDo_partialACIs = flags.do_fig4 || flags.do_fig4a || flags.do_fig4b;

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
        flags_for_input = {TF_type, ...
                   glmfct, ... % 'dir_noise',dir_noise, 'dir_target',dir_target, ...
                   'trialtype_analysis', 'total', ...
                   'add_signal',0, ...
                   'apply_SNR',0, ...
                   'skip_if_on_disk',1, ...
                   'expvar_after_reversal', 4, ...
                   'lambda', Lambdas, ...
                   'no_permutation', 'no_bias'}; 
        if bPlot_ACIs == 0
            flags_for_input{end+1} = 'no_plot';
        end

        dir_subj = [fastACI_dir_data experiment filesep subj filesep];
        dir_res = [dir_subj 'Results' filesep];
        fname_results = Get_filenames(dir_res,['savegame*' noise_type '.mat']);
        if length(fname_results) ~= 1
            error('More than one file was found...')
        end
        fname_results = [dir_res fname_results{1}];
      
        % %%% Extra flags if needed:
        % if ~isempty(dir_noise)
        %     flags_for_input{end+1} = 'dir_noise';
        %     flags_for_input{end+1} = dir_noise;
        % end
        % if ~isempty(dir_target)
        %     flags_for_input{end+1} = 'dir_target';
        %     flags_for_input{end+1} = dir_target;
        % end
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
            
            % if isempty(Data_matrix)
            %     % Nothing to do
            % else
            %     flags_here{end+1} = 'Data_matrix';
            %     flags_here{end+1} = Data_matrix;
            % end
            
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

if nargout == 0
    dir_out = fastACI_paths('dir_output');

    for i = 1:length(h)
        opts = [];
        opts.format = 'epsc';
        Saveas(h(i),[dir_out hname{i}],opts);

        opts.format = 'fig';
        Saveas(h(i),[dir_out hname{i}],opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_Post_figure

pause_len = 1;
fprintf('Pausing for %.1f seconds\n',pause_len);
pause(pause_len)