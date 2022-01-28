function publ_varnet2013_FrontHumNeurosci_figs(varargin)
% function [h,hname] = publ_varnet2013_FrontHumNeurosci_figs(varargin)
%
% Generates the figures

% % To display Fig. 1 of Varnet et al. (2013, Front. Hum. Neurosci.) use :::
%     publ_varnet2013_FrontHumNeurosci_figs('fig1');
%
% % To display Fig. 4a of Varnet et al. (2013, Front. Hum. Neurosci.) use :::
%     publ_varnet2013_FrontHumNeurosci_figs('fig4a');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_varnet2013_FrontHumNeurosci_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig1','fig4a'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

f_limits = [0 4050]; % Hz, the exact bin is at 4048 Hz (see the paper, p. 4, top)
t_limits = [0 0.3425]; % s, first 0.34 s were accounted for (see the paper, p. 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig1
    [bStored,dir_subj] = publ_varnet2013_FrontHumNeurosci_0_checkdata;
    if bStored == 0
        fprintf('The experimental data were not found on disk\n')
        fprintf('No further actions will be performed...\n')
        return;
    end
    dir_target = Get_filenames(dir_subj,'*Signal');
    if length(dir_target) ~= 1
        error('More than one folder found...')
    end
    dir_target = [dir_subj dir_target{1} filesep];
    %%%
    opts = [];
    opts.window = 'hamming';
    files = {'Aba.wav','Ada.wav'};
    [T_dB,f_spec,t_spec] = Time_frequency_converter(dir_target,files,length(files),opts);
    %%%
    tf_ms = t_limits(end); 
    ff    = f_limits(end); 
    
    idx_t = find(round(100*t_spec)/100 <= tf_ms/1000);
    idx_f = find(round(f_spec) <= ff);
    f_spec = f_spec(idx_f);
    t_spec = t_spec(idx_t);
    T_dB   = T_dB(idx_f,idx_t,:);
    %%%

    max_dB = max(max(max(T_dB)));

    figure;
    subplot(1,3,1)
    plot_stft(t_spec,f_spec,T_dB(:,:,1)-max_dB);
    title('/aba/');
    
    % figure;
    subplot(1,3,2)
    plot_stft(t_spec,f_spec,T_dB(:,:,2)-max_dB);
    title('/ada/');
    
    % figure;
    subplot(1,3,3)
    plot_stft(t_spec,f_spec,T_dB(:,:,1)-T_dB(:,:,2));
    title('/aba/-/ada/');
    
    Pos = get(gcf,'Position');
    Pos(3) = 2.5*Pos(3); % widening the figure
    set(gcf,'Position',Pos);
    
    h(end+1) = gcf;
    hname{end+1} = 'fig1-aba-ada-diff';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig4a
    [bStored,dir_subj] = publ_varnet2013_FrontHumNeurosci_0_checkdata;
    if bStored == 0
        fprintf('The experimental data were not found on disk\n')
        fprintf('No further actions will be performed...\n')
        return;
    end
    %%% Updating dir_target:
    dir_target = Get_filenames(dir_subj,'*Signal');
    if length(dir_target) ~= 1
        error('More than one folder found...')
    end
    dir_target = [dir_subj dir_target{1} filesep];
    %%% Updating dir_noise:
    dir_noise = Get_filenames(dir_subj,'*Bruit');
    if length(dir_noise) ~= 1
        error('More than one folder found...')
    end
    dir_noise = [dir_subj dir_noise{1} filesep];
    %%% The expected results should be something like:
    % dir_noise  = '/home/alejandro/Documents/Databases/data/fastACI_data/data_varnet2013/Sujet_Leo_S1/ListeBruit/';
    % dir_target = '/home/alejandro/Documents/Databases/data/fastACI_data/data_varnet2013/Sujet_Leo_S1/ListeSignal/';
    
    fname_results = [dir_subj 'savegame_final.mat'];

    flags_for_input = {'dir_target',dir_target, ...
                       'dir_noise' ,dir_noise, ...
                       'varnet2013'}; % loads the group flags
    % The use of the flag 'varnet2013' is equivalent to specify:
    % TF_type = 'spect';
    % glmfct = 'glmfitqp';
    % flags_for_input = {TF_type,glmfct,'trialtype_analysis', 'total',...
    %                 'N_folds', 10, 'add_signal',0, ...
    %                 'apply_SNR',0, 'skip_if_on_disk',1, ...
    %                 'permutation', ...
    %                 'f_limits', f_limits,'t_limits', t_limits};
    [ACI,cfg_ACI,results] = fastACI_getACI(fname_results,flags_for_input{:});
    h(end+1) = gcf;
    hname{end+1} = 'fig4a-middle-SLV';
end