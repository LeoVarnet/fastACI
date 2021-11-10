function outs_from_Praat = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI,outs_from_Praat, Styles, Colours, LineWidth)
% function outs = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI,outs_from_Praat, Styles, Colours, LineWidth)
% 
% dir_where - it is the directory where the waveforms are located
%
% Test this script by running (for instance):
%     publ_osses2021c_DAGA_2_figs('fig1a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    outs_from_Praat = [];
end
if nargin < 4
    Styles = {'-','-.'};
end
if nargin < 5
    Colours = {[0.5 0.5 0.5],'k'};
end
if nargin < 6
    LineWidth = 1; % default as used in osses2021c
end

if isempty(outs_from_Praat)
    warning('Praat results not found on disk, Praat will be run again using default values')
    
    par_formants.timestep = 0.01; % positive timestep 0.01
    par_formants.nformants = 5; % positive nformants 5
    par_formants.maxformant = 5500; % positive maxformant 5500
    par_formants.windowlength = 0.025; % positive windowlength 0.025
    par_formants.dynamicrange = 30; % positive dynamic range 20

    par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
    par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
    par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

    % Before 4/11/2021, I_min set to 40 dB:
    par_formants.I_min = 59;%75; %, arbitrary value

    outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
else
    % Nothing to do: just using the information in outs_from_Praat
end

Nsounds = length(outs_from_Praat.f0);
if Nsounds ~= length(Styles)
    warning('Metrics of more than two waveforms are being plotted: Repeating the format used in the second sound...');
    disp('')
end

labels2add = [];
for kk = 1:Nsounds
    
    if isfield(cfg_ACI,'target_names')
    %    if strfind(outs_from_Praat.filesF{kk},cfg_ACI.target_names{kk})
            labels2add{end+1} = cfg_ACI.target_names{kk};
    %    end
    end

    f2plot = affichage_get_freq_resolution(outs_from_Praat.f0{kk},cfg_ACI); % figure; plot(outs.t_f0{1},outs.f0{1},'k--');
    pl(kk) = plot(outs_from_Praat.t_f0{kk},f2plot,'LineStyle',Styles{kk},'Color',Colours{kk},'LineWidth',LineWidth);
    
    for ii = 1:size(outs_from_Praat.F{kk},2)
        % Adding each formant:
        f2plot = affichage_get_freq_resolution(outs_from_Praat.F{kk}(:,ii),cfg_ACI);
        hold on; % figure; plot(t_F{kk},F{kk}(:,ii),Style{kk});
        plot(outs_from_Praat.t_F{kk},f2plot,'LineStyle',Styles{kk},'Color',Colours{kk},'LineWidth',LineWidth);
    end
end

if ~isempty(labels2add)
    hl = legend(pl,labels2add);
    
    outs_from_Praat.hl = hl;
    outs_from_Praat.hl_description = 'handle of the legend';
end