function outs_from_Praat = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI,outs_from_Praat, Styles, Colours, LineWidth)
% function outs = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI,outs_from_Praat)
% 
% dir_where - it is the directory where the waveforms are located

if nargin < 3
    outs_from_Praat = [];
end
if nargin < 4
    Styles = {'-','-.'};
end
if nargin < 5
    Colours = {[0.5 0.5 0.5],[0 0 0]};
end
if nargin < 6
    LineWidth = 2;
end

if isempty(outs_from_Praat)
    par_formants.timestep = 0.01; % positive timestep 0.01
    par_formants.nformants = 5; % positive nformants 5
    par_formants.maxformant = 5500; % positive maxformant 5500
    par_formants.windowlength = 0.025; % positive windowlength 0.025
    par_formants.dynamicrange = 30; % positive dynamic range 20

    par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
    par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
    par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

    par_formants.I_min = 59;%75;% 40; %, arbitrary value

    outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
else
    % Nothing to do: just using the information in outs_from_Praat
end


Nsounds = length(outs_from_Praat.f0);

labels2add = [];
for kk = 1:Nsounds
    
    if isfield(cfg_ACI,'target_names')
    %    if strfind(outs_from_Praat.filesF{kk},cfg_ACI.target_names{kk})
            labels2add{end+1} = cfg_ACI.target_names{kk};
    %    end
    end

    f2plot = affichage_get_freq_resolution(outs_from_Praat.f0{kk},cfg_ACI); % figure; plot(outs.t_f0{1},outs.f0{1},'k--');
    pl(kk) = plot(outs_from_Praat.t_f0{kk},f2plot,'LineStyle',Styles{1},'Color',Colours{kk},'LineWidth',LineWidth);
    
    for ii = 1:size(outs_from_Praat.F{kk},2)
        f2plot = affichage_get_freq_resolution(outs_from_Praat.F{kk}(:,ii),cfg_ACI);
        hold on; % figure; plot(t_F{kk},F{kk}(:,ii),Style{kk});
        plot(outs_from_Praat.t_F{kk},f2plot,'LineStyle',Styles{2},'Color',Colours{kk},'LineWidth',LineWidth);
    end
end

if ~isempty(labels2add)
    hl = legend(pl,labels2add);
    
    outs_from_Praat.hl = hl;
    outs_from_Praat.hl_description = 'handle of the legend';
end