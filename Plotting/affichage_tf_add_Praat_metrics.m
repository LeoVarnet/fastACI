function outs = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI)
% function outs = affichage_tf_add_Praat_metrics(dir_where,cfg_ACI)
% 
% dir_where - it is the directory where the waveforms are located

par_formants.timestep = 0.01; % positive timestep 0.01
par_formants.nformants = 5; % positive nformants 5
par_formants.maxformant = 5500; % positive maxformant 5500
par_formants.windowlength = 0.025; % positive windowlength 0.025
par_formants.dynamicrange = 30; % positive dynamic range 20

par_formants.minpitch = 200; % positive minimum pitch 50 (for intensity)
par_formants.pitchfloor = 100; % positive pitch floor 100 (for f0)
par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)

par_formants.I_min = 60; % arbitrary value

outs = Get_all_metrics_from_Praat(dir_where,par_formants);

%%%
Styles = {'-','--'};
Colours = {[0.5 0.5 0.5],'k'};

Nsounds = length(outs.f0);

labels2add = [];
for kk = 1:Nsounds
    
    if isfield(cfg_ACI,'target_names')
        if strfind(outs.filesF{kk},cfg_ACI.target_names{kk})
            labels2add{end+1} = cfg_ACI.target_names{kk};
        end
    end

    f2plot = il_get_freq_affichage(outs.f0{kk},cfg_ACI); % figure; plot(outs.t_f0{1},outs.f0{1},'k--');
    pl(kk) = plot(outs.t_f0{kk},f2plot,Styles{kk},'Color',Colours{kk},'LineWidth',2);
    
    for ii = 1:size(outs.F{kk},2)
        f2plot = il_get_freq_affichage(outs.F{kk}(:,ii),cfg_ACI);
        hold on; % figure; plot(t_F{kk},F{kk}(:,ii),Style{kk});
        plot(outs.t_F{kk},f2plot,Styles{kk},'Color',Colours{kk});
    end
end

if ~isempty(labels2add)
    hl = legend(pl,labels2add);
    
    outs.hl = hl;
    outs.hl_description = 'handle of the legend';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2plot = il_get_freq_affichage(f,cfg_ACI)

spacing_type = [];
%%% Detect the frequency spacing:
% 1. Detecting if ERB-spaced:
ferb = freqtoaud(cfg_ACI.f);

df(1) = ferb(2)-ferb(1);
df(2) = ferb(end)-ferb(end-1);

if round(10*df(1)) == round(10*df(2))
    fprintf('\t%s: ERB-spacing detected...\n',upper(mfilename));
    spacing_type = 'erb';
    delta_f = df(1);
end

switch spacing_type
    case 'erb'
        step = (max(cfg_ACI.f)-min(cfg_ACI.f))/(length(cfg_ACI.f)-1); % this is a linear step equal to 0.5 ERB_N
        f2plot = step*(freqtoaud(f)-min(ferb))/delta_f;
        
    otherwise
        warning('Frequency spacing not detected, returning an empty frequency vector')
        f2plot = nan(size(f));
end