function outs_from_Praat = affichage_tf_add_Praat_metrics_one_sound(fname_full,cfg_ACI,outs_from_Praat, Style, Colour, LineWidth,bPlot)
% function outs = affichage_tf_add_Praat_metrics_one_sound(fname_full,cfg_ACI,outs_from_Praat, Style, Colour, LineWidth,bPlot)
% 
% fname_full should be the name of a sound to be processed. As a difference
%   to affichage_tf_add_Praat_metrics.m, only one waveform will be plotted
%   and no legend will be added to the plot.
%
% See also: affichage_tf_add_Praat_metrics.m
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    outs_from_Praat = [];
end
if nargin < 4
    Style = {'-','-.'};
end
if nargin < 5
    Colour = {[0.5 0.5 0.5],'k'};
end
if nargin < 6
    LineWidth = 1; % default as used in osses2021c
end
if nargin < 7
    bPlot = 1;
end
[dir_where,fname,ext] = fileparts(fname_full); % fname will not have an extension
dir_where = [dir_where filesep];
% fname = [fname ext];

if isempty(outs_from_Praat)
    warning('Praat results not found on disk, Praat will be run again using default values')
    
    par_formants.timestep = 0.005; % positive timestep 0.01
    par_formants.nformants = 5; % positive nformants 5
    
    % Formants
    par_formants.maxformant = 5250; % positive maxformant 5500
    par_formants.windowlength = 0.025;% 0.025 % positive windowlength 0.025
    par_formants.dynamicrange = 20;%30; % positive dynamic range 20
    
    % F0
    %par_formants.minpitch = 200; % previous parameter value (14/10/2022)
    par_formants.minpitch = 250; % positive minimum pitch 50 (for intensity)
    %par_formants.pitchfloor = 100; % previous parameter value (14/10/2022)
    par_formants.pitchfloor = 50; % positive pitch floor 100 (for f0)
    par_formants.pitchceiling = 500; % positive pitch ceiling 500 (for f0)
    
    % Before 4/11/2021, I_min set to 40 dB:
    par_formants.I_min = 0;%59;%75; %, arbitrary value
    outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
else
    % Nothing to do: just using the information in outs_from_Praat
    if ~isfield(outs_from_Praat,'f0')
        % This means that maybe you have not yet the Praat parameters
        par_formants = outs_from_Praat;
        outs_from_Praat = Get_all_metrics_from_Praat(dir_where,par_formants);
    end
end

if isfield(outs_from_Praat,'f0')
    Nsounds = length(outs_from_Praat.f0);
else
    disp('');
end
    
% labels2add = [];
bAdd_traces = zeros([1 Nsounds]);


for kk = 1:Nsounds

    % if isfield(cfg_ACI,'target_names')
    try
       if strfind(outs_from_Praat.filesF{kk},fname)
            bAdd_traces(kk) = 1;
       end
    end
    % end

    delay = cfg_ACI.t(1); %compensation for the gammatone delay 
    if bPlot
        if bAdd_traces(kk)
            [f2plot,cfg_ACI] = affichage_get_freq_resolution(outs_from_Praat.f0{kk},cfg_ACI); % figure; plot(outs.t_f0{1},outs.f0{1},'k--');
            pl(kk) = plot(outs_from_Praat.t_f0{kk}+delay,f2plot,'LineStyle',Style,'Color',Colour,'LineWidth',LineWidth);

            for ii = 1:size(outs_from_Praat.F{kk},2)
                % Adding each formant:
                f2plot = affichage_get_freq_resolution(outs_from_Praat.F{kk}(:,ii),cfg_ACI);
                hold on; % figure; plot(t_F{kk},F{kk}(:,ii),Style{kk});
                plot(outs_from_Praat.t_F{kk}+delay,f2plot,'LineStyle',Style,'Color',Colour,'LineWidth',LineWidth);
            end
        end
    end
end

if nargout ~= 0
    for kk = Nsounds:-1:1
        if bAdd_traces(kk) == 0
            if isfield(outs_from_Praat,'t_f0')
                outs_from_Praat.t_f0(kk) = [];
            end
            if isfield(outs_from_Praat,'t_F')
                outs_from_Praat.t_F(kk) = [];
            end
            if isfield(outs_from_Praat,'t_I')
                outs_from_Praat.t_I(kk) = [];
            end
            if isfield(outs_from_Praat,'f0')
                outs_from_Praat.f0(kk) = [];
            end
            if isfield(outs_from_Praat,'F')
                outs_from_Praat.F(kk) = [];
            end
            if isfield(outs_from_Praat,'t_f0')
                outs_from_Praat.I(kk) = [];
            end
        end
    end
end
