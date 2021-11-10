function outs = Get_all_metrics_from_Praat(dir_where,params)
% function outs = Get_all_metrics_from_Praat(dir_where,params)
%
% Based on l20210713_AnalyseStims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    if ~isunix
        dir_where = 'C:\Users\R_user\Logatome\';
    else
        % dir_where = '/home/alejandro/Documents/Databases/data/fastACI/speechACI_Logatome-abda-S41F/SAO/speech-samples-TEST/'; 
        dir_where = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Stimuli/Logatome/';
    end
end

params = Ensure_field(params,'timestep',0.01); % positive timestep 0.01
params = Ensure_field(params,'nformants',5); % positive nformants 5
params = Ensure_field(params,'maxformant',6000); % positive maxformant 5500
params = Ensure_field(params,'windowlength',0.01); % positive windowlength 0.025
params = Ensure_field(params,'dynamicrange',20); % positive dynamic range 20
params = Ensure_field(params,'minpitch',200); % positive minimum pitch 50 (for intensity)
params = Ensure_field(params,'pitchfloor',50); % positive pitch floor 100 (for f0)
params = Ensure_field(params,'pitchceiling',500); % positive pitch ceiling 500 (for f0)

params = Ensure_field(params,'I_min',0); % if not specified the minimum intensity is set to 0 dB

filesF  = Get_formants_from_dir(dir_where,params);
filesI  = Get_intensity_from_dir(dir_where,params);
filesf0 = Get_f0_from_dir(dir_where,params);

Nsounds = length(filesF);
% minIforF = 75; % Minimum intensity to be added to the plot (relative value)

if Nsounds ~= 0
    for i = 1:Nsounds

        [t_f0{i},f0{i}] = Get_f0_from_txt([dir_where filesf0{i}]);
        [t_I{i} ,I{i}]  = Get_intensity_from_txt([dir_where filesI{i}]);
        [t_F{i} ,F{i}]  = Get_formants_from_txt([dir_where filesF{i}]);

        try
            minIforF = max(max(I{i})-20, params.I_min); % In case the user requests
                % an I_min value that is too low, then the limit is set to 20 dB
                % below the maximum assessed intensity value
        catch
            minIforF = params.I_min;
        end

        idxs = find(I{i}<minIforF | isnan(I{i}));
        F{i}(idxs,:) = nan;
        f0{i}(idxs) = nan;

    end

    outs.t_f0 = t_f0;
    outs.f0 = f0;
    outs.t_I = t_I;
    outs.I = I;
    outs.t_F = t_F;
    outs.F = F;
end

outs.filesf0 = filesf0;
outs.filesF = filesF;
outs.filesI = filesI;
outs.params = params;