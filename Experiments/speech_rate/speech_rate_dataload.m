function [Data_matrix,cfg_ACI] = segmentation_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
% function [Data_matrix,cfg_ACI] = segmentation_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
%
% In this dataload function, the output array Data_matrix is built from the
%   stored f0 vector (cfg_game.f0vec) and time vector (cfg_game.timevec).
%   This means that this dataload version, exceptionally, does not use 
%   directly the list of stimuli (ListStim) or the data_passation variable.
%
% Author: Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_stim  = cfg_game.stim_order; % data_passation.n_stim; % should be the same, remove n_stim
numerase = 0; % to refresh the screen after fprintf (see below)
    
N_trialselect = length(n_stim);

%undersampling = round((undersampling_ms*1e-3)*fs);
%fprintf('%s: Undersampling every %.0f samples will be applied\n',upper(mfilename),undersampling);

for i_trial=1:N_trialselect
    
    fprintf(repmat('\b',1,numerase));
    msg=sprintf('%s: Loading stim no %.0f of %.0f\n',upper(mfilename),i_trial,N_trialselect);
    fprintf(msg);
    numerase=numel(msg);
    
    randomvectors = [cfg_game.scalevec(:,n_stim(i_trial))];

    Data_matrix(i_trial,:,:) = randomvectors;%[noise, fs] = audioread([dir_noise ListStim(n_stim(i_trial)).name ]);

end

cfg_ACI.t = [1];
cfg_ACI.t_description = 'Dimension';
cfg_ACI.f = 1:size(cfg_game.scalevec,1);
cfg_ACI.f_description = 'Segment edge number';