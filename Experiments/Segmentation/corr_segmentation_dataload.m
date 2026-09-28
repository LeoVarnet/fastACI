function [Data_matrix,cfg_ACI] = segmentation_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
% function [Data_matrix,cfg_ACI] = modulationACI_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
%
% data_passation is an input parameter to keep the same function structure as
%   fastACI_getACI_dataload.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_target = cfg_game.dir_target;
dir_noise  = cfg_game.dir_noise;

n_stim  = cfg_game.stim_order; % data_passation.n_stim; % should be the same, remove n_stim
% fcut_noiseE = 30; % Hz
% undersampling_ms = 10; % ms

numerase = 0; % to refresh the screen after fprintf (see below)
% 
% if isfield(cfg_game,'SNR')
%     SNR = cfg_game.SNR; % dB
%     lvl = cfg_game.SPL; % dB SPL
% % else
%     cfg_tmp = modulationACI_set(cfg_game);
%     SNR = cfg_tmp.SNR;
%     lvl = cfg_tmp.SPL; 
% end
% dBFS =  93.6139; % based on cal signal which is 72.6 dB using Sennheiser HD650

% switch cfg_ACI.flags.TF_type
%     case 'gammatone'
%         method = 'Gammatone_proc';
%         fcut = [40 8000]; % Hz, fixed parameter, related to the processing 'gammatone'
% end
    
Nchannel = 2;%length(fcut)-1;

% % ATTENTION BRICOLAGE %
% if ~isempty(dir_target)
%     [tone, fs] = audioread([dir_target 'nontarget.wav']);
% else
%     tone = [];
% end

N_trialselect = length(n_stim);

%undersampling = round((undersampling_ms*1e-3)*fs);
%fprintf('%s: Undersampling every %.0f samples will be applied\n',upper(mfilename),undersampling);

for i_trial=1:N_trialselect
    
    fprintf(repmat('\b',1,numerase));
    msg=sprintf('%s: Loading stim no %.0f of %.0f\n',upper(mfilename),i_trial,N_trialselect);
    fprintf(msg);
    numerase=numel(msg);
    
    randomvectors = [cfg_game.f0vec(:,n_stim(i_trial))*cfg_game.timevec(:,n_stim(i_trial))'];

    Data_matrix(i_trial,:,:) = randomvectors;%[noise, fs] = audioread([dir_noise ListStim(n_stim(i_trial)).name ]);

end

cfg_ACI.t = 1:size(cfg_game.timevec,1);
cfg_ACI.t_description = 'f0';
cfg_ACI.f = 1:size(cfg_game.timevec,1);
cfg_ACI.f_description = 'timing';

%Data_matrix = permute(noise_E, [3 2 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [somme, A] = il_Addition_RSB(S, B, SNR)
% function [somme, A] = il_Addition_RSB(S, B, SNR)
% 
% It adds a sound S to a noise B with a given SNR (in dB), and updates the
% factor A such that the sum is somme = A*S + B;

Ps = mean(S.^2);
Pb = mean(B.^2);

A=sqrt((Pb/Ps)*10^(SNR/10));
somme = A*S + B;
