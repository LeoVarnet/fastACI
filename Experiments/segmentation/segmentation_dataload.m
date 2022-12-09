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
    
    randomvectors = [cfg_game.f0vec(:,n_stim(i_trial)) cfg_game.timevec(:,n_stim(i_trial))];

    Data_matrix(i_trial,:,:) = randomvectors;%[noise, fs] = audioread([dir_noise ListStim(n_stim(i_trial)).name ]);

end
%     
%     % ATTENTION BRICOLAGE %
%     if ~isempty(tone)
%         %   Creating the trial (but always using the non-modulated sound)
%         S = il_Addition_RSB(tone,noise,SNR);
%         S = scaletodbspl(S,lvl,dBFS); % same as dBlvl(S,lvl) with 72.6 as lvl_ref
%     else
%         error('Alejandro: Not validated')
%     end
%     
%     for i_channel = 1:Nchannel
%         switch method
%             case 'butter'
%                 % Filter design:
%                 [B, A]   = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
%                 % Applying the filter:
%                 Sout = abs(filterfilter(B, A, S)); % absolute value is full wave rectification
%                 % Sout = max(filterfilter(B, A, S),0); % This would be half-wave rectification
%                 if ~isempty(tone)
%                      % Regular processing:
%                      [B2, A2] = butter(1, fcut_noiseE/(fs/2));
%                      E = filterfilter(B2, A2, Sout);
%                 else
%                     if i_trial == 1 && i_channel == 1
%                         warning('Alejandro: I do not agree')
%                         E = Sout;
%                     end
%                 end
%                 E = E(1:undersampling:end);
%                 noise_E(i_channel, :, i_trial) = E;
% 
%                 t=(1:size(E,2))/(fs/undersampling);
%                 fc = fcut;
%             case 'Gammatone_proc'
%                  binwidth = undersampling/fs; % 'undersampling' samples
%                  flags_gamma = {'basef',fcut(i_channel+1),'flow',fcut(i_channel),'fhigh',fcut(i_channel+1), ...
%                      'bwmul',.5,'dboffset',dBFS,'no_adt','binwidth',binwidth, ...
%                     'no_outerear','no_middleear', 'hilbert'};
%                  [E, fc, t, outs] = Gammatone_proc(S,fs,flags_gamma{:});
%                  noise_E(:, :, i_trial) = E;
%         end
%     end
% end
cfg_ACI.t = [1,2];
cfg_ACI.t_description = 'Dimension';
cfg_ACI.f = 1:size(cfg_game.timevec,1);
cfg_ACI.f_description = 'Segment edge number';

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
