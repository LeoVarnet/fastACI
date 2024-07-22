function [Data_matrix,cfg_ACI] = replication_ahumada1975_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
% function [Data_matrix,cfg_ACI] = modulationACI_dataload(cfg_ACI, ListStim, cfg_game, data_passation)
%
% data_passation is an input parameter to keep the same function structure as
%   fastACI_getACI_dataload.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dir_target = cfg_game.dir_target;
dir_noise  = cfg_game.dir_noise;

n_stim  = cfg_game.stim_order; % data_passation.n_stim; % should be the same, remove n_stim


numerase = 0; % to refresh the screen after fprintf (see below)

if isfield(cfg_game,'SNR')
    SNR = cfg_game.SNR; % dB
    lvl = cfg_game.SPL; % dB SPL
else
    cfg_tmp = modulationACI_set(cfg_game);
    SNR = cfg_tmp.SNR;
    lvl = cfg_tmp.SPL; 
end
dBFS =  93.6139; % based on cal signal which is 72.6 dB using Sennheiser HD650
% 
% switch cfg_ACI.flags.TF_type
%     case 'spect'
%         
%     case 'gammatone'
%         method = 'Gammatone_proc';
%         fcut = [40 8000]; % Hz, fixed parameter, related to the processing 'gammatone'
% end
%     
%Nchannel = length(fcut)-1;

N_trialselect = length(n_stim);

for i_trial=1:N_trialselect

    fprintf(repmat('\b',1,numerase));
    msg=sprintf('%s: Loading stim no %.0f of %.0f\n',upper(mfilename),i_trial,N_trialselect);
    fprintf(msg);
    numerase=numel(msg);

    [noise, fs] = audioread([dir_noise ListStim(n_stim(i_trial)).name ]);

    Nsamples = length(noise);
    Nsample_seg = 0.1*fs;
    Nseg = floor(Nsamples/Nsample_seg);
    f = (0:Nsample_seg/2-1)*fs/Nsample_seg;
    %%% TODO %%%
    %     %   Creating the trial (but always using the non-modulated sound)
    %     S = il_Addition_RSB(tone,noise,SNR);
    %     S = scaletodbspl(S,lvl,dBFS); % same as dBlvl(S,lvl) with 72.6 as lvl_ref

    for i_seg = 1:Nseg
        noise_seg = noise((i_seg-1)*Nsample_seg+1:i_seg*Nsample_seg);
        noise_seg_fft = fft(noise_seg); noise_seg_fft = noise_seg_fft(1:end/2);
        noise_spect(1,i_seg,i_trial) = sum(abs(noise_seg_fft(f>=375 & f<425)).^2);
        noise_spect(2,i_seg,i_trial) = sum(abs(noise_seg_fft(f>=425 & f<475)).^2);
        noise_spect(3,i_seg,i_trial) = sum(abs(noise_seg_fft(f>=475 & f<525)).^2);
        noise_spect(4,i_seg,i_trial) = sum(abs(noise_seg_fft(f>=525 & f<575)).^2);
        noise_spect(5,i_seg,i_trial) = sum(abs(noise_seg_fft(f>=575 & f<625)).^2);

        %figure; plot(f,abs(noise_seg_fft));
    end

end
cfg_ACI.t = [0.5:4.5]*0.1;
cfg_ACI.t_description = 'Time (s)';
cfg_ACI.f = [400 450 500 550 600];
cfg_ACI.f_description = 'Frequency (Hz)';

Data_matrix = permute(noise_spect, [3 1 2]);

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
