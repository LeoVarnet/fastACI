function [undersmplE] = noisetone_converter(dir_where,foldername,ListStim,n_stim,fcut,undersampling,fEcut,lvl,SNR)

if nargin<7
    fEcut=20;
end
Nchannel = length(fcut)-1;

% ATTENTION BRICOLAGE %
[tone, fs] = audioread([dir_where 'TargetStims' filesep 'nontarget.wav']);

for i_trial=1:length(n_stim)
    fprintf(['stim #' num2str(i_trial) '\n'])
    [noise, fs] = audioread([dir_where foldername filesep ListStim(n_stim(i_trial)).name]);
    
% ATTENTION BRICOLAGE %
    S = generate_stim( tone, noise, SNR, 0, 'white');
    S = dBlvl(S,lvl);
    
    %S=S/rms(S);
    for i_channel = 1:Nchannel
         [B, A] = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
         % %option hilbert envelope
         %E = abs(hilbert(filtfilt(B, A, S)));
         %option filtering
          [B2, A2] = butter(1, fEcut/(fs/2));
          E = filterfilter(B2, A2, abs(filterfilter(B, A, S)));
          % %
         undersmplE(i_channel, :, i_trial) = E(1:undersampling:end);
    end
end
