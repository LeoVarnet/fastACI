function [undersmplE] = noise_converter(dir_where,foldername,ListStim,n_stim,fcut,undersampling,fEcut)

if nargin<7
    fEcut=20;
end
Nchannel = length(fcut)-1;

for i_trial=1:length(n_stim)
    fprintf(['stim #' num2str(i_trial) '\n'])
    [S, fs] = audioread([dir_where foldername filesep ListStim(n_stim(i_trial)).name ]);
    
    S=S/rms(S);
    for i_channel = 1:Nchannel
         [B, A] = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
         % %option hilbert envelope
         %E = abs(hilbert(filtfilt(B, A, S)));
         % %option filtering
          [B2, A2] = butter(1, fEcut/(fs/2));
          E = filterfilter(B2, A2, abs(filterfilter(B, A, S)));
          % %
         undersmplE(i_channel, :, i_trial) = E(1:undersampling:end);
    end
end