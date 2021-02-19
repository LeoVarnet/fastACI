function S_envelope_undersampled = noisetone_converter(dir_where,foldername,ListStim,n_stim,fcut,undersampling,fEcut,lvl,SNR)
% function S_envelope_undersampled = noisetone_converter(dir_where,foldername,ListStim,n_stim,fcut,undersampling,fEcut,lvl,SNR)

if nargin<7
    fEcut=20;
end
Nchannel = length(fcut)-1;

% 1. Reads the stored pure tone:
[tone, fs] = audioread([dir_where 'TargetStims' filesep 'nontarget.wav']);

for i_trial=1:length(n_stim)
    fprintf(['stim #' num2str(i_trial) '\n'])
    % 2. Reads each stored noise:
    [noise, fs] = audioread([dir_where foldername filesep ListStim(n_stim(i_trial)).name]);
        
    % 3. Generates a noise+tone at a given SNR, no ramps
    S = generate_stim( tone, noise, SNR, 0, 'white');
    
    % 4. Sets the level of the interval (T+N):
    S = dBlvl(S,lvl);
    
    if Nchannel > 1
        error('%s, AO: Continue validating here...',upper(mfilename));
    end
    for i_channel = 1:Nchannel
        [B, A] = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
        S_bandpass = filterfilter(B, A, S); % bandpassed signal
        
        [B, A] = butter(1, fEcut/(fs/2));
        S_envelope = filterfilter(B, A, abs(S_bandpass));
        
        S_envelope_undersampled(i_channel, :, i_trial) = S_envelope(1:undersampling:end);
    end
end

if nargout == 0
    figure; 
    plot(S); hold on;
    plot(S_bandpass,'m--')
    plot(S_envelope,'r--','LineWidth',2)
    
    xlabel('Time (Samples)')
    ylabel('Amplitude')
    
    legend('Broadband signal','Bandpassed signal','Envelope of the bandpassed signal');
end


%     [S, fs]   = audioread([dir_where foldername filesep ListStim(n_stim(i_trial)).name ]);  
%     S=S/rms(S);

%     for i_channel = 1:Nchannel
%          [B, A] = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
%          % %option hilbert envelope
%          %E = abs(hilbert(filtfilt(B, A, S)));
%          % %option filtering
%           [B2, A2] = butter(1, fEcut/(fs/2));
%           E = filterfilter(B2, A2, abs(filterfilter(B, A, S)));
%           % %
%          undersmplE(i_channel, :, i_trial) = E(1:undersampling:end);
%     end