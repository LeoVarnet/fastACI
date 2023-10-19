function [AmplifiedS, gain, CF, percentClipped] = AmplifyNALR(S, fs, HL2, HLF, doplot)
% [AmplifiedS, gain, CF, percentClipped] = AmplifyNALR(S, fs, HL2, HLF, doplot)
% Applies NAL-R amplification to the signal
%
% Input:
% - S: input sound
% - fs: sampling freq
% - HL: hearing thresholds (dB HL)
% - HLF: frequencies (Hz) at which hearing thresholds are given 
% - doplot: 1 to plot the figures
%
% Output:
% - AmplifiedS: amplified sound
% - gain: gains applied to the signal (at frequencies CF)
% - CF: center frequencies of the processing channels
%
% Author: C. Micheyl, Hearing Sciences Dept, Starkey France
% Copyright: Starkey Hearing Technologies
% Date: 5 Dec. 2017
% Modified by Leo Varnet Spet. 2023

if nargin < 3
error('AmplifyNALR_Leo requires hearing thresholds HL2')
end
if nargin < 4 
    HLF = [250 500 750 1000 1500 2000 4000 6000 8000]; 
end
if nargin < 5
    doplot = 0;
end

% PTA
HL = mean(HL2); % Mean across the left and right ears
HLat3F = interp1(HLF, HL, [500 1000 2000]); % HL at 500, 1000, and 2000 Hz
PTA = mean(HLat3F); % PTA over 500, 1000, and 2000 Hz (dB HL)

% Filterbank for signal processing
% Frequency-band cutoff frequencies (Hz)
hiloF = [  50 [1:15]*500;    
    [1:15]*500 7950]; % Leo: I changed this so that the amplification applies up to 8000Hz
nBands = size(hiloF, 2);
filterOrder = 512; % Desired filter order (including filtfilt)
for bandN = 1:nBands
    B(bandN,:) = fir1(filterOrder/2, [hiloF(1,bandN) hiloF(2,bandN)]/(fs/2)); % Divide filter order by 2 to account for filtfilt
end

if doplot
    for bandN = 1:nBands
        [H,F] = freqz(B(bandN,:), 1, 512, fs);
        figure(1); semilogx(F, 20*log10(abs(H))); hold on; title('processing filters'); xlabel('frequency [Hz]'); ylabel('gain [dB]')
        figure(2); plot(F, 20*log10(abs(H))); hold on; title('processing filters'); xlabel('frequency [Hz]'); ylabel('gain [dB]')
    end
    figure(1); hold off;
    figure(2); hold off;
end

% NAL-R gain calculation

% IGdB = 0.05 * PTAdB + 0.31 * dBHL + C
% Source: 
% Sandlin, RE (2000) Textbook of hearing aid amplification. Singular Publishing. p.75
% Palmer CV, Lindley GA Overview and Rationale for Prescriptive Formulas
% for Linear and Nonlinear Hearing Aids.
% Rajkumar S, Muttan S, Jaya V, Vignesh SS (2013) Comparative  Analysis of Different Prescriptive Formulae Used in the Evaluation of Real Ear Insertion Gain 
% for Digital Hearing Aid. Universal Journal of Biomedical Engineering
% 1(2): 32-41 (note: this paper contains an error in a coefficient, 0.15 instead of 0.05!)
NALCF = [  0  250 500 750 1000 1500 2000 3000 4000 6000 8000]; % Frequencies for correction values (Hz)
NALC  = [-17  -17  -8  -3    1   1    -1   -2   -2   -2   -2]; % Correction value (dB)
% Note: -17 value at 0 Hz added by Leo (original value by Christophe: -30)
% IGdB = 0.05 * PTAdB + 0.31 * dBHL + C
CF = geomean(hiloF);
gain = 0.05* PTA + 0.31 * interp1([0 HLF], [HL(1) HL], CF) + interp1(NALCF, NALC, CF); % 'Insertion gain' per band

nAudioChannels = size(S, 2);
T =[1:length(S)]/fs;

% Amplification
for audioChannelN = 1:nAudioChannels 
    signal = S(:,audioChannelN ); % Mono signal
    
    for bandN = 1:nBands
        signalBand(bandN,:) = filtfilt(B(bandN,:), 1, signal);
    end
    
    AmplifiedS(:, audioChannelN) = 10.^(gain/20) * signalBand ;
    max(abs(AmplifiedS(:, audioChannelN)));
    percentClipped = length(find(abs(AmplifiedS(:, audioChannelN)) > 1))/length(AmplifiedS(:, audioChannelN))*100;
end

if doplot
    % plot original and amplified signals
   figure(3); plot(T,AmplifiedS(:,1),'r',T,S(:,1),'k'); title('signal before/after'); xlabel('time [sec]'); ylabel('amplitude')
   Sfft = fft(S(:,1));
   AmplifiedSfft = fft(AmplifiedS(:,1));
   f = ([1:length(Sfft)]/length(Sfft))*fs;
   figure(4); plot(f, 10*log10(abs(AmplifiedSfft)), 'r', f, 10*log10(abs(Sfft)), 'k'); xlim([0, fs/2]); title('signal before/after'); xlabel('frequency [Hz]'); ylabel('amplitude [dB]')
   figure(5); plot(f, 20*log10(abs(AmplifiedSfft)./abs(Sfft))); xlim([0, fs/2]); title('actual gain'); xlabel('frequency [Hz]'); ylabel('gain [dB]')
end