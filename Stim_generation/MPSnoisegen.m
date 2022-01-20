function [fhat] = MPSnoisegen(Ns, fs, cutoff_t, cutoff_f)
% [fhat] = MPSnoisegen_debug(Ns, fs, cutoff_t, cutoff_f)
% MPSNOISEGEN_DEBUG generates a MPS noise of length Ns sample and sampling
% rate fs which is low-pas in the MPS domain with cutoffs cutoff_t (in Hz)
% in the rate domain and cutoff_f (in cyc/Hz) in the scale domain. 
% This function uses the LTFAT+phaseret instead of stft-istft because the
% latter does not produces reliable results on filtered spectrograms
%
% Author: Leo Varnet - 2021
% Original name (fastACI_sim): MPSnoisegen_debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPhaseret_toolbox = exist('rtpghi.m','file');
if bPhaseret_toolbox == 0
    if ~exist('fastACI_dir_phaseret.m','file')
        fastACI_set_phaseret;
    end
    dir_phaseret = fastACI_dir_phaseret;
    addpath(dir_phaseret);
    phaseretstart;
    
    phaseretmex_unix(dir_phaseret);
    
    error('%s.m: The phaseret toolbox was not found in the path. Please find the toolbox and initialise it (phaseretstart)',upper(mfilename));
end
bPlot = 0;

flow = 0; % Hz, fixed parameter
fhigh = min(10000,fs/2); % Hz, fixed parameter
noise = il_bpnoise(Ns,flow,fhigh,fs);
noise = scaletodbspl(noise,65,100);

if bPlot
    figure; sgram(noise, fs)
end

% built-in parameters for the spectrograms
fr =  16.3*8; %16.3; % frequency resolution (Hz)
hop_rate = 0.002; % 0.001 % hop rate (re: fs) temporal resolution

% other parameters for spectrogram
M = 2*floor((fs/2)/fr)-1; % number of frequency channels ('fr' Hz resolution)
a = floor(hop_rate*fs); % hop size in samps (1 ms)
L = dgtlength(Ns,a,M); % length of gabor system
[gnum,~]=gabwin({'gauss','width',2*fs*(1/(2*pi*fr)),'atheight',0.6065},a,M,L); % ensure frequency resolution (sigmaf) of 'fr' for gaussian window

% % built-in parameters for phase retrieval
gam = 2.414042867677931e+03;%

% compute the log spectrogram
[sabs_carrier, t_TF, f_TF] = il_logspec(noise,fs,fr,gnum,a,M,L);%,DBNOISE);
sabs = sabs_carrier;
% cut spectro to original length (to avoid introducing an artificial component corresponding to the )
tmax = Ns/fs;
size_sabs = size(sabs);
sabs = sabs(:,t_TF<=tmax);
t_TF = t_TF(t_TF<=tmax);

if bPlot % Display the spectrogram
    figure; 
    h=pcolor(t_TF, f_TF, sabs); set(h, 'EdgeColor','none')
    title('Spectrogram');
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    colorbar
end

% compute the MPS
[amp_fabs, phase_fabs, mt, mf] = il_spec2MPS(sabs, t_TF, f_TF);

if bPlot  % Display the amplitude and phase spectrum
    figure
    h=pcolor(mt, mf, 20*log10(amp_fabs)); set(h, 'EdgeColor','none')
    title('Amplitude Spectrum');
    axis xy;
    xlabel('Hz');ylabel('cyc/Hz')
    colormap(jet);colorbar
    
    % figure;
    % h=pcolor(mt, mf, phase_fabs); set(h, 'EdgeColor','none')
    % title('Phase Spectrum');
    % axis xy;
    % cmap = colormap('HSV');
    % cmap(1,:) = [1 1 1];
    % colormap(cmap);
    % colorbar()
end

% Filtering in the MPS domain
mt_idx = abs(mt)<cutoff_t;
mf_idx = abs(mf)<cutoff_f;
filterMPS = mf_idx'*mt_idx;

filtered_amp_fabs = amp_fabs;
filtered_amp_fabs(find(~filterMPS)) = 0;

if bPlot  % Display the filtered amplitude spectrum
    figure
    h=pcolor(mt, mf, 20*log10(filtered_amp_fabs)); set(h, 'EdgeColor','none')
    title('Filtered Amplitude Spectrum');
    axis xy;
    xlabel('Hz');ylabel('cyc/Hz')
    colormap(jet);colorbar
end

filtered_amp_fabs = ifftshift(filtered_amp_fabs);
filtered_phase_fabs = ifftshift(phase_fabs);
filtered_fabs = filtered_amp_fabs.*exp(complex(0,filtered_phase_fabs));

filtered_sabs = il_MPS2spec(filtered_fabs);

if bPlot
    figure; 
    h=pcolor(t_TF, f_TF, filtered_sabs); set(h, 'EdgeColor','none')
    title('Filtered spectrogram');
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    colorbar
end

% Restore silent part of spectrogram and normalize between 0 and 1 to get
% the modulation spectrogram
temp = zeros(size_sabs);
min_filtered_sabs = min(min(filtered_sabs));
max_filtered_sabs = max(max(filtered_sabs));

temp(1:size(filtered_sabs,1),1:size(filtered_sabs,2)) = (filtered_sabs-min_filtered_sabs)/(max_filtered_sabs-min_filtered_sabs);
filtered_sabs = temp;
clear temp

% Invert the spectrogram to get new synthetic sound using phase-gradient
% heap integration (phaseret)

lin_sabs = realpow(10.0, (sabs_carrier)./20.0);% -DBNOISE *maxsabs   % First transform back to linear scale
lin_sabs = filtered_sabs.*lin_sabs;

%gam = pghi_findgamma(gnum,a,M,L);
chatint = rtpghi(lin_sabs,gam,a,M,'freqinv'); % here is the phaseret part
fhat = idgtreal(chatint,{'dual',gnum},a,M,Ns); % inverse discrete gabor transform

% checking the output
if bPlot
    [sabs_fhat, t_TF, f_TF] = il_logspec(fhat,fs,fr,gnum,a,M,L);%,DBNOISE);

    % cut spectro to original length (to avoid introducing an artificial component corresponding to the )
    tmax = Ns/fs;
    sabs_fhat = sabs_fhat(:,t_TF<=tmax);
    t_TF = t_TF(t_TF<=tmax);

    figure; 
    h=pcolor(t_TF, f_TF, sabs_fhat); set(h, 'EdgeColor','none')
    title('Spectrogram');
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    colorbar

    % compute the MPS
    [amp_fabs_fhat, phase_fabs_fhat, mt, mf] = il_spec2MPS(sabs_fhat, t_TF, f_TF);

    figure
    h=pcolor(mt, mf, 20*log10(amp_fabs_fhat)); set(h, 'EdgeColor','none')
    title('Amplitude Spectrum');
    axis xy;
    xlabel('Hz');ylabel('cyc/Hz')
    colormap(jet);colorbar
    
    figure;
    h=pcolor(mt, mf, phase_fabs_fhat); set(h, 'EdgeColor','none')
    title('Phase Spectrum');
    axis xy;
    cmap = colormap('HSV');
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    colorbar()
end
% End of the script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inline functions, used by the main code in this script:
function [sabs, t_TF, f_TF] = il_logspec(noise,fs,fr,gnum,a,M,L)%,DBNOISE)

% get spectrogram
[c,~,gnum] = dgtreal(noise,gnum,a,M,L); % dgt for real-valued signal
sabs = abs(c);
% [c_carrier] = dgtreal(carrier,gnum,a,M,L); % dgt for real-valued signal
% sabs_carrier = abs(c_carrier);

t_TF = (1:size(c,2))*a/fs;
f_TF = (1:size(c,1))/(1/fr);

% convert to log spectrogram
maxsabs = max(max(sabs));
sabsDisp = 20*log10(sabs./maxsabs);%+DBNOISE;
%sabsDisp(sabsDisp<0.0) = 0.0;
sabs = sabsDisp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = il_bpnoise(len,flow,fhigh,fs)
% bpnoise.m - generates spectrally rectangular shaped band-pass noise.
%
% Usage: out = bpnoise(len,flow,fhigh,fs)
%
% len    = output length in samples
% flow   = lower cutoff frequency in Hz
% fhigh  = upper cutoff frequency in Hz
% fs     = sampling rate in Hz
%
% out    = output vector
%
% Function taken from AFC toolbox, programmed by Stephan Ewert and colleagues

out = real(ifft(il_scut(fft(randn(len,1)),flow,fhigh,fs)));

out = out/(norm(out,2)/sqrt(len));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cut = il_scut(in,flow,fhigh,fs)
% function cut = il_scut(in,flow,fhigh,fs)
%
% Function taken from AFC toolbox, programmed by Stephan Ewert and colleagues

len = length(in);
flow = round(flow*len/fs);
fhigh = round(fhigh*len/fs);
cut = zeros(len,1);
cut(flow+1:fhigh+1) = in(flow+1:fhigh+1);

% HACK: if lowpass ( flow = 0) index would be greater than len (len +1)
if flow == 0
	flow = 1;
end

cut(len-fhigh+1:len-flow+1) = in(len-fhigh+1:len-flow+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [amp_fabs, phase_fabs, mt, mf] = il_spec2MPS(sabs, t_TF, f_TF)
% calculate the 2D fft
fabs = fftshift(fft2(sabs));

% calculate amplitude and phase
amp_fabs = abs(fabs);
phase_fabs = angle(fabs);

% Calculate the labels for temporal and spectral frequencies in physical units
mt = linspace(-(1/t_TF(1))/2, (1/t_TF(1))/2, length(t_TF));
mf = linspace(-(1/f_TF(1))/2, (1/f_TF(1))/2, length(f_TF));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_sabs = il_MPS2spec(new_fabs)

% recover filtered spectrogram by inverse 2D FFT
new_sabs = real(ifft2(new_fabs));

% The amplitude must also stay positive
%new_sabs(new_sabs < 0.0) = 0.0;