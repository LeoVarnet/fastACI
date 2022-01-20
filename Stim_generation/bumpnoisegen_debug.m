function [S, opts] = bumpnoisegen_debug(Ns, fs, Nb, sigma_t, sigma_ERB, A, lvl, dBFS,loc_t,loc_f,version)
%BUMPNOISEGEN_DEBUG Summary of this function goes here
%   Generates a bump noise of length Ns sample and sampling rate fs. The
%   bump noise is composed of Nb gaussian bumps randomly placed in the T-F
%   domain (with an uniform distribution in the ERB scale). Each bump has a
%   temporal width sigma_t, a spectral width sigma_ERB and an amplitude A
%   (in dB). 
% ref: Varnet, Langlet, Lorenzi, Lazard, Micheyl, "High-Frequency
% Sensorineural Hearing Loss Alters Cue-Weighting Strategies for
% Discriminating Stop Consonants in Noise", Trends in Hearing (2019)
%
% This script uses the stft-istft toolbox instead of LTFAT+phaseret
% because the latter uses mex files which do not compile on Mac.
%
% Léo Varnet - 2021
% Versions:
% 1: As adapted by Alejandro and used in July 2021 for a pilot. This version
%    was used in g20210915_characterising_bump_speech.m
% 1.1: - latest bump at dur-2*sigma_t
%      - one bump per ERB
%      - lmax = lvl + 15 remopved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 11
    version = 1;
end

if nargin < 7
    lvl = 50;
    warning('Default lvl of 50 dB will be used')
end
if nargin < 8
    dBFS = 100; % dBFS convention, a waveform amplitude of 1 will equal this value
    warning('Default full-scale convention of 100 dBFS will be used')
end

if fs ~= 16000
    error('Bump noise so far validated for a sampling frequency of 16 kHz')
end

switch version
    case 1
        offset_t = 0;
    case {1.1, 1.2, 1.3}
        offset_t = 50e-3;
end

if offset_t > 0
    Ns_extra = round(offset_t*fs);
    Ns = Ns + 2*Ns_extra;
end

% built-in parameters
nF = 1024; % floor(0.064*fs); % 512; % number of points in the FFT
nW = gausswin(nF) ; % window
nH = nF/8; % 512/8; % offset between successive windows (in samples)
n_it = 5; % number of iterations for the phase reconstruction algorithm.

t = (0:nH:Ns-nF)/fs;
f = linspace(0,(fs/2),nF/2+1);

%[max(t)/2, max(t)*20, 5000/2, max(f)*8, +1;...
% locations and scale of the blobs (t, sigma_t, f, sigma_f, sign)
[Tgrid,Fgrid] = meshgrid(t,f);

BdB = zeros(size(Tgrid));
%A = 1-exp(-(Tgrid).^2/0.05^2)-exp(-(Tgrid-max(max(Tgrid))).^2/0.05^2);

%%% Random locations of the bumps:
% In time:
dur = Ns/fs;
switch version
    case 1
        time_first_possible_bump = 0;
        time_last_possible_bump = dur;
        
    case {1.1, 1.2, 1.3}
        time_first_possible_bump = 0+.5*offset_t;
        time_last_possible_bump = dur-2*sigma_t-offset_t;
end
per_time = time_last_possible_bump-time_first_possible_bump; % period of time

Nb_old = Nb;
switch version 
    case 1
        if nargin < 10
            % In frequency:
            loc_f = ERB2f( rand(Nb,1)*f2ERB(fs/2) ); % in ERB_N
        end
    case 1.1
        loc_f = audtofreq(freqtoaud(80):1:freqtoaud(fs/2)-sigma_ERB*2);
        Nb = length(loc_f);
        loc_f = loc_f(:);
        
    case {1.2, 1.3}
        fmin = freqtoaud(80);
        fmax = freqtoaud(fs/2)-sigma_ERB*2;
        loc_f_idle = audtofreq(fmin:1:fmax);
        Nb = length(loc_f_idle);
        loc_f = audtofreq( rand(Nb,1)*(fmax-fmin)+fmin ); 
        loc_f = loc_f(:);
end

if nargin < 9
    loc_t = rand(Nb,1)*time_last_possible_bump;
end

if length(loc_t) ~= length(loc_f)
    loc_t = time_first_possible_bump + rand(Nb,1)*per_time;
end

if max(loc_t) > dur
    error('loc_t at a time later than the duration of the noise...')
end

% 'BW' of the bubble, sigma_ERB/2 above and below the random location
sigma_f = ERB2f(f2ERB(loc_f)+sigma_ERB/2)-ERB2f(f2ERB(loc_f)-sigma_ERB/2);

for i=1:Nb
    % Summing all the Gaussian-distributed bumps:
    BdB = BdB + A * exp(-(Tgrid-loc_t(i)).^2/sigma_t^2-(Fgrid-loc_f(i)).^2/sigma_f(i)^2);
end

switch version
    case 1
        idx = find(BdB>A);
        BdB(idx) = A;
end
% BdB(BdB>A) = A; % In case of overlapping bumps, ensure that the max height will not be higher than A 
B = 10.^(BdB/20);

%%% Creating the white noise that will be 'bumped':
% noise = randn(Ns,1);
% noise = scaletodbspl(noise,lvl,dBFS);

flow = 0; % Hz, fixed parameter
switch version
    case 1
        fhigh = 10000; % Hz, fixed parameter
    case {1.1, 1.2, 1.3}
        fhigh = min(10000,fs/2); % Hz, fixed parameter
end
noise = il_bpnoise(Ns,flow,fhigh,fs);
% lvl = lvl-10*log10(10000)+10*log10(fhigh-flow);
noise = scaletodbspl(noise,lvl,dBFS);

%%%
noiseSpec = il_stft(noise, nF, nW, nH);
noiseSpec = noiseSpec.*B; % introducing the bumps

% adding a last time sample to the spectrogram to ensure that the
% reconstructed noise will have the desired length
noiseSpec = [noiseSpec noiseSpec(:,end)];

S = il_phaseRecon(noiseSpec, noiseSpec, n_it, nF, nH); % randomises the phase after adding the bumps
S = S(1:Ns);
S = S(:); % Converting the waveform back to a column array (istft produces
          % waveforms as raw arrays...

switch version
    case 1
        lvl_max = lvl+15;
        if rmsdb(S)+dBFS > lvl_max
            S = scaletodbspl(S,lvl_max,dBFS);
            fprintf('\tLimited maximum level...');
        end
end

if Ns_extra ~= 0
    S = S(Ns_extra:end-Ns_extra-1);
    
    ramp_dur_ms = 20; %
    S = Do_cos_ramp(S,fs,ramp_dur_ms,ramp_dur_ms);
    
    loc_t = loc_t - offset_t;
end

if nargout >= 2
    opts = [];
    opts.BdB = BdB;
    opts.loc_t = loc_t; % location(s) of the bubbles in the time domain
    opts.loc_f = loc_f; % location(s) of the bubbles in the frequency domain
    opts.f_t   = fs/nH; % sampling frequency of the time grid
    opts.t     = t;
    opts.f     = f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = il_stft(x, f, w, h)
% D = stft(X, F, W, H)                            Short-time Fourier transform.
%	Returns some frames of short-term Fourier transform of x.  Each 
%	column of the result is one F-point fft; each successive frame is 
%	offset by H points until X is exhausted.  Data is hamm-windowed 
%	at W pts, or rectangular if W=0, or with W if it is a vector.
%	See also 'istft.m'.
% dpwe 1994may05.  Uses built-in 'fft'
% $Header: /homes/dpwe/public_html/resources/matlab/pvoc/RCS/stft.m,v 1.2 2009/01/07 04:32:42 dpwe Exp $

s = length(x);

if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,f);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = f/2;   % midpoint of win
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, f);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    %win = win(:);
  end
else
  win = w;
  w = length(w);
end

c = 1;

% pre-allocate output array
d = zeros((1+f/2),1+fix((s-f)/h));

for b = 0:h:(s-f)
  u = win.*x((b+1):(b+f));
  t = fft(u);
  d(:,c) = t(1:(1+f/2))';
  c = c+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,spec] = il_phaseRecon(targetMag, initialPhase, nIter, nFft, hop)

% Reconstruct a signal from a target magnitude spectrogram and initial
% phase.  Does several reconstructions, keeping the phase from the last
% reconstruction with the target magnitude.

targetMag = abs(targetMag);
phase = initialPhase ./ abs(initialPhase);

for i = 1:nIter
    x = il_istft(targetMag .* phase, nFft, nFft, hop);
    spec = il_stft(x, nFft, nFft, hop);
    phase = spec ./ abs(spec); % real and imaginary part of the phase are constrained to -1 to 1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = il_istft(d, ftsize, w, h)
% X = istft(D, F, W, H)                   Inverse short-time Fourier transform.
%	Performs overlap-add resynthesis from the short-time Fourier transform 
%	data in D.  Each column of D is taken as the result of an F-point 
%	fft; each successive frame was offset by H points. Data is 
%	hamm-windowed at W pts.  
%       W = 0 gives a rectangular window; W as a vector uses that as window.
% dpwe 1994may24.  Uses built-in 'ifft' etc.
% $Header: /homes/dpwe/public_html/resources/matlab/pvoc/RCS/istft.m,v 1.4 2009/01/07 04:20:00 dpwe Exp $

s = size(d);
 
cols = s(2);
xlen = ftsize + (cols-1)*h;
x = zeros(1,xlen);

if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,ftsize);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = ftsize/2;
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, ftsize);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    % 2009-01-06: Make stft-istft loop be identity
    win = 2/3*win;
  end
else
  win = w;
  w = length(win);
end
   
for b = 0:h:(h*(cols-1))
  ft = d(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x((b+1):(b+ftsize)) = x((b+1):(b+ftsize))+px.*win;
end;

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

% eof
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