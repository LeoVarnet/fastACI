function [ AM, fs ] = create_AM( fc, AM_fm, AM_depth, duration, fs, A, carrier_waveform, AM_waveform)
%[ AM, fs ] = CREATE_AM( fc, AM_fm, AM_depth, duration, fs, A, carrier_waveform, AM_waveform )
% fc: carrier frequency (Hz)
% AM_fm: rate of amplitude modulation (Hz)
% AM_depth: m (between 0 and 1)
% duration: time (in s) (default: 1)
% fs: sampling frequency (Hz) (default: 44100)
% A: amplitude of the signal (default = 1)
% AM_waveform: 'sin' (default), 'triangle', 'square'
% carrier_waveform: 'sin' (default), 'triangle', 'square', 'noise' (in the
% latter case, white noise will be used unless you specify fc = [fmin, fmax])
%
% The stimulus is described by the following equation:
% s(t)=[1+AM_depth*sin(2*pi*AM_fm*t+3*pi/2)]*sin(2*pi*fc*t)
% Raised-cosine ramps (20 ms) are applied at the onset and at the offset. 
%
% Leo Varnet 2016 - last modified 2018

%defaults
if nargin<=7
    AM_waveform = 'sin';
end
if nargin<=6
    carrier_waveform = 'sin';
end
if nargin<=5
    A = 1;
end
if nargin<=4
    fs = 44100;
end
if nargin<=3
    duration = 1;
end
if nargin<2
    error('Not enough input arguments');
end

if strcmp(carrier_waveform, 'sin')
    function_carrier = @sin;
elseif strcmp(carrier_waveform, 'triangle')
    function_carrier = @triangle;
elseif strcmp(carrier_waveform, 'square')
    function_carrier = @square;
end
if strcmp(AM_waveform, 'sin')
    function_AM = @sin;
elseif strcmp(AM_waveform, 'triangle')
    function_AM = @triangle;
elseif strcmp(AM_waveform, 'square')
    function_AM = @square;
end

% generate t
t = 1/fs:(1/fs):duration;

% stim
if ~strcmp(carrier_waveform, 'noise')
    stim = (1+AM_depth*function_AM(2*pi*AM_fm*t+3*pi/2)).*function_carrier(2*pi*fc*t);
else
    stim = (1+AM_depth*function_AM(2*pi*AM_fm*t+3*pi/2)).*randn(size(t));
    % bandpass noise
    if length(fc)==2
        [Bfilt, Afilt] = butter(2,fc/(fs/2));
        stim=filtfilt(Bfilt, Afilt, stim);
    end
end

% add ramps
ramp_duration = 0.02;%(1/AM_fm)/2;%duration*0.05;
t_ramp = 0:(1/fs):ramp_duration;
ramp_up = cos((2*pi*t_ramp(end:-1:1))/(2*max(t_ramp)))/2+0.5;
ramp_down = cos((2*pi*t_ramp)/(2*max(t_ramp)))/2+0.5;
weighting = ones(1,length(t));
weighting(1:length(t_ramp)) = ramp_up;
weighting(end-length(t_ramp)+1:end) = ramp_down;

AM = A*stim.*weighting;

end
