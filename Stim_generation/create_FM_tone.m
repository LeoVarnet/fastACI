function [outsig, fs] = create_FM_tone( fc, FM_fm, fdev, duration, fs, A)
% [outsig, fs ] = CREATE_FM_TONE( fc, FM_fm, fdev, duration, fs, A)
%
% fc : carrier frequency (Hz)
% FM_fm : rate of frequency modulation (Hz)
% fdev     : frequency deviation
% duration : time (in s) (default: 1)
% fs       : sampling frequency (Hz) (default: 44100)
% A        : amplitude of the signal (default = 1)
% FM_waveform : 'sin' (default), 'triangle', 'square'
% carrier_waveform : 'sin' (default), 'triangle', 'square'
%
% The equation s(t)=1*sin[2*pi*fc*t+beta*sin(2*pi*FM_fm*t)] describes the
% stimulus.
% 
% Author: Leo Varnet 2016, 
% Author: Alejandro Osses 2022, modifications

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

function_carrier = @sin;
function_FM      = @sin;

% generate t
t = 1/fs:(1/fs):duration;

factor = fdev/FM_fm;
% stim
stim = function_carrier(2*pi*fc*t+factor*function_FM(2*pi*FM_fm*t));

outsig = A*stim;