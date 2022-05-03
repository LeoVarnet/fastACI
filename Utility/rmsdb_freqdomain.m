function y = rmsdb_freqdomain(f,X)
% function y = rmsdb_freqdomain(f,X)
%
% 1. Description:
%       Root-Mean-Square value of the spectrum X. If X comes from an N-point
%       FFT, then X has K = N/2 elements.
%       N_time corresponds to the length in samples of the time series x used
%       (at some point) as input to the FFT (or DFT).
%
% 2.1 Example 1:
%   [x, Fs] = wavread('Choice.wav'); 
%   rmsdb(x)
% 
% 2.2 Example 2:
%   y = rmsdb('Choice.wav');
% 
% 2.3 Example 3, rms value between 0.1 and 0.2 seconds:
%   [x, fs] = Wavread('Choice.wav'); 
%   ti = 0.1;
%   tf = 0.2;
%   y = rmsdb(x,fs,ti,tf);
%
% Programmed by ExpORL, KU Leuven, Belgium
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 28/10/2014 % Update this date manually
% Last use on   : 31/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r,c]=size(X);

if r == 1
    K = c;
else 
    K = r;
end

N = 2*K;

if nargin < 3
    N_time = N; % by default we assume that an N-length x waveform was used
                % to obtain an N-point FFT
end

df = f(2)-f(1); % frequency resolution

X2 = X.^2;

if c == 1
    y = 10*log10( df*X'*X );
elseif r == 1
    Nf = c;
    y = 10*log10( df*X*X' );
else % Generic case:
    y = 10*log10( df*sum(X.*X) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end