function [yrms,yrss] = rmsdb(x,fs,ti,tf)
% function [yrms,yrss] = rmsdb(x,fs,ti,tf)
%
% 1. Description:
%       Root-Mean-Square value of x, in dB
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
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2018
% Last update on: 28/10/2014 
% Last use on   : 14/11/2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(x)
    try
        x = Wavread(x);
    catch
        error('variable x interpreted as char, but no wav file with such a name was found')
    end
end

if nargin < 4
    Nf = length(x);
else
    Nf = round(tf*fs);
end

if nargin < 3
    Ni = 1;
else
    Ni = ceil(ti*fs + 1e-6); % to make min idx equal to 1
end

[r,c]=size(x);

if c == 1
    N = length(x(Ni:Nf));
    yrms = 10*log10( x(Ni:Nf)'*x(Ni:Nf)/N );
elseif r == 1
    Nf = c;
    N = length(x(Ni:Nf));
    yrms = 10*log10( x(Ni:Nf)*x(Ni:Nf)'/N );
else % Generic case:
    N = length(x(Ni:Nf,:));
    yrms = 10*log10( sum(x(Ni:Nf,:).*x(Ni:Nf,:))/N );
end
yrss = yrms + 10*log10(N); % Non-normalised sum of squares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
