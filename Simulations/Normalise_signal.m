function [insig_norm,c,idx] = Normalise_signal(insig,fs,method)
% function [insig_norm,c,idx] = Normalise_signal(insig,fs,method)
%
% 1. Description:
%       Normalisation of a signal x. If x has more than 1 column, each 
%       column is interpreted as a different time series.
%       This normalisation is approximately independent respect to the 
%       sampling frequency.
% 
%       This script makes the multiplication sum(y^2) equal to (N/fs), where
%       N/fs is the duration of the digital signal in seconds.
% 
%       If x (at fs) has a duration of 3.65 [s] then sum(y.*y) will be 3.65
%       For the same signal x, optimaldetector(y4optdet,y4optdet) = 1
% 
% 2. Stand-alone example:
%       x = wgn(1,10,1); % 10 element white noise
%       fs = 44100;
%       y = Normalise_signal(x,fs);
% 
%       fs = 44100;
%       x = wgn(1,3.65*fs,1); % 3.65-s length white noise
%       method = 0;
%       y = Normalise_signal(x,fs,method);
%       disp(num2str(sum(y.*y))) % this should give 3.65 of total energy
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also optimaldetector.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 20/10/2014
% Last update on: 27/12/2016 
% Last use on   : 15/05/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 2;
end

idx = find(isnan(insig));
if length(idx)~=0
    insig(idx) = [];
    warning('NaN elements were found and deleted');
end

insig_norm = nan(size(insig)); % just allocation

N = length(insig);
if size(insig,2) ~= 1
    error('Have a look at this...')
    for i = 1:size(insig,2)
        
        switch method
            case 0
                c = sqrt( N/(fs* sum(insig(:,i).^2)) );
            case 1
                c = sqrt( 1/sum(insig(:,i).^2) );
            case 2
                c = sqrt(fs/sum(insig(:,i).^2) );
            case 3
                c = fs/sum(insig(:,i));
        end
        insig_norm(:,i) = insig(:,i) * c;
    end
else
    
    switch method
        case 0
            c = sqrt( N/(fs * sum(insig.^2)));
        case 1
            c = sqrt( 1/sum(insig.*insig) );
        case 2
            c = sqrt(fs/sum(insig.*insig) );
        case 3
            c = fs/sum(insig);
        case 4 % AMT default for templates
            c = sqrt( length(insig)/sum(insig.*insig) ); % 1/'rms'
    end
    insig_norm = c*insig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
