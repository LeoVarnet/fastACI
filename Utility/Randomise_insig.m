function outsig = Randomise_insig(insig,N)
% function outsig = Randomise_insig(insig,N)
%
% 1. Description:
%       outsig containes the same samples (same signal length) than insig1 
%       but with a starting sample that is randomly defined (following a 
%       uniform distribution). 
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/08/2015
% Last update on: 21/12/2016 
% Last use on   : 21/12/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = round( length(insig)*random('unif',0,1) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];

if nargin >= 2
    outsig = outsig(1:N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
