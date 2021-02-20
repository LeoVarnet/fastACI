function bIs_the_same = Still_same_seed_status(seed1)
% function bIs_the_same = Still_same_seed_status(seed1)
%
% 1. Description:
% 
% Programmed by Alejandro Osses, ENS PSL University, Paris, France, 2021-
% Created on     : 20/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bIs_the_same = 1; % starts assuming that they are the same

s_now = rng;

if ~strcmp(seed1.Type,s_now.Type)
    bIs_the_same = 0;
end

if seed1.Seed ~= s_now.Seed
    bIs_the_same = 0;
end

N = length(s_now.State);
if sum(seed1.State == s_now.State) ~= N % all elements of 'State' should be the same
    bIs_the_same = 0;
end
