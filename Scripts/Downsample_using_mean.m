function outsig = Downsample_using_mean(insig,binwidth_samples)
% function outsig = Downsample_using_mean(insig,binwidth_samples)
%
% Author: Alejandro Osses
% See also: Gammatone_proc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_blocks = size(insig,1)/binwidth_samples;

for i = 1:N_blocks
    idxi = binwidth_samples*(i-1)+1;
    idxf = binwidth_samples*i;
    
    outsig(i,:) = mean(insig(idxi:idxf,:));
end