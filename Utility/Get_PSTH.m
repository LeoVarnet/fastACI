function [psth_avg,t_psth] = Get_PSTH(insig,fs,psth_binwidth)
% function psth_avg = Get_PSTH(insig,fs,psth_binwidth)
%
% Author: Alejandro Osses
% Date: 4/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psthbins = round(psth_binwidth*fs);  % number of psth_ft bins per psth bin
N_samples= size(insig,1);
% N_bins   = length(insig)/psthbins;

t = (1:N_samples)/fs;

count_i = 1;
for i = 1:psthbins:N_samples
    idxi = i;
    idxf = min(i+psthbins-1,N_samples);
    psth_avg(count_i,:) = sum(insig(idxi:idxf,:));
    t_psth(count_i) = t(idxi);
    count_i = count_i + 1;
end

% psth_avg = psth_avg/psth_binwidth;
psth_avg = mean(psth_avg,2)/psth_binwidth;
