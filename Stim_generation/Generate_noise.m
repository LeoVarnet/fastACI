function [insig] = Generate_noise(N_samples,noise_type)
% function [insig] = Generate_noise(N_samples,noise_type)
%
% 1. Description:
%       Noise generation, so far only validated for white and pink noise.
%       The level adjustment and sampling frequency control occurs outside
%       this script.
% 2. Additional info:
%       See also speechACI_Logatome_init.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(noise_type)
case 'white'
    insig = randn(N_samples,1); % This N_samples already includes the ramp times
case 'pink'
    insig = noise(N_samples,'pink'); % function from LTFAT toolbox
end

disp('')

