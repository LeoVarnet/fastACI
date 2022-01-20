function [insig] = Generate_noise(N_samples,noise_type,fs)
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
    case 'smps'
        error('Discarded noise algorithm. See previous versions of this script (%s) to still retrieve this algorithm',mfilename);
    case 'bumpv1p2_10db'
        % Needs creation
        sigma_t = 0.02; % temporal width of the bumps (in s)
        sigma_ERB = 0.5; % spectral width of the bumps (in ERB)
        
        A = 10;
        lvl_target = 65;
        % 10 is  6.4 dB (- 3.6 dB)
        % 20 is  9   dB (-11   dB)
        % 30 is 14.3 dB (-15.7 dB)
        switch A
            case 5
                lvl = lvl_target;
            case 10
                lvl = lvl_target - 3.6;
            case 20
                lvl = lvl_target - 9;
            case 30
                lvl = lvl_target - 11.3;
            otherwise
                error('No level empirically obtained...')
        end

        dBFS = 100;
        version = 1.2;
        insig = bumpnoisegen(N_samples, fs, [], sigma_t, sigma_ERB, A, lvl, dBFS,[],[],version);
            
    case 'smpsv1p3'
        cutoff_t = 35;
        cutoff_f = 10/1000;
        
        insig = MPSnoisegen_debug(N_samples, fs, cutoff_t, cutoff_f);
end

disp('')