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
        % temporary built-in parameters
        var = load('mAMPS','mAMPS');
        refAMPS = var.mAMPS; % = [];
        
        gam = 239766.68634682; 
        fr = 16.3;
        
        tempnoise = randn(N_samples,1);%noise(N_samples,'pink'); % function from LTFAT toolbox
        [~,PMPSN,tmodN,smodN] = MPSpec(tempnoise,fs);
        if any(size(PMPSN)~=size(refAMPS))
            error(['Incompatible MPS sizes. Please check that your reference MPS has been generated from sounds with sampling frequency = ' num2str(fs) ' and number of samples = ' num2str(N_samples) '.\n'])
        end
        insig = MPS_gen(refAMPS, PMPSN, tmodN, smodN, N_samples, fs, gam, fr);
end

disp('')

