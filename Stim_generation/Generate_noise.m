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
        % N_samples = 10432
        N_samples_temp = 15664;
        

        tmod_lim = [-25 25];%[-200 200];%2*[-200 200];%
        smod_lim = [-2 2]/1000;%[-2.5 2.5]/1000;%5*[-3 3]/1000;%
        
        % var = load('mAMPS','refAMPS');
        % refAMPS = var.refAMPS; % = [];
        
        gam = 239766.68634682; 
        fr = 16.3;
        
        N = randn(N_samples_temp,1);%N = noise(N_samples_temp,'pink');%
        [AMPSN,PMPSN,tmodN,smodN] = MPSpec(N,fs);
        
        tmodN_idx = tmodN>tmod_lim(1) & tmodN<tmod_lim(2);
        smodN_idx = smodN>smod_lim(1) & smodN<smod_lim(2);
        filterMPS = smodN_idx'*tmodN_idx;
        
        % generate noise with the filtered AMPS
        insig = MPS_gen(AMPSN.*filterMPS, PMPSN, tmodN, smodN, length(N), fs);%AMPSN.*(refAMPS/max(refAMPS(:)))
        
        % Cut stim to desired length        
        insig = insig(1:N_samples);
end

disp('')

