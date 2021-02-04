function [CI, CIrand, CIboot, ResponseMatrix, CI_targetpresent, CIrand_targetpresent, CIboot_targetpresent, CI_targetabsent, CIrand_targetabsent, CIboot_targetabsent] = computeCI(list_signal,list_response,noise_E,n_rand,n_boot,is_norm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n_trials = length(list_response);
if nargin < 4
    n_rand = 0;
end
if nargin < 5
    n_boot = 0;
end
if nargin < 6
    is_norm = 'yes';
end

Enorm = nan(size(noise_E));
for i_channel = 1:size(noise_E,1)
    if isyes(is_norm)
        allsamples = noise_E(i_channel,:,:);
        Enorm(i_channel,:,:) = (noise_E(i_channel,:,:) - mean(allsamples(:)))/std(allsamples(:));%squeeze(noise_E(i_channel,:,:) - mean(allsamples(:)))/std(allsamples(:));
    else
        Enorm(i_channel,:,:) = (noise_E(i_channel,:,:));%squeeze(noise_E(i_channel,:,:));
    end

end

ResponseMatrix(1,1) = sum(list_signal==1 & list_response==1);
ResponseMatrix(1,2) = sum(list_signal==1 & list_response==2);
ResponseMatrix(2,1) = sum(list_signal==2 & list_response==1);
ResponseMatrix(2,2) = sum(list_signal==2 & list_response==2);

CI_H = mean(Enorm(:,:,list_signal==2 & list_response==2),3);
CI_M = mean(Enorm(:,:,list_signal==2 & list_response==1),3);
CI_CR = mean(Enorm(:,:,list_signal==1 & list_response==1),3);
CI_FA = mean(Enorm(:,:,list_signal==1 & list_response==2),3);
CI = CI_H + CI_FA - CI_M - CI_CR;
CI_targetpresent = CI_H - CI_M;
CI_targetabsent = CI_FA - CI_CR;

CIrand_H = [];
CIrand_M = [];
CIrand_CR = [];
CIrand_FA = [];
CIrand = [];
CIrand_targetpresent = [];
CIrand_targetabsent = [];
CIboot_H = [];
CIboot_M = [];
CIboot_CR = [];
CIboot_FA = [];
CIboot = [];
CIboot_targetpresent = [];
CIboot_targetabsent = [];

if n_rand>0
    for i_rand = 1:n_rand
        rand_response=list_response(randperm(n_trials));
        CIrand_H = mean(Enorm(:,:,list_signal==2 & rand_response==2),3);
        CIrand_M = mean(Enorm(:,:,list_signal==2 & rand_response==1),3);
        CIrand_CR = mean(Enorm(:,:,list_signal==1 & rand_response==1),3);
        CIrand_FA = mean(Enorm(:,:,list_signal==1 & rand_response==2),3);
        CIrand(:,:,i_rand) = CIrand_H + CIrand_FA - CIrand_M - CIrand_CR;%CIrand(i_channel,:,i_rand) + Enorm(i_channel,:,i_trial)*correlator/n_trials;
        CIrand_targetpresent(:,:,i_rand) = CIrand_H - CIrand_M;
        CIrand_targetabsent(:,:,i_rand) = CIrand_FA - CIrand_CR;
    end
end

if n_boot>0
    for i_boot = 1:n_boot
        list_trial = randi(n_trials, 1, n_trials);
        boot_response=list_response(list_trial);
        boot_signal=list_signal(list_trial);
        boot_Enorm=Enorm(:,:,list_trial);
        CIboot_H = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==2),3);
        CIboot_M = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==1),3);
        CIboot_CR = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==1),3);
        CIboot_FA = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==2),3);
        CIboot(:,:,i_boot) = CIboot_H + CIboot_FA - CIboot_M - CIboot_CR;%CIrand(i_channel,:,i_rand) + Enorm(i_channel,:,i_trial)*correlator/n_trials;
        CIboot_targetpresent(:,:,i_boot) = CIboot_H - CIboot_M;
        CIboot_targetabsent(:,:,i_boot) = CIboot_FA - CIboot_CR;
    end
end
end

