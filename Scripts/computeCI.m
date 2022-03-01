function [CI, CIrand, CIboot, ResponseMatrix, CI_targetpresent, CIrand_targetpresent, CIboot_targetpresent, CI_targetabsent, CIrand_targetabsent, CIboot_targetabsent] = computeCI(list_signal,list_response,noise_E,n_rand,n_boot,is_norm)
% function [CI, CIrand, CIboot, ResponseMatrix, CI_targetpresent, CIrand_targetpresent, CIboot_targetpresent, CI_targetabsent, CIrand_targetabsent, CIboot_targetabsent] = ...
%                          computeCI(list_signal,list_response,noise_E,n_rand,n_boot,is_norm)
%
% Summary of this function goes here
%
%   Input parameters:
%       list_signal: list of size (1 x N_trials), if 1 then the signal was a 'reference'
%           if 2 the signal was the 'target'
%       list_response: list of size (1 x N_trials) with the participants' responses
%           1 or 2 if the participant indicated the sound to be 'reference' or 'target'.
%       noise_E: Envelope of the noises (dim: N_channels x N_samples x N_noises)
%       n_rand,
%       n_boot  (default = 0, i.e., bootstrapping deactivated)
%       is_norm (default = 1, i.e., normalisation activated)
%
%   Output parameters:
%       CI: Classification image, in agreement with Eq. (1) of Murray (2011)
%       
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_trials = length(list_response);

resp_if_tar = 2; % 'if target' (AM)       - 2 modulated tone
resp_if_ref = 1; % 'if reference' (no AM) - 1 pure tone

if nargin < 4
    n_rand = 0;
end
if nargin < 5
    n_boot = 0;
end
if nargin < 6
    is_norm = 1; 
end

Enorm = nan(size(noise_E)); % memory allocation
for i_channel = 1:size(noise_E,1)
    if is_norm == 1
        error('Continue validating here...')
        allsamples = noise_E(i_channel,:,:);
        Enorm(i_channel,:,:) = (noise_E(i_channel,:,:) - mean(allsamples(:)))/std(allsamples(:));%squeeze(noise_E(i_channel,:,:) - mean(allsamples(:)))/std(allsamples(:));
    else
        Enorm(i_channel,:,:) = (noise_E(i_channel,:,:));
    end
end
 
idxs_H  = find(list_signal==resp_if_tar & list_response==resp_if_tar);
idxs_M  = find(list_signal==resp_if_tar & list_response==resp_if_ref); 
idxs_CR = find(list_signal==resp_if_ref & list_response==resp_if_ref);
idxs_FA = find(list_signal==resp_if_ref & list_response==resp_if_tar);

ResponseMatrix(1,1) = length(idxs_CR); % Correct rejection
ResponseMatrix(1,2) = length(idxs_FA); % False alarm
ResponseMatrix(2,1) = length(idxs_M); % Miss
ResponseMatrix(2,2) = length(idxs_H); % Hit rate

%%% Before (Leo's code), it should be the same:
% ResponseMatrix(1,1) = sum(list_signal==1 & list_response==1);
% ResponseMatrix(1,2) = sum(list_signal==1 & list_response==2);
% ResponseMatrix(2,1) = sum(list_signal==2 & list_response==1);
% ResponseMatrix(2,2) = sum(list_signal==2 & list_response==2);

% Classification Image (CI) assessed from the mean of the envelope (Enorm) 
%     across all stimuli (dim=3). One value for each time sample is obtained:
dim_for_mean = 3;
if size(Enorm,dim_for_mean) ~= n_trials
    error('Dimension %s is not the dimension containing the different trial sounds, change manually this in the script.',dim_for_mean);
end

CI_H  = mean(Enorm(:,:,idxs_H) ,dim_for_mean);
CI_M  = mean(Enorm(:,:,idxs_M) ,dim_for_mean);
CI_CR = mean(Enorm(:,:,idxs_CR),dim_for_mean);
CI_FA = mean(Enorm(:,:,idxs_FA),dim_for_mean);

CI_targetpresent = CI_H  - CI_M;
CI_targetabsent  = CI_FA - CI_CR;

%%% Before (Leo's code):
% CI_H = mean(Enorm(:,:,list_signal==2 & list_response==2),3);
% CI_M = mean(Enorm(:,:,list_signal==2 & list_response==1),3);
% CI_CR = mean(Enorm(:,:,list_signal==1 & list_response==1),3);
% CI_FA = mean(Enorm(:,:,list_signal==1 & list_response==2),3);
% CI = 0.5 * (CI_H + CI_FA - CI_M - CI_CR);
% CI_targetpresent = CI_H - CI_M;
% CI_targetabsent = CI_FA - CI_CR;

CI = (CI_H + CI_FA) - (CI_M + CI_CR); % see Eq. 1 in Murray 2011, equivalent to
                                      %     c = (n12+n22)-(n11+n21)
% CIrand_H = [];
% CIrand_M = [];
% CIrand_CR = [];
% CIrand_FA = [];
% CIrand = [];
% CIrand_targetpresent = [];
% CIrand_targetabsent = [];
% CIboot_H = [];
% CIboot_M = [];
% CIboot_CR = [];
% CIboot_FA = [];
% CIboot = [];
% CIboot_targetpresent = [];
% CIboot_targetabsent = [];
 
% 
% if n_boot>0
%     for i_boot = 1:n_boot
%         list_trial = randi(n_trials, 1, n_trials);
%         boot_response=list_response(list_trial);
%         boot_signal=list_signal(list_trial);
%         boot_Enorm=Enorm(:,:,list_trial);
%         CIboot_H = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==2),3);
%         CIboot_M = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==1),3);
%         CIboot_CR = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==1),3);
%         CIboot_FA = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==2),3);
%         CIboot(:,:,i_boot) = 0.5 * (CIboot_H + CIboot_FA - CIboot_M - CIboot_CR);%CIrand(i_channel,:,i_rand) + Enorm(i_channel,:,i_trial)*correlator/n_trials;
%         CIboot_targetpresent(:,:,i_boot) = CIboot_H - CIboot_M;
%         CIboot_targetabsent(:,:,i_boot) = CIboot_FA - CIboot_CR;
%     end
% end
% end
% 

if n_rand>0
    for i_rand = 1:n_rand
        rand_response = list_response(randperm(n_trials));
        CIrand_H  = mean(Enorm(:,:,list_signal==resp_if_tar & rand_response==resp_if_tar),dim_for_mean);
        CIrand_M  = mean(Enorm(:,:,list_signal==resp_if_tar & rand_response==resp_if_ref),dim_for_mean);
        CIrand_CR = mean(Enorm(:,:,list_signal==resp_if_ref & rand_response==resp_if_ref),dim_for_mean);
        CIrand_FA = mean(Enorm(:,:,list_signal==resp_if_ref & rand_response==resp_if_tar),dim_for_mean);
        CIrand_targetpresent(:,:,i_rand) = CIrand_H - CIrand_M;
        CIrand_targetabsent(:,:,i_rand) = CIrand_FA - CIrand_CR;
        % CIrand(:,:,i_rand) = 0.5 * (CIrand_H + CIrand_FA - CIrand_M - CIrand_CR); %%% From Leo's
        
        CIrand(:,:,i_rand) = CIrand_H + CIrand_FA - CIrand_M - CIrand_CR; %%% As Alejandro thinks
        % CIrand(i_channel,:,i_rand) + Enorm(i_channel,:,i_trial)*correlator/n_trials;
    end
end

if n_boot>0
    for i_boot = 1:n_boot
        list_trial    = randi(n_trials, 1, n_trials);
        boot_response = list_response(list_trial);
        boot_signal   = list_signal(list_trial);
        boot_Enorm    = Enorm(:,:,list_trial);
        CIboot_H  = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==2),dim_for_mean);
        CIboot_M  = mean(boot_Enorm(:,:,boot_signal==2 & boot_response==1),dim_for_mean);
        CIboot_CR = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==1),dim_for_mean);
        CIboot_FA = mean(boot_Enorm(:,:,boot_signal==1 & boot_response==2),dim_for_mean);
        CIboot(:,:,i_boot) = 0.5 * (CIboot_H + CIboot_FA - CIboot_M - CIboot_CR); %%% From Leo's
        % CIboot(:,:,i_boot) = CIboot_H + CIboot_FA - CIboot_M - CIboot_CR; %%% As Alejandro thinks
        
        CIboot_targetpresent(:,:,i_boot) = CIboot_H - CIboot_M;
        CIboot_targetabsent(:,:,i_boot) = CIboot_FA - CIboot_CR;
    end
end