function multnoise = wallaert2017_internalnoisemult(intrep,fs,intnoise_multratio_dB)
% function multnoise = wallaert2017_internalnoisemult(intrep,fs,intnoise_multratio)
%
% 1. Description:
%      Multiplicative internal noise. The multiplicative internal noise is 
%      proportional to the energy in each audio-modulation band of intrep.
%      In this implementation, intnoise_multratio is expressed in 'dB'.
%      
%      The sampling frequency fs is not used in this implementation but it
%      was kept as input parameter for consistency with the other internal
%      noise scripts.
%      
% Author: Andrew King, Leo Varnet, Nicolas Wallaert, Stephan Ewert, and
%         Christian Lorenzi: Original implementation
% Author: Alejandro Osses: Integration in fastACI, extension to intrep 
%         contained as cell arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multratio = 10^(intnoise_multratio_dB/20)-1;

if ~iscell(intrep)
    rms_Emod = sqrt(mean(intrep.^2,1));

    [N_samples,N_chann,N_mod_chann] = size(intrep);
    multnoise = randn([N_samples,N_chann,N_mod_chann]);

    for i_chan = 1:N_chann
        for i_mod = 1:N_mod_chann
            mean_multnoise = sqrt(mean(multnoise(:,i_chan,i_mod).^2,1));
            
            multnoise_scaled = multnoise(:,i_chan,i_mod)/mean_multnoise; % rms = 1
            multnoise_scaled = multnoise_scaled*rms_Emod(1,i_chan,i_mod);
            multnoise(:,i_chan,i_mod) = multratio*multnoise_scaled;
        end
    end
else
    % As often in AMT, with internal representations as cell arrays
    [N_samples,N_chann,N_mod_chann] = size(intrep{end}); % last filter should be the largest
    
    error('Not validated for cell arrays yet')
end

% [idx_value, ~, ~] = ndgrid(1:N_samples,1:N_chann,1:N_mod_chann);
% decay_curve = exp(t(flip(idx_value,1))/decay_tau);
% 
% if ~iscell(intrep)
%     memnoise = intnoise_std*randn(size(intrep));
%     memnoise = memnoise.*decay_curve;
% else
%     memnoise = cell([length(intrep) 1]); % initialisation
%     for i = 1:length(intrep)
%         [N,M] = size(intrep{i});
%         memnoise{i} = intnoise_std*randn([N M]);
%         memnoise{i} = memnoise{i}.*decay_curve(1:N,1:M);
%     end
% end
