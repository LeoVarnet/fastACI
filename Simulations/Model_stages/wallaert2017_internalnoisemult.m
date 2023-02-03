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
            
            if ~isnan(rms_Emod(1,i_chan,i_mod))
                multnoise_scaled = multnoise(:,i_chan,i_mod)/mean_multnoise; % rms = 1
                multnoise_scaled = multnoise_scaled*rms_Emod(1,i_chan,i_mod);
            else
                multnoise_scaled = multnoise(:,i_chan,i_mod)*0;
            end
            multnoise(:,i_chan,i_mod) = multratio*multnoise_scaled;
        end
    end
else
    % As often in AMT, with internal representations as cell arrays
    multnoise = cell(size(intrep)); % initialisation
    for i = 1:length(intrep)
        [N_samples,N_chann,N_mod_chann] = size(intrep{i});
        
        rms_Emod = sqrt(mean(intrep{i}.^2,1));
        
        multnoise{i} = randn(size(intrep{i}));
        for i_chan = 1:N_chann
            for i_mod = 1:N_mod_chann
                mean_multnoise = sqrt(mean(multnoise{i}(:,i_chan,i_mod).^2,1));

                if ~isnan(rms_Emod(1,i_chan,i_mod))
                    multnoise_scaled = multnoise{i}(:,i_chan,i_mod)/mean_multnoise; % rms = 1
                    multnoise_scaled = multnoise_scaled*rms_Emod(1,i_chan,i_mod);
                else
                    multnoise_scaled = multnoise{i}(:,i_chan,i_mod)*0;
                end
                multnoise{i}(:,i_chan,i_mod) = multratio*multnoise_scaled;
            end
        end
        
    end
end
