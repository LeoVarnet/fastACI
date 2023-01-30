function addnoise = wallaert2017_internalnoiseadd(intrep,fs,intnoise_addstd)
% function addnoise = wallaert2017_internalnoiseadd(intrep,fs,intnoise_addstd)
%
% 1. Description:
%      Additive internal noise. 
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

if ~iscell(intrep)
    addnoise  = randn(size(intrep))*intnoise_addstd;
else
    % As often in AMT, with internal representations as cell arrays
    % [N_samples,N_chann,N_mod_chann] = size(intrep{end}); % last filter should be the largest
    
    error('Not validated for cell arrays yet')
end