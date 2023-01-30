function memnoise = wallaert2017_internalnoisemem(intrep,fs,intnoise_std,decay_tau)
% function memnoise = wallaert2017_internalnoisemem(intrep,fs,intnoise_std,memorydecay_tau)
%
% 1. Description:
%      Memory internal noise. The memory internal noise is exponentially 
%      decaying with a (default) time constant decay_tau of 1.2 s (half life).
%      
%      fs is the sampling frequency of the internal representation 'intrep'.
%      intnoise_std is the standard deviation of the internal noise.
%      decay_tau is the time constant of the exponentially-decreasing curve
%
% Author: Andrew King, Leo Varnet, Nicolas Wallaert, Stephan Ewert, and
%         Christian Lorenzi: Original implementation
% Author: Alejandro Osses: Integration in fastACI, extension to intrep 
%         contained as cell arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(intrep)
    % Regular case, if the internal representation is a double array
    [N_samples,N_chann,N_mod_chann] = size(intrep);
    t = (1:size(intrep,1))'/fs;    
else
    % As often in AMT, with internal representations as cell arrays
    [N_samples,N_chann,N_mod_chann] = size(intrep{end}); % last filter should be the largest
    t = (1:size(intrep{1},1))'/fs;
end

[idx_value, ~, ~] = ndgrid(1:N_samples,1:N_chann,1:N_mod_chann);
decay_curve = exp(t(flip(idx_value,1))/decay_tau);

if ~iscell(intrep)
    memnoise = intnoise_std*randn(size(intrep));
    memnoise = memnoise.*decay_curve;
else
    memnoise = cell([length(intrep) 1]); % initialisation
    for i = 1:length(intrep)
        [N,M] = size(intrep{i});
        memnoise{i} = intnoise_std*randn([N M]);
        memnoise{i} = memnoise{i}.*decay_curve(1:N,1:M);
    end
end
