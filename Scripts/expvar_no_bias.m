function [idx_analysis, label_expvar_no_bias, N_excluded] = expvar_no_bias(data_passation,idx_analysis)
% function [select_after_reversal, label_expvar_after_reversal, N_excluded] = expvar_no_bias(data_passation,idx_analysis)
%
%     TRIAL EQUALISATION: if equal to 1, this processing discards trials so
%     that the number of responses (1 or 2) are equal, starting with the 
%     more extreme expvar values.
%
% Author: Leo Varnet, adapted by Alejandro Osses
% Created on: 21/07/2022 - before: this code was hardcoded in fastACI_getACI_preprocess.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expvar      = data_passation.expvar;
n_responses = data_passation.n_responses;
if nargin < 2
    idx_analysis = ones(size(expvar));
end

N_r1 = sum(n_responses(idx_analysis)==1);
N_r2 = sum(n_responses(idx_analysis)==2);
[~,sorted_idx] = sort(abs(expvar(idx_analysis)-median(expvar(idx_analysis)))); % sorting the trials according to the distance to the median expvar
sorted_response = n_responses(idx_analysis(sorted_idx));
sorted_idx_r1 = sorted_idx(sorted_response==1);
sorted_idx_r2 = sorted_idx(sorted_response==2);
if N_r1>N_r2
    idxf = length(sorted_idx_r1);
    idxi = idxf-(N_r1-N_r2);
    trials2exclude = idx_analysis(sorted_idx_r1(idxi:idxf));
elseif N_r2>N_r1
    idxf = length(sorted_idx_r2);
    idxi = idxf-(N_r2-N_r1);
    trials2exclude = idx_analysis(sorted_idx_r2(idxi:idxf));
else
    idxf = 0;
    idxi = 1; % so that idxf-idxi+1 == 0
    trials2exclude = []; % already without a bias
end
N_excluded = idxf-idxi+1;

idx_analysis = setdiff(idx_analysis,trials2exclude);
label_expvar_no_bias = sprintf('\t\t\t(expvar_no_bias: %.0f extra trials are being excluded)\n',N_excluded);
