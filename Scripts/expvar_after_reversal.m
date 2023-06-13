function [select_after_reversal, label_expvar_after_reversal, N_excluded] = expvar_after_reversal(data_passation,expvar_after_reversal,idx_trialselect)
% function [select_after_reversal, label_expvar_after_reversal, N_excluded] = expvar_after_reversal(data_passation,expvar_after_reversal,idx_trialselect)
%
% Author: Alejandro Osses
% Created on: 21/07/2022 - before: this code was hardcoded in fastACI_getACI_preprocess.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_excluded = [];

if nargin < 3
    idx_trialselect = ones(size(data_passation.expvar));
end

select_after_reversal = ones(size(idx_trialselect));

expvar = data_passation.expvar;
% n_responses = data_passation.n_responses;

if expvar_after_reversal > 0
    N_reversals_i = expvar_after_reversal;
    idx_sessions = find(data_passation.resume_trial < length(idx_trialselect));

    L_session = length(data_passation.expvar)-1;

    % L_session = median(diff(data_passation.resume_trial));
    for i = 1:length(idx_sessions)
        idxi = max(data_passation.resume_trial(i),1);
        if i~= length(idx_sessions)
            % L_session_here: It will be different than L_session if the
            %     participant paused the measurements during the sessions:
            L_session_here = data_passation.resume_trial(i+1) - data_passation.resume_trial(i);
            idxf = data_passation.resume_trial(i) + L_session_here-1;
        else
            idxf = data_passation.resume_trial(i) + L_session;
        end
        idxf = min(idxf, length(idx_trialselect)); % limited by the selected number of trials
        [~,idx2check] = Get_mAFC_reversals(expvar(idxi:idxf));
        if ~isempty(idx2check) & length(idx2check)>=N_reversals_i
            idxs4null = (1:idx2check(N_reversals_i)-1) + idxi-1;
        else
            fprintf('%s: Less than %.0f reversals were found, %.0f trials (trials between %.0f and %.0f) will be excluded anyway\n',upper(mfilename),N_reversals_i,idxf-idxi+1,idxi,idxf);
            idxs4null = idxi:idxf;
        end

        select_after_reversal(idxs4null) = 0;
    end
    N_excluded(end+1) = sum(select_after_reversal == 0);
    label_expvar_after_reversal = sprintf('\t\t\t(expvar_after_reversal=%.0f: %.0f extra trials are being excluded)\n', ...
        N_reversals_i, N_excluded(end));
else
    % Nothing to do
end