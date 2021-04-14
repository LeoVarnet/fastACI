function [version,version_description] = Check_ACI_data_type(file2load)
% function [version,version_description] = Check_ACI_data_type(file2load)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var = load(file2load);

version = 'current'; % starting version, if no other version is detected, then
               % the 'current' version is assumed
version_description = '2021: Current fastACI format';

if isfield(var,'stim_order')
    % This most likely is the Zenodo format...
    
    version = 2015;
    version_description = '2015: as in Varnet et al. 2015';
    %         stim_order: [1×10000 double]
    %           n_signal: [1×10000 double]
    %                SNR: [1×10000 double]
    %               date: [6×10000 double]
    %      response_time: [1×10000 double]
    %     correct_answer: [1×10000 logical]
elseif isfield(var,'i') && isfield(var,'ListStim')
    if isfield(var,'ordre_aleatoire')
        version = 2013;
        version_description = '2013: as in Varnet et al. 2013';
    elseif isfield(var,'cfg_game') && isfield(var,'data_passation')
        version = 2021.1;
        version_description = '2021: Leo''s version as used in his AMrevcorr paper';
    end
end