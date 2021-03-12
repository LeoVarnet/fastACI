function version = Check_ACI_data_type(file2load)
% function version = Check_ACI_data_type(file2load)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var = load(file2load);

version = 'current'; % starting version, if no other version is detected, then
               % the 'current' version is assumed

if isfield(var,'stim_order')
    % This most likely is the Zenodo format...
    
    version = 2015;
    %         stim_order: [1×10000 double]
    %           n_signal: [1×10000 double]
    %                SNR: [1×10000 double]
    %               date: [6×10000 double]
    %      response_time: [1×10000 double]
    %     correct_answer: [1×10000 logical]
elseif isfield(var,'i') && isfield(var,'ordre_aleatoire') && isfield(var,'ListStim')
    version = 2013;
end