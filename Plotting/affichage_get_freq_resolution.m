function [f2plot,cfg_ACI] = affichage_get_freq_resolution(f,cfg_ACI)
% function [f2plot,cfg_ACI] = affichage_get_freq_resolution(f,cfg_ACI)
%
%%% Examples:
%   publ_osses2021c_DAGA_figs('fig1a');
%   publ_osses2021c_DAGA_figs('fig1b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~isfield(cfg_ACI,'f')
%     cfg_ACI.f = audtofreq(3:.5:33);
%     fprintf('%s.m: Using ERB frequencies by default...\n',mfilename);
% end
if isfield(cfg_ACI,'f')
    %%% Detect the frequency spacing:
    % 1. Detecting if ERB-spaced:
    ferb = freqtoaud(cfg_ACI.f);

    df(1) = ferb(2)-ferb(1);
    df(2) = ferb(end)-ferb(end-1);

    if round(10*df(1)) == round(10*df(2))
        % This should be a constant for at least the first and second element of df
        fprintf('\t%s: ERB-spacing detected...\n',upper(mfilename));
        spacing_type = 'erb';
        delta_f = df(1);
    else
        df_lin = [min(diff(cfg_ACI.f)) min(diff(cfg_ACI.f))];
        if df_lin(1) == df_lin(2) && df_lin(1) == 1
            spacing_type = 'bin';
        else
            spacing_type = 'Hz';
        end

        fprintf('%.m: Assuming a regular spacing in %s\n',mfilename,spacing_type);
    end

    switch spacing_type
        case 'erb'
            step = (max(cfg_ACI.f)-min(cfg_ACI.f))/(length(cfg_ACI.f)-1); % this is a linear step equal to 0.5 ERB_N
            for i = 1:size(f,2)
                % For each column of f2plot
                f2plot(:,i) = step*(freqtoaud(f(:,i))-min(ferb))/delta_f;
            end

        otherwise
            % warning('Frequency spacing not detected, returning an empty frequency vector')
            f2plot = f;
    end
    
else
    spacing_type = 'Hz';
    f2plot = f;
end