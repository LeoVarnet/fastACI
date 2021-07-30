function f2plot = affichage_get_freq_resolution(f,cfg_ACI)
% function f2plot = affichage_get_freq_resolution(f,cfg_ACI)

spacing_type = [];
%%% Detect the frequency spacing:
% 1. Detecting if ERB-spaced:
ferb = freqtoaud(cfg_ACI.f);

df(1) = ferb(2)-ferb(1);
df(2) = ferb(end)-ferb(end-1);

if round(10*df(1)) == round(10*df(2))
    fprintf('\t%s: ERB-spacing detected...\n',upper(mfilename));
    spacing_type = 'erb';
    delta_f = df(1);
end

switch spacing_type
    case 'erb'
        step = (max(cfg_ACI.f)-min(cfg_ACI.f))/(length(cfg_ACI.f)-1); % this is a linear step equal to 0.5 ERB_N
        for i = 1:size(f,2)
            f2plot(:,i) = step*(freqtoaud(f(:,i))-min(ferb))/delta_f;
        end

    otherwise
        warning('Frequency spacing not detected, returning an empty frequency vector')
        f2plot = nan(size(f));
end