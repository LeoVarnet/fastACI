function pres_osses2023_09_FA(varargin)
% function pres_osses2023_09_FA(varargin)
%
% % To display Fig. 3 of the Forum Acusticum presentation by Osses and 
% %   Varnet, (2023, Forum Acusticum) use :::
%     pres_osses2023_09_FA('fig3');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

if nargin == 0
    help pres_osses2023_09_FA;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig3'}; % Ahumada's data
    
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_out = [pwd filesep];

if flags.do_fig3
    % Weights in pixels:
    Weights_RM = [  18  -9  7 -8 18; ... % 600
                    0    0  0  0  0; ... % 550
                   -14   0 67 22 15; ... % 500
                   -20   9  0  0  0; ... % 450 Hz
                     0  14  0  0  0];    % 400

    Weights_KL = [   0   0 -16  9 13; ... % 600
                     9   0 -11  0 10; ... % 550
                     0   9  57 18 25; ... % 500
                     0   0 -12  0  0; ... % 450 Hz
                    -4  -3   0  0 -12];   % 400

    val_max = max([Weights_KL(:); Weights_RM(:)]);
    val_max = 1.2*val_max;

    Weights_RM = Weights_RM/val_max;
    Weights_KL = Weights_KL/val_max;

    YT  = [0 1 2 3 4];
    YTL = [400 450 500 550 600];

    offsety = repmat([4; 3; 2; 1; 0],1,5);

    ti = [0:.1:.4];
    figure;
    for i = 1:size(Weights_RM,1)
        y_var = Weights_RM(i,:)+offsety(i,:);
        plot(ti,y_var,'b-'); hold on; grid on;
        idx = find(y_var~=offsety(i,1));
        if ~isempty(idx)
            plot(ti(idx),y_var(idx),'bo','MarkerFaceColor','b'); 
        end
    end
    ylim([-.5 4.5]);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YTL);
    set(gca,'XTick',ti);
    title('Participant RM');

    ylabel('Frequency (Hz)');
    xlabel('Starting time of the segment (s)');

    Pos = get(gcf,'Position');
    Pos(3:4) = [300 425];
    set(gcf,'Position',Pos);

    opts = [];
    if iswindows
        opts.format = 'emf';
    else
        opts.format = 'eps';
    end
    h(end+1) = gcf;

    hname{end+1} = 'fig3-Ahumada-data-1';
    Saveas(h(end),[dir_out hname{end}],opts);
    %%%
    for i = 1:size(Weights_KL,1)
        y_var = Weights_KL(i,:)+offsety(i,:);
        plot(ti,y_var,'k--'); 
        idx = find(y_var~=offsety(i,1));
        if ~isempty(idx)
            plot(ti(idx),y_var(idx),'ks','MarkerFaceColor','w');
        end
    end
    ylim([-.5 4.5]);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',YTL);
    set(gca,'XTick',ti);
    title('Participants RM+KL');

    h(end+1) = gcf;
    hname{end+1} = 'fig3-Ahumada-data-2';
    Saveas(h(end),[dir_out hname{end}],opts);
    
end
