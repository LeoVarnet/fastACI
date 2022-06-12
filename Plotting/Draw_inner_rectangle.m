function files_out = Draw_inner_rectangle(XL,YL,Colour,Style,LineWidth)
% function files_out = Draw_inner_rectangle(XL,YL,Colour,Style,LineWidth)
%
% 1. Description:
%
% 2. Additional info:
%
% 3. Stand-alone example:
%   % Example in Alejandro's computer:
%   folder = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Stimuli/Logatome/';
%   Get_formants_from_dir(folder);
% 
% Programmed by Alejandro Osses, ENS, France, 2022
% Created on    : 10/06/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    XL = get(gca,'XLim');
end
if nargin < 2
    YL = get(gca,'YLim');
end
if nargin < 3
    Colour = 'k';
    fprintf('%s: Using the default colour for the rectangle\n',upper(mfilename));
end
if nargin < 4
    Style = '-'; % continuous line
end
if nargin < 5
    LineWidth = 2;
end
offx = (XL(2)-XL(1))*.01; % 1% of the X-axis length
offy = (YL(2)-YL(1))*.01; % 1% of the Y-axis length

hold on;
plot_options = {Style,'Color',Colour,'LineWidth',LineWidth};
plot(XL,(YL(1)+offy)*[1 1],plot_options{:}); % lower side
plot(XL,(YL(2)-offy)*[1 1],plot_options{:}); % upper side
plot((XL(1)+offx)*[1 1],YL,plot_options{:}); % left side
plot((XL(2)-offx)*[1 1],YL,plot_options{:}); % left side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
