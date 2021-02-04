function Saveas(h, filename, option)
% function Saveas(h, filename, option)
%
% 1. Description:
%       Function to save figures in vectorial formats. EPS is more suitable 
%       for LaTeX documents and EMF is more suitable if figures are going 
%       to be included in a microsoft PPT presentation.
%           - option.format = 'epsc'; % for eps in color
%           - option.format = 'emf'; for emf format
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium 2014
% Created in    : January-March 2014-2015
% Last update on: 20/08/2014 
% Last use on   : 11/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    option = [];
else
    if ischar(option)
        option.format = option;
    else
        option = Ensure_field(option,'format','epsc');
    end
end
option = Ensure_field(option,'format','epsc');
if ~isfield(option,'Colour')
    option.Colour = 'rgb';
end

str = fileparts(filename);

if strcmp(str,'')
    str = [pwd filesep];
    filename = [str filename];
end

option = Ensure_field(option,'bPrint',1);
option = Ensure_field(option,'bScale',0);

try
    set(h,'PaperType', 'A4')
end

if option.bScale == 1
    try
        Check_figure_size(h); % Constrains figure width to A4 Paper size
    end
end

try
    set(h,'PaperPositionMode', 'auto')
end

if option.bPrint
%     figure(h);
%     Print_date_on_figure;
end

switch option.Colour
    case 'rgb'
        saveas(h, filename, option.format);
    case 'gray'
        style = hgexport('factorystyle');
        style.Color = option.Colour;
        hgexport(h,filename,style);
    otherwise
        saveas(h, filename, option.format);
end

disp(['Figure saved as: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
