function var2latex(matrix, filename, varargin)
% function var2latex(matrix, filename, varargin)
%
% 1. Description:
%       Converts input 'matrix' (MATLAB variable) to a LaTeX table
% 
% 2. Additional info:
%       MATLAB Command line tested: YES
%       Generation in a file tested: NOT YET
% 
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       x = rand(5,4);
%       var2latex(x); % Make sure 'x' does exist
% 
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file: latex.m, programmed by  M. Koehler (koehler@in.tum.de)
% Last update on: 19/06/2014
% Last use on   : 31/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rowLabels   = [];
colLabels   = [];
alignment   = 'l';
hlines      = 0;
vlines      = 0;
format      = [];
textsize    = [];

if nargin < 2
    filename = '__MCL'; % default printed on MATLAB command line
end

okargs = {'rowlabels','columnlabels', 'alignment', 'format', 'size', 'vlines', 'hlines'};
for j=1:2:(nargin-2)
    pname = varargin{j};
    pval = varargin{j+1};
    k = strmatch(lower(pname), okargs);
    if isempty(k)
        error('matrix2latex: ', 'Unknown parameter name: %s.', pname);
    elseif length(k)>1
        error('matrix2latex: ', 'Ambiguous parameter name: %s.', pname);
    else
        switch(k)
            case 1  % rowlabels
                rowLabels = pval;
                if isnumeric(rowLabels)
                    rowLabels = cellstr(num2str(rowLabels(:)));
                end
            case 2  % column labels
                colLabels = pval;
                if isnumeric(colLabels)
                    colLabels = cellstr(num2str(colLabels(:)));
                end
            case 3  % alignment
                alignment = lower(pval);
                if alignment == 'right'
                    alignment = 'r';
                end
                if alignment == 'left'
                    alignment = 'l';
                end
                if alignment == 'center'
                    alignment = 'c';
                end
                if alignment ~= 'l' && alignment ~= 'c' && alignment ~= 'r'
                    alignment = 'l';
                    warning('matrix2latex: ', 'Unkown alignment. (Set it to \''left\''.)');
                end
            case 4  % format
                format = lower(pval);
            case 5  % Font size
                textsize = pval;
            case 6  % Vertical lines or not
                vlines = pval;
            case 7  % Horizontal lines or not
                hlines = pval;
        end
    end
end


if isequal(filename, '__MCL')

    outputstr = '';

    width = size(matrix, 2);
    height = size(matrix, 1);

    if isnumeric(matrix)
        matrix = num2cell(matrix);
        for h=1:height
            for w=1:width
                if(~isempty(format))
                    matrix{h, w} = num2str(matrix{h, w}, format);
                else
                    matrix{h, w} = num2str(matrix{h, w});
                end
            end
        end
    end

    if(~isempty(textsize))
        outputstr = strcat(outputstr, sprintf('\\begin{%s}', textsize) );
    end
    
    tabletype = 'tabular';
   
    outputstr = strcat(outputstr,sprintf(['\r\\n\b\b\r\\n\b\b\r\n\\begin{table}\\n\b\b\r\b\b\r\n\\begin{' tabletype '}{']) );

    if(~isempty(rowLabels))
        outputstr = strcat(outputstr, sprintf( 'l|') );
    end
    for i=1:width
        
        % IF vlines is enabled then add them in the tabular.
        if isequal(vlines, 1)
            outputstr = strcat(outputstr, sprintf('%c|', alignment) );
        else
            % Do not add them
                outputstr = strcat(outputstr, sprintf('%c', alignment) );            
        end
    end
    outputstr = strcat(outputstr, sprintf('}\r\\n\b\b') );

    
    if(~isempty(colLabels))
        if(~isempty(rowLabels))
            outputstr = strcat(outputstr, sprintf('&') );
        end
        for w=1:width-1
            outputstr = strcat(outputstr, sprintf('\\textbf{%s}&', colLabels{w}) );
        end
        outputstr = strcat(outputstr, sprintf('\\textbf{%s}\\\\\\hline\r\\n\b\b', colLabels{width}) );
    end

    for h=1:height
        if(~isempty(rowLabels))
            outputstr = strcat(outputstr, sprintf('\\textbf{%s}&', rowLabels{h}) );
        end
        for w=1:width-1
            outputstr = strcat(outputstr, sprintf('%s&', matrix{h, w}) );
        end
        
        % If the horizontal lines are enabled, add them.
        if isequal(hlines, 1)
            outputstr = strcat(outputstr, sprintf('%s\\\\\r\\n\b\b\\hline\r\\n\b\b', matrix{h, width}) );
        else
            if ~isequal(h, height)
                outputstr = strcat(outputstr, sprintf('%s\\\\\r\\n\b\b', matrix{h, width}) );
            else
                outputstr = strcat(outputstr, sprintf('%s\\\\\r\\n\b\b', matrix{h, width}) );
            end
        end
    end

    outputstr = strcat(outputstr, sprintf(['\\end{tabular}\n\\end{table}\n']) );

    if(~isempty(textsize))
        outputstr = strcat(outputstr, sprintf('\\end{%s}', textsize) );
    end
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('% Start copying here')
    disp( sprintf('%s',outputstr) )
    disp('% Finish copying here')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

else

    fid = fopen(filename, 'w');

    width = size(matrix, 2);
    height = size(matrix, 1);

    if isnumeric(matrix)
        matrix = num2cell(matrix);
        for h=1:height
            for w=1:width
                if(~isempty(format))
                    matrix{h, w} = num2str(matrix{h, w}, format);
                else
                    matrix{h, w} = num2str(matrix{h, w});
                end
            end
        end
    end

    if(~isempty(textsize))
        fprintf(fid, '\\begin{%s}', textsize);
    end

    fprintf(fid, '\\begin{tabular}{|');

    if(~isempty(rowLabels))
        fprintf(fid, 'l|');
    end
    for i=1:width
        fprintf(fid, '%c|', alignment);
    end
    fprintf(fid, '}\r\n');

    fprintf(fid, '\\hline\r\n');

    if(~isempty(colLabels))
        if(~isempty(rowLabels))
            fprintf(fid, '&');
        end
        for w=1:width-1
            fprintf(fid, '\\textbf{%s}&', colLabels{w});
        end
        fprintf(fid, '\\textbf{%s}\\\\\\hline\r\n', colLabels{width});
    end

    for h=1:height
        if(~isempty(rowLabels))
            fprintf(fid, '\\textbf{%s}&', rowLabels{h});
        end
        for w=1:width-1
            fprintf(fid, '%s&', matrix{h, w});
        end
        fprintf(fid, '%s\\\\\\hline\r\n', matrix{h, width});
    end

    fprintf(fid, '\\end{tabular}\r\n');

    if(~isempty(textsize))
        fprintf(fid, '\\end{%s}', textsize);
    end

    fclose(fid);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
