function [t,formants] = Get_formants_from_txt(filename)
% function [t,formants] = Get_formants_from_txt(filename)
%
% Reads Fundamental frequency F0 from text file generated in Praat:
%   In Praat:
%       1. Open wav file
%       2. Praat objects window: press View and edit
%       3. On plots, select the whole file
%       4. Go to: Pitch -> Pitch listing
%       5. On the same window select Save as... 
%
% Desirable: to set the Pitch Range from 75 to 300 Hz (Pitch Settings...)
%
% % Example:
%   [t f0] = Get_F0_praat_from_txt('/home/alejandro/Documenten/MATLAB/MATLAB_svn/new_audio/wdz6.txt');
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, 2014
% Created in    : 2014
% Last update on: 19/05/2014
% Last use on   : 19/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

literal         = 'name'; % Character when no F0 is detected
count           = 1;

fid             = fopen(filename);
y               = 0;
tline           = fgetl(fid);

while ischar(tline) && length(tline)>2
    matches = strfind(tline, literal);
    
    if isempty(matches)
        
        tlinesplit = il_strsplit_line(tline);
        % Then it is a non-null F0 value
        disp('')
        match_title     = strfind(tline, 'time');
        if isempty(match_title) % then it is not a header
            t(count)    = str2double(tlinesplit{1}); 
            for i = 2:length(tlinesplit)
                formants(count,i-1)   = str2double(tlinesplit{i}); 
            end
            count       = count+1;
        end
    else
        % disp('')
        num = str2double(tline(1:end-length(literal)-1));
        if ~isnan(num)
            t(count) = str2double(tline(1:end-length(literal)-1));
            formants(count) = NaN;
            count = count+1;
        end
    end
     
    tline = fgetl(fid); % Get next line
end
fclose(fid);

t = t(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end

function tlinesplit = il_strsplit_line(tline)

tlinesplit = strsplit(tline,' ');

for i = length(tlinesplit):-1:1
    if strcmp(tlinesplit{i},sprintf('\t')) || isempty(tlinesplit{i})
            tlinesplit(i) = [];
    end
end