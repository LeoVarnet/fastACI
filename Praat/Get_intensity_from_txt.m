function [t,I] = Get_intensity_from_txt(filename)
% function [t,I] = Get_intensity_from_txt(filename)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

literal         = '--undefined--'; % Character when no F0 is detected
count           = 1;

fid             = fopen(filename);
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
            I(count)   = str2double(tlinesplit{2}); 
            count       = count+1;
        end
    else
        disp('')
        t(count) = str2double(tline(1:end-length(literal)-1));
        I(count) = NaN;
        count = count+1;
    end
     
    tline = fgetl(fid); % Get next line
end
fclose(fid);

t = t(:);
I = I(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end

function tlinesplit = il_strsplit_line(tline)

tlinesplit = strsplit(tline,' ');

for i = length(tlinesplit):-1:1
    if strcmp(tlinesplit{i},sprintf('\t')) || isempty(tlinesplit{i})
            tlinesplit(i) = [];
    end
end