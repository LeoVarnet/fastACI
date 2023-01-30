function r=strsplit(string_in,delimiter)
% function r=strsplit(s,d)
%
% 1. Description:
%     Splits string s into cell of strings using delimiter d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

r={}; % initialisation
prev=1;
count=1;
for i=1:length(string_in)
    if (string_in(i)==delimiter)
        if (prev<=i-1)
            r{count}=string_in(prev:i-1);
        end
        count=count+1;
        prev=i+1;
    end
end

r{count}=string_in(prev:end);
