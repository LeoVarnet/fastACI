function r=strsplit(s,d)
% function r=strsplit(s,d)
%
% split string s into cell of strings using delimiter d
 
% r={};
% count=1;
% while (1)
%     [P,s] = strtok(s,d);
%     if (~isempty(s) )
%         r{count}=P;
%         count=count+1;
%     else
%         if (~isempty(P))
%             r{count}=P;
%             count=count+1;
%         end
%         break;
%     end
% end

r={};
prev=1;
count=1;
for i=1:length(s)
    if (s(i)==d)
        if (prev<=i-1)
            r{count}=s(prev:i-1);
        end
        count=count+1;
        prev=i+1;
    end
end

r{count}=s(prev:end);
%
