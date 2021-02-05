function Xlabel(string,FontSize)
% function Xlabel(string,FontSize)

if nargin == 1
    FontSize = 14;
end

h = xlabel(string);
set(h,'FontSize',FontSize)

end
