function Ylabel(string,FontSize)
% function Ylabel(string,FontSize)

if nargin == 1
    FontSize = 14;
end

h = ylabel(string);
set(h,'FontSize',FontSize)

end
