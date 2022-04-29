function Pyra_size = fastACI_impyramid_get_size(Pyramid,Nlevelmin,Nlevel)
% FASTACI_IMPYRAMID_GET_SIZE
%
% It assumes that dimension 1 is 'samples', dimension 2 is frequency and 
%   dimension 3 is time.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pyra_size = zeros([Nlevel 2]); % memory allocation

for i_level = Nlevelmin:Nlevel
    Pyra_here = Pyramid{i_level};
    Pyra_size(i_level,:) = [size(Pyra_here,2) size(Pyra_here,3)];
end