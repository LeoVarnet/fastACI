function [version,version_str] = fastACI_version
% function version = fastACI_version
%
% Versions:
% 1.3 on 03/05/2023
% 1.2 on 11/11/2022
% 1.1 on 09/05/2022
% 1.0 on 10/09/2021
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = 1.3;
if nargout >= 2
    version_str=sprintf('version %.1f',version);
end