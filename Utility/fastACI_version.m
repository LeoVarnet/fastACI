function [version,version_str] = fastACI_version
% function version = fastACI_version

version = 1.0;
if nargout >= 2
    version_str=sprintf('version %.1f',version);
end