function [mue,corrmue] = optimal_detector(ir_stim,template,fs)
% function [mue,corrmue] = optimal_detector(ir_stim,template,fs)
% 
%   1. Description:
%       OPTIMALDETECTOR  Generic optimal detector for the CASP and Breebaart models
%
%       This is a correlation-based optimal detector for a signal known 
%       exactly see Green & Swets (1966) for more details.
%
%       Url: http://amtoolbox.sourceforge.net/doc/modelstages/optimaldetector.php
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 % if fs is not specified, then AMT default is used
    normmethod = 1; % default in AMT
else
    normmethod = 2;
end

corrmue = ir_stim.*template;

N = numel(corrmue);
switch normmethod
    case 1
        optfactor = sqrt(N); % default in AMT
        % Take mean over all dimensions of internal representation and correct for
        % optimalityfactor.
        mue = optfactor*mean(corrmue(:));
    case 2 % According to exact Equation as shown in Green & Swets 1966, pp. 163
        optfactor = 1/fs;
        mue = optfactor*sum(corrmue(:));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
