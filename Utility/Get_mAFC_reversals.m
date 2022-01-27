function [r,idx,reve,staircase] = Get_mAFC_reversals(staircase)
% function [r idx reve staircase] = Get_mAFC_reversals(staircase)
%
% 1. Description:
%
% 2. Stand-alone example:
%       r20151120_update_MA;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/11/2015
% Last update on: 20/11/2015 
% Last use on   : 12/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
    plot(staircase); hold on;
end

reve = [];
initial = -1; % no direction defined at the beginning
for i = 2:length(staircase)
    diffe = staircase(i)-staircase(i-1);
    
    if diffe > 0
        tmp(i-1) = 1;
        if initial == -1; initial = 1; end; % staircase starts increasing
        if initial ~= 1
            reve(end+1) = i-1;
            initial = ~initial;
        end
    elseif diffe < 0
        tmp(i-1) = 0;
        if initial == -1; initial = 0; end; % staircase starts decreasing
        if initial ~= 0
            reve(end+1) = i-1;
            initial = ~initial;
        end
    elseif diffe == 0 % no reversal
        tmp(i-1) = initial;
    end
    
end

if nargout == 0
    plot(reve, staircase(reve), 'o','LineWidth',2)
end

r =  staircase(reve);
idx = reve;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
