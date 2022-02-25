function [NewLevel, NewDirection, NewSpeedFactor] = LevelControlVlLX(LastLevel, LastResult, TargetIntelligibility, LastDirection, LastSpeedFactor)
% function [NewLevel, NewDirection, NewSpeedFactor] = LevelControlVlLX(LastLevel, LastResult, TargetIntelligibility, LastDirection, LastSpeedFactor)
%
% Level control OldenburgerL50 without slope according to: 
%   Brand, T. and Kollmeier, B. (2002). "Efficient adaptive procedures for 
%       threshold and concurrent slope estimations for psychophysics and 
%       speech intelligibility tests", J. Acoust. Soc. Am. 111(6), 2801-2810
%       Implemented by Daniel Berg , HoerTEch gGmbH November 2009.
%
% Extended for target intellgibillities other than 0.5 with additional
% 'stepsizer-multiplier':
% level steps are multiplied by StepMultiplier if 
%    - current SpeedFactor > StepMultiplierMinSpeedFactor  
%    AND one of the next conditions is true: 
%    + direction is 'down' AND TargetIntelligibility >= DownStepTargetMultiplier OR 
%    + direction is 'up' AND TargetIntelligibility <= UpStepTargetMultiplier
%
% Author: Daniel Berg 2012
% Version history
% 10.5.2012 First version
% 1.6.2012  Bugfix if trial intelligibility of first trial is identical to
%           target intelligibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default parameters from VlMatrix levfit-configuration:
Slope = 0.15; % This parameter is different from OLSA default
MinSpeedFactor = 0.25;
Divider = 1.41;
SpeedFactor = 1.5;
StepMultiplier = 2.0;
StepMultiplierMinSpeedFactor = 0.5;
DownStepTargetMultiplier = 0.8;
UpStepTargetMultiplier = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LocalLastDirection = [];

if (nargin == 5)
    SpeedFactor = LastSpeedFactor;
    LocalLastDirection = LastDirection;
elseif  (nargin ~= 3)
    error('wrong number of arguments');
end

% special if result is exactly TargetIntelligibility: return unchanged values 
if (LastResult == TargetIntelligibility)
   % if we are within the track (i.e. a 'lastdirection' already exists then
   % we keep this direction
   if (nargin == 5)
       NewDirection = LastDirection;
   % otherwise (first trial) we return a direction of 0
   % to indicate that no turning point if direction changes in next trial
   else
       NewDirection = 0;
   end
   NewSpeedFactor = SpeedFactor;
   NewLevel = LastLevel;
   return;
end

% check, if StepMultiplier is to be skipped due to current SpeedFactor
% BEFORE calculating the new SpeedFactor!
ApplyMultiplier = 1;
if SpeedFactor <= StepMultiplierMinSpeedFactor
    ApplyMultiplier = 0;
end
    
% calculate new 'direction' (level up (1) or down (-1)) or identical
NewDirection = 1;
if (LastResult > TargetIntelligibility)
    NewDirection = -1;
end

% if there is already a 'LastDirection', then this is not the first trial
% and we have to check for a reversal point and calculate new speed factor
if ~isempty(LocalLastDirection)
    % direction to be changed?
    if LocalLastDirection ~= 0 & LocalLastDirection ~= NewDirection
        SpeedFactor = SpeedFactor / Divider;
        % apply constraint for minimum speed factor
        if (SpeedFactor < MinSpeedFactor)
            SpeedFactor = MinSpeedFactor;
        end
        
    end
end

% calculate step to apply
Step = (LastResult - TargetIntelligibility) * (SpeedFactor / Slope);
 
% no check if to apply multiplier at all, and if one of the other
% constraints (see above) is fulfilled
if ApplyMultiplier  
    if (     (NewDirection == -1 && TargetIntelligibility >= DownStepTargetMultiplier) ...
        ||  (NewDirection == 1  && TargetIntelligibility <= UpStepTargetMultiplier))
        Step = Step * StepMultiplier;
    end
end

% return direction and speed factor (direction set above)
NewSpeedFactor = SpeedFactor;

% return new level
NewLevel = LastLevel - Step;
