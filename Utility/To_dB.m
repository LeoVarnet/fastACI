function gain_dB = To_dB(gain)
% function gain_dB = To_dB(gain)
% 
% To_dB: Convert voltage gain to decibels.
% gain_dB = To_dB(gain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright     : Cochlear Ltd
% $Change       : 86418 $
% $Revision     : #1 $
% $DateTime     : 2008/03/04 14:27:13 $
% Authors       : Brett Swanson
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 30/07/2014 % Update this date manually
% Last use on   : 30/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gain_dB = 20 * log10(gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
