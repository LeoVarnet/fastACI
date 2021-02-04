function [ ERB ] = f2ERB( f )
% [ ERB ] = F2ERB( f )
%   Calculates the value on the ERB scale corresponding to a frequency f
%   (see Hohmann, 2002) 
%
% Leo Varnet 2016

l = 24.7;
q = 9.265;
ERB = q*log(1+f/(l*q));

end

