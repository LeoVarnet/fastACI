function [y,env] = ch_am(sig,fmod,x,option,fs,start_phase)
% function [y,env] = ch_am(sig,fmod,x,option,start_phase);
% 
% 1. Description:
% DE:   am.m amplitudenmoduliert beliebiges Signal 'sig' mit Modulationsfrequenz 
%       fmod. 'option' ist entweder m oder d, je nachdem, ob Modulationsgrad 
%       'm' (in %) oder Modulationsmass d (dB) in 'x' angegeben wird. Falls 
%       sig eine Matrix ist, arbeitet am spaltenweise, y ist dann eine matrix
% 
% EN:   am.m modulates in amplitude signal 'sig' with modulation frequency
%       f0mod. 'option' is either 'm' or 'd', depending on whether modulation
%       factor 'm' or modulation depth 'd' is given. If 'sig' is a matrix, 
%       working on columns, 'y' is then a also a matrix.
% 
% optional: Abtastrate (default 44,1 kHz)
% optional: bei start_phase = -pi/2 beginnt Modulation im Minimum (default)
%
% 2. Additional info:
%       Tested cross-platform: Yes
% 
% 3. Example:
%       fc = 2000; % carrier frequency
%       T = 1/fc; % period
%       fmod = 200; % modulation frequency
%       Tmod = 3/fmod;
%       dur = 1*Tmod; 
%       fs = 44100; % sampling frequency
%       sig = Create_sin(fc,dur,fs,0);
%       m = 100-11;
%       option = 'm';
%       start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
%       ch_am(sig,fmod,m,option,fs,start_phase);
% 
% References: Schoene1979
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    fs=44100;
end    

if nargin < 6
   start_phase=-pi/2;
end

if nargin < 7
   korr=0;
else
   % EN: Calculate correction so that maximum level remains the same (as unmodulated signal)
   % DE: Korrektur berechnen, damit Maximalpegel gleich bleibt (wie unmoduliertes Signal)
   if strcmp(option,'d')
      d=x;
   else
      d=20*log10((1+x/100)/(1-x/100));
   end
   korr=20.*log10(2.*10.^(d./20)./(1+10.^(d./20))); %vgl. Schone1979
end   
   
if strcmp(option,'m')
    if x <= 1 & x ~= 0
        disp('Assuming that the value a percentage is...')
    end
    m=x/100;
   
elseif strcmp(option,'d')
   m=(10^(x/20)-1)/(10^(x/20)+1);
else
   error('falsche Parameter! y=am(sig,fmod,x,option); option: ''d'' oder ''m')
end
   
[r,c]=size(sig);
if r*c == 0,
    y = []; return  % falls leere Matrix: austeigen
end
if (r==1),   % convert row vector to column
    sig = sig(:);  len = c;
else
    len = r;
end
   
t = (0:1/fs:((len-1)/fs))';   % Zeitvektor
t = t(:,ones(1,size(sig,2))); % scalar expansion, falls sig Matrix ist

% env=amp((1 + m * sin(2*pi*fmod*t+start_phase)),-korr);
env = (1 + m * sin(2*pi*fmod*t+start_phase)       );
env = env/max(abs(env));
try
    y = sig .* env; % Modulation beginnt in Minimum
catch
    y = sig .* env'; % Modulation beginnt in Minimum
end

% Anhang: einige Formeln zur Amplitudenmodulation
% Pegelabstand zwischen Tr�ger und Seitenlinien
% deltal=20*log10(2/m);
% d als Funktion von m
% d=20*log10((1+m)/(1-m));
% d bei QFM
% d=10*log10(m^2+1); 

if nargout == 0
    figure;
    subplot(2,1,1)
    plot( t,sig,t,y )
    legend('no AM','with AM');
    xlabel('time [s]')
    grid on
    
    subplot(2,1,2)
    plot( t,env,t,y )
    xlabel('time [s]')
    grid on
    
    legend('Envelope', 'AM-signal')
end