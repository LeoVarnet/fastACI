function outs_from_Praat = Add_formants(fname,cfg_ACI,par_formants,Colour)
% function outs_from_Praat = Add_formants(fname,cfg_ACI,par_formants,Colour)
%
% 1. Description:
%   This function assesses the Praat parameters for the sound 'fname'. The
%   Praat parameters are stored on a txt file in the same folder as fname.
%   Subsequently, the txt files are read and are added to a preexisting figure
%   that was generated using affichage_tf function.
%
%   This function requires that Praat is configured in your system. Please
%   configure the location of Praat such that fastACI_paths('praat') returns
%   a valid location in your computer.
%
% 2. Example:
%   publ_osses2022b_JASA_figs('fig1');
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    Colour = 'w';
end
if nargin < 3
    par_formants = [];
end
LW = 1;

if nargout == 0
    bPlot = 1;
else
    bPlot = 0;
end

par_formants = Ensure_field(par_formants,'timestep',0.01); % positive timestep 0.01
par_formants = Ensure_field(par_formants,'nformants',5); % positive nformants 5

% Formants
par_formants = Ensure_field(par_formants,'maxformant',5500); % positive maxformant 5500
par_formants = Ensure_field(par_formants,'windowlength',0.025);% 0.025 % positive windowlength 0.025
par_formants = Ensure_field(par_formants,'dynamicrange',30); % positive dynamic range 20

% F0
par_formants = Ensure_field(par_formants,'minpitch',200); % positive minimum pitch 50 (for intensity)
%par_formants.pitchfloor = 100; % previous parameter value (14/10/2022)
par_formants = Ensure_field(par_formants,'pitchfloor',50); % positive pitch floor 100 (for f0)
par_formants = Ensure_field(par_formants,'pitchceiling',500); % positive pitch ceiling 500 (for f0)

% Before 4/11/2021, I_min set to 40 dB:
par_formants = Ensure_field(par_formants,'I_min',59);%75; %, arbitrary value

outs_from_Praat = affichage_tf_add_Praat_metrics_one_sound(fname,cfg_ACI,par_formants, '-', Colour, LW, bPlot);