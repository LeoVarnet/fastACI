function Generate_stimuli_from_cfg_crea(cfg_crea_file)
% function Generate_stimuli_from_cfg_crea(cfg_crea_file)
%
% 1. Description:
%   Make sure that the cfg_crea_file is visible to MATLAB. The best way is
%   to specify the full location of the file.
%
% 2. Stand-alone example:
%       dir_where = [fastACI_paths('dir_data') 'speechACI_Logatome-abda-S43M' filesep 'SLV' filesep 'Results' filesep];
%       % Local dir_data for Alejandro: /home/alejandro/Documents/Databases/data/fastACI_data/ (unix)
%       crea_file = 'cfgcrea_2021_12_06_09_54_SLV_speechACI_Logatome-abda-S43M_sMPSv1p3.mat';
%       cfg_crea_file = [dir_where crea_file]; % full location
%       Generate_stimuli_from_cfg_crea(cfg_crea_file);
%
% See also:
%       g20211209_recreating_waveforms_sMPS.m (fastACI_sim repo)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help Generate_stimuli_from_cfg_crea;
    return;
end

cfg_crea = [];
load(cfg_crea_file,'cfg_crea');

% Updates the cfg_crea directory to the local settings, if needed:
cfg_crea = Check_cfg_crea_dirs(cfg_crea);

% Looks for the init file for the current experiment and runs it:
script_init = [cfg_crea.experiment '_init'];
exp2eval = sprintf('%s(cfg_crea)',script_init);
eval(exp2eval);