function files_out = Get_formants_from_dir(folder, params)
% function files_out = Get_formants_from_dir(folder, params)
%
% 1. Description:
%
% 2. Additional info:
%
% 3. Stand-alone example:
%   % Example in Alejandro's computer:
%   folder = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Stimuli/Logatome/';
%   Get_formants_from_dir(folder);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/08/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files_out = [];

% bDiary = 1;
% Diary(mfilename,bDiary,folder);

if nargin < 2
    params = [];
end

params = Ensure_field(params,'timestep'    , 0.01); % positive timestep 0.01
params = Ensure_field(params,'nformants'   ,    5); % positive nformants 5
params = Ensure_field(params,'maxformant'  , 5500); % positive maxformant 5500
params = Ensure_field(params,'windowlength',0.025); % positive windowlength 0.025
params = Ensure_field(params,'dynamicrange',   20); % positive dynamic range 20

%%%
bRun_Praat = 0;

files_wav = Get_filenames(folder,'*.wav');
for i = 1:length(files_wav)
	files_out{i} = [files_wav{i}(1:end-4) '_F.txt'];
    if ~exist([folder files_out{i}],'file') % Check if they are on disk:
        bRun_Praat = 1;
    end
end
%%%

if bRun_Praat
    try
        dir_praat = fastACI_paths('praat');
    catch 
        if isunix
            % Typical location for praat on a unix system:
            dir_praat = '/usr/bin/praat';
        end
    end
    local_praat     = dir_praat;
    local_praat_sc  = [fastACIsim_basepath 'MATLAB' filesep 'Praat' filesep]; % Get_TUe_paths('praat_scripts');
    script_name     = 'Get_formants_from_dir.praat';
    script          = [local_praat_sc script_name];

    if isunix
        command4system = [local_praat ' ' script ' ' folder ' ' ...
                                                  '"' num2str(params.timestep)     '" ' ...
                                                  '"' num2str(params.nformants)    '" ' ...
                                                  '"' num2str(params.maxformant)   '" ' ...
                                                  '"' num2str(params.windowlength) '" ' ...
                                                  '"' num2str(params.dynamicrange) '"'];
    else
         command4system = [local_praat ' ' script ' ' folder ' ' ...
                                                   num2str(params.timestep)     ' ' ...
                                                   num2str(params.nformants)    ' ' ...
                                                   num2str(params.maxformant)   ' ' ...
                                                   num2str(params.windowlength) ' ' ...
                                                   num2str(params.dynamicrange)];
    end

    disp([mfilename '.m: ' command4system])
    [s,r] = system( command4system );

    if s == -1
        disp([mfilename '.m: problem running Praat, please check that praat is installed and the scripts do exist...'])
    else
        %%% Then it was sucessful
        fprintf('Praat script %s successfully completed\n',script_name);

        for i = 1:length(files_wav)

            if isunix
                % In unix (Ubuntu), the formant file gets a backslash, here this is corrected:
                if exist([folder '\' files_out{i}],'file')
                    movefile([folder '\' files_out{i}],[folder files_out{i}]);
                end 
            end
        end

        % if bDiary
        % 	diary off
        % end
    end
else
    fprintf('\t%s: All requested Praat files are already on disk. Praat was not run, and the stored results are read instead\n',upper(mfilename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
