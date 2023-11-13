function paths = fastACI_paths(type)
% function paths = fastACI_paths(type)
%
% 1. Description: It gets all the paths related to the fastACI toolbox or only
%   one directory (if 'type' is defined, see the examples below).
%   The directory dir_data is expected to contain one folder for each fastACI
%   experiment, within which each participant has one folder containing target
%   and noise stimuli. 
%   The current list of directories  that are loaded in this script is:
%
% dir_data:   Directory where the stimuli and other 'binary' files of each
%             experiment will be stored
% dir_output: Directory where output variables will be stored 
% dir_fastACI: Directory of the fastACI repository
% dir_fastACI_sim: Directory of the fastACI_sim (internal) repository is located
%
% 2. Example:
% %%% 2.1. Getting all the paths in the current computer:
%    paths = fastACI_paths;
%    fprintf('The current fastACI folder is: %s\n',paths.dir_fastACI);
%    fprintf('The current dir_output is: %s\n',paths.dir_output);
%
% %%% 2.2. Getting one specific folder:
%    dir_fastACI = fastACI_paths('dir_fastACI');
%    dir_output  = fastACI_paths('dir_output');
%    fprintf('The current fastACI folder is: %s\n',dir_fastACI);
%    fprintf('The current dir_output is: %s\n',dir_output);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paths.dir_data = fastACI_dir_data;
if ispc
    % Directories for Leo:
    switch version
        case {'9.7.0.1261785 (R2019b) Update 3', '9.11.0.1873467 (R2021b) Update 3'} % Leo's old office laptop
            paths.praat = 'C:\Praat\praatcon.exe';%'C:\Program files\Praat\praatcon.exe';
            paths.dir_output = [];
            paths.dir_output_fastACI2021_JASA = 'C:\Users\LeoVarnet\ownCloud\Professionnel\Publications et communications\2020-x - ANR fastACI\';
            paths.dir_output_fastACI2021_JASA_eps = 'C:\Users\LeoVarnet\ownCloud\Professionnel\Publications et communications\2020-x - ANR fastACI\';
        case '9.3.0.713579 (R2017b)' % Leo's home laptop
            paths.praat = 'C:\Users\Varnet L�o\Desktop\Programmes\praat.exe';%'C:\Program files\Praat\praatcon.exe';
            paths.dir_output = [];
            paths.dir_output_fastACI2021_JASA =  'C:\Users\Varnet L�o\Dropbox\Professionnel\Publications et communications\2020-x - ANR fastACI\2021 - JASA\';
            paths.dir_output_fastACI2021_JASA_eps = 'C:\Users\Varnet L�o\Dropbox\Professionnel\Publications et communications\2020-x - ANR fastACI\2021 - JASA\';
        case '9.10.0.1684407 (R2021a) Update 3' % Leo's new PC
            paths.praat = 'C:\Users\LSP005\Desktop\Programmes\praatcon5353_win64\praatcon.exe';
            paths.dir_output = [];
            warning('Path not set up.\n')
        otherwise
            paths.dir_output = [pwd filesep];
    end
    
elseif ismac
    
    % Directories in lab's computer:
    paths.dir_output = '/Users/leovarnet/ownCloud/Data/fastACI_data/outputs/';
    if ~exist(paths.dir_output,'dir')
        dir_output = [fastACI_basepath 'outputs' filesep];
        if ~exist(dir_output,'dir')
            mkdir(dir_output)
        end
        paths.dir_output = dir_output;
    end
    % warning('Leo: define and update the following folders...')
    paths.dir_output_fastACI2021_JASA = '';
    paths.dir_output_fastACI2021_JASA_eps = '';
    
elseif isunix
    % Directories for Alejandro:
    paths.praat       = '/usr/bin/praat';
    paths.dir_amtoolbox = fastACI_dir_amtoolbox;
    dir_up = userpath;
    if isempty(dir_up)
        % If empty:
        dir_up = [pwd filesep];
    end
    if strcmp(dir_up(end),':')
        dir_up = dir_up(1:end-1);
    end
    paths.dir_output  = [dir_up filesep 'outputs' filesep]; % '/home/alejandro/Documents/MATLAB/outputs/'
    
    paths.dir_output_fastACI2021_JASA     = '/home/alejandro/Documents/Databases/data/Osses-Varnet-2021-JASA/';
    paths.dir_output_fastACI2021_JASA_eps = '/home/alejandro/Documents/Documenten-ENS/01-Text/05-Doc/pr2021-05-20-ideas-4-paper/Figures-new/';
end

paths.dir_fastACI     = fastACI_basepath;
if exist('fastACIsim_basepath','file')
    paths.dir_fastACI_sim = fastACIsim_basepath;
end

if (nargin==1)
    if (~isfield(paths,type))
        error('Invalid type');
    end
    paths=paths.(type);
end