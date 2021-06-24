function noise_waveforms = modulationACI_seeds_getstimuli(file, dir_subject) %(cfg_in)
% function noise_waveforms = modulationACI_seeds_getstimuli(file, dir_where)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ==0
    Subject_ID = 'SAOtest';
    dir_where = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Interim_results/';
    file = [dir_where 'cfgcrea_2021_03_03_16_19_' Subject_ID '_modulationACI_seeds.mat'];
    
    dir_subject = [dir_where Subject_ID filesep];
end

if nargin < 2
    % dir_where = fileparts(file);
    % dir_where = [dir_where filesep];
    % % [path,name,ext]=fileparts(path);
    % disp('')
    error('Specify a dir_subject')
end

dir_audio   = [dir_subject 'NoiseStims' filesep];

if ~exist(dir_subject,'dir')
    mkdir(dir_subject);
end
if ~exist(dir_audio,'dir')
    mkdir(dir_audio);
end

var_type = strsplit(file,filesep);
var_type = var_type{end};
var_type = strsplit(var_type,'_');
var_type = var_type{1};

switch var_type
    case 'cfg_crea'
        var = load(file,'cfg_crea');
        
        cfg = var.cfg_crea;
    case 'savegame'
        var = load(file,'cfg_game');
        
        cfg = var.cfg_game;
end

N = cfg.N;
N = input(['Enter the number of trial noises that you want to generate (max=' num2str(N) '): ']);

[cfg.n_targets_sorted, cfg.n_response_correct_target] = Get_n_signals(cfg);
 
fs = cfg.fs;

noise_waveforms = [];

for i_current = 1:N
    
    stimnumber=num2str(i_current);
    if i_current<1000; stimnumber=['0' stimnumber]; end
    if i_current<100;  stimnumber=['0' stimnumber]; end
    if i_current<10;   stimnumber=['0' stimnumber]; end
    
    fname = [dir_audio 'Noise_' stimnumber '.wav'];
    if nargout == 0
        bDo = ~exist(fname,'file');
    end
    if nargout >= 0
        bDo = 1;
    end
    
    if bDo
        n_stim = cfg.stim_order(i_current); % copied from Script2, as of 3 March 2021
    
        data_passation.i_current = i_current;
        data_passation.n_stim(i_current) = n_stim;
        data_passation.expvar(i_current) = 1; % Idle variable, not required
        signals = modulationACI_seeds_user(cfg,data_passation);

        noise = signals.stim_noise_alone;
        if nargout >= 1
            noise_waveforms(:,end+1) = noise;
        end
        if nargout == 0
            audiowrite(fname,noise,fs); % ,'BitsPerSample',24);
        end
    else
        fprintf('\t%s. file %s found on disk\n',stimnumber,fname)
        % error('Directory %s is not empty, remove the wave files and re-run this script',dir_audio);
    end
    disp('')
end
% save([dir_audio 'Noises_all.mat'],'noise_all','fs');

disp('')