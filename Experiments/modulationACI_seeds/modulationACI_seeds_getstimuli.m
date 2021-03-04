function cfg_out = modulationACI_seeds_getstimuli %(cfg_in)
% function str_stim = modulationACI_seeds_getstimuli(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subject_ID = 'SAOtest';
dir_where = '/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Interim_results/';
file = [dir_where 'cfgcrea_2021_03_03_16_19_' Subject_ID '_modulationACI_seeds.mat'];

dir_subject = [dir_where Subject_ID filesep];
dir_audio   = [dir_subject 'NoiseStims' filesep];

if ~isfolder(dir_subject)
    mkdir(dir_subject);
end
if ~isfolder(dir_audio)
    mkdir(dir_audio);
end

cfg_crea = [];
load(file,'cfg_crea');

N = cfg_crea.N;
N = input(['Enter the number of trial noises that you want to generate (max=' num2str(N) '): ']);

[cfg_crea.n_signals, cfg_crea.n_response_correct_target] = Get_n_signals(cfg_crea);
 
fs = cfg_crea.fs;

noise_all = [];

for i_current = 1:N
    
    n_stim = cfg_crea.stim_order(i_current); % copied from Script2, as of 3 March 2021
    
    data_passation.i_current = i_current;
    data_passation.n_stim(i_current) = n_stim;
    data_passation.expvar(i_current) = 1; % Idle variable, not required
    signals = modulationACI_seeds_user(cfg_crea,data_passation);
    
    noise = signals.stim_noise_alone;
    noise_all(:,end+1) = noise;
    
    stimnumber=num2str(i_current);
    if i_current<1000; stimnumber=['0' stimnumber]; end
    if i_current<100;  stimnumber=['0' stimnumber]; end
    if i_current<10;   stimnumber=['0' stimnumber]; end
    
    fname = [dir_audio 'Noise_' stimnumber '.wav'];
    if ~exist(fname,'file')
        audiowrite(fname,noise,fs,'BitsPerSample',24);
    else
        error('Directory %s is not empty, remove the wave files and re-run this script',dir_audio);
    end
    disp('')
end
save([dir_audio 'Noises_all.mat'],'noise_all','fs');

disp('')

% if nargin == 0
%     cfg_in = [];
% end
% 
% cfg = [];
% cfg_out = cfg_in; % copying input to output struct
% 
% cfg.response_names = {'pure tone', 'modulated tone'}; 
% cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation
% 
% bDebug = 0;
% cfg.displayN       = bDebug; % 'oui'
% cfg.feedback       = 1;
% 
% %cfg_game.end_sessions   = [500 1000 1500 2000 2500]; 
% cfg.sessionsN      = 500; % CAUTION: Overwritten in the case of simulation
% cfg.adapt          = 1; % 'out';%
% cfg.randorder      = 1;
% 
% cfg.startvar = -8;  % old name 'm_start'
% cfg.expvar_description = 'modulation depth (dB)';
% 
% cfg.maxvar = 0;
% 
% % Staircase algorithm parameters
% if cfg.adapt == 1
% 	cfg.start_stepsize     = 4;
% 	cfg.min_stepsize       = 1;
%     cfg.adapt_stepsize     = 90/100;
% else
%     error('Not validated yet...')
% end
% 
% cfg_out = Merge_structs(cfg,cfg_out);