function cfg_inout = modulationACI_init(cfg_inout)
% function cfg_crea = modulationACI_init(cfg_crea)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_subject = [cfg_inout.dir_stim cfg_inout.Subject_ID filesep];

if ~exist(dir_subject,'dir')
    mkdir(dir_subject)
    fprintf('Directory %s was created\n',dir_subject);
end
dir_where  = [dir_subject cfg_inout.folder_name filesep]; 
% dir_where = 'C:\Users\Varnet Leo\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM\'; % dir at Leo's
        
fs         = cfg_inout.fs;
dur        = cfg_inout.stim_dur;
noise_type = cfg_inout.noise_type;
N_presentation  = cfg_inout.N_presentation;
N_target   = cfg_inout.N_target;
N_total    = N_presentation*N_target; 

disp('The new stimuli will be generated using the following parameters:')
fprintf('\tNoise type=%s\n',noise_type);
fprintf('\tDuration=%.2f s\n',dur);
fprintf('\tSampling frequency=%.2f Hz\n',fs);
fprintf('\tA total of %.0f stimuli will be generated: %.0f ''noises'', %.0f ''conditions''\n',N_total,N_presentation,N_target);
fprintf('\tTarget folder: %s\n',dir_where);
disp('Press any button to continue, press ctrl+C to abort')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
pause()

if ~isdir(dir_where)
    mkdir(dir_where);
    bGenerate_stimuli = 1;
elseif ~isempty(dir(dir_where))
    warning('%s: Folder %s is not empty.',upper(mfilename),dir_where);
    bInit_participant_only = input('Do you to initialise a participant without generating new sounds? (1=yes; 0=no): ');
    bGenerate_stimuli = ~bInit_participant_only;

    if bGenerate_stimuli
        error('%s: Folder %s is not empty. Choose a new ''folder_name'' or remove manually the existing content in the target folder',upper(mfilename),dir_where);
    end
end

% Initialises the seed numbers on the fly:
seed_number_max = 4*cfg_inout.N; % arbitrary number, seeds will range from 0 to 4*N

method = 'perm';
% method = 'unif';

s_current = rng; % gets current seed

type_seed = 'suffle';
try
    rng(type_seed); % seed based on the computer's clock
catch me
    rng('default');
    rng(type_seed); % seed based on the computer's clock
end

switch method
    case 'unif'
        error('Not tested recently')
    case 'perm'
        numbers = randperm(seed_number_max);
        seed_numbers = numbers(1:cfg_inout.N); % takes only the first 'N' numbers
end

cfg_inout.seeds_order = seed_numbers; % to be used sequentially
cfg_inout.seeds_order_method = method;

if cfg_inout.randorder == 1
    cfg_inout.stim_order = randperm(cfg_inout.N); 
else
    error('Not tested recently')
end

%%% Seed set back
rng(s_current);
%%% 

%% Noise samples generation
for i=1:N_total
    clc
    if bGenerate_stimuli
        fprintf('Creating noise stimulus # %.0f of %.0f\n',i,N_total);

        %%% Fixing the seed
        str_stim = [];
        str_inout = []; 
        % str_inout.istarget = 0;
        str_inout.expvar(i)   = 0; % idle
        str_inout.i_current = i; % idle
        str_inout.n_stim(i) = i;
        if ~isfield(cfg_inout,'n_targets_sorted');
            bRemove_after = 1;
            cfg_inout.n_targets_sorted(i) = 1; % we generate noises only, not important whether is actually target or not
        else
            bRemove_after = 0;
        end
        
        s_current = rng; % gets current seed
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number,type_seed); % N_samples = round(cfg_inout.stim_dur * fs); noise = randn(N_samples,1); rng(seed_number,type_seed)
        
        str2eval = sprintf('str_stim=%s_user(cfg_inout,str_inout);',cfg_inout.experiment);
        eval(str2eval);
        if bRemove_after
            cfg_inout = rmfield(cfg_inout,'n_targets_sorted');
        end
        
        noise_cal = str_stim.stim_noise_alone;
    end
    % Ensuring that the noise iteration number has four characters:
    stimnumber=num2str(i);
    if i<1000; stimnumber=['0' stimnumber]; end
    if i<100;  stimnumber=['0' stimnumber]; end
    if i<10;   stimnumber=['0' stimnumber]; end
    % Ensuring 'End'

    fname = ['Noise_' stimnumber '.wav'];
    ListStim(i).name = fname;
    if bGenerate_stimuli
        audiowrite([dir_where fname], noise_cal, fs);
    end
end
cfg_inout.ListStim = ListStim;