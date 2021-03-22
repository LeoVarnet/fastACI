function cfg_inout = modulationACI_seeds_init(cfg_inout)
% function cfg_inout = modulationACI_seeds_init(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
    cfg_inout.N = 3000;
    warning('Test only...')
end

if ~isfield(cfg_inout,'randorder')
    warning('%s: the field randorder has not been given. It is assumed that the stimuli will be presented randomly...',upper(mfilename));
    cfg_inout.randorder = 1;
end
    
% Initialises the seed numbers on the fly:
seed_number_max = 4*cfg_inout.N; % arbitrary number, seeds will range from 0 to 4*N

method = 'perm';
% method = 'unif';

s_current = rng; % gets current seed
rng('default');
rng('shuffle'); % seed based on the computer's clock

switch method
    case 'unif'
        warning('This method has not been tested recently. It may be removed soon from this list of options...')
        seed_numbers = round(seed_number_max*random('unif',0,1,[1,cfg_inout.N]));
                                        % (allows repeted seed numbers)
    case 'perm'
        numbers = randperm(seed_number_max);
        seed_numbers = numbers(1:cfg_inout.N); % takes only the first 'N' numbers
end

cfg_inout.seeds_order = seed_numbers; % to be used sequentially
cfg_inout.seeds_order_method = method;

if cfg_inout.randorder == 1
    cfg_inout.stim_order = randperm(cfg_inout.N); 
else
    cfg_inout.stim_order = 1:cfg_inout.N; 
end

%%% Trials for validation:
% In each session, the waveforms of 1 trial are stored for the later validation
% of the reconstructed waveforms, after the experimental session.
N_sessions = cfg_inout.N/cfg_inout.sessionsN;
% values between 1 and sessionsN:
sessionN_trials_validation = round((cfg_inout.sessionsN-1)*random('unif',0,1,[1,N_sessions]))+1;
cfg_inout.sessionN_trials_validation = cfg_inout.sessionsN*[0:N_sessions-1] + sessionN_trials_validation;

%%% Seed set back
rng(s_current);
%%% 
