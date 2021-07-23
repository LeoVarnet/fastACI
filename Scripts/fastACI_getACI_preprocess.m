function [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix)
% function [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix)
%
% Changing the parameter names:
% New name      Old name        Changed on:
%
% Old name: Script4_Calcul_ACI_preprocess.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_trialselect = cfg_ACI.idx_trialselect;
if length(idx_trialselect) > length(data_passation.n_responses)
    fprintf('\t%s: incomplete data_passation structure, not all trials seem to have been completed\n',upper(mfilename));
    cfg_ACI.idx_trialselect = cfg_ACI.idx_trialselect(1:length(data_passation.n_responses));
    idx_trialselect = cfg_ACI.idx_trialselect;
    
    cfg_ACI.N_trialselect = length(data_passation.n_responses);
end
n_responses = data_passation.n_responses(idx_trialselect); % responses given by the participant (trial order)
n_targets   = data_passation.n_targets(idx_trialselect); % expected responses (trial order)
expvar      = data_passation.expvar(idx_trialselect); % value of the experimental variable (trial order)
is_correct  = data_passation.is_correct(idx_trialselect);
% Selection of the condition
% cfg_ACI.NameCond = Condition;
%cfg_ACI.n_signal_analysis = []; % perform analysis only on trials corresponding to specific signals numbers
if ~isempty(cfg_ACI.keyvals.expvar_limits)
    cfg_ACI.SNR_analysis  = cfg_ACI.keyvals.expvar_limits; % select a range of SNR
else
    cfg_ACI.SNR_analysis  = [min(expvar) max(expvar)];
end
%cfg_ACI.n_trials_analysis = []; % select a range of trial numbers

do_permutation = cfg_ACI.flags.do_permutation; % By default the permutation test is 'on'
if do_permutation
    N_perm = cfg_ACI.keyvals.N_perm; 
end 

%cfg_ACI.trialtype_analysis = [];
switch cfg_ACI.keyvals.trialtype_analysis
    case 'incorrect'
        select_trialtype = ~is_correct;
    case 'correct'
        select_trialtype = is_correct;
    case 'total'
        select_trialtype = ones(size(is_correct));
    otherwise
        if cfg_ACI.keyvals.trialtype_analysis(1) == 't' 
            if str2double(cfg_ACI.keyvals.trialtype_analysis(2:end))<=cfg_ACI.N_target
                select_trialtype = (n_targets == str2double(cfg_ACI.keyvals.trialtype_analysis(2:end)));
            else
                error(['Trialtype condition unrecognised: ' cfg_ACI.trialtype_analysis ' (but there are only ' num2str(cfg_ACI.N_target) ' targets)\n'])
            end
        else
            error(['Trialtype condition unrecognised: ' cfg_ACI.trialtype_analysis '\n'])
        end
end

cfg_ACI = set_default_cfg(cfg_ACI, ...
    'n_signal_analysis', 1:cfg_ACI.N_target, 'n_trials_analysis', [1 cfg_ACI.N_trialselect]);

N_trialselect = cfg_ACI.N_trialselect;
select_n_trials = (1:N_trialselect>=cfg_ACI.n_trials_analysis(1) & 1:N_trialselect<=cfg_ACI.n_trials_analysis(2));
 
expvar_trialselect = (expvar>=cfg_ACI.SNR_analysis(1) & expvar<=cfg_ACI.SNR_analysis(2));

%%% select_n_signal: to process target 1 or 2, or 1 and 2
% select_n_signal = zeros(1,N_trialselect);
% 
% for i=cfg_ACI.n_signal_analysis
%     select_n_signal = ((n_targets==i) | select_n_signal);
% end
idx_analysis = find(select_n_trials & expvar_trialselect & select_trialtype); % & select_n_signal);

if length(idx_analysis) ~= size(Data_matrix,1)
    fprintf('\t%s: Selecting a subset of the experimental trials:\n',upper(mfilename));
    fprintf('\t\t %.0f trials out of %.0f are being processed (expvar_limits between %.1f and %.1f)\n',length(idx_analysis),size(Data_matrix,1),cfg_ACI.SNR_analysis);
    
    N_trialselect = length(idx_analysis);
    cfg_ACI.N_trialselect = N_trialselect;
end

cfg_ACI.N_trials = length(idx_analysis);
y_all     = double((n_responses==1)'); % all trials that for which target 1 has been chosen
y         = y_all(idx_analysis);

n_targets_select = n_targets(idx_analysis); % idx = find(n_targets_select==0); n_targets_select(idx) = 1;
try
    y_correct = double((cfg_ACI.response_correct_target(n_targets_select)==1)');
catch
    idx = find(n_targets_select==0); 
    n_targets_select(idx) = 1;
    fprintf('\t%s: %.0f case(s) solved\n',upper(mfilename),length(idx));
    y_correct = double((cfg_ACI.response_correct_target(n_targets_select)==1)');
end
% y_correct indicates the trials where target '1' has been expected to 
%     be chosen or 'target 1 present' trials
    
%%% Start: Testing Alejandro on 27/04/2021
if do_permutation
    cfg_perm = [];
    cfg_perm.N_perm = N_perm;
    for i = 1:N_perm
        cfg_perm.idxs_perm(:,i) = transpose(randperm(cfg_ACI.N_trials));
        cfg_perm.y_perm(:,i)    = y(cfg_perm.idxs_perm(:,i)); % random answers
        n_targets_perm = n_targets_select(cfg_perm.idxs_perm(:,i));
        cfg_perm.y_correct_perm(:,i) = double((cfg_ACI.response_correct_target(n_targets_perm)==1)');

        n_responses_perm = data_passation.n_responses(idx_analysis);
        n_responses_perm = n_responses_perm(cfg_perm.idxs_perm(:,i));
        cfg_perm.perc_correct_if_perm(i) = sum(n_responses_perm == data_passation.n_response_correct_target(idx_analysis))/length(idx_analysis);
        if i == 1
            cfg_perm.perc_correct = sum(data_passation.is_correct(idx_analysis))/length(idx_analysis);
        end
    end
end
    
% y_correct_perm(:,i), X, U_perm_here

%%% End: Testing Alejandro on 27/04/2021
    
switch cfg_ACI.flags.glmfct
    case {'lassoglm','lasso'}
        
        preX = Data_matrix(idx_analysis,:,:);
        %% Create Gaussian pyramid (modified_impyramid)
        %computes a Gaussian pyramid reduction of preX. 

        Nlevel = 5; % number of levels (= degrees of filtering) in the Gaussian pyramid
        Nlevelmin = 2; % level minimum considered in the analysis
 
        % The following step ensures that Nf and Nt have M*2^(Nlevel) elements (with M
        % integer) by discarding samples (or adding dummy samples if needed). This
        % is mandatory for accurate reconstruction of the pyramid.

        Nt_input = size(Data_matrix,3);
        
        Nt_X = 2*2^5; %256; % 
        if Nt_input < 256
            % Nothing to do
        else
            while Nt_X < Nt_input
                Nt_X = Nt_X*2;
            end
        end
        
        Nf_X = 2*2^5;%= 128;%
        t_X = cfg_ACI.t;
        f_X = cfg_ACI.f;

        N_t = cfg_ACI.N_t;
        N_f = cfg_ACI.N_f;
        
        if Nt_X>N_t % Along time dimension, zero padding
            [N,M,P] = size(preX);
            nullMatrix = zeros(N,M,Nt_X-P);
            preX = cat(3,preX,nullMatrix);
            dt = diff(cfg_ACI.t(1:2));
            t_X = (1:Nt_X)*dt;
        end
        
        if Nt_X<N_t % Along time dimension, truncating
            preX = preX(:,:,1:Nt_X);
            t_X = cfg_ACI.t(1:Nt_X);
        end
        if Nf_X<N_f % Along frequency dimension, truncating
            preX = preX(:,1:Nf_X,:);
            f_X = cfg_ACI.f(1:Nf_X);
        else
            warning('Choose a higher value for NFFT')
        end
        
        % Gaussian pyramid reduction

        Pyramid = {preX}; % structure containing the different levels of the pyramid

        % Filters the matrix level by level into a Nlevel Gaussian pyramid

        for i_level = Nlevelmin:Nlevel
            Pyramid{i_level} = Script4_Calcul_ACI_modified_impyramid(Pyramid{i_level-1}, 'reduce'); 
        end

        %% Reconstructing the Pyramid 
        % Refer to Script4_Calcul_ACI_calculate.m ['%%% 2. Getting the 
        %    'expanded' ACIs per level (Expand the weight matrix) %%%] to
        %    see how to reconstruct the 'Pyramid'. In the _calculate.m 
        %    script, the resulting ACI is obtained from the expansion of the
        %    resulting pyramid.
           
        X = [];
        for i_level = Nlevelmin:Nlevel
            Pyra_here = squeeze(Pyramid{i_level});
            Pyra_size(i_level,:) = [size(Pyra_here,2) size(Pyra_here,3)];
            X = cat(2, X, reshape(Pyra_here,N_trialselect,[])); % along Dim=2
        end

        cfg_ACI.lasso_Nlevelmin = Nlevelmin;
        cfg_ACI.lasso_Nlevel    = Nlevel;
        cfg_ACI.lasso_Pyra_size = Pyra_size;
        
        cfg_ACI.t_X = t_X;
        cfg_ACI.f_X = f_X;
        
    otherwise 
        % Nothing to do
        X = Data_matrix(idx_analysis,:); % 
end

switch cfg_ACI.withU
    case {1,'yes','oui'}
        U = [y_correct ones(cfg_ACI.N_trials, 1)];
        if do_permutation
            for i = 1:N_perm
                cfg_perm.U_perm(:,:,i) = [cfg_perm.y_correct_perm(:,i) ones(cfg_ACI.N_trials, 1)]; %%% Testing Alejandro
            end
        end

    otherwise
        U = [];
        if do_permutation
            for i = 1:N_perm
                cfg_perm.U_perm(:,:,i) = []; %%% Testing Alejandro
            end
        end
end
cfg_ACI.N_trials_true  = sum(y==1);
cfg_ACI.N_trials_false = sum(y==0);

if do_permutation
    cfg_ACI.cfg_perm = cfg_perm;
end
    
% warning('Temporal')
%%% Commented by Leo on 28/04/2021
% if cfg_ACI.zscore
%     X=(X-mean(X(:)))/std(X(:)); % Again a normalisation...
% end
%%% End commented by Leo on 28/04/2021

if cfg_ACI.zscore
    %%%% Added by Leo on 28/04/2021:
    % Here each row of X (along dimension 1, each trial) is 
    %     normalised to have a mean of 0 and a std of 1.
    X = zscore(X,[],1);
    % U = zscore(U,[],1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inline functions: