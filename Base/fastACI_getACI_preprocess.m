function [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix)
% function [y, y_correct, X, U, cfg_ACI] = fastACI_getACI_preprocess(cfg_ACI, data_passation, Data_matrix)
%
% Changing the parameter names:
% New name      Old name        Changed on:
%
% Old name: Script4_Calcul_ACI_preprocess.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kv = cfg_ACI.keyvals; % List of configuration parameters defined in arg_fastACI_getACI.m
fg = cfg_ACI.flags; % List of modules that are activated or bypassed defined in arg_fastACI_getACI.m

do_permutation = fg.do_permutation; % By default the permutation test is 'on'
do_no_bias     = fg.do_no_bias;

idx_trialselect = cfg_ACI.idx_trialselect;
if length(idx_trialselect) > length(data_passation.n_responses)
    fprintf('\t%s: incomplete data_passation structure, not all trials seem to have been completed\n',upper(mfilename));
    cfg_ACI.idx_trialselect = cfg_ACI.idx_trialselect(1:length(data_passation.n_responses));
    idx_trialselect = cfg_ACI.idx_trialselect;
    
    cfg_ACI.N_trialselect = length(data_passation.n_responses);
end
if ~isfield(cfg_ACI,'N_trialselect')
    cfg_ACI.N_trialselect = length(data_passation.n_responses);
end
n_responses = data_passation.n_responses(idx_trialselect); % responses given by the participant (trial order)
n_targets   = data_passation.n_targets(idx_trialselect); % expected responses (trial order)
expvar      = data_passation.expvar(idx_trialselect); % value of the experimental variable (trial order)
is_correct  = data_passation.is_correct(idx_trialselect);

if ~isempty(kv.expvar_limits)
    cfg_ACI.SNR_analysis  = kv.expvar_limits; % select a range of SNR
    label_expvar_limits = sprintf('\t\t\t(expvar_limits between %.1f and %.1f)\n',cfg_ACI.SNR_analysis);
else
    cfg_ACI.SNR_analysis  = [min(expvar) max(expvar)];
    label_expvar_limits = [];
end

if do_permutation
    N_perm = kv.N_perm; 
end 

%cfg_ACI.trialtype_analysis = [];
switch kv.trialtype_analysis
    case 'incorrect'
        select_trialtype = ~is_correct;
    case 'correct'
        select_trialtype = is_correct;
    case 'total'
        select_trialtype = ones(size(is_correct));
    otherwise
        if kv.trialtype_analysis(1) == 't' 
            if str2double(kv.trialtype_analysis(2:end))<=cfg_ACI.N_target
                select_trialtype = (n_targets == str2double(kv.trialtype_analysis(2:end)));
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

% if ~isempty(kv.idx_trialselect)
%     idx_set_zero = setxor(1:length(select_n_trials), kv.idx_trialselect);
%     select_n_trials(idx_set_zero) = 0;
% end

expvar_trialselect = (expvar>=cfg_ACI.SNR_analysis(1) & expvar<=cfg_ACI.SNR_analysis(2));

select_after_reversal = ones(size(idx_trialselect)); % It will be overwritten if kv.expvar_after_reversal
if isfield(kv,'expvar_after_reversal')
    if kv.expvar_after_reversal > 0
        [select_after_reversal,label_expvar_after_reversal] = expvar_after_reversal(data_passation,kv.expvar_after_reversal,idx_trialselect);
    else
        label_expvar_after_reversal = ''; % empty label
    end
end

%%% Here is where all the filters are applied:
idx_analysis = find(select_n_trials & expvar_trialselect & select_trialtype & ...
    select_after_reversal); % & select_n_signal);

if do_no_bias
    %%% do_no_bias=TRIAL EQUALISATION
    [idx_analysis,label_expvar_no_bias] = expvar_no_bias(data_passation,idx_analysis);
else
    label_expvar_no_bias = ''; % empty label
end
%%% END do_no_bias

if length(idx_analysis) ~= size(Data_matrix,1)
    
    label_actual_expvar = sprintf('\t\t\t(Actual expvar values between %.1f and %.1f)\n', ...
        min(expvar(idx_analysis)),max(expvar(idx_analysis)));
    
    fprintf('\t%s: Selecting a subset of the experimental trials:\n',upper(mfilename));
    fprintf('\t\t %.0f trials out of %.0f are being processed:\n%s%s%s%s\n', ...
        length(idx_analysis),size(Data_matrix,1), ...
        label_expvar_limits, ...
        label_expvar_after_reversal, ...
        label_expvar_no_bias, ...
        label_actual_expvar);
    
    N_trialselect = length(idx_analysis);
    cfg_ACI.N_trialselect = N_trialselect;
end
cfg_ACI.idx_analysis = idx_analysis;

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
    
switch fg.glmfct
    case {'lassoglm','lasso', 'l1lm','l1glm'}
        
        %%% Loading defaults for 'lassoglm' or 'lasso', if not previously loaded
        cfg_ACI   = arg_glmfct(cfg_ACI,fg,kv);
        %%%
        Nlevel    = cfg_ACI.lasso_Nlevel;
        Nlevelmin = cfg_ACI.lasso_Nlevelmin;
        
        preX = Data_matrix(idx_analysis,:,:);
        %% Create Gaussian pyramid (modified_impyramid)
        %computes a Gaussian pyramid reduction of preX. 

        % The following step ensures that Nf and Nt have M*2^(Nlevel) elements (with M
        % integer) by discarding samples (or adding dummy samples if needed). This
        % is mandatory for accurate reconstruction of the pyramid.

        Nt_input = size(Data_matrix,3);
        
        Nt_X = 2*2^5; %256; % 2^7;%
        if Nt_input < Nt_X
            % Nothing to do
        else
            while Nt_X < Nt_input
                Nt_X = Nt_X*2;
            end
        end
        
        Nf_X = 2*2^5;%= 128;%2^7;%
        t_X = cfg_ACI.t;
        f_X = cfg_ACI.f;

        if ~isfield(cfg_ACI,'N_t')
            cfg_ACI.N_t = length(cfg_ACI.t);
        end
        if ~isfield(cfg_ACI,'N_f')
            cfg_ACI.N_f = length(cfg_ACI.f);
        end
        N_t = cfg_ACI.N_t;    
        N_f = cfg_ACI.N_f;
        
        switch kv.pyramid_script
            case 'imgaussfilt'
                    % Nothing to do
            case 'imresize'
                % The use of this function requires that the time and frequency
                % dimensions are multiples of a poiwer of 2, so, zero padding
                % is needed:
                
                %%% Time dimension:
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

                %%% Frequency dimension:
                if Nf_X<N_f % Along frequency dimension, truncating
                    preX = preX(:,1:Nf_X,:);
                    f_X = cfg_ACI.f(1:Nf_X);
                elseif Nf_X == N_f
                    % Nothing to do
                else
                    warning('Choose a higher value for NFFT')
                end
        end
                
        %%% Gaussian pyramid reduction
        Pyramid = {preX}; % structure containing the different levels of the pyramid

        % Filters the matrix level by level into a Nlevel Gaussian pyramid

        if Nlevelmin <= 2
            i_level_start = max(Nlevelmin,2);
        else
            warning('Nlevelmin > 2, this option has not been comprehensively tested yet');
            i_level_start = 2;
            
            % TODO:
            %     Nlevelmin >= 3 could be useful in the long run if we want 
            %     to use ACIs with higher resolution but the same size of 
            %     gaussian elements in the Pyramid basis. But with the 
            %     current resolution Nlevelmin should be 1 or 2
        end
        
        for i_level = i_level_start:Nlevel
            switch kv.pyramid_script
                case 'imgaussfilt'
                    kv.i_level = i_level; % only used if kv.pyramid_script is 'imgaussfilt'
                    Pyramid{i_level} = fastACI_impyramid(Pyramid{1}, 'reduce',kv);
                otherwise
                    Pyramid{i_level} = fastACI_impyramid(Pyramid{i_level-1}, 'reduce',kv);
            end
        end

        %% Reconstructing the Pyramid 
        % Refer to Script4_Calcul_ACI_calculate.m ['%%% 2. Getting the 
        %    'expanded' ACIs per level (Expand the weight matrix) %%%] to
        %    see how to reconstruct the 'Pyramid'. In the _calculate.m 
        %    script, the resulting ACI is obtained from the expansion of the
        %    resulting pyramid.
           
        X = [];
        Pyra_size = fastACI_impyramid_get_size(Pyramid,Nlevelmin,Nlevel);
        
        for i_level = Nlevelmin:Nlevel
            Pyra_here = squeeze(Pyramid{i_level});
            X = cat(2, X, reshape(Pyra_here,N_trialselect,[])); % along Dim=2
        end
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
    
if cfg_ACI.zscore
    %%%% Added by Leo on 28/04/2021:
    % Here each row of X (along dimension 1, each trial) is 
    %     normalised to have a mean of 0 and a std of 1.
    X = zscore(X,[],1);
    % U = zscore(U,[],1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inline functions:
