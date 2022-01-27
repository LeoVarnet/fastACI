function [CI, results, cfg] = CI_glmqpoptim_fct( cfg, y, y_correct_NOTUSED, X, U )
% function [CI, results, cfg] = CI_glmqpoptim_fct( cfg, y, y_correct_NOTUSED, X, U )
%
% Calculating GLM with smooth prior.
%
% Input parameters:
%       cfg: Configuration file
%       y
%       ~: the third input is not used for compatibility with other fitting functions
%       X: Contains the T-F samples of the test waveforms: dimensions N_intervals x N_t x N_f
%       U: Two columns vector. y_correct (0 or 1) is contained in the first column,
%          only ones are contained in column 2. The explanation of this variable
%          may be: Given that the distribution here is 'binomial', U is a
%          binary vector indicating success/failure. The second column contains
%          the total number of trials per observation (so always 1).
%
% Output parameters:
%       results: is the output of 'cvglmfitqp'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%% Calculating GLM with smoothness prior %%\n');

if ~isfield(cfg,'N_folds')
    cfg.N_folds=5;
    disp('Assigning the field ''N_folds'' to the cfg structure...')
end

prior_type = cfg.prior;
N_folds    = cfg.N_folds;

bUse_as_received = 0; % 0 => is using Alejandro''s recoding

switch prior_type
    case 'weight decay'
        warning('%s: Not debugged by Alejandro yet...',upper(mfilename))
        qf = speye(size(X,2)+size(U,2)); % sparse identity
        Q = tailmdx(qf,cfg.lambda_1);
    case 'none'
        warning('%s: Not debugged by Alejandro yet...',upper(mfilename))
        qf = 0*speye(size(X,2)+size(U,2));
        Q=qf;
    case 'smoothness'
        N_t = length(cfg.t); % Number of time samples
        N_f = length(cfg.f); % Number of frequencies
        
        % This is some sort of initialisation, independent of the number of
        % trials, but dependent on the number of T-F samples. It takes about 30 s
        if bUse_as_received == 1
            tic
            qf_ref = il_qfsmooth_as_received(cfg,X,U);
            toc
            
            qf = il_qfsmooth(N_t,N_f,U);

            figure; plot(qf_ref(:)-qf(:));
        end
        if bUse_as_received == 0
            % tic
            qf = il_qfsmooth(N_t,N_f,U); 
            % toc    
        end
    otherwise
        error('The specified prior does not exist');
end

folds = getcvfolds(length(y),N_folds);
cfg.folds = folds;
% % if N_folds = 10, and length(y) = 1000:
%  su = sum(folds,2); % is 9    for all y(1,:)
%  su = sum(folds,1); % is 900  for all y(:,1)
% %  It seems that if y(i,j) == 0 then it sample the index y(i,j) is passed 
% %     to the fold 'j'.

opts.family = 'binomlogit';
opts.getH = 0;
if isfield(cfg,'lambda0')
    opts.lambda0 = cfg.lambda0;
end
if isfield(cfg,'stepsize')
    opts.stepsize = cfg.stepsize;
end
if isfield(cfg,'lambdamax')
    opts.lambdamax = cfg.lambdamax;
end
if isfield(cfg,'precision')
    opts.precision = cfg.precision;
end
if isfield(cfg,'maxiter')
    opts.maxiter = cfg.maxiter;
end

[results,evaluation] = cvglmfitqp(y,[X,U],qf,folds, opts);
% My understanding:
%   qf: is the priors
%   y: should be the vector with responses equal 1 or 0
%   X: should be the T-F vector 
%   U: Number of trials for each X_i (so it is always equal to 1, as X does not represent an histogram)
results.evaluation = evaluation;

CI = reshape( results.w(1:end-size(U,2)), length(cfg.f), length(cfg.t));

if isfield(cfg,'filtrage')
    switch cfg.filtrage % Filtrage si demande
        case {1,'yes','oui'}
            error('Not validated yet...')
            % [ CI ] = lowpass_CI( CI, cfg.t, cfg.f, cfg.flt_freqcoup, cfg.flt_quefrcoup );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qf = il_qfsmooth(N_t,N_f,U)

[~,~,qf1,qf2] = qfsmooth(N_t, N_f);
N_qfs = size(qf1)+size(U,2); % AO: seems to be something like degrees of freedom?
qf = zeros(N_qfs); % memory allocation

N = N_t*N_f;
idxs = 1:N;
% qf0 = 0*speye(size(qf1)+size(U,2));
qf(idxs,idxs) = qf1 + qf2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qf = il_qfsmooth_as_received(cfg,X,U)

[~,~,qf1,qf2] = qfsmooth(length(cfg.t), length(cfg.f));
qfs = zeros([size(qf1)+size(U,2),2]);

qfs(1:size(X,2),1:size(X,2),1) = qf1;
qfs(1:size(X,2),1:size(X,2),2) = qf2; 
qf0 = 0*speye(size(qf1)+size(U,2));
% clear qf1 qf2
qf = qfs(:,:,1) + qfs(:,:,2);
% clear qfs
