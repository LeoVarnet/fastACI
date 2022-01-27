function results = glmfitqp(y,X,qf,opts)
% function [thefit] = glmfitqp(y,X,qf,opts)
%
% Optimizes a GLM model through MAP where the internal variable is generated through:
% r = X*w
%
% Where there is a quadratic penalty on the model parameters:
%
% p(w) = 1/2*w'*qf*w
%
% opts: a struct with the following props:
%  family: The GLM model type. Acceptable values:
%          normid     - normal w/identity inverse link (least-squares)
%          binomlogit - binomial w/ logistic inverse link (default)
%          poissexp   - Poisson w/ exponential inverse link
%  familyextra: extra parameter depending on family type
%      for normid: equals sigma, the standard deviation of the noise
%                  (default 1)
%      for binomlogit: equals ntrials, the number of trials
%                  (default 1)
%  w0: if left blank will be initialized to random values
%  baseline: adds a constant offset vector to r
%  getH: returns the hessian H matrix if = 1 (default 0)
%        Note that the Hessian matrix is the inverse of the covariance
%        matrix of the estimated w; error bars can be derived from it
%  algo: Variant of fitting algorithm to use.
%        choices: 'medium' - Medium-scale fminunc
%                 'large'  - Large-scale fminunc
%                 'lbfgs' - Limited memory BFGS - requires minFunc.m
%                  http://www.di.ens.fr/~mschmidt/Software/minFunc.html
%                 
%  Display: how much information to spit out, either off, final or iter
%           (iter is default)
%  weights: weights to give each observation. Default 1 for each obs.
% 
% returns:
%  thefit: a struct with the optimal w, total log likelihood, log penalty,
%          H if requested, log likelihood of each observation
%
% Example use:
%
% %Smoothness prior
% qf = blkdiag(qfsmooth1D(16),.01);
% rg = (-7.5:7.5)';
% w = exp(-rg.^2/3^2).*sin(rg*2*pi/6);
% 
% %Simulate some data w/ gaussian noise 
% nobs = 30;
% X = [randn(nobs,length(w)),ones(nobs,1)];
% y0 = X*[w;.01] + .1*randn(nobs,1);
% 
% opts.family = 'normid';
% results = glmfitqp(y0,X,.5*qf,opts);
%
% Algorith: uses minFunc by Mark Schmidt (2005)
%           patched to avoid Cholesky decomposition of Hessian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    opts = struct;
end

if size(X,2) ~= size(qf,2)
    error('The size of X is not compatible with the size of qf');
end

defaults.family   = 'binomlogit';
defaults.familyextra = 1;
defaults.w0     = [];
defaults.baseline = zeros(length(y),1,class(y));
defaults.getH = 0;
defaults.algo = 'large';
defaults.Display = 'Iter';
defaults.weights = 1;

opts = setdefaults(opts,defaults);

%Check if link is acceptable
goodlink = any(strcmp(opts.family,{'normid','binomlogit','poissexp'}));

if ~goodlink
    error('Family %s is not supported', opts.family);
end

%Check inputs
if strcmp(opts.family,'binomlogit')
    if any(y ~= round(y))
        error('y can only take integer values for family binomlogit; change opts.familyextra for multiple trials');
    elseif any(y < 0 | y > opts.familyextra)
        error('y must be >= 0 and <= opts.familyextra for family binomlogit');
    end
elseif strcmp(opts.family,'poissexp')
    if any(y ~= round(y))
        error('y can only take integer values for family poissexp');
    elseif any(y < 0)
        error('y must take positive values for family poissexp');
    end
end

if numel(opts.weights) ~= 1 && any(size(opts.weights)~=size(y))
    error('Weights and y incompatible sizes');
end

if numel(opts.baseline) ~= 1 && any(size(opts.baseline)~=size(y))
    error('Baseline and y incompatible sizes');
end

if isempty(opts.w0) || strcmp(opts.family,'normid')
    %Get an initial estimate of w
    yd = y - opts.baseline;
    try
        w = [X;qf*opts.familyextra^2]\[opts.weights.*yd;zeros(size(X,2),1)];
    catch me
        if strcmp(opts.family,'normid')
            error('Could not solve system through \');
        end
        w = zeros(size(X,2),1);
    end
else
    w = opts.w0;
end

if ~strcmp(opts.family,'normid')
    % Refine, adjusting w:
    w = il_irls(y,X,w,qf,opts.baseline,opts.family,opts.familyextra,opts.algo,opts.Display,opts.weights);
end

l0 = evalMaxGlmLikelihood(y,opts.family,opts.familyextra,opts.weights);

if opts.getH
    [ll,~,Hd,~,lleach] = evalGlmLikelihood(y,X,w,opts.baseline,opts.family,opts.familyextra,opts.weights);
    Y = bsxfun(@times,sqrt(Hd),X);
    H = Y'*Y;
    results.H = H + qf;
else
    [ll,~,~,~,lleach] = evalGlmLikelihood(y,X,w,opts.baseline,opts.family,opts.familyextra,opts.weights);
end

ll = ll - l0;

lp = .5*w'*qf*w;

results.w = w;
results.loglikelihood = -ll;
results.logpenalty    = -lp;
results.loglikelihoodeach = -lleach;
results.opts = opts;
results.eval_script_call = '[ll,~,~,~,lleach] = evalGlmLikelihood(y,X,w,opts.baseline,opts.family,opts.familyextra,opts.weights);';

%%% End of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = il_irls(y,X,w,qf,b,family,familyextra,algo,Display,weights)

mat_ver = version('-release');
mat_ver = str2double(mat_ver(1:end-1)); % Added by AO, for compatibility

opts.Display = Display;
switch algo
    case 'medium'
        opts.LargeScale = 'off';
        opts.GradObj = 'on';
        w = fminunc(@(w) il_evalLplusp(y,X,w,b,qf,family,familyextra,weights),double(w),opts);
    case 'large'
        opts.LargeScale = 'on';
        opts.GradObj = 'on';
        opts.Hessian = 'on'; 
    
        if mat_ver >= 2019
            w = fminunc_local(@(w) il_evalLplusp(y,X,w,b,qf,family,familyextra,weights),double(w),opts);
        else
            w = fminunc(@(w) il_evalLplusp(y,X,w,b,qf,family,familyextra,weights),double(w),opts);
        end
        
    case 'lbfgs'
        opts.Method = 'lbfgs';
        w = minFunc(@(w) il_evalLplusp(y,X,w,b,qf,family,familyextra,weights),double(w),opts);
    case 'newton'
        opts.Method = 'mnewton';
        opts.HessianIter = 3;
        opts.HessianModify = 6; %Requires patching minFunc.m
        w = minFunc(@(w) il_evalLplusp(y,X,w,b,qf,family,familyextra,weights),double(w),opts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,g,H] = il_evalLplusp(y,X,w,b,qf,family,familyextra,weights)

if nargout > 2
    [f0,g0,Hd] = evalGlmLikelihood(y,X,w,b,family,familyextra,weights);
    dw = qf*w;
    f = double(f0+.5*w'*dw);
    g = double(g0+ dw);

    Y = bsxfun(@times,sqrt(Hd),X);
    H = Y'*Y + qf;
elseif nargout > 1
    [f0,g0] = evalGlmLikelihood(y,X,w,b,family,familyextra,weights);
    dw = qf*w;
    f = double(f0+.5*w'*dw);
    g = double(g0+ dw);
else
    error('Not validated yet...')
    f0 = evalGlmLikelihood(y,X,w,b,family,familyextra);
    f = f0+.5*w'*qf*w;
end
