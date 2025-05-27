function results = crossValidate(y,y_idxs,lambda0,fitfun,valfun,folds,opts)
% function results = crossValidate(y,y_idxs,lambda0,fitfun,valfun,folds,opts)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults.precision = .2;
defaults.stepsize = 10;
defaults.minDiffSecondRound = 10;
defaults.maxiter  = 15;
defaults.parallel = 0;
defaults.oldfit = [];
defaults.Display = 'iter';

kfolds  = size(folds,2);
opts    = setdefaults(opts,defaults,true);
lambdas = lambda0*[1/sqrt(opts.stepsize),sqrt(opts.stepsize)];
fits    = cell(kfolds,opts.maxiter);
maxiter = opts.maxiter;

if strcmp(opts.Display,'off')
    printfun = @il_nullfun;
else
    printfun = @fprintf;
end

printfun('\nStarting cross validation\n\n');
printfun('                         CV error (smaller is better)\n');
printfun('    lambda        ');

for ii = 1:kfolds
    printfun('fold %d   ', ii);
end

printfun('    total\n');

cvgofs = zeros(kfolds,maxiter);
extrainfos = fits;

jj = 0;
while jj < maxiter
    jj = jj + 1;

    lambda = lambdas(jj);
    % LEO: I removed this arbitrary parameter
%     if lambda > 1e6 || lambda < 1e-6
%         jj = jj - 1;
%         break; %too big | too small
%     end

    printfun('  %10.2g    ',lambda);

    if ~opts.parallel
        for ii = 1:kfolds
            if jj == 1
                oldfit = opts.oldfit;
            else
                oldfit = fits{ii,jj-1};
            end
            [fits{ii,jj},gof,extrainfo] = il_fit_a_fold(y,y_idxs,ii,folds,fitfun,valfun,oldfit,lambda);
            cvgofs(ii,jj) = gof;
            extrainfos{ii,jj} = extrainfo;

            %Dispay during%
            printfun('%8.1f ',cvgofs(ii,jj));
        end
    else
        if jj == 1
            oldfits = cell(kfolds,1);
            for ii = 1:kfolds
                oldfits{ii} = opts.oldfit;
            end
        else
            oldfits = fits(:,jj-1);
        end

        parfor ii = 1:kfolds
            [fits{ii,jj},gof,extrainfo] = il_fit_a_fold(y,y_idxs,ii,folds,fitfun,valfun,oldfits{ii},lambda);
            cvgofs(ii,jj) = gof;
            extrainfos{ii,jj} = extrainfo;
        end

        %Display after if parallel
        for ii = 1:kfolds
            printfun('%8.1f ',cvgofs(ii,jj));
        end
    end

    [sortedl,sortidx] = sort(lambdas(1:jj));
    sortedcvs = sum(cvgofs(:,sortidx),1);
    [themax,maxidx] = min(sortedcvs);

    printfun('    %8.2f\n',sum(cvgofs(:,jj),1));

    if jj > 1
        if maxidx == jj
            nextlambda = max(sortedl)*opts.stepsize;
        elseif maxidx == 1
            nextlambda = min(sortedl)/opts.stepsize;
        else

            xs = sortedl(  maxidx-1:maxidx+1);
            ys = sortedcvs(maxidx-1:maxidx+1);

            if(xs(3)/xs(1) > opts.minDiffSecondRound)
                if diff(log(xs(1:2))) > diff(log(xs(2:3)))
                    nextlambda = exp(mean(log(xs(1:2))));
                else
                    nextlambda = exp(mean(log(xs(2:3))));
                end
            else
                nextlambda = il_giveGuess(xs(:),ys(:));
            end

            if nextlambda < xs(1)
                nextlambda = exp(mean(log(xs(1:2))));
            elseif nextlambda > xs(3)
                nextlambda = exp(mean(log(xs(2:3))));
            end

            if any(abs(nextlambda./sortedl-1) < opts.precision) 
                lambda = sortedl(maxidx);
                break;
            end
        end
        lambdas(jj+1) = nextlambda;
        lambda = nextlambda;
    end


end

results.finallambda = lambda;
results.lambdas = lambdas(1:jj);
results.fits = fits(:,1:jj);
results.cvgofs = cvgofs(:,1:jj);
results.cvgof  = themax;
results.extrainfos = extrainfos;


printfun('\nPerforming final fit with lambda = %8.2f\n',lambda);

results.finalfit = fitfun(y,y_idxs,fits{1,jj},lambda,true);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_nullfun(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g] = il_giveGuess(x,y)
%Assume that y ~ a*x^2 + b*x + c, solve for a, b, c, return -b/2a
x = log(x);
H = [x.^2,x,ones(size(x))];
v = [H;1e-6*eye(3)]\[y;zeros(3,1)];
g = exp(-v(2)/2/v(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thefit,gof,extrainfo] = il_fit_a_fold(y,X,ii,folds,fitfun,valfun,oldfit,lambda)

fitset = folds(:,ii);
valset = ~fitset;

thefit = fitfun(y(fitset),X(fitset,:),oldfit,lambda,false);
try
    [gof,extrainfo]    = valfun(y(valset),X(valset,:),thefit);
catch me
    error('Error here...')
    gof = valfun(y(valset),X(valset,:),thefit);
    extrainfo = [];
end
