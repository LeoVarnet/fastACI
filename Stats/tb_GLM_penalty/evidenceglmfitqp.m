function r = evidenceglmfitqp(y,X,qfs,opts)
    %[thefit] = evidenceglmfitqp(y,X,qfs,opts)
    %
    % Optimizes a GLM model through MAP where the internal variable is generated through:
    % r = X*w
    %
    % Where there is a quadratic penalty on the model parameters:
    %
    % p(w) = 1/2*opts.qf0 + sum_i lambda(i)*w'*qfs(:,:,i)*w
    %
    % lambda are optimized by evidence optimization with a Laplace 
    % approximation 
    %
    % Reference: Bishop (2007), Pattern Recognition and Machine Learning, sections 3.5.1, 4.4.1 & 5.7.2
    %
    % Example use:
    % %%
    % %Figure out optimal strength of prior through evidence optimization
    % %Here the weight parameter is a 2d Gabor and the smoothness varies along
    % %the two dims
    % [xi,yi] = meshgrid(-7.5:7.5,-7.5:7.5);
    % 
    % %Filter is a 2D Gabor
    % w = sin(xi).*exp(-(xi.^2+yi.^2)/2/3^2);
    % w = w(:);
    % 
    % [~,~,qf1,qf2] = qfsmooth(16,16);
    % Q = zeros([size(qf1)+1,2]);
    % Q(1:end-1,1:end-1,1) = qf1;
    % Q(1:end-1,1:end-1,2) = qf2;
    % 
    % %Simulate a model with w = Gabor function
    % nobs = 800;
    % X = [randn(nobs,length(w)),ones(nobs,1)];
    % r = .5*X*[w;.01];
    % 
    % %output is binary -> logistic regression
    % r = binornd(1,1./(1+exp(-r)));
    % 
    % %Fit the data
    % clear opts
    % opts.family = 'binomlogit';
    % 
    % 
    % %Because there are two hyperparameters to optimize, use evidence framework
    % results = evidenceglmfitqp(r,X,Q,opts);
    % results0 = glmfitqp(r,X,speye(257)*.1,opts);
    % 
    % clf;
    % subplot(1,3,1);imagesc(reshape(w,16,16));title('Actual filter');
    % subplot(1,3,2);imagesc(reshape(results0.w(1:end-1),16,16));title('Maximum likelihood');
    % subplot(1,3,3);imagesc(reshape(results.w(1:end-1),16,16));title('Maximum A Posteriori w/ smoothness prior');
    %
    % See also: glmfitqp
    defaults.maxdelta = 2;
    defaults.evidenceeps = 1e-2;
    defaults.lambdaeps = 1e-2;
    defaults.weps = 1e-2;
    defaults.maxlambda = 1e10;
    defaults.minlambda = 1e-3;
    defaults.burnin = 1;
    defaults.maxiter  = 30;
    defaults.familyextra = 1;
    defaults.family = 'binomlogit';
    defaults.baseline = zeros(length(y),1);
    defaults.lambda0 = mean(sum(X.^2))*.01*ones(size(qfs,3),1); %finds a reasonable range for lambda0
    defaults.qf0 = 1e-3*speye(size(X,2));
    defaults.getH = 0;
    defaults.weights = ones(size(y));
    defaults.Display = 'iter';
    defaults.maxregressions = 5;
    
    opts = setdefaults(opts,defaults);
    
    subopts = opts;
    subopts.getH = true;
    
    subopts = rmfield(subopts,{'minlambda','maxregressions','evidenceeps','lambdaeps','weps','maxlambda','burnin','maxiter','qf0','lambda0','maxdelta'});
    subopts.Display = 'off';
    
    iter = 1;
    lambda = opts.lambda0;
    
    evidence = zeros(opts.maxiter,1);
    lambdah = zeros(length(lambda),opts.maxiter);
    
    w0 = zeros(size(X,2),1);
    w = w0;
    
    %Main loop
    nregressions = 0;
    results0 = [];
    
    %Initialize with reasonable value
    
    while iter <= opts.maxiter
        
        %Fit assuming lambdas are known
        subopts.w0 = w;
        Q = opts.qf0 + tailmdx(qfs,lambda);
        results = glmfitqp(y,X,Q,subopts);
        
        if iter == 1
            if strcmp(opts.family,'normid')
                %Estimate variance
                subopts.familyextra = std(y-X*results.w);
            end
            results = glmfitqp(y,X,Q,subopts);
        end
        
        %Compute evidence
        %Equation 4.137, Bishop
        evidence(iter) = results.loglikelihood + results.logpenalty - .5*mylogdet(results.H) + .5*mylogdet(Q);
        
        if strcmp(opts.family,'normid')
            evidence(iter) = evidence(iter) + length(y)/2*log(1/subopts.familyextra^2);
        end
        
        %Optimize evidence wrt lambda
        H0 = results.H - Q;
        
        Es = zeros(length(lambda),1);
        for ii = 1:length(lambda)
            Es(ii) = results.w'*squeeze(qfs(:,:,ii))*results.w;
        end
        
        active = lambda ~= opts.maxlambda;
        if all(~active)
            warning('blogreghp:divergentlambda','All lambdas diverged to Inf');
            break;
        end
        
        lambdah(:,iter) = lambda;
        
        if iter == 1 && size(qfs,3) > 1
            %Scale all items at once on first try
            if ~strcmp(opts.family,'normid')
                %Use trust-region Newton to get next guess
                [~,g,H] = nevidence(0,H0,opts.qf0,tailmdx(qfs,lambda),Es);
                lambda = lambda*exp(bound(0,-H\g,opts.maxdelta));
            else
                logprecision = log(1/subopts.familyextra^2);
                E0 = -2*results.loglikelihood*subopts.familyextra^2;
                [~,g,H] = nevidencegauss([logprecision,0],H0*subopts.familyextra^2,opts.qf0,tailmdx(qfs,lambda),Es(active),E0,size(X,1));

                %Check results
                %g2 = gradest(@(x) nevidencegauss( x,H0*subopts.familyextra^2,opts.qf0,qfs(:,:,active),Es(active),E0,size(X,1)),[logprecision;log(lambda(active))]);
                %H2 = hessian(@(x) nevidencegauss( x,H0*subopts.familyextra^2,opts.qf0,qfs(:,:,active),Es(active),E0,size(X,1)),[logprecision;log(lambda(active))]);

                loglambda = bound([logprecision;0],getDir(H,g),opts.maxdelta);
                lambda = lambda*exp(loglambda(2));
                logprecision = loglambda(1);
                
                subopts.familyextra = exp(-1/2*logprecision);
            end
        else
            if ~strcmp(opts.family,'normid')
                %Use trust-region Newton to get next guess
                [~,g,H] = nevidence(log(lambda(active)),H0,opts.qf0,qfs(:,:,active),Es(active));
                loglambda = bound(log(lambda(active)),getDir(H,g),opts.maxdelta);
            else

                logprecision = log(1/subopts.familyextra^2);
                E0 = -2*results.loglikelihood*subopts.familyextra^2;
                [~,g,H] = nevidencegauss([logprecision;log(lambda(active))],H0*subopts.familyextra^2,opts.qf0,qfs(:,:,active),Es(active),E0,size(X,1));

                %Check results
                %g2 = gradest(@(x) nevidencegauss( x,H0*subopts.familyextra^2,opts.qf0,qfs(:,:,active),Es(active),E0,size(X,1)),[logprecision;log(lambda(active))]);
                %H2 = hessian(@(x) nevidencegauss( x,H0*subopts.familyextra^2,opts.qf0,qfs(:,:,active),Es(active),E0,size(X,1)),[logprecision;log(lambda(active))]);

                loglambda = bound([logprecision;log(lambda(active))],getDir(H,g),opts.maxdelta);
                logprecision = loglambda(1);
                subopts.familyextra = exp(-1/2*logprecision);
                loglambda = loglambda(2:end);
            end
            lambda(active) = exp(loglambda);
        end

        if iter > opts.burnin
            lambda = max(min(lambda,opts.maxlambda),opts.minlambda);
        end
        
        %Verify everything has converged
        w = results.w;
        if iter > opts.burnin && std(w - w0) < opts.weps*std(w) && ...
                       all(abs(lambda - lambdah(:,iter-1)) < opts.lambdaeps*lambda) && ...
                       abs(evidence(iter) - evidence(iter - 1)) < opts.evidenceeps;
            fprintf('\n\n Converged in %d iterations\n',iter);
            break;
        end
        
        if iter > opts.burnin && evidence(iter) < max(evidence(1:iter - 1))
            %Use heuristics to figure out some new lambda
            if iter == opts.burnin + 1  %Early failure
                if norm(lambdah(:,1)) > norm(lambdah(:,2))
                    lambda = 10*lambdah(:,1);
                else
                    lambda = .1*lambdah(:,1);
                end
                iter = iter - 1;
                continue;
            elseif nregressions < opts.maxregressions
                fprintf(' * Regression at ');
                [~,maxiter] = max(evidence(1:iter-1));
                lambda = exp(-.5*(log(lambdah(:,iter)) - log(lambdah(:,maxiter))) + log(lambdah(:,maxiter)));
                nregressions = nregressions + 1;
                w = w0;
                results = results0;
            else
                a = sprintf('%5.2e ',lambdah(:,iter));
                fprintf(' * Regression at iteration %03d, evidence %9.2e, lambdas = [%s]; giving up\n',iter,evidence(iter),a);
                w = w0;
                results = results0;
                break;
            end 
        end
            
        a = sprintf('%5.2e ',lambdah(:,iter));
        fprintf('Iteration %03d, evidence %9.2e, lambdas = [%s]\n',iter,evidence(iter),a);
        iter = iter+1;
        
        w0 = w;
        results0 = results;
    end
    
    iter = min(iter,opts.maxiter);
    [~,maxiter] = max(evidence(1:iter));
    
    r.w = w;
    r.lambda = lambdah(:,maxiter);
    
    
    
    r.loglikelihood = results.loglikelihood;
    r.logpenalty = results.logpenalty;
    
    r.logevidence = evidence(maxiter);
    r.logevidences = evidence(1:iter);
    r.H = results.H;
    
    if strcmp(opts.family,'normid')
        r.sigmanoise = subopts.familyextra;
    end
end

function d = getDir(H,g)
    [R,posDef] = chol(H);

    % If the Cholesky factorization was successful, then the Hessian is
    % positive definite, solve the system
    if posDef == 0
        d = -R\(R'\g);
    else
        H = diag(diag(H));
        % otherwise, adjust the Hessian to be positive definite based on the
        % minimum eigenvalue, and solve with QR
        % (expensive, we don't want to do this very much)
        %H = H + eye(length(g)) * max(0,1e-12 - min(real(eig(H))));
        d = -H\g;
    end
end

function x = bound(x0,delta,maxdelta)
    x = x0 + delta/norm(delta)*min(norm(delta),maxdelta);
end

function d = mylogdet(H)
    L = cholcov(H);
    d = 2*sum(log(diag(L)));
end

function [f,gl,H] = nevidencegauss(loglambdas,H0,Q0,qfs,Es,E0,Nobs)
    logprecision = loglambdas(1);
    loglambda = loglambdas(2:end);
    precision = exp(logprecision);
    lambda = exp(loglambda);
    At = tailmdx(qfs,lambda)+Q0;
    H = precision*H0 + At;
    
    if nargout > 1
        %Do through eigenvalue decomposition
        [VH,DH] = eig(H);
        [VA,DA] = eig(At);

        logdetA = sum(log(diag(DA)));
        logdetH = sum(log(diag(DH)));
    else
        %Do through Choleski decomp
        logdetA = mylogdet(At);
        logdetH = mylogdet(H);
    end

    f = -(1/2*logdetA - 1/2*logdetH + 1/2*Nobs*logprecision - 1/2*lambda'*Es-1/2*precision*E0);
    
    if nargout > 1
        Ainvp = (VA*bsxfun(@times,VA',1./diag(DA)))';
        Hinvp = (VH*bsxfun(@times,VH',1./diag(DH)))';

        gl = zeros(length(lambda)+1,1);
        gl(1) = .5*precision*(sum(sum(Hinvp.*H0))+E0) - 1/2*Nobs;
        
        H  = zeros(length(lambda));
        
        H(1,1) = (gl(1) + 1/2*Nobs + .5*precision^2*sum(sum(-(Hinvp*H0*Hinvp).*H0)))*.5;
        
        for ii = 1:size(qfs,3)
            Ap = squeeze(qfs(:,:,ii));
            gl(ii+1) = -.5*lambda(ii)*(sum(sum((Ainvp-Hinvp).*Ap))-Es(ii));
            
            %Compute del g del lambda
            H(ii+1,ii+1) = (gl(ii+1) - 1/2*lambda(ii)^2*sum(sum( ( -Ainvp*Ap*Ainvp + Hinvp*Ap*Hinvp ).*Ap)))*.5;
            
            H(1,ii+1)    = .5*precision*lambda(ii)*sum(sum( -(Hinvp*Ap*Hinvp).*H0));
            
            %for jj = ii+1:size(qfs,3)
            %    Ap2 = squeeze(qfs(:,:,jj));
            %    H(ii+1,jj+1) = -1/2*lambda(ii)*lambda(jj)*sum(sum( (-Ainvp*Ap2*Ainvp + Hinvp*Ap2*Hinvp).*Ap) );
            %end
        end
        
        H = H + H';
    end
end

%Computes the negative log-evidence and its derivatives
function [f,gl,H] = nevidence(loglambda,H0,Q0,qfs,Es)
    lambda = exp(loglambda);
    At = tailmdx(qfs,lambda)+Q0;
    H = H0 + At;
    
    if nargout > 1
        %Do through eigenvalue decomposition
        [VH,DH] = eig(H);
        [VA,DA] = eig(At);

        logdetA = sum(log(diag(DA)));
        logdetH = sum(log(diag(DH)));
    else
        %Do through Choleski decomp
        logdetA = mylogdet(At);
        logdetH = mylogdet(H);
    end

    f = -(1/2*logdetA - 1/2*logdetH - 1/2*lambda'*Es);
    
    if nargout > 1
        Ainvp = (VA*bsxfun(@times,VA',1./diag(DA)))';
        Hinvp = (VH*bsxfun(@times,VH',1./diag(DH)))';

        gl = zeros(length(lambda),1);
        H  = zeros(length(lambda));
        for ii = 1:size(qfs,3)
            Ap = squeeze(qfs(:,:,ii));
            gl(ii) = -.5*lambda(ii)*(sum(sum((Ainvp-Hinvp).*Ap))-Es(ii));
            
            %Compute del g del lambda
            H(ii,ii) = (gl(ii) - 1/2*lambda(ii)^2*sum(sum( ( -Ainvp*Ap*Ainvp + Hinvp*Ap*Hinvp ).*Ap)))*.5;
            %for jj = ii+1:size(qfs,3)
            %    Ap2 = squeeze(qfs(:,:,jj));
            %    H(ii,jj) = -1/2*lambda(ii)*lambda(jj)*sum(sum( (-Ainvp*Ap2*Ainvp + Hinvp*Ap2*Hinvp).*Ap) );
            %end
        end
        H = H + H';
    end
end

function [M] = tailmdx(M,v)
    if numel(v) == 1
        M = M*v;
        return;
    end
    sizes = size(M);
    M = reshape(M,prod(sizes(1:end-1)),sizes(end));
    %if nnz(M)/numel(M) < .1
    %    M = reshape(full(sparse(M)*v(:)),sizes(1:end-1));
    %else
        M = reshape(M*v(:),sizes(1:end-1));
    %end
end