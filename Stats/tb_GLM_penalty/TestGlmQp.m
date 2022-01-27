%Smoothness prior
qf = blkdiag(qfsmooth1D(16),.01);
rg = (-7.5:7.5)';
w = exp(-rg.^2/3^2).*sin(rg*2*pi/6);

%Simulate some data w/ gaussian noise 
nobs = 30;
X = [randn(nobs,length(w)),ones(nobs,1)];
y0 = X*[w;.01] + .1*randn(nobs,1);

clear opts
opts.family = 'normid';
results = evidenceglmfitqp(y0,X,.5*qf,opts);

plot(results.w(1:end-1))


%%
%Same, but with logistic regression
qf = blkdiag(qfsmooth1D(16),.01);
rg = (-7.5:7.5)';
w = exp(-rg.^2/3^2).*sin(rg*2*pi/6);

%Simulate some data w/ gaussian noise 
nobs = 150;
X = [randn(nobs,length(w)),ones(nobs,1)];
r = 3*X*[w;.01];
r = binornd(1,1./(1+exp(-r)));

clear opts
opts.family = 'poissexp';
results = glmfitqp(r,X,2*qf,opts);

plot(results.w(1:end-1))

%%
%Figure out optimal strength of prior through cross validation
%or through evidence optimization
%Assume smoothness of the model parameters
rg = (-15.5:15.5)';
qf = blkdiag(qfsmooth1D(length(rg)),.01);


%Simulate a model with w = Gabor function
w = exp(-rg.^2/4^2).*sin(rg*2*pi/15);
nobs = 150;
X = [randn(nobs,length(w)),ones(nobs,1)];
r = 3*X*[w;.01];

%output is binary -> logistic regression
r = binornd(1,1./(1+exp(-r)));

%Set up 5-fold CV
folds = getcvfolds(length(r),5,1001);

%Fit the data
clear opts
opts.family = 'binomlogit';

results = cvglmfitqp(r,X,qf,folds,opts);
results0 = glmfitqp(r,X,qf*0,opts);
results2 = evidenceglmfitqp(r,X,qf,opts);

clf;
subplot(1,4,1);plot(w);title('Actual filter');axis([.5,length(w)+.5,-.7,.7]);
subplot(1,4,2);plot(results0.w(1:end-1));title('Maximum likelihood');axis([.5,length(w)+.5,[-1.1,1.1]*max(abs(results0.w(1:end-1)))]);
subplot(1,4,3);plot(results.w(1:end-1));title('MAP w/ cross-val');axis([.5,length(w)+.5,[-1.1,1.1]*max(abs(results2.w(1:end-1)))]);
subplot(1,4,4);plot(results2.w(1:end-1));title('MAP w/ evidence');axis([.5,length(w)+.5,[-1.1,1.1]*max(abs(results2.w(1:end-1)))]);

%%
%Figure out optimal strength of prior through evidence optimization
%Here the weight parameter is a 2d Gabor and the smoothness varies along
%the two dims
[xi,yi] = meshgrid(-7.5:7.5,-7.5:7.5);

%Filter is a 2D Gabor
w = sin(xi).*exp(-(xi.^2+yi.^2)/2/3^2);
w = w(:);

[~,~,qf1,qf2] = qfsmooth(16,16);
Q = zeros([size(qf1)+1,2]);
Q(1:end-1,1:end-1,1) = qf1;
Q(1:end-1,1:end-1,2) = qf2;

%Simulate a model with w = Gabor function
nobs = 800;
X = [randn(nobs,length(w)),ones(nobs,1)];
r = .5*X*[w;.01];

%output is binary -> logistic regression
r = binornd(1,1./(1+exp(-r)));

%Fit the data
clear opts
opts.family = 'binomlogit';


%Because there are two hyperparameters to optimize, use evidence framework
results = evidenceglmfitqp(r,X,Q,opts);
results0 = glmfitqp(r,X,speye(257)*.1,opts);

clf;
subplot(1,3,1);imagesc(reshape(w,16,16));title('Actual filter');
subplot(1,3,2);imagesc(reshape(results0.w(1:end-1),16,16));title('Maximum likelihood');
subplot(1,3,3);imagesc(reshape(results.w(1:end-1),16,16));title('Maximum A Posteriori w/ smoothness prior');

%%

%Parallel processing
%Assuming you have the distributed computing toolbox, you can parallelize
%cross-validation. This works best when k is equal to a multiple of the
%number of cores in your CPU, so for 4-core computers, 4 or 8-fold
%validation is optimal. 
%
%Note that this may require a lot of RAM

%Call this once:
matlabpool open

%%
%Figure out optimal strength of prior through cross validation
%Assume smoothness of the model parameters
rg = (-15.5:15.5)';
qf = blkdiag(qfsmooth1D(length(rg)),.01);


%Simulate a model with w = Gabor function
w = exp(-rg.^2/4^2).*sin(rg*2*pi/15);
nobs = 1500;
X = [randn(nobs,length(w)),ones(nobs,1)];
r = 3*X*[w;.01];

%output is binary -> logistic regression
r = binornd(1,1./(1+exp(-r)));

%Set up 8-fold CV
folds = getcvfolds(length(r),8,1001);

%Fit the data
clear opts
opts.family = 'binomlogit';

%Regular
tic;results = cvglmfitqp(r,X,qf,folds,opts);toc;

%Parallel
opts.parallel = 1;
tic;results = cvglmfitqp(r,X,qf,folds,opts);toc;
%Same results, much faster (2.5x here)
