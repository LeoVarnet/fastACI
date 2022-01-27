function qf = qfsmooth1D(numx)
% function qf = qfsmooth1D(numx)
%
% qf = qfsmooth1D(numx)
% Create a quadratic form for smoothness regularization based on
% second-order derivative operator, for a one-dimensional signal

D = zeros(numx+1,numx);
for ii = 1:numx-1
    D(ii,ii:ii+1) = [-1 1];
end

qf = D'*D;