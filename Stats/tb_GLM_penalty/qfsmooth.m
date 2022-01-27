function [qf,D,qfx,qfy] = qfsmooth(numx,numy)
% function [qf,D,qfx,qfy] = qfsmooth(numx,numy)
% qf = qfsmooth(numx,numy)
%
% Creates a quadratic form for smoothness regularisation, with D a first-
%   order partial derivative operator and qf = D'*D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = zeros((numx-1)*numy + (numy-1)*numx,numx*numy);

for jj = 1:numy
    for ii = 1:numx-1
        [xi,yi] = meshgrid(1:numx,1:numy);
        dd = (xi == ii & yi == jj) - (xi == ii + 1 & yi == jj);
        D(ii + (jj-1)*(numx-1),:) = dd(:);
    end
end

for jj = 1:numy-1
    for ii = 1:numx
        [xi,yi] = meshgrid(1:numx,1:numy);
        dd = (xi == ii & yi == jj) - (xi == ii & yi == jj + 1);
        D((numx -1)*numy + ii + (jj-1)*numx,:) = dd(:);
    end
end

qf  = D'*D;
idxf = floor(size(D,1)/2);
D1  = D(1:idxf,:);
qfx = D1'*D1;
D1  = D(idxf+1:end,:);
qfy = D1'*D1;