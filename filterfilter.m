function [Y] = filterfilter(B,A,X)
%FILTERFILTER my own simple version of filtfilt with a 'symmetric' option
%on correcting for edge effects (see imfilter)
%
% Leo Varnet - 2019

nsamples = floor(length(X)/2);
X = X(:);
Xsym = [flipud(X(1:nsamples)); X; flipud(X(nsamples:end))];
Ysym = flipud(filter(B,A,flipud(filter(B,A,Xsym))));
Y = Ysym(nsamples+1:nsamples+length(X));

end

