function [ll,g,Hd,r,lleach,y_est] = evalGlmLikelihood(y,X,w,b,family,familyextra,weights)
% function [ll,g,Hd,r,lleach,y_est] = evalGlmLikelihood(y,X,w,b,family,familyextra,weights)
%
% [ll,g] = evalGlmLikelihood(y,X,w,b,family,familyextra,weights)
% Evaluates the negative log-likelihood of the data assuming a GLM with the 
% specified properties
%
% Comment by Alejandro on 14/01/2022:
% In this script, ll is not referenced to a saturated model (see https://en.wikipedia.org/wiki/Deviance_(statistics)). 
%     Furthermore, this estimate of deviance is not normalised by the number
%     of samples, meaning that the estimate will be a higher number if the
%     number of samples within a fold is higher, even if the same dataset
%     is fitted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if nargin < 7
    weights = 1;
end

r = X*w + b;
switch family
    case 'normid'
        y_est = r;
        gain = 1/familyextra^2;
        lleach = .5/familyextra^2*(y-y_est).^2;
    case 'binomlogit'
        gain = familyextra; % (gain = 1 if familyextra == 1)
        y_est = 1./(1+exp(-r));
        if familyextra == 1
            % 'Can skip a few computations in this case'. In fact, it should
            %    be the same calculation as for familyextra ~= 1, but here 
            %    it is visually clearer how the assessment is done.
            lleach = zeros(size(y));
            lleach(y==1) = -log(  y_est(y==1)+eps);
            lleach(y==0) = -log(1-y_est(y==0)+eps);
        else
            y = y/familyextra;
            lleach = -familyextra*(y.*log(y_est+eps) + (1-y).*log(1-y_est+eps));
        end
    case 'poissexp'
        y_est = exp(r);
        gain = 1;
        lleach = -y.*r + y_est;
    otherwise 
        error('Unsupported family');
end

lleach = lleach.*weights;
ll = sum(lleach);

if nargout > 1
    %These are all canonical links
    r = weights.*(y_est-y);
    g = gain*(X'*r);
end
if nargout > 2
    %Compute the diagonal term that occurs in the Hessian (canonical
    %links only)
    switch family
        case 'normid'
            Hd = 1/familyextra^2*ones(size(X,1),1).*weights;
        case 'binomlogit'
            Hd = familyextra*y_est.*(1-y_est).*weights;
        case 'poissexp'
            Hd = y_est.*weights;
    end
end