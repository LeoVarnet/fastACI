function [ll] = evalMaxGlmLikelihood(y,family,familyextra,weights)
    %Evaluate the negative of the maximum attainable likelihood for a given observation
    if nargin < 4
        weights = 1;
    end
    
    switch family
        case 'normid'
            ll = 0;
        case 'binomlogit'
            y = y/familyextra;
            ll = -familyextra*((weights.*y)'*log(y+eps) + (weights.*(1-y))'*log(1-y+eps));
        case 'poissexp'
            ll = -(weights.*y)'*log(y+eps) + sum(weights.*y);
    end
end