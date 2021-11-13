function [ACI,cfg_ACI,sumReWeight] = Convert_lasso_B2ACI(B, cfg_ACI, results)
% function [ACI,cfg_ACI,sumReWeight] = Convert_lasso_B2ACI(B, cfg_ACI, results)
%
% See Script_LassoPyramid_21042021.m by Leo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(B,2) ~= 1
    % Then multiple lambda values have been obtained
    idxlambda = results.idxlambda;
else
    idxlambda = 1;
end

temp = B;
Nlevelmin = cfg_ACI.lasso_Nlevelmin;
Nlevel    = cfg_ACI.lasso_Nlevel;
Pyra_size = cfg_ACI.lasso_Pyra_size;

for i_level = Nlevelmin:Nlevel
    %%% 1. Getting the 'reduced' ACIs per level
    idxi = 1;
    idxf = Pyra_size(i_level,1)*Pyra_size(i_level,2);
    N_iterations = size(B,2);
    Pyra_here = temp(idxi:idxf,:);
    Size_reshape = [Pyra_size(i_level,1) Pyra_size(i_level,2) N_iterations];
    Pyra_here = reshape(Pyra_here,Size_reshape);
    Pyra_here = permute(Pyra_here,[3 1 2]); % time and frequency dimensions in the proper location

    WeightPyramid{i_level} = Pyra_here;

    temp = temp(idxf+1:end,:); % The assigned bins are removed

    %%% 2. Getting the 'expanded' ACIs per level (Expand the weight matrix)
    % Interpolate so that dimension 1 has length Nt_X in each level of the
    % pyramid -- this is the same procedure as before with RePyramid 

    % Pyra_here = reshape(WeightPyramid{i_level},[Pyra_size(i_level,1),Pyra_size(i_level,2), N_iterations]);
    for j_level = 1:i_level-1
        Pyra_here = Script4_Calcul_ACI_modified_impyramid(Pyra_here, 'expand');
    end
    ReWeightPyramid{i_level} = squeeze(Pyra_here(:,:,:));
end

%%% Plot fit corresponding to the best lambda
sumReWeight = zeros(size(ReWeightPyramid{Nlevel}));
for i_level = Nlevelmin:Nlevel
    sumReWeight = sumReWeight + ReWeightPyramid{i_level};
end

if length(cfg_ACI.t_limits_idx) < size(sumReWeight,3) % Dim 3 is time 
    % time was padded for the pyramid calculation, so we truncate it back
    sumReWeight = sumReWeight(:,:,cfg_ACI.t_limits_idx);
    if cfg_ACI.t_X(1) == cfg_ACI.t(1) && cfg_ACI.t_X(cfg_ACI.N_t) == cfg_ACI.t(cfg_ACI.N_t)
        % then t_X is equal (but longer than) t
        cfg_ACI = rmfield(cfg_ACI,'t_X');
    end
else
    if isfield(cfg_ACI,'t_X')
        if length(cfg_ACI.t_X) ~= length(cfg_ACI.t)

            if length(cfg_ACI.t_X) < length(cfg_ACI.t)
                cfg_ACI.t = cfg_ACI.t(1:length(cfg_ACI.t_X));
            end

            cfg_ACI = rmfield(cfg_ACI,'t_X');
        end
    end
end

if length(cfg_ACI.f_limits_idx) > size(sumReWeight,2) % Dim 2 is frequency
    % we need to truncate the frequency .f

    if cfg_ACI.f(1)==cfg_ACI.f_X(1) && cfg_ACI.f(length(cfg_ACI.f_X))==cfg_ACI.f_X(end)
        cfg_ACI.f = cfg_ACI.f_X;
        cfg_ACI.f_limits_idx = cfg_ACI.f_limits_idx(1:length(cfg_ACI.f_X));

        cfg_ACI = rmfield(cfg_ACI,'f_X');
    end
end

if length(size(sumReWeight)) == 3
    ACI = squeeze( sumReWeight(idxlambda,:,:) );
else
    ACI = squeeze( sumReWeight(:,:) );
end

% ACI = ACI(:,cfg_ACI.t_limits_idx);