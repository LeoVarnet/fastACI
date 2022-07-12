function [ACI,cfg_ACI,B_times_w] = Convert_lasso_B2ACI(w, cfg_ACI, idxlambda,keyvals) 
% function [ACI,cfg_ACI,sumReWeight] = Convert_lasso_B2ACI(w, cfg_ACI, results,keyvals)
%
% See Script_LassoPyramid_21042021.m by Leo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    idxlambda = 1;
end

temp = w;
Nlevelmin = cfg_ACI.lasso_Nlevelmin;
Nlevel    = cfg_ACI.lasso_Nlevel;
Pyra_size = cfg_ACI.lasso_Pyra_size;

N_iterations = size(w,2);

for i_level = Nlevelmin:Nlevel
    %%% 1. Getting the 'reduced' ACIs per level
    idxi = 1;
    idxf = Pyra_size(i_level,1)*Pyra_size(i_level,2);
    Pyra_here = temp(idxi:idxf,:);
    Size_reshape = [Pyra_size(i_level,1) Pyra_size(i_level,2) N_iterations];
    Pyra_here = reshape(Pyra_here,Size_reshape);
    
    Pyra_here = permute(Pyra_here,[3 1 2]); % time and frequency dimensions in the proper location
    
    WeightPyramid{i_level} = Pyra_here;

    temp = temp(idxf+1:end,:); % The assigned bins are removed

    %%% 2. Getting the 'expanded' ACIs per level (Expand the weight matrix)
    % Interpolate so that dimension 1 has length Nt_X in each level of the
    % pyramid -- this is the same procedure as before with RePyramid 

    switch keyvals.pyramid_script
        case 'imgaussfilt'
            keyvals.i_level   = i_level;
            keyvals.Pyra_size = Pyra_size;
            Pyra_here = fastACI_impyramid(Pyra_here,'expand',keyvals);
            
        case 'imresize'
            % Pyra_here = reshape(WeightPyramid{i_level},[Pyra_size(i_level,1),Pyra_size(i_level,2), N_iterations]);
            for j_level = 1:i_level-1
                Pyra_here = fastACI_impyramid(Pyra_here, 'expand',keyvals);
            end
    end
    Blevel_times_w{i_level} = squeeze(Pyra_here(:,:,:));
end

%%% Plot fit corresponding to the best lambda
B_times_w = zeros(size(Blevel_times_w{Nlevel})); % here we will sum up each level of the pyramid
for i_level = Nlevelmin:Nlevel
    B_times_w = B_times_w + Blevel_times_w{i_level};
end

if size(B_times_w,3) ~= 1
    dim_time = 3; % Dim 3 is time when multiple lambdas are tested
else
    dim_time = 2;
end

if ~isfield(cfg_ACI,'t_limits_idx')
    cfg_ACI.t_limits_idx = 1:length(cfg_ACI.t);
end

if length(cfg_ACI.t_limits_idx) < size(B_times_w,dim_time) 
    % time was padded for the pyramid calculation, so we truncate it back
    switch dim_time
        case 3 % default
            B_times_w = B_times_w(:,:,cfg_ACI.t_limits_idx);
        case 2
            % when only one lambda is tested
            B_times_w = B_times_w(:,cfg_ACI.t_limits_idx);
    end
    if isfield(cfg_ACI,'t_X')
        if cfg_ACI.t_X(1) == cfg_ACI.t(1) && cfg_ACI.t_X(cfg_ACI.N_t) == cfg_ACI.t(cfg_ACI.N_t)
            % then t_X is equal (but longer than) t
            cfg_ACI = rmfield(cfg_ACI,'t_X');
        end
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

if ~isfield(cfg_ACI,'f_limits_idx')
    cfg_ACI.f_limits_idx = 1:length(cfg_ACI.f);
end
if length(cfg_ACI.f_limits_idx) > size(B_times_w,2) % Dim 2 is frequency
    % we need to truncate the frequency .f

    if isfield(cfg_ACI,'f_X')
        if cfg_ACI.f(1)==cfg_ACI.f_X(1) && cfg_ACI.f(length(cfg_ACI.f_X))==cfg_ACI.f_X(end)
            cfg_ACI.f = cfg_ACI.f_X;
            cfg_ACI.f_limits_idx = cfg_ACI.f_limits_idx(1:length(cfg_ACI.f_X));

            cfg_ACI = rmfield(cfg_ACI,'f_X');
        end
    end
end

if length(size(B_times_w)) == 3
    ACI = squeeze( B_times_w(idxlambda,:,:) ); % ACI equation
else
    ACI = squeeze( B_times_w(:,:) );
end

% ACI = ACI(:,cfg_ACI.t_limits_idx);