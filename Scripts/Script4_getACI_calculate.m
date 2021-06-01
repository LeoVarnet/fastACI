function [ACI,results,cfg_ACI] = Script4_getACI_calculate(cfg_ACI, y, y_correct, X, U)
% function [ACI,results,cfg_ACI] = Script4_getACI_calculate(cfg_ACI)
%
% 1. Description (FR):
%       Se placer dans le dossier contenant le dossier contenant les donnees 
%       du participant a analyser
%
%    For shortening the calculation process in this script, the following 
%       simplifications can be done (please run with NO simplification when
%       preparing data for publication):
%         - IdxTrialsLoad set to a lower range (e.g. IdxTrialsLoad = 1:1000)
%         - cfg_ACI.lambda0   = 85; % Valeur initiale de lambda
%
% To set:
%       DimCI: ('tf' or 'lyon') - Choice of auditory model: Lyon or something else (see Varnet2015: 'Cochleograms')
%       opts_ACI.glmfct (default: glm)
%
% 1. Reading data from a MAT file and check compatibility (il_convert_tb_ACI_data)
% 2. Reading/setting options for calculation
% 3. Reading the sound waveforms, getting T-F representations (data_load)
% 4. Preprocessing before the ACI assessment
%
% Scripts from where 'Script4_Calcul_ACI_debug' is called from:
%   1. g20210301_recreating_varnet2013.m
%   2. g20210413_SAO5000_varnet2013.m
%
% Changing the parameter names:
% New name      Old name        Changed on:
%
% Old name: Script4_Calcul_ACI_calculate.m (changed on 21/05/2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmfct = cfg_ACI.glmfct;

%% 5. Calculation of the ACI

N_f = length(cfg_ACI.f);
N_t = length(cfg_ACI.t);

fprintf('Starting the ACI assessment\n');% ,Analysis_condition);
tic

do_permutation = cfg_ACI.flags.do_permutation;
if do_permutation
    cfg_perm = cfg_ACI.cfg_perm;
    N_perm = cfg_perm.N_perm;
end

switch glmfct
    case 'glmfitqp'
        [ACI, results, cfg_ACI] = CI_glmqpoptim_fct(cfg_ACI, y, y_correct, X, U); 

        if do_permutation
            cfg_perm = Merge_structs(cfg_perm,cfg_ACI);
            cfg_perm.maxiter = 1;
            cfg_perm.lambda0 = results.finallambda * sqrt(cfg_ACI.stepsize); % crossValidate.m: stepsize will be 'compensated'

            for i = 1:N_perm
                fprintf('\t Assessing permuted ACI: %.0f of %.0f\n',i,N_perm);

                U_perm_here = squeeze(cfg_perm.U_perm(:,:,i));
                [ACI_perm_here, results_perm, cfg_perm] = CI_glmqpoptim_fct(cfg_perm, cfg_perm.y_perm(:,i), [], X, U_perm_here); 
                ACI_perm(:,:,i) = reshape(ACI_perm_here, [N_f,N_t]);
            end

            ACI_perm_CI_low  = prctile(ACI_perm,5,3);
            ACI_perm_CI_high = prctile(ACI_perm,95,3);

            results.idxs_perm = cfg_perm.idxs_perm;
            results.ACI_perm  = ACI_perm;
            results.ACI_perm_CI_low   = ACI_perm_CI_low;
            results.ACI_perm_CI_high  = ACI_perm_CI_high;
        end

        % The following warning appears when using the Optimization
        %   toolbox of MATLAB in version R2020a and it seems to
        %   have originated as of version R2018a (according to
        %   some forums on the internet).
        [msg,warnID] = lastwarn;
        msg_ref='The quasi-newton algorithm does not use analytic Hessian. Hessian flag in options will be ignored (supplied Hessian will not be used).';
        warnID_ref = 'optim:fminunc:HessIgnored';
        if strcmp(msg,msg_ref) || strcmp(warnID,warnID_ref)
            warning('Turn on manually the use of the old optimisation toolbox of MATLAB (go to glmfitqp.m->function il_irls)')
        end

        % results.evaluation.opts

    case 'classic_revcorr'
        [ACI, pval_ACI] = corr(y,X);
        ACI = reshape(ACI, [N_f,N_t]);
        pval_ACI = reshape(pval_ACI, [N_f,N_t]);
        results.ACI = ACI;
        results.pval_ACI = pval_ACI;

        if do_permutation
            for i = 1:N_perm
                fprintf('\t Assessing permuted ACI: %.0f of %.0f\n',i,N_perm);

                U_perm_here = squeeze(cfg_perm.U_perm(:,:,i));
                ACI_perm_here = corr(cfg_perm.y_perm(:,i), X); 
                ACI_perm(:,:,i) = reshape(ACI_perm_here, [N_f,N_t]);
            end

            ACI_perm_CI_low  = prctile(ACI_perm,5,3);
            ACI_perm_CI_high = prctile(ACI_perm,95,3);

            results.idxs_perm = cfg_perm.idxs_perm;
            results.ACI_perm  = ACI_perm;
            results.ACI_perm_CI_low   = ACI_perm_CI_low;
            results.ACI_perm_CI_high  = ACI_perm_CI_high;
        end

    case {'lassoglm','lasso'}
        N_folds = cfg_ACI.N_folds; % k_folds validation
        % cfg_ACI = Ensure_field(cfg_ACI,'lambda0',[]);
        lambda0 = cfg_ACI.lambda0;

        tic
        switch glmfct
            case 'lassoglm'
                [B,FitInfo] = lassoglm(X,y,'binomial','CV',N_folds,'lambda',lambda0,'link','logit');
                results.cvgofs  = FitInfo.Deviance; % this is already the mean value
                idxlambda = FitInfo.IndexMinDeviance;
                
            case 'lasso'
                [B,FitInfo] = lasso(X,y,'CV',N_folds);
                idxlambda = FitInfo.IndexMinMSE;
                
        end
        toc

        results.lambdas = FitInfo.Lambda; 
        % results.cvgofs  = FitInfo.Deviance; % this is already the mean value
        results.FitInfo = FitInfo;
        results.B       = B;
        % results.w       = B(:,FitInfo.IndexMinDeviance);
        % results.finalfit.w = results.w;
                
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
        
        disp('')
        % %%
        % figure; 
        % for i_level = Nlevelmin:Nlevel
        %     subplot(2,4,i_level)
        %     h=pcolor(t_X, f_X, ReWeightPyramid{i_level}(:,:,1)); set(h,'EdgeColor','none'); xlabel('time (s)'); ylabel('freq (Hz)');
        % end
        % 
        % %%
        % figure;
        % %subplot(2,3,1)
        % plot(FitInfo.Lambda, FitInfo.MSE); hold on
        % plot(FitInfo.Lambda(FitInfo.IndexMinMSE), FitInfo.MSE(FitInfo.IndexMinMSE),'ro');
        % title('MSE')
        %     %%
        %     figure
        % for i_level = Nlevelmin:Nlevel
        %     subplot(2,3,i_level)
        %     plot(FitInfo.Lambda, reshape(WeightPyramid{i_level},[],size(B,2)))
        %     title(['Weights, level ' num2str(i_level)])
        % end
         
        %%% Plot fit corresponding to the best lambda
        sumReWeight = zeros(size(ReWeightPyramid{Nlevel}));
        for i_level = Nlevelmin:Nlevel
            sumReWeight = sumReWeight + ReWeightPyramid{i_level};
        end
        
        if length(cfg_ACI.t_limits_idx) < size(sumReWeight,3) % Dim 3 is time  
            sumReWeight = sumReWeight(:,:,cfg_ACI.t_limits_idx);    
            cfg_ACI.t_X = cfg_ACI.t_X(cfg_ACI.t_limits_idx);
        end
        
        ACI = squeeze( sumReWeight(idxlambda,:,:) );
        % ACI = ACI(:,cfg_ACI.t_limits_idx);
        
        results.ACI = ACI;
        results.ACIs = sumReWeight;
        
        bPlot = 1; warning('This is temporal here...')
        if bPlot
            figure; 
            % h=pcolor(t_X, f_X, -sumReWeight(:,:,idxlambda)); set(h, 'EdgeColor', 'none'); xlabel('time (s)'); ylabel('freq (Hz)'); colorbar; title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); colorbar; caxis([-1 1]*max(abs(caxis)));
            h=pcolor(cfg_ACI.t_X, cfg_ACI.f_X, ACI); 
            set(h, 'EdgeColor', 'none'); 
            xlabel('time (s)'); 
            ylabel('freq (Hz)'); colorbar; 
            title('betaSmooth: ACI obtained with a lasso regression on smooth basis'); 
            colorbar; % caxis([-1 1]*max(abs(caxis)));
            % %set(gca, 'YScale', 'log');
            
            Nt = length(cfg_ACI.t_X);
            Nf = length(cfg_ACI.f_X);
            
            t_spec = cfg_ACI.t_X;
            f_spec = cfg_ACI.f_X;
            figure;
            if N_iterations >= 50
                betalasso50 = reshape(sumReWeight(50,:,:), [Nf,Nt]);
                
                subplot(5,1,1)
                hp = pcolor(t_spec, f_spec, betalasso50); 
                xlabel('time (s)'); 
                ylabel('freq (Hz)');colorbar
                set(hp, 'EdgeColor', 'none'); 
                title(['ACI obtained with a lasso regression, \lambda = ' num2str(FitInfo.Lambda(50))])
            end
            if N_iterations >= 60
                betalasso60 = reshape(sumReWeight(60,:,:), [Nf,Nt]);
                
                subplot(5,1,2)
                hp = pcolor(t_spec, f_spec, betalasso60); xlabel('time (s)'); ylabel('freq (Hz)');colorbar
                set(hp, 'EdgeColor', 'none'); 
                title(['ACI obtained with a lasso regression, \lambda = ' num2str(FitInfo.Lambda(60))])
            end
            if N_iterations >= 70
                betalasso70 = reshape(sumReWeight(70,:,:), [Nf,Nt]);
                
                subplot(5,1,3)
                hp = pcolor(t_spec, f_spec, betalasso70); xlabel('time (s)'); ylabel('freq (Hz)');colorbar
                set(hp, 'EdgeColor', 'none'); 
                title(['ACI obtained with a lasso regression, \lambda = ' num2str(FitInfo.Lambda(70))])
            end
            if N_iterations >= 80
                betalasso80 = reshape(sumReWeight(80,:,:), [Nf,Nt]);
                
                subplot(5,1,4)
                hp = pcolor(t_spec, f_spec, betalasso80); xlabel('time (s)'); ylabel('freq (Hz)');colorbar
                set(hp, 'EdgeColor', 'none'); 
                title(['ACI obtained with a lasso regression, \lambda = ' num2str(FitInfo.Lambda(80))])
            end
            if N_iterations >= 90
                betalasso90 = reshape(sumReWeight(90,:), [Nf,Nt]);
                
                subplot(5,1,5)
                hp = pcolor(t_spec, f_spec, betalasso90); xlabel('time (s)'); ylabel('freq (Hz)');colorbar
                set(hp, 'EdgeColor', 'none'); 
                title(['ACI obtained with a lasso regression, \lambda = ' num2str(FitInfo.Lambda(90))])
            end
            % % [~,bestlambda] = min(FitInfo.Deviance);
            % betalassobest = reshape(B(:,idxlambda), [Nf,Nt]);
            % 
            % figure;
            % plot(FitInfo.Lambda,FitInfo.Deviance,'-');hold on;plot(FitInfo.Lambda([50,60,70,80,90]),FitInfo.Deviance([50,60,70,80,90]),'ro');plot(FitInfo.Lambda(bestlambda),FitInfo.Deviance(bestlambda),'b*')
            % 
            % figure;
            % pcolor(t_spec, f_spec, betalassobest); xlabel('time (s)'); ylabel('freq (Hz)');colorbar
            % title(['ACI obtained with a lasso regression, best \lambda = ' num2str(FitInfo.Lambda(bestlambda))])
        end

    otherwise
        error('Function %s not debugged yet...',glmfct);
        % [ACI, results, cfg_ACI] = feval(cfg_ACI.glmfct, cfg_ACI, y, y_correct, X, U );
end
fprintf('Completed condition \n');
toc
