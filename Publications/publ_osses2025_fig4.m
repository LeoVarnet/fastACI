function [] = publ_osses2025_fig4()
% function publ_osses2025_fig4()
%%%%% Figure 4 from fastACI paper, hyperparameter selection figure %%%%%
%
% % To display Fig. 4 of Osses et al use
%     publ_osses2025_fig4();
% You should plot figure 7 once first to compute the ACI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

path_to_ACI = [fastACI_paths('dir_data') 'speechACI_Logatome-abda-S43M' filesep 'S01' filesep 'Results/Results_ACI/'];

%path_to_ACI = '/home/leovarnet/Nextcloud/Data/fastACI_data/speechACI_Logatome-abda-S43M/S01/Results/Results_ACI/';
ACI_name = 'ACI-S01-speechACI_Logatome-sMPSv1p3-nob-gt-l1glm+pyrga-rev4.mat';

load([path_to_ACI ACI_name])

lambdas_to_plot = [2 8 results.idxlambda 19];
N_lambdas_to_plot = length(lambdas_to_plot);

%%

figure('Position',[100 100 800 500]);
loglog(results.lambdas,mean(results.FitInfo.Dev_test,2),'LineWidth',1.5); hold on
loglog(results.lambdas(lambdas_to_plot),mean(results.FitInfo.Dev_test(lambdas_to_plot,:),2),'*','MarkerSize',10); hold on

xlabel('\lambda')
ylabel('CV deviance')
ylim([500 800])
%loglog(results.lambdas,mean(results.FitInfo.Dev_train,2)); hold on

%%
figure('Position',[100 100 800 200]);
t = tiledlayout(1,N_lambdas_to_plot,'TileSpacing', 'compact');

for i_lambdas = 1:N_lambdas_to_plot

nexttile(i_lambdas); %(i_subject,i_masker)
affichage_tf(squeeze(results.ACIs(lambdas_to_plot(i_lambdas),:,:)),'CI', 'cfg', cfg_ACI,'NfrequencyTicks',5); %hold on
    %title([Maskers_titles(i_masker)], 'interpreter','none');
    if i_lambdas > 1
        ylabel('');set(gca,'YTickLabels',{});
    end
    if i_lambdas == N_lambdas_to_plot
        c_lim = caxis;
        colorbar('Ticks',c_lim,'TickLabels',{'ada','aba'});% change colorbar
    end
    if i_lambdas < N_lambdas_to_plot
        colorbar off
    end
    set(gca, 'FontSize', 12)
end

end