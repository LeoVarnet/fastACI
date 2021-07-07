%% Closing message

clc
if i_current==cfg_game.N
%     fprintf('\n   *** EXPERIMENT COMPLETE! *** \n\n');
%     fprintf('   Thank you for participating! \n');
    fprintf('\n   *** FIN DE L''EXPERIENCE ! *** \n\n');
    fprintf('   Merci de votre participation ! \n');
else
%     fprintf('\n   *** BREAK *** \n\n');
%     fprintf(['   You have already categorized ' num2str(i_current-1) ' of ' num2str(cfg_game.N) ' stimuli.\n']);
    fprintf('\n   *** PAUSE *** \n\n');
    fprintf(['   Vous avez ecoute ' num2str(i_current-1) ' sons sur un total de ' num2str(cfg_game.N) '.\n']);
end
%fprintf(['   Saving game to ' savename '.\n']);
fprintf(['   Sauvegarde dans le fichier ' savename '.\n']);