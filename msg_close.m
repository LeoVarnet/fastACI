%% Closing message

clc
if i_current==cfg_game.N
    switch cfg_game.Language
        case 'EN'
            fprintf('\n   *** EXPERIMENT COMPLETE! *** \n\n');
            fprintf('   Thank you for participating! \n');
            
        case 'FR'
            fprintf('\n   *** FIN DE L''EXPERIENCE ! *** \n\n');
            fprintf('   Merci de votre participation ! \n');
    end
else
    switch cfg_game.Language
        case 'EN'
            fprintf('\n   *** BREAK *** \n\n');
            fprintf('   You have already categorised %.0f of %.0f stimuli.\n',i_current-1,cfg_game.N);
        case 'FR'
            fprintf('\n   *** PAUSE *** \n\n');
            fprintf('   Vous avez ecoute %.0f sons sur un total de %.0f.\n',i_current-1,cfg_game.N);
    end
end
switch cfg_game.Language
    case 'EN'
        fprintf('   Saving game to %s.\n',savename);
    case 'FR'
        fprintf('   Sauvegarde dans le fichier %s.\n',savename);
end