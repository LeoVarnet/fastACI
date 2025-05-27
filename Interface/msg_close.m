%% Closing message

clc
if i_current==cfg_game.N_trials
    switch cfg_game.Language
        case 'EN'
            fprintf('\n   *** EXPERIMENT COMPLETE! *** \n\n');
            fprintf('   Thank you for participating! \n');
            
        case 'FR'
            fprintf('\n   *** FIN DE L''EXP\311RIENCE ! *** \n\n');
            fprintf('   Merci pour votre participation ! \n');
    end
else
    switch cfg_game.Language
        case 'EN'
            fprintf('\n   *** BREAK *** \n\n');
            fprintf('   You have already categorised %.0f of %.0f stimuli.\n',i_current-1,cfg_game.N_trials);
        case 'FR'
            fprintf('\n   *** PAUSE *** \n\n');
            fprintf('   Vous avez \351cout\351 %.0f sons sur un total de %.0f.\n',i_current-1,cfg_game.N_trials);
    end
end
switch cfg_game.Language
    case 'EN'
        fprintf('   Saving game to %s.\n',savename);
    case 'FR'
        fprintf('   Sauvegarde dans le fichier %s.\n',savename);
end