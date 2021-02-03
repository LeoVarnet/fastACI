%% Closing message

clc
if i==cfg_game.N
    fprintf('\n   *** EXPERIMENT COMPLETE! *** \n\n');
    fprintf('   Thank you for participating! \n');
else
    fprintf('\n   *** BREAK *** \n\n');
    fprintf(['   You have already categorized ' num2str(i-1) ' of ' num2str(cfg_game.N) ' stimuli.\n']);
end
fprintf(['   Saving game to ' savename '.\n']);