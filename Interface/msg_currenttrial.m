if isfield(cfg_game,'expvar_description')
    expvar_description = [', ' cfg_game.expvar_description];
else
    expvar_description = '';
end

switch cfg_game.Language
    case 'EN'
        
    case 'FR'
        
end
switch cfg_game.Language
    case 'EN'
        fprintf('\n\t*** MAIN EXPERIMENT ***\n\n');
        fprintf('\tPlaying stimulus # %.0f of %.0f (Next session stop in %.0f trials)\n',i_current,cfg_game.N,N_for_next_stop);
        fprintf('\tDependent variable: expvar = %.2f%s \n',expvar,expvar_description);
        fprintf('\n');
    case 'FR'
        fprintf('\n\t*** EXP\311RIENCE PRINCIPALE ***\n\n');
        fprintf('\t\311coute num\351ro %.0f sur %.0f (prochaine pause dans %.0f \351coutes)\n',i_current,cfg_game.N,N_for_next_stop);
        fprintf('\tVolume relatif de la voix : %.2f dB\n',expvar);
        
        if i_current>100
            bias_r1 = 100*sum(data_passation.n_responses(end-100+1:end)==1)/100;
            fprintf('\tPourcentage de r\351ponses "%s" : %.1f %%\n',cfg_game.response_names{1},bias_r1);
            if cfg_game.is_simulation == 1
                bias_r2 = 100*sum(data_passation.n_responses(end-100+1:end)==2)/100;
                % Extra info for the simulations:
                fprintf('\tPourcentage de r\351ponses "%s" : %.1f %%\n',cfg_game.response_names{2},bias_r2);
                fprintf('\tPourcentage de bonnes r\351ponses: %.1f %%\n',100*sum(data_passation.is_correct(end-100+1:end)==1)/100);
            else
                if bias_r1>60
                    fprintf('\tPourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{1});
                elseif bias_r1<40
                    fprintf('\tPourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{2});
                else
                    fprintf('\tPourcentage de r\351ponses "%s" normal \n',cfg_game.response_names{1});
                end
            end
        end
        fprintf('\n');
end