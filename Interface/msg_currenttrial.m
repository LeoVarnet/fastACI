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
        fprintf('\n    *** MAIN EXPERIMENT ***\n\n');
        fprintf('    Playing stimulus # %.0f of %.0f (Next session stop in %.0f trials)\n',i_current,cfg_game.N,N_for_next_stop);
        fprintf('    Dependent variable: expvar = %.2f%s \n',expvar,expvar_description);
        fprintf('\n');
    case 'FR'
        fprintf('\n    *** EXP\311RIENCE PRINCIPALE ***\n\n');
        fprintf('    \311coute num\351ro %.0f sur %.0f (prochaine pause dans %.0f \351coutes)\n',i_current,cfg_game.N,N_for_next_stop);
        fprintf('    Volume relatif de la voix : %.2f dB\n',expvar);
        if i_current>100
            bias_r1 = 100*sum(data_passation.n_responses(end-100+1:end)==1)/100;
            if bias_r1>60
                fprintf('    Pourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{1});
            elseif bias_r1<40
                fprintf('    Pourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{2});
            else
                fprintf('    Pourcentage de r\351ponses "%s" normal \n',cfg_game.response_names{1});
            end
        end
        fprintf('\n');
end