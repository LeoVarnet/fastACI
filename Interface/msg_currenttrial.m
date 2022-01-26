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
        
        N_here = min(i_current-1,100);
        if i_current>10
            if N_here < 100 && cfg_game.is_simulation
                fprintf('Averaging first %.0f trials\n',N_here);
            end
            bias_r1 = 100*sum(data_passation.n_responses(end-N_here+1:end)==1)/N_here;
            fprintf('\tPercentage of "%s" responses: %.1f %%\n',cfg_game.response_names{1},bias_r1);
            if cfg_game.is_simulation == 1
                bias_r2 = 100*sum(data_passation.n_responses(end-N_here+1:end)==2)/N_here;
                % Extra info for the simulations:
                fprintf('\tPercentage of "%s" responses: %.1f %%\n',cfg_game.response_names{2},bias_r2);
                fprintf('\tPercentage of correct responses: %.1f %%\n',100*sum(data_passation.is_correct(end-N_here+1:end)==1)/N_here);
            else
                if i_current > 100
                    if bias_r1>60
                        fprintf('\tPercentage of "%s" responses = %.0f %% (too much "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{1});
                    elseif bias_r1<40
                        fprintf('\tPercentage of "%s" responses = %.0f %% (too much "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{2});
                    else
                        fprintf('\tThe percentages of responses "%s" et "%s" are normally balanced \n',cfg_game.response_names{1},cfg_game.response_names{2});
                    end
                end
            end
        end
        fprintf('\n');
        
    case 'FR'
        fprintf('\n\t*** EXP\311RIENCE PRINCIPALE ***\n\n');
        fprintf('\t\311coute num\351ro %.0f sur %.0f -- Prochaine pause dans %.0f \351coutes\n',i_current,cfg_game.N,N_for_next_stop);
        fprintf('\tVolume relatif de la voix : %.2f dB\n',expvar);
        
        if i_current>100
            bias_r1 = 100*sum(data_passation.n_responses(end-100+1:end)==1)/100;
            if cfg_game.is_simulation == 1
                bias_r2 = 100*sum(data_passation.n_responses(end-100+1:end)==2)/100;
                % Extra info for the simulations:
                fprintf('\tPourcentage de r\351ponses "%s" : %.0f %%\n',cfg_game.response_names{1},bias_r1);
                fprintf('\tPourcentage de r\351ponses "%s" : %.0f %%\n',cfg_game.response_names{2},bias_r2);
                fprintf('\tPourcentage de bonnes r\351ponses: %.0f %%\n',100*sum(data_passation.is_correct(end-100+1:end)==1)/100);
            else
                if bias_r1>60
                    fprintf('\tPourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{1});
                elseif bias_r1<40
                    fprintf('\tPourcentage de r\351ponses "%s" = %.0f %% (trop de "%s") \n',cfg_game.response_names{1},bias_r1,cfg_game.response_names{2});
                else
                    fprintf('\tPourcentages de r\351ponses "%s" et "%s" normaux \n',cfg_game.response_names{1},cfg_game.response_names{2});
                end
            end
        end
        fprintf('\n');
end
