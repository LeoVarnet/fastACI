function Response = Response_keyboard( Liste_Reponses, cfg, varargin )
% function Reponse = Response_keyboard( Liste_Reponses, cfg, varargin )
%
% Response_keyboard waits for the response of the user (the participant) 
% according to the specifications and parameters of the experiment; 
% possibility of specifying the responses 'cheat code' and additional 
% argument, for instance [5,45]
%
% Old name: Reponse_clavier.m

Liste_Cheatcode = [];
if ~isempty(varargin)
    if length(varargin)==1
        Liste_Cheatcode = varargin{1};
    else
        error('Too many additional arguments');
    end
end

test_boucle=1;
while test_boucle
    Response=[];
    
    switch cfg.Language
        case 'EN'
            fprintf('    Press:\n');
        case 'FR'
            fprintf('    Appuyez sur :\n');
    end
    for idx=1:length(Liste_Reponses)
        switch cfg.Language
            case 'EN'
                fprintf(['     - ' num2str(idx) ' ' Liste_Reponses{idx} '\n']);
            case 'FR'
                fprintf(['     - ' num2str(idx) ' pour ' Liste_Reponses{idx} '\n']);
        end
        if idx == length(Liste_Reponses)
            switch cfg.Language
                case 'EN'
                    Response = input('    then press Enter\n','s');
                case 'FR'
                    fprintf('    puis appuyez sur Entr\351e');
                    Response = input('\n','s');
            end
        end
    end
    Response=str2double(Response);
    test_boucle = isempty(Response);
    if ~test_boucle
        test_touche=0;
        for i=1:length(Liste_Reponses)
            test_touche=test_touche||Response==i;
        end
        for i=1:length(Liste_Cheatcode)
            test_touche=test_touche||Response==Liste_Cheatcode(i);
        end
        test_boucle=~test_touche;
    end
end