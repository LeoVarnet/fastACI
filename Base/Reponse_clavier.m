function [ Reponse ] = Reponse_clavier( Liste_Reponses, varargin )
%REPONSE_CLAVIER Attend une r�ponse de l'utilisateur parmis celles
%sp�cifi�es en param�tre; possibilit� de sp�cifier des r�ponses 'cheat
%code' en argument additionnel, par ex. [5,45]

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
    Reponse=[];
    
    %fprintf('    Press:\n');
    fprintf('    Appuyez sur :\n');
    for idx=1:length(Liste_Reponses)
        fprintf(['     - ' num2str(idx) ' ' Liste_Reponses{idx} '\n']);
        if idx == length(Liste_Reponses)
            %Reponse = input('    then press Enter\n','s');
            Reponse = input('    puis appuyez sur Entree\n','s');
        end
    end
    Reponse=str2double(Reponse);
    test_boucle = isempty(Reponse);
    if ~test_boucle
        test_touche=0;
        for i=1:length(Liste_Reponses)
            test_touche=test_touche||Reponse==i;
        end
        for i=1:length(Liste_Cheatcode)
            test_touche=test_touche||Reponse==Liste_Cheatcode(i);
        end
        test_boucle=~test_touche;
    end
end