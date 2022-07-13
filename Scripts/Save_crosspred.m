function fname = Save_crosspred(fnameACI,crosspred,script_str,suffix_crosspred)
% function fname = Save_crossprediction(fnameACI,crosspred,script_str,suffix_crosspred)
%
% Stores the crossprediction results into a folder named after fnameACI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    suffix_crosspred = '';
else
    if ~isempty(suffix_crosspred)
        if ~strcmp(suffix_crosspred,'-')
            suffix_crosspred = ['-' suffix_crosspred];
        end
    end
end
info_toolbox = Get_toolbox_info(script_str);

[dir_where,fname] = fileparts(fnameACI);
dir_where = [dir_where filesep fname filesep]; mkdir(dir_where); 

fname = [dir_where 'Crosspred' suffix_crosspred];
if nargout == 0
    if ~exist(fname,'file')
        save(fname,'crosspred','info_toolbox');
    else
        error('Program here the function to append ACIs...')
    end
else
    fprintf('%s: Only getting the cross-prediction filename\n',upper(mfilename));
end