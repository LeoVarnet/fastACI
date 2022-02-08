function fname = Save_crosspred(fnameACI,crosspred,script_str)
% function fname = Save_crossprediction(fnameACI,crosspred,script_str)
%
% Stores the crossprediction results into a folder named after fnameACI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info_toolbox = Get_toolbox_info(script_str);

[dir_where,fname] = fileparts(fnameACI);
dir_where = [dir_where filesep fname filesep]; mkdir(dir_where); 

fname = [dir_where 'Crosspred'];
if nargout == 0
    if ~exist(fname,'file')
        save(fname,'crosspred','info_toolbox');
    else
        error('Program here the function to append ACIs...')
    end
else
    fprintf('%s: Only getting the cross-prediction filename\n',upper(mfilename));
end