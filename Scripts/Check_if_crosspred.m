function [crosspred, fname, bComplete] = Check_if_crosspred(fnameACI,file_crosspred)
% function crosspred = Check_if_crosspred(fnameACI,file_crosspred)
%
% Used in l20211217_crossprediction_version_Alejandro_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crosspred = [];

% info_toolbox = Get_toolbox_info(script_str);
if ischar(file_crosspred)
    file_crosspred = {file_crosspred};
end
N_crosspred = length(file_crosspred); % test ACIs for the cross prediction
idx = nan([1 N_crosspred]);

bComplete = 0;

[dir_where,fname] = fileparts(fnameACI); % reference ACI from which the waveforms/data are taken
dir_where = [dir_where filesep fname filesep]; 
if exist(dir_where,'dir')
    try
        cfg_file = [dir_where(1:end-1) '.mat'];
        var = load(cfg_file);
        idxlambda = var.results.idxlambda;
    end
    fname = [dir_where 'Crosspred.mat'];
    if exist(fname,'file')
        var = load(fname);
        N_existing = length(var.crosspred); % Number of fields in the stored cross prediction
        if N_existing ~= N_crosspred
            error('%s: Mismatch between number of requested cross predictions and number of fields in the existing cross prediction file',upper(mfilename));
        end
        for i = 1:N_crosspred % each requested file_crosspred
            for j = 1:N_existing % each file name in the crossprediction
                file1 = file_crosspred{i};
                file2 = var.crosspred(j).ACI_crosspred;
                if strcmp(file1,file2)
                    idx(i) = j;
                    % idx(i) = j;
                    break; % ends the for loop for 'j'
                else
                    try
                        % Then the files are 'different' but maybe because they
                        %   differ in their absolute path, or were obtained using
                        %   different operating systems.

                        % file1 is formatted using the current OS:
                        [xx,fname1,ext] = fileparts(file1); % fileparts will work
                        fname1 = [fname1 ext];
                        fname2 = file2(end+1-length(fname1):end);

                        fprintf('%s: Comparing the file names (%s) but not the absolute paths...\n',upper(mfilename),fname1);

                        if strcmp(fname1,fname2)
                            idx(i) = j;
                            % idx(j) = j;
                            break;
                        end
                    catch
                        if isempty(file2)
                            warning('No cross prediction found...')
                            idx(i) = nan;
                        end
                    end
                
                end
            end
        end
        
        crosspred = var.crosspred;
        
        if sum(isnan(idx))==0 
            bComplete = 1;
        else
            warning('This cross prediction file seems to be incomplete...')
        end
    end
end