function cfg_final = import_cfg( cfg_init, varargin )
% cfg_final = import_cfg ( cfg_init, parameters )
% Imports fields with names specified in the cell structure parameters
% from the structure array cfg_init to cfg_final. Used as import_cfg (
% cfg_init, cfg_basis, parameters ), fields from cfg_init are added to
% cfg_basis

cfg_final = [];
if isstruct(varargin{1})
    cfg_final = varargin{1};
    varargin=varargin(2:end);
end

while length(varargin)>=1
    if ischar(varargin{1})
        if isfield(cfg_init, varargin{1})
            cfg_final = setfield(cfg_final, varargin{1}, getfield(cfg_init, varargin{1}));
        else
            error(['importation error : ' varargin{1} ' is not the name of a field in the structure array ' inputname(1) '\n']);
        end
    else
        error(['importation error : field name specified in parameters must be strings \n']);
    end
    varargin=varargin(2:end);
end

end