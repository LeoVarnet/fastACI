function cfg_final = set_default_cfg( cfg_init, varargin )
% cfg_final = set_default_cfg( cfg_init, varargin )
% Completes the structure array cfg_init with field names and associated
% values specified in varargin (fills empty fields and creates non-existing fields)

cfg_final = cfg_init;
while length(varargin)>=2
    if ~isfield(cfg_final, varargin{1}) || isempty(getfield(cfg_final, varargin{1}))
        cfg_final = setfield(cfg_final, varargin{1}, varargin{2});
    end
    varargin = varargin(3:end);
end

end