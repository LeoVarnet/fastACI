function [] = check_cfg( cfg, varargin )
%[] = check_cfg( cfg )
%checks in the structure array cfg for fields with names specified in the
%cell structure parameters that do not exist or are empty. If so, the
%function displays an error message.

while length(varargin)>=1
    if ~isfield(cfg, varargin{1})
        error(['error : field ' varargin{1} ' is missing in the structure array ' inputname(1)])
    elseif isempty(getfield(cfg, varargin{1}))
        error(['error : field ' varargin{1} ' is empty in the structure array ' inputname(1)])
    end
    varargin=varargin(2:end);
end

end