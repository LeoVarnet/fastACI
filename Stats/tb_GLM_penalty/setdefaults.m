function [opts] = setdefaults(opts,defaults,nowarning)
    fn = fieldnames(defaults);
    for ii = 1:length(fn)
        if ~isfield(opts,fn{ii}) || isempty(opts.(fn{ii}))
            opts.(fn{ii}) = defaults.(fn{ii});
        end
    end
    
    if nargin < 3 || nowarning == 0
        fn = fieldnames(opts);
        for ii = 1:length(fn)
            if ~isfield(defaults,fn{ii});
                warning('setdefaults:unknownarg','Argument %s is unsupported',fn{ii});
            end
        end
    end
end