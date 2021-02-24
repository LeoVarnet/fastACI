function cfg_inout = modulationACI_seeds_init(cfg_inout)
% function cfg_inout = modulationACI_seeds_init(cfg_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end

% Initialises the seed numbers on the fly:
seed_number_max = 4*cfg_inout.N; % arbitrary number, seeds will range from 0 to 4*N
seed_numbers = round(seed_number_max*random('unif',0,1,[1,cfg_inout.N]));
                                % (allows repeted seed numbers)
cfg_inout.seeds_order = seed_numbers; % to be used sequentially