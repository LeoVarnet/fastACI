function definput = arg_Gammatone_proc(definput)

definput.flags.fc_nextpow2 = {'no_fc_nextpow2','fc_nextpow2'}; % added on 30/04/2023 for lasso GLMs
definput.keyvals.bwmul    = .5; % spaced at 0.5 ERB
definput.keyvals.binwidth = 1e-3; % s
