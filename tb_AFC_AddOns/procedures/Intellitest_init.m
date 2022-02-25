% Intellitest_init - init script for VlMatrix procedure

% change save function to method specific default if default save function is requested
if ( strcmp( def.savefcn, 'default' ) )
    def.savefcn = 'Intellitest_srt';
end

% sanity check
if ( def.allowpredict ~= 0 )
	def.allowpredict = 0; % not implemented for OLSA
	warning('def.allowpredict must equal 0 for OLSA or other speech matrials');
end

% prepare for interleaving
if ~iscell(def.TargetIntelligibility)
    expvarTmp = def.TargetIntelligibility;
    def=rmfield(def,'TargetIntelligibility');

    if ( def.interleaved )
        if size( expvarTmp, 1 ) == 1
            expvarTmp = repmat(expvarTmp, def.interleavenum, 1);
        elseif size( expvarTmp, 1 ) ~= def.interleavenum
            error('afc_main: def.expvar dimensions mismatch def.interleavenum');
        end
    end

    for ( i=1:def.interleavenum )
        def.TargetIntelligibility{i} = expvarTmp(i);
    end
end

for (i=1:def.interleavenum)

    % initialize starting value for all interleaved tracks
    work.expvarnext{i} = def.startvar{i};				% pre buffer (cell array) for expvaract, copied to expvaract in afc_interleave

    % add procedure specific entries to structure work
    work.OLSA_Direction{i} = [];
    work.OLSA_PercentCorrect{i} = 0;
    work.OLSA_SelectedWords{i} = [];
    work.OLSA_PresentedWords{i} = [];

end

% eof
