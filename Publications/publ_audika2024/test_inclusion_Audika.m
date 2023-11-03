function [SNRthres,is_included] = test_inclusion_Audika(participantname)

expename = 'speechACI_Audika-adga-S43M';
dir_data = fastACI_dir_data();
dir_savegame = [dir_data expename filesep participantname filesep 'Results' filesep];
D = dir([dir_savegame 'savegame*.mat']);
if isempty(D)
    error(['No savegame for this participant in experiment ' expename])
end
load([dir_savegame D(1).name])
r = Get_mAFC_reversals(data_passation.expvar(1:400));
SNRthres = median(r(5:end));
is_included = (SNRthres<-12);
fprintf(['Participant ' participantname ' obtained a SNR threshold of ' num2str(SNRthres) ' dB in the first session.\n'])
end