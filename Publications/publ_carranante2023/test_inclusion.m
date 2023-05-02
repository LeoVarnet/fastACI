function [SNRthres,is_included] = test_inclusion(expename, participantname)

dir_data = fastACI_dir_data();
dir_savegame = [dir_data expename filesep participantname filesep 'Results' filesep];
D = dir([dir_savegame 'savegame*.mat']);
load([dir_savegame D(1).name])
r = Get_mAFC_reversals(data_passation.expvar(1:400));
SNRthres = median(r(5:end));
is_included = (SNRthres<-12);
fprintf(['Participant ' participantname ' obtained a SNR threshold of ' num2str(SNRthres) ' dB in the first session.\n'])
if is_included
    fprintf(['This participant can therefore be included in the experiment.\n'])
else
    fprintf(['This participant should therefore be rejected.\n'])
end
end