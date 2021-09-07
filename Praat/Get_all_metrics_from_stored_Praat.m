function outs = Get_all_metrics_from_stored_Praat(dir_where,params)
% function outs = Get_all_metrics_from_stored_Praat(dir_where,params)
%
% Based on l20210713_AnalyseStims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filesF  = Get_filenames(dir_where,'*_F*.txt');
filesI  = Get_filenames(dir_where,'*_I*.txt');
filesf0 = Get_filenames(dir_where,'*_f0*.txt');

Nsounds = length(filesF);

for i = 1:Nsounds
    
    [t_f0{i},f0{i}] = Get_f0_from_txt([dir_where filesf0{i}]);
    [t_I{i} ,I{i}]  = Get_intensity_from_txt([dir_where filesI{i}]);
    [t_F{i} ,F{i}]  = Get_formants_from_txt([dir_where filesF{i}]);
    
    minIforF = max(max(I{i})-20, params.I_min);
    
    idxs = find(I{i}<minIforF | isnan(I{i}));
    F{i}(idxs,:) = nan;
    f0{i}(idxs) = nan;
    
end
    
outs.t_f0 = t_f0;
outs.f0 = f0;
outs.t_I = t_I;
outs.I = I;
outs.t_F = t_F;
outs.F = F;

outs.filesf0 = filesf0;
outs.filesF = filesF;
outs.filesI = filesI;