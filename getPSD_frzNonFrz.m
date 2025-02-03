function getPSD_frzNonFrz(file1, path1)
tic
load(fullfile(path1, [file1, '_goodlfp_select.mat']));
disp(['Loaded: ', file1]);
load(fullfile(path1, 'evtTime.mat'));
load(fullfile(path1, 'frzTime.mat'));

fs_lfp = 1250;
frzTime_ctxA = frzTime.frz_ctxA*fs_lfp;
nonFrzTime_ctxA = frzTime.nonFrz_ctxA*fs_lfp;
frzTime_ctxB = frzTime.frz_ctxB*fs_lfp;
nonFrzTime_ctxB = frzTime.nonFrz_ctxB*fs_lfp;

winCut = 10; 
frzTime_select_ctxA = round(frzTime_ctxA((frzTime_ctxA(:,2)-frzTime_ctxA(:,1))>=winCut*fs_lfp,:));
nonFrzTime_select_ctxA = round(nonFrzTime_ctxA((nonFrzTime_ctxA(:,2)-nonFrzTime_ctxA(:,1))>=winCut*fs_lfp,:));
frzTime_select_ctxB = round(frzTime_ctxB((frzTime_ctxB(:,2)-frzTime_ctxB(:,1))>=winCut*fs_lfp,:));
nonFrzTime_select_ctxB = round(nonFrzTime_ctxB((nonFrzTime_ctxB(:,2)-nonFrzTime_ctxB(:,1))>=winCut*fs_lfp,:));


if isempty(frzTime_select_ctxA)
    winCut = 1; frzTime_select_ctxA = round(frzTime_ctxA((frzTime_ctxA(:,2)-frzTime_ctxA(:,1))>=winCut*fs_lfp,:));
end
if isempty(nonFrzTime_select_ctxA)
    winCut = 1; nonFrzTime_select_ctxA = round(nonFrzTime_ctxA((nonFrzTime_ctxA(:,2)-nonFrzTime_ctxA(:,1))>=winCut*fs_lfp,:));
end
if isempty(frzTime_select_ctxB)
    winCut = 1; frzTime_select_ctxB = round(frzTime_ctxB((frzTime_ctxB(:,2)-frzTime_ctxB(:,1))>=winCut*fs_lfp,:));
end
if isempty(nonFrzTime_select_ctxB)
    winCut = 1; nonFrzTime_select_ctxB = round(nonFrzTime_ctxB((nonFrzTime_ctxB(:,2)-nonFrzTime_ctxB(:,1))>=winCut*fs_lfp,:));
end

rezPSD_nonFrz_ctxA = getPSDseg_singleFile_unequal(lfpData, nonFrzTime_select_ctxA);
rezPSD_nonFrz_ctxA.nonFrzTime_select_ctxA = nonFrzTime_select_ctxA;
rezPSD_frz_ctxA = getPSDseg_singleFile_unequal(lfpData, frzTime_select_ctxA);
rezPSD_frz_ctxA.frzTime_select_ctxA = frzTime_select_ctxA;

rezPSD_nonFrz_ctxB = getPSDseg_singleFile_unequal(lfpData, nonFrzTime_select_ctxB);
rezPSD_nonFrz_ctxB.nonFrzTime_select_ctxB = nonFrzTime_select_ctxB;
rezPSD_frz_ctxB = getPSDseg_singleFile_unequal(lfpData, frzTime_select_ctxB);
rezPSD_frz_ctxB.frzTime_select_ctxB = frzTime_select_ctxB;

PSDforband_frznonFrz_ctxAB = [mean(rezPSD_frz_ctxA.PSDforband,2),mean(rezPSD_nonFrz_ctxA.PSDforband,2),...
                              mean(rezPSD_frz_ctxB.PSDforband,2), mean(rezPSD_nonFrz_ctxB.PSDforband,2)];
relative_PSD_frznonFrz_ctxAB = [mean(rezPSD_frz_ctxA.relative_PSD,2),mean(rezPSD_nonFrz_ctxA.relative_PSD,2),...
                                mean(rezPSD_frz_ctxB.relative_PSD,2), mean(rezPSD_nonFrz_ctxB.relative_PSD,2)];

PSDforband_frznonFrz_ctxAB = array2table(PSDforband_frznonFrz_ctxAB, 'RowNames', {'delta (2-6 Hz)', 'theta (6-12 Hz)','alpha (8-13 Hz)', 'beta (13-30 Hz)', 'gamma (30-90 Hz)', 'high gamma (50-100 Hz)'}, ...
                            'VariableNames', {'rezPSD_frz_ctxA', 'rezPSD_nonFrz_ctxA', 'rezPSD_frz_ctxB', 'rezPSD_nonFrz_ctxB'});
save(fullfile(path1, [file1, '_PSDforband_frznonFrz_ctxAB.mat']), "PSDforband_frznonFrz_ctxAB");
disp(['Saved: ', file1, '_PSDforband_frznonFrz_ctxAB.mat']);

relative_PSD_frznonFrz_ctxAB = array2table(relative_PSD_frznonFrz_ctxAB, 'RowNames', {'delta (2-6 Hz)', 'theta (6-12 Hz)','alpha (8-13 Hz)', 'beta (13-30 Hz)', 'gamma (30-90 Hz)', 'high gamma (50-100 Hz)'}, ...
                            'VariableNames', {'rezPSD_frz_ctxA', 'rezPSD_nonFrz_ctxA', 'rezPSD_frz_ctxB', 'rezPSD_nonFrz_ctxB'});
save(fullfile(path1, [file1, '_relative_PSD_frznonFrz_ctxAB.mat']), "relative_PSD_frznonFrz_ctxAB");
disp(['Saved: ', file1, '_relative_PSD_frznonFrz_ctxAB.mat']);

rezPSD = struct;
rezPSD.rezPSD_nonFrz_ctxA = rezPSD_nonFrz_ctxA;
rezPSD.rezPSD_frz_ctxA = rezPSD_frz_ctxA;
rezPSD.rezPSD_nonFrz_ctxB = rezPSD_nonFrz_ctxB;
rezPSD.rezPSD_frz_ctxB = rezPSD_frz_ctxB;
save(fullfile(path1, [file1, '_rezPSD_frznonFrz_ctxAB.mat']), "rezPSD");
disp(['Saved: ', file1, '_rezPSD_frznonFrz_ctxAB.mat']);
toc
end

function rezPSD = getPSDseg_singleFile_unequal(lfpData, t)
params.tapers=[3,5];
params.pad=1;
params.Fs=1250;
params.fpass=[1,100];
params.err=0;
params.trialave=1;
frequencyband={[2 6],[6 12],[8 13], [13 30],[30 90],[60 90]}; 
rezPSD = struct('f',{}, 'PSD',{}, 'relative_PSD',{}, 'PSDforband',{});
f = cell(size(t,1),1);
PSD = cell(size(t,1),1);
PSDforband = ones(length(frequencyband), size(t,1),1); 
RelativePSDforband =  ones(length(frequencyband), size(t,1),1);
for i = 1:size(t,1)
    [PSDtmp,ftmp]=mtspectrumc(lfpData(t(i,1):t(i,2),:),params);
    f{i} = ftmp;
    PSD{i} = PSDtmp;
    index = cellfun(@(x) find(ftmp>x(1)&ftmp<x(2)), frequencyband, 'UniformOutput',0);
    PSDindex = cellfun(@(x) trapz_xy(ftmp(x), PSDtmp(x)), index, 'UniformOutput',0);
    RelativePSDindex = cellfun(@(x) x/trapz_xy(ftmp,PSDtmp), PSDindex, 'UniformOutput',0);
    PSDforband(:,i) = cell2mat(PSDindex)';
    RelativePSDforband(:,i) = cell2mat(RelativePSDindex)';
end
rezPSD(1).f=f;
rezPSD(1).PSD=PSD;
rezPSD(1).PSDforband=PSDforband;
rezPSD(1).relative_PSD=RelativePSDforband;
end