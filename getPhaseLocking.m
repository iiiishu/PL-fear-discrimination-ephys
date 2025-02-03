function  getPhaseLocking(file1, path1)
tic
load(fullfile(path1, [file1, '_goodlfp_select.mat']));
disp(['Loaded: ', file1]);
load(fullfile(path1, 'evtTime.mat'));
load(fullfile(path1, 'frzTime.mat'));
load(fullfile(path1, [file1, '_myCellMetrics.mat']));
spkData = myCellMetrics.spkData;

clear('band_def','band_name');
cnt=1;
band_def{cnt}=[2,6];band_name{cnt}='4Hz';cnt=cnt+1;
band_def{cnt}=[6,12];band_name{cnt}='theta';cnt=cnt+1;
band_def{cnt}=[30,90];band_name{cnt}='gamma';
bandsel = 3;
pass_band=band_def{bandsel};

ctxField = 'ctxB';
win_ctxB =[evtTime.(ctxField), evtTime.(ctxField)+duration];
[spikephase_ctxB, pval_ctxB] = getPhaseLockingBatch_whole(spkData, lfpData, win_ctxB, 'freqband', pass_band);

ctxField = 'ctxA';
win_ctxA =[evtTime.(ctxField), evtTime.(ctxField)+duration]; 
[spikephase_ctxA, pval_ctxA] = getPhaseLockingBatch_whole(spkData, lfpData, win_ctxA, 'freqband', pass_band);

[spikephase_all, pval_all] = getPhaseLockingBatch_whole(spkData, lfpData, [win_ctxA; win_ctxB], 'freqband', pass_band);

frzTime_ctxA = frzTime.frz_ctxA;
nonFrzTime_ctxA = frzTime.nonFrz_ctxA;
frzTime_ctxB = frzTime.frz_ctxB;
nonFrzTime_ctxB = frzTime.nonFrz_ctxB;

frzTime_all = [frzTime_ctxA; frzTime_ctxB];
nonFrzTime_all = [nonFrzTime_ctxA; nonFrzTime_ctxB];

winCut = 1;
d_nonFrz = nonFrzTime_all(:,2)-nonFrzTime_all(:,1);
d_frz = frzTime_all(:,2)-frzTime_all(:,1);


nonFrzTime_sel = nonFrzTime_all(d_nonFrz>winCut,:);
[spikephase_nonFrz, pval_nonFrz] = getPhaseLockingBatch_whole(spkData, lfpData, nonFrzTime_sel, 'freqband', pass_band);
frzTime_sel = frzTime_all(d_frz>winCut,:);
[spikephase_frz, pval_frz] = getPhaseLockingBatch_whole(spkData, lfpData, frzTime_sel, 'freqband', pass_band);

spikephase = struct('all', spikephase_all, 'pval_all', pval_all,...
                    'ctxA', spikephase_ctxA, 'pval_ctxA', pval_ctxA, 'ctxB', spikephase_ctxB, 'pval_ctxB', pval_ctxB, ...
                    'frz', spikephase_frz, 'pval_frz', pval_frz, 'nonFrz', spikephase_nonFrz, 'pval_nonFrz', pval_nonFrz);
fieldName = join(['spikephase_', band_name{bandsel}]);
myCellMetrics.(fieldName) = spikephase;
save(fullfile('G:\LST\preprocData\DCtra_PhaseLocking\gamma\', [file1, '_myCellMetrics.mat']), "myCellMetrics");
toc
end