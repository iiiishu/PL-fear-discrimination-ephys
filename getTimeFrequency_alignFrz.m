function getTimeFrequency_alignFrz(file1, path1, f_lim, savePath)
% Time-frequency analysis of LFPs during freezing and non-freezing periods
% was performed using Wavelet Transform. 

load(fullfile(path1, [file1, '_goodlfp_select.mat']));
disp(['Loaded: ', file1]);
load(fullfile(path1, 'evtTime.mat'));
load(fullfile(path1, 'frzTime.mat'));

frzTime_ctxA = frzTime.frz_ctxA-evtTime.ctxA;
nonFrzTime_ctxA = frzTime.nonFrz_ctxA-evtTime.ctxA;
frzTime_ctxB = frzTime.frz_ctxB-evtTime.ctxB;
nonFrzTime_ctxB = frzTime.nonFrz_ctxB-evtTime.ctxB;

fs_lfp = 1250;
timescale = linspace(0, size(lfpData,1)/fs_lfp, size(lfpData,1)); %改数
tic
lfp_spectro0 = [];
for i = 1: size(lfpData,2)
    lfpDataSelect = lfpData(:,i);    
    lfp_spectro0(:,:,i) = abs(awt_freqlist(lfpDataSelect,fs_lfp,f_lim(1):f_lim(2)));
end

lfp_spectro = squeeze(mean(lfp_spectro0,3));
baseWin = [min(evtTime.hcA, evtTime.hcB), min(evtTime.hcA, evtTime.hcB)+180];
lfp_spectro_basecorr = basecorrect(lfp_spectro, timescale, baseWin(1), baseWin(2),'zscore');  % 设定在hc1的起始2 min范围作基线校准
y0 = smoothdata(lfp_spectro_basecorr, 1, 'gaussian', 20*1250);
y = smoothdata(y0, 2, 'gaussian', 2);
fprintf('Cal spectral: done!\n')
toc

colorLim = [-0.5,1.5];
t_idx = ((timescale>=evtTime.ctxA) & (timescale<=evtTime.ctxA+300));
t_plot = timescale(t_idx)-evtTime.ctxA;
figure();
imagesc(t_plot,f_lim(1):f_lim(2),y(t_idx,:)');
axis xy; axis tight;
xlabel('Time from context switch (s)');
ylabel('Frequency (Hz)');
colorbar;
caxis(colorLim);
hold on
plotSpectralFrz(frzTime_ctxA, nonFrzTime_ctxA, f_lim);
hold off
set(gca, 'Position', [0.1,0.25,0.7,0.45]); 
title(['timeFreq_ctxA_' file1], 'Interpreter','none');
annotation('textbox', [0.957142857142857 0.393217685638299 0.123469391352234 0.0659863959943588],...
    'String',{'Z-score'},...
    'Rotation',90,...
    'FontSize',10.45,...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8);
print(gcf, '-dpdf', '-painters', [savePath file1 '_timeFreq_ctxA_' '.pdf']);

t_idx2 = ((timescale>=evtTime.ctxB) & (timescale<=evtTime.ctxB+300));
t_plot2 = timescale(t_idx2)-evtTime.ctxB;
figure();
imagesc(t_plot2,f_lim(1):f_lim(2),y(t_idx2,:)');
axis xy; axis tight;
xlabel('Time from context switch (s)');
ylabel('Frequency (Hz)');
colorbar;
caxis(colorLim);

hold on
plotSpectralFrz(frzTime_ctxB, nonFrzTime_ctxB, f_lim);
hold off
set(gca, 'Position', [0.1,0.25,0.7,0.45]); 
title(['timeFreq_ctxB_' file1], 'Interpreter','none');
annotation('textbox', [0.957142857142857 0.393217685638299 0.123469391352234 0.0659863959943588],...
    'String',{'Z-score'},...
    'Rotation',90,...
    'FontSize',10.45,...
    'FitBoxToText','off',...
    'EdgeColor','none');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8);
print(gcf, '-dpdf', '-painters', [savePath file1 '_timeFreq_ctxB_' '.pdf']);


end

