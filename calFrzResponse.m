function myCellMetrics = calFrzResponse(myCellMetrics)
% To classify neuronal responses during freezing and non-freezing periods,
% a permutation test approach was employed.
total_startT = tic;
frzResp = string;
for  i = 1:length(myCellMetrics.unitID)
    spkData = myCellMetrics.spkData{i};
    onsetTime = myCellMetrics.onsetTime(i);
    frzTime = myCellMetrics.frzTime(i);
    try
        frzResp(i,:) = frzResponse_singleUnit(spkData, onsetTime, frzTime , 0);
    catch ME
        detListTable = detUnit(myCellMetrics, i);
        figure
        histogram(spkData(spkData>onsetTime.hcA & spkData<onsetTime.ctxA+300)-onsetTime.hcA+onsetTime.ctxB+300, 'BinWidth',5); hold on
        histogram(spkData(spkData>onsetTime.hcB & spkData<onsetTime.ctxB+300)-onsetTime.hcB, 'BinWidth',5); hold off
        legend('ctxB', 'ctxA')
        title(detListTable.dayInfo, 'Interpreter','none');
        subtitle(['UnitIdx: ',num2str(detListTable.unitIdx), '. cluID: ',num2str(detListTable.cluID), '. unitID: ',num2str(detListTable.unitID)])
        fprintf('An error occurred: %s\n', ME.message);
        frzResp(i,:) = "invalid";
    end
end
myCellMetrics.frzResponse = frzResp;

total_endT = toc(total_startT);
disp(['CalFrzResp: ', char(myCellMetrics.dayInfo(1)) ,'. Total time: ' num2str(total_endT), ' s.']);
end

%% ----------------------------Function------------------------------- %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resp = frzResponse_singleUnit(spike, evt, frzTime, plotChoice)
nboot = 1e3;
win = 300; % in sec

spike_ctxA = spike(spike>evt.ctxA & spike<evt.ctxA+win) - evt.ctxA;
spike_ctxB = spike(spike>evt.ctxB & spike<evt.ctxB+win) - evt.ctxB + win; 
spike_ctx = [spike_ctxA; spike_ctxB];

isi = [0; diff(spike_ctx)];
isi_resmp = zeros(length(isi), nboot);
for i = 1:nboot
    isi_resmp(:,i) = randsample(isi, length(isi)); 
end
spk_resamp = cumsum(isi_resmp,1); 

t_frz_A = frzTime.frz_ctxA - evt.ctxA;
t_nonFrz_A = frzTime.nonFrz_ctxA - evt.ctxA;
t_frz_B = frzTime.frz_ctxB - evt.ctxB + win;
t_nonFrz_B = frzTime.nonFrz_ctxB - evt.ctxB + win;

t_frz = [t_frz_A; t_frz_B];
t_nonFrz = [t_nonFrz_A; t_nonFrz_B];

frz_resamp = zeros(size(t_frz,1),nboot); 
frz_ori = zeros(size(t_frz,1),1);

for j = 1:size(t_frz,1)
    frz_be = t_frz(j, 1);
    frz_en = t_frz(j, 2);
    frz_resamp(j,:) = sum(spk_resamp>=frz_be & spk_resamp<=frz_en, 1); % 对重排的spike序列，计算给定时间窗的spike count
    frz_ori(j,:) =  sum(spike_ctx>=frz_be & spike_ctx<=frz_en, 1);
end
t_frz_sum = sum((t_frz(:,2)-t_frz(:,1)),1);
frz_resamp_mean = sum(frz_resamp,1)/t_frz_sum;
frz_ori_mean = sum(frz_ori,1)/t_frz_sum;

nonFrz_resamp = zeros(size(t_nonFrz,1),nboot); 
nonFrz_ori = zeros(size(t_nonFrz,1),1);
for j = 1:size(t_nonFrz,1)
    nonFrz_be = t_nonFrz(j, 1);
    nonFrz_en = t_nonFrz(j, 2);
    nonFrz_resamp(j,:) = sum(spk_resamp>=nonFrz_be & spk_resamp<=nonFrz_en, 1); 
    nonFrz_ori(j,:) =  sum(spike_ctx>=nonFrz_be & spike_ctx<=nonFrz_en, 1);
end
t_nonFrz_sum = sum((t_nonFrz(:,2)-t_nonFrz(:,1)),1);
nonFrz_resamp_mean = sum(nonFrz_resamp,1)/t_nonFrz_sum;
nonFrz_ori_mean = sum(nonFrz_ori,1)/t_nonFrz_sum;

delta_resamp = nonFrz_resamp_mean-frz_resamp_mean;
delta_orig = nonFrz_ori_mean-frz_ori_mean;

alpha = 0.05; 
pct_low = alpha*100/2; 
pct_high = 100-alpha*100/2;
rank_low = prctile(delta_resamp, pct_low); 
rank_high = prctile(delta_resamp, pct_high);
if delta_orig < rank_low
    resp = "frz_Activation";
elseif delta_orig > rank_high
    resp = "nonFrz_Activation";
else 
    resp = "none";
end

% Plot histogram
if plotChoice
    figure;
    subplot(2,1,1) 
    histogram(delta_resamp); 
    h1 = histfit(delta_resamp); 
    h1(2).Color = [.2 .2 .2];
    hold on
    plot([delta_orig delta_orig], ylim, 'Color', 'r', 'LineWidth', 0.8); hold on
    plot([rank_low rank_low], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    plot([rank_high rank_high], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    box off
    hold off
end
end

%%
