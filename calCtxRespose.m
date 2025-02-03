function myCellMetrics = calCtxRespose(myCellMetrics)
% To classify neuronal responses to contexts, a permutation test approach
% was employed, given that rats were exposed to each context only once
% during each test session

total_startT = tic;
ctxRespA = string;
ctxRespB = string;
for  i = 1:length(myCellMetrics.unitID)
     [ctxRespA(i,:), ctxRespB(i,:)] = ctxResponse_singleUnit(myCellMetrics.spkData{i}, myCellMetrics.onsetTime(i), 'randsample', 0);
end
FR_zscoreA = cell2mat({myCellMetrics.psth_zscore.ctxA}');
ctxRespA = zscoreCutoff(FR_zscoreA(:,31:60), ctxRespA);
FR_zscoreB = cell2mat({myCellMetrics.psth_zscore.ctxB}');
ctxRespB = zscoreCutoff(FR_zscoreB(:,31:60), ctxRespB);
myCellMetrics.ctxResponseA = ctxRespA;
myCellMetrics.ctxResponseB = ctxRespB;

total_endT = toc(total_startT);
disp(['CalCtxResp: ', char(myCellMetrics.dayInfo(1)) ,'. Total time: ' num2str(total_endT), ' s.']);

ctxType = string;
ctxResponseA = myCellMetrics.ctxResponseA;
ctxResponseB = myCellMetrics.ctxResponseB;

for j = 1:length(myCellMetrics.unitID)
    if strcmp(ctxResponseA(j), "activation") && ~strcmp(ctxResponseB(j), "activation")
        ctxType(j,1) = "ctxA activated";
    elseif strcmp(ctxResponseB(j), "activation") && ~strcmp(ctxResponseA(j), "activation")
        ctxType(j,1) = "ctxB activated";
    elseif strcmp(ctxResponseA(j), "activation") && strcmp(ctxResponseB(j), "activation")
        ctxType(j,1) = "Both activated";
    else
        ctxType(j,1) = "Others";
    end
end
myCellMetrics.ctxType = ctxType;
end


%% ----------------------------Function------------------------------- %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [respA, respB] = ctxResponse_singleUnit(spike, evt, permutationChoice, plotChoice)
nboot = 1e3;
win = 300;
switch permutationChoice
    case 'randsample' 
        transTime = 0; 
        spike_hcA = spike(spike>evt.hcA & spike<evt.hcA+win-transTime) -evt.hcA;
        spike_ctxA = spike(spike>evt.ctxA & spike<evt.ctxA+win) -evt.hcA-transTime;
        spikeA = [spike_hcA; spike_ctxA];
        isiA = [0; diff(spikeA)];
        evtA = [0 evt.ctxA-evt.hcA-transTime];
        isi_resmpA = zeros(length(isiA), nboot);
        for i1 = 1:nboot
            isi_resmpA(:,i1) = randsample(isiA, length(isiA)); 
        end
        spk_resamp1 = cumsum(isi_resmpA,1); 
        count_resampA = zeros(length(evtA),nboot); 
        count_oriA = zeros(length(evtA),1);
        for n1 = 1:2 
            count_resampA(n1,:) = sum(spk_resamp1>evtA(n1) & spk_resamp1<(evtA(n1)+win),1);
            count_oriA(n1,:) =  sum(spikeA>evtA(n1) & spikeA<evtA(n1)+win,1);
        end

        spike_hcB = spike(spike>evt.hcB & spike<evt.hcB+win-transTime) -evt.hcB;
        spike_ctxB = spike(spike>evt.ctxB & spike<evt.ctxB+win) -evt.hcB-transTime;
        spikeB = [spike_hcB; spike_ctxB];
        isiB = [0; diff(spikeB)];
        evtB = [0 evt.ctxB-evt.hcB-transTime];
        isi_resmpB = zeros(length(isiB), nboot);
        for i2 = 1:nboot
            isi_resmpB(:,i2) = randsample(isiB, length(isiB)); 
        end
        spk_resampB = cumsum(isi_resmpB,1);
        count_resampB = zeros(length(evtB),nboot); 
        count_oriB = zeros(length(evtB),1);
        for n2 = 1:2 
            count_resampB(n2,:) = sum(spk_resampB>evtB(n2) & spk_resampB<(evtB(n2)+win),1);
            count_oriB(n2,:) =  sum(spikeB>evtB(n2) & spikeB<evtB(n2)+win,1);
        end
end

alpha = 0.05; 
pct_low = alpha*100/2; 
pct_high = 100-alpha*100/2;
A_low = prctile(count_resampA(2,:), pct_low); 
A_high = prctile(count_resampA(2,:), pct_high);
if count_oriA(2) < A_low
    respA = "inhibition";
elseif count_oriA(2) > A_high
    respA = "activation";
else 
    respA = "none";
end

B_low = prctile(count_resampB(2,:), pct_low); 
B_high = prctile(count_resampB(2,:), pct_high);
if count_oriB(2) < B_low
    respB = "inhibition";
elseif count_oriB(2) > B_high
    respB = "activation";
else
    respB = "none";
end

% Plot histogram
if plotChoice
    figure;
    subplot(2,1,1) 
    histogram(count_resampA(2,:)); 
    h1 = histfit(count_resampA(2,:)); h1(2).Color = [.2 .2 .2];
    hold on
    plot([count_oriA(2) count_oriA(2)], ylim, 'Color', 'r', 'LineWidth', 0.8); hold on
    plot([A_low A_low], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    plot([A_high A_high], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    box off
    hold off
    subplot(2,1,2)
    histogram(count_resampB(2,:)); 
    h2 = histfit(count_resampB(2,:)); h2(2).Color = [.2 .2 .2];
    hold on
    plot([count_oriB(2) count_oriB(2)], ylim, 'Color', 'r', 'LineWidth', 0.8); hold on
    plot([B_low B_low], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    plot([B_high B_high], ylim, 'Color', 'k', 'LineWidth', 0.8); 
    box off
    hold off
end
end

function spk_resamp = spkResamp(spike, nboot, isi_index)
for ii = 1:nboot
    isi = [0; diff(spike)];
    isi_resamp = isi(isi_index(:,ii));
    spk_resamp(:,ii) = cumsum(isi_resamp);
end
end

function ctxResp = zscoreCutoff(FR_zscore, ctxResp)
for i = 1:size(FR_zscore, 1)
    if strcmp(ctxResp(i), "activation")
        above_threshold = FR_zscore(i,:) > 1.96;
        moreThan_count = length(find(above_threshold)) >= 3;
        consecutive_count = conv(above_threshold, ones(1,2), 'valid');
        has_consecutive_2 = any(consecutive_count >= 2);
        if  ~moreThan_count && ~has_consecutive_2
            ctxResp(i) = "none";
        end
    elseif strcmp(ctxResp(i), "inhibition")
        below_threshold = FR_zscore(i,:) < -1.96;
        lessThan_count = length(find(below_threshold)) >= 3;
        consecutive_count = conv(below_threshold, ones(1,2), 'valid');
        has_consecutive_2 = any(consecutive_count >= 2);
        if ~lessThan_count && ~has_consecutive_2
            ctxResp(i) = "none";
        end
    end
end
end