function myCellMetrics = calCtxSelectivity(myCellMetrics)
% Neuronal context-selectivity was determined using a permutation test to
% compare firing rate Z-scores between threatening and safe contexts.

ctxSelect = string;
total_startT = tic;
for  i = 1:length(myCellMetrics.unitID)
    psth_zscore = myCellMetrics.psth_zscore(i);
    ctxResponseA = myCellMetrics.ctxResponseA(i);
    ctxResponseB = myCellMetrics.ctxResponseB(i);
    try
        ctxSelect(i,:) = ctxSelect_singleUnit(psth_zscore.ctxA, psth_zscore.ctxB, 0);
        if strcmpi(ctxSelect(i,:), "ctxA") && strcmpi(ctxResponseA, "activation")
            ctxSelect(i,:)= "ctxA";
        elseif strcmpi(ctxSelect(i,:), "ctxB") && strcmpi(ctxResponseB, "activation")
            ctxSelect(i,:)= "ctxB";
        else
            ctxSelect(i,:)= "none";
        end
    catch ME
        fprintf('An error occurred: %s\n', ME.message);
        keyboard
    end
end
myCellMetrics.ctxSelect = ctxSelect;
total_endT = toc(total_startT);
disp(['CalCtxSelectivity: ', char(myCellMetrics.dayInfo(1)) ,'. Total time: ' num2str(total_endT), ' s.']);

%% ----------------------------Function------------------------------- %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctxSelect = ctxSelect_singleUnit(FR_ctxA, FR_ctxB, plotChoice)
FR_ctxA = FR_ctxA(31:end); 
FR_ctxB = FR_ctxB(31:end);
FR_ctx = [FR_ctxA, FR_ctxB];
nboot = 1e3;
FR_A_perm = zeros(30,nboot);
FR_B_perm = zeros(30,nboot);
for j = 1:nboot
    ind_A = randsample(1:60, 30);
    ind_B = setdiff(1:60, ind_A);
    FR_A_perm(:,j) = FR_ctx(ind_A);
    FR_B_perm(:,j) = FR_ctx(ind_B);
end
delta_orig = mean(FR_ctxA)-mean(FR_ctxB); % A-B
delta_perm = mean(FR_A_perm,1)-mean(FR_B_perm,1);
alpha = 0.05; 
pct_low = alpha*100/2; 
pct_high = 100-alpha*100/2;
rank_low = prctile(delta_perm, pct_low);
rank_high = prctile(delta_perm, pct_high);
if delta_orig < rank_low 
    ctxSelect = "ctxB";
elseif delta_orig > rank_high
    ctxSelect = "ctxA";
else 
    ctxSelect = "none";
end

% Plot histogram
if plotChoice
    figure
    histogram(delta_perm); 
    h1 = histfit(delta_perm); 
    h1(2).Color = [.2 .2 .2];
    hold on
    plot([delta_orig delta_orig], ylim, 'Color', 'r', 'LineWidth', 0.8, 'LineStyle','-'); hold on
    plot([rank_low rank_low], ylim, 'Color', 'k', 'LineWidth', 0.8, 'LineStyle','-'); 
    plot([rank_high rank_high], ylim, 'Color', 'k', 'LineWidth', 0.8, 'LineStyle','-'); 
    box off
    hold off
end

end
end

