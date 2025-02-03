function selecRatio = calSelectivityRatio(myCellMetrics, filter)
% To quantify the degree of context preference, a selectivity index was
% calculated for each neuron
duration = 300;
bin = 10; 
idx_ctx = round(duration/bin)+1:round(duration*2/bin);
idx = propertyFilter(myCellMetrics, filter);
choice = 'psth';
FR_all_ctxA = cell2mat({myCellMetrics.(choice).ctxA}');
FR_all_ctxB = cell2mat({myCellMetrics.(choice).ctxB}');
FR_ctxA = FR_all_ctxA(idx,idx_ctx);
FR_ctxB = FR_all_ctxB(idx,idx_ctx);
clear FR_all*

selecRatio = (mean(FR_ctxB,2)-mean(FR_ctxA,2))./(mean(FR_ctxB,2)+mean(FR_ctxA,2));
end