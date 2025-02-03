function rez_decode = calCtxDecoding(myCellMetrics,filter)
nrepeat = 20;
idx_unit = propertyFilter(myCellMetrics, filter);
bin = 1; 
[~, psth_zscore, psth] = calFiringRate(myCellMetrics, 'Bin', bin);
FR_ctxB = [psth_zscore.ctxB];
FR_ctxA = [psth_zscore.ctxA];

duration = 300; 
idx_hc = 1: round(duration/bin);
idx_ctx = round(duration/bin)+1:round(duration*2/bin);
data_hc1 = FR_ctxB(idx_unit, idx_hc);
data_N = FR_ctxB(idx_unit, idx_ctx);
data_hc2 = FR_ctxA(idx_unit, idx_hc);
data_T = FR_ctxA(idx_unit, idx_ctx);
idx_ctx_cmb = round(duration/bin)-30/bin+1:round(duration*2/bin);
data_N_cmb = FR_ctxB(idx_unit, idx_ctx_cmb);
data_T_cmb = FR_ctxA(idx_unit, idx_ctx_cmb);


%% GPFA
FR1 = [ data_N];
FR2 = [ data_T];
trialId = [1 2];
trials_FR1 = struct('spikes', FR1, 'trialId', 1);
trials_FR2 = struct('spikes', FR2, 'trialId', 2);
all_trials = [trials_FR1; trials_FR2];

method = 'gpfa';  
xDim = 8;  
kernSD = 30;
binWidth = 10;   

result = neuralTraj(1, all_trials, 'method', method, ...
    'xDim', xDim, 'kernSDList', kernSD, 'useSqrt', false, 'binWidth',binWidth);
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
xspec = "xorth";
dimsToPlot = 1:3;
rez_gpfa = extractFields(seqTrain, ["trialId" , xspec]);

id1 = [rez_gpfa.trialId]==1;
id2 = [rez_gpfa.trialId]==2;
xorth = [rez_gpfa(id1).xorth, ...
         rez_gpfa(id2).xorth]';
data_decode = xorth(:,dimsToPlot);
lgh = size(data_decode,1);
label = trialId;

%% Create the label vector for your two stimuli
labels = reshape(repmat(label, [lgh/length(label),1]),[], 1); 

%% Randomly split the data into training and testing sets
% rng('default') % for reproducibility
roc_decoding = cell(nrepeat,1);
acc_rates = [];
predict_prob = [];
train_fraction = 0.7;
for i = 1:nrepeat
    idx_train = sort(randperm(lgh, round(train_fraction*lgh)));
    idx_test = setdiff(1:lgh, idx_train);
    data_train = data_decode(idx_train, :);
    data_test = data_decode(idx_test, :);
    label_train = labels(idx_train);
    label_test = labels(idx_test);

    % Train a SVM classifier on the training data
    mdl = fitcsvm(data_train, label_train, ...
        'KernelFunction', 'linear', ...
        'ClassNames', label, ...
        'BoxConstraint',0.3, 'KernelScale', 100);
    [label_predicted, scores] = predict(mdl, data_test);
    ac = [label_predicted, label_test];
    acc_rates(i) = sum(string(label_predicted) == string(label_test))/length(label_predicted)*100;
    fprintf('acc_rates =  %.2f\n', acc_rates(i));

    confu_matrix = confusionmat(label_test, label_predicted, 'Order', label);
    predict_prob(:,:,i) = confu_matrix./sum(confu_matrix,2)*100;

    labels_sfl = labels(randperm(length(labels)));
    label_tra_sfl = labels_sfl(idx_train);
    label_test_sfl = labels_sfl(idx_test);
    
    mdl_sfl = fitcsvm(data_train, label_tra_sfl, ...
        'KernelFunction', 'linear', ...
        'ClassNames', label, ...
        'BoxConstraint',0.5, 'KernelScale',100);
    [label_pred_sfl, scores_sfl] = predict(mdl_sfl, data_test);
    acc_rates_sfl(i) = sum(string(label_pred_sfl) == string(label_test_sfl))/length(label_pred_sfl)*100;
    fprintf('acc_rates_sfl =  %.2f\n', acc_rates_sfl(i));
    
    confu_matrix_shuffle = confusionmat(label_test_sfl, label_pred_sfl, 'Order', label);
    predict_prob_sfl(:,:,i) = confu_matrix_shuffle./sum(confu_matrix_shuffle,2)*100;     
    
end
rez_decode.acc_rates = acc_rates;
rez_decode.predict_prob = predict_prob;
rez_decode.acc_rates_sfl = acc_rates_sfl;
rez_decode.predict_prob_sfl = predict_prob_sfl;

predict_prob_sfl_mean = mean(predict_prob_sfl,3);
predict_prob_mean = mean(predict_prob,3);
labels_show = ["CtxD", "CtxA"];
figure
subplot(2,2,1)
h1 = heatmap(predict_prob_mean);
xlabel('Predicted label'); ylabel('True label');
h1.XDisplayLabels = labels_show; 
h1.YDisplayLabels = labels_show;
h1.Colormap = colormap(gradColorMap(256, [[60 58 200]; [50 124 252]; [9 189 189]; [255 190 54]; [249 249 22]]'/255)); 
h1.ColorLimits = [0 100]; 


subplot(2,2,2)
h2 = heatmap(predict_prob_sfl_mean);
xlabel('Predicted label'); ylabel('True label');
h2.XDisplayLabels = labels_show; 
h2.YDisplayLabels = labels_show;
h2.Colormap = h1.Colormap;
h2.ColorLimits = [0 100]; 

subplot(2,2,3); 
boxplot([acc_rates', acc_rates_sfl']);

end



