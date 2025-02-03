function [myCellMetrics, optoTag] = optoTaggingIdentification(myCellMetrics, varargin)
% The Stimulus-associated spike latency test (SALT) was employed to assess
% the light-following response of each neuron.
p = inputParser;
addParameter(p, 'calChoice', true, @islogical);   
addParameter(p, 'plotChoice', false, @islogical);
parse(p, varargin{:});
calChoice = p.Results.calChoice;
plotChoice = p.Results.plotChoice;

fs_spk = 2e4;
idx = propertyFilter(myCellMetrics); 
optoON_all = {myCellMetrics.optoTime.optoON}';
numTrials_1hz = 60;                
numTrials_10hz = 60;
numTrials_20hz = 60;
numTrials = numTrials_1hz+numTrials_10hz+numTrials_20hz;
optoON = cell(length(idx),1);
for i = 1:length(idx)
    optoON_initial = optoON_all{idx(i)}*fs_spk;
    optoON_1hz = find(diff(optoON_initial)==1/1*fs_spk,1);
    if isempty(optoON_1hz)
        optoON_1hz = 1;
    end
    optoON_10hz = find(diff(optoON_initial)==1/10*fs_spk,1); 
    if isempty(optoON_10hz)
        optoON_10hz = 61;
    end
    optoON_20hz = find(diff(optoON_initial)==1/20*fs_spk,1);
    if isempty(optoON_20hz)
        optoON_20hz = 661;
    end
    optoON_sel = optoON_initial([optoON_1hz:(optoON_1hz+numTrials_1hz-1), ...
                optoON_10hz:(optoON_10hz+numTrials_10hz-1)...
                optoON_20hz:(optoON_20hz+numTrials_20hz-1)])/fs_spk;
    optoON{i} = optoON_sel;
end

spike = myCellMetrics.spkData(idx);
plot_pre = 0.04;                         
plot_post = 0.04;
plot_bin = 0.001;                       
plot_preIdx = 1:plot_pre/plot_bin;
plot_postIdx = plot_pre/plot_bin+1:(plot_pre+plot_post)/plot_bin;
plot_edges = -plot_pre:plot_bin:plot_post;
plot_nbin = length(plot_edges)-1;

statsAlpha = 0.001;                
statsWin = 0.010;                   
statsBin = 0.001;                   
statsStep = 0.005;                 
stats_edgesBase = -statsWin:statsBin:0;
stats_edgesTest = 0:statsBin:statsWin;
stats_preIdx = 1:statsWin/statsBin;
stats_postIdx = 1:statsWin/statsBin;

tagType = string;
pValue = ones(length(idx),1);
spikeFR = cell(length(idx),1);
relativeSpikes =  cell(length(idx),1);
for i = 1:length(idx)            
    spike_now = spike{i};
    optoON_now = optoON{i};
    relativeSpikes_byTrial = cell(numTrials,1);
    spikeFR_byTrial = zeros(numTrials, plot_nbin);
    spikeBinary_base = zeros(numTrials, length(stats_edgesBase)-1);
    spikeBinary_test = zeros(numTrials, length(stats_edgesTest)-1);
    for t = 1:numTrials           
        trialStart = optoON_now(t)-plot_pre;
        trialEnd = optoON_now(t)+plot_post;
        trialSpikes = spike_now(spike_now >= trialStart & spike_now <= trialEnd); 
        relativeSpikes_now = trialSpikes - optoON_now(t);                         
        relativeSpikes_byTrial{t} = relativeSpikes_now;
        spikeFR_byTrial(t,:) = histcounts(relativeSpikes_now,plot_edges)/plot_bin;
        spikeBinary_base(t,:) = histcounts(relativeSpikes_now, stats_edgesBase);
        spikeBinary_test(t,:) = histcounts(relativeSpikes_now, stats_edgesTest);        
    end
    spikeFR{i} = spikeFR_byTrial;
    relativeSpikes{i} = relativeSpikes_byTrial;

    if calChoice
        [pValue(i), ~, spkLatency] = saltAddhlsi(spikeBinary_base,spikeBinary_test,statsBin,statsStep);
        if pValue(i)<statsAlpha 
            tagType(i) = "optoAct";
        else
            tagType(i) = "none";                               
        end
        optoTag(i).tagType = tagType(i);
        optoTag(i).pValue = pValue(i);
        optoTag(i).spkLatency = spkLatency;
        if i == length(idx) 
            myCellMetrics.optoTag = optoTag;
        end
    end
end
optoTag = myCellMetrics.optoTag;
tagType = cat(1,optoTag.tagType);
pValue = cat(1,optoTag.pValue);
spkLatency = {optoTag.spkLatency}';

%% plot by population
idxTag_optoAct = find(strcmp(tagType, "optoAct"));
idxTag_none = find(strcmp(tagType, "none"));
if plotChoice
    basecorrectWin = statsWin;             
    if ~isempty(idxTag_optoAct)
        spikeFR_sel_optoAct =spikeFR(idxTag_optoAct);
        spikeFR_popu_optoAct = cellfun(@(x) mean(x,1), spikeFR_sel_optoAct, 'UniformOutput',false);
        spikeFR_popuCor_optoAct = basecorrect(cell2mat(spikeFR_popu_optoAct)',1:plot_nbin,plot_preIdx(end)-basecorrectWin/plot_bin+1,plot_preIdx(end),'zscore')';
        spikeFR_popuCor_Sort_optoAct = sortFR(spikeFR_popuCor_optoAct, spikeFR_popuCor_optoAct(:,plot_postIdx(1):plot_postIdx(1)+round(0.005/plot_bin)-1));
        spikeFR_popuCor_smooth_optoAct = smoothdata(spikeFR_popuCor_Sort_optoAct,2, 'gaussian',5);
    else
        spikeFR_popuCor_smooth_optoAct = [];
    end
    if ~isempty(idxTag_none)
        spikeFR_sel_none =spikeFR(idxTag_none);
        spikeFR_popu_none = cellfun(@(x) mean(x,1), spikeFR_sel_none, 'UniformOutput',false);
        spikeFR_popuCor_none = basecorrect(cell2mat(spikeFR_popu_none)',1:plot_nbin,plot_preIdx(end)-basecorrectWin/plot_bin+1,plot_preIdx(end),'zscore')';
        spikeFR_popuCor_Sort_none = sortFR(spikeFR_popuCor_none, spikeFR_popuCor_none(:,plot_postIdx(1):plot_postIdx(1)+round(0.005/plot_bin)-1));
        spikeFR_popuCor_smooth_none = smoothdata(spikeFR_popuCor_Sort_none,2, 'gaussian',5);
    else
        spikeFR_popuCor_smooth_none = [];
    end
    
    figure; 
    ax01 = axes;
    ax01.OuterPosition = [0 0 1 1];
    ax01.Position = [0.1, 0.3, 0.35, 0.5];
    x = linspace(-plot_pre,plot_post,plot_nbin);
    plotShadedError(spikeFR_popuCor_smooth_optoAct,'x',x, 'color', [1.0000    0.6431    0.3451]); hold on
    plotShadedError(spikeFR_popuCor_smooth_none,'x',x, 'color', [0.5020 0.5020 0.5020]); hold on
    plot([0 0], ax01.YLim, '--', 'LineWidth', 1, 'Color',[0.5 0.5 0.5]); hold off
    box off
    ylabel('Z-score');
    xlabel('Time from light ON (s)');
    titleTXT = [num2str(myCellMetrics.dayInfo(idx(idxTag_none(1))))];
    sgtitle(titleTXT, 'Interpreter','none', 'Fontsize', 10);
    
    width = 0.35; 
    left = 0.52;
    figGap = 0.005;
    height1 = (0.6-2*figGap)*length(idxTag_optoAct)/(length(idxTag_none)+length(idxTag_optoAct)); % 第一个子图高度，optoAct
    height2 = height1*length(idxTag_none)/length(idxTag_optoAct);
    y1 = 0.85-height1;
    y2 = y1-height2-figGap;
    colorLim = [-1 3];
    
    if ~isempty(idxTag_optoAct)
        h1 = heatmap(spikeFR_popuCor_smooth_optoAct, 'Position', [left, y1, width, height1]);
        h1.GridVisible = 'off';
        h1.Colormap = colormap(flip(othercolor('RdYlBu4')));
        h1.ColorLimits = colorLim;
        x_idx = find(str2double(h1.XDisplayLabels)); h1.XDisplayLabels(x_idx) = {''}; % 删除heatmap的x，y轴标注
        y_idx = find(str2double(h1.YDisplayLabels)); h1.YDisplayLabels(y_idx) = {''};
        ylabel('optoAct');
        colorbar off
        ax1 = axes;
        ax1.Color = 'none'; 
        ax1.OuterPosition = [0 0 1 1];
        ax1.Position = [left, y1, width, height1];
        ax1.XLim = [0, plot_nbin]; 
        ax1.XTick = []; ax1.YTick = []; 
        hold on
        plot(ax1, [plot_preIdx(end) plot_preIdx(end)], ylim, '--', 'LineWidth', 1, 'Color', 'w');
        hold off
    end
    
    if ~isempty(idxTag_none)
        h2 = heatmap(spikeFR_popuCor_smooth_none, 'Position', [left, y2, width, height2]);
        h2.GridVisible = 'off';
        h2.Colormap = h1.Colormap;
        h2.ColorLimits = colorLim;
        x_idx = find(str2double(h2.XDisplayLabels)); h2.XDisplayLabels(x_idx) = {''}; % 删除heatmap的x，y轴标注
        y_idx = find(str2double(h2.YDisplayLabels)); h2.YDisplayLabels(y_idx) = {''};
        h2.XDisplayLabels(x_idx([1, plot_postIdx(1), plot_postIdx(end)])) = {['-' num2str(plot_pre)], '0', num2str(plot_post)};
        set(struct(h2).NodeChildren(3), 'XTickLabelRotation', 0); 
        xlabel('Time from light ON (s)');
        ylabel('No response');
        colorbar off
        ax2 = axes;
        ax2.Color = 'none'; 
        ax2.OuterPosition = [0 0 1 1];
        ax2.Position = [left, y1-height2-figGap, width, height2];
        ax2.XLim = [0, plot_nbin]; 
        ax2.XTick = []; ax2.YTick = []; 
        hold on
        plot(ax2, [plot_preIdx(end) plot_preIdx(end)], ylim, '--', 'LineWidth', 1, 'Color', 'w');
        hold off
    end
    
    c = colorbar;
    set(c, 'Position', [0.9 0.4 0.03 0.4]);
    clim(colorLim);
    
    ax4 = axes;
    ax4.Color = 'none'; 
    ax4.OuterPosition = [0 0 1 1];
    ax4.Position = [left, y1+height1+figGap, width, 0.015];
    optoMarkWidth = 0.005/plot_bin; 
    fill([plot_preIdx(end), plot_preIdx(end)+optoMarkWidth, plot_preIdx(end)+optoMarkWidth, plot_preIdx(end)], ...  
         [0, 0, 1, 1], ...
         [0.0745 0.6235 1.0000], 'FaceAlpha', 1.0, 'EdgeColor', 'none');
    ax4.XLim = [1 plot_nbin];
    ax4.YLim = [0 1];
    ax4.XTick = []; ax4.YTick = []; 
    axis off
    titleTXT = [num2str(myCellMetrics.dayInfo(idx(idxTag_none(1))))];
    sgtitle(titleTXT, 'Interpreter','none', 'Fontsize', 10);
    set(gca,'FontSize',8.5);
    
    ax5 = axes;
    ax5.Color = 'none'; 
    ax5.OuterPosition = [0 0 1 1]; 
    ax5.Position = [0.1, 0.3+0.5+figGap, 0.35, 0.015];
    optoMarkWidth = 0.005; 
    fill([0, optoMarkWidth, optoMarkWidth, 0], ...  
         [0, 0, 1, 1], ...
         [0.0745 0.6235 1.0000], 'FaceAlpha', 1.0, 'EdgeColor', 'none');
    ax5.XLim = [-plot_pre plot_post];
    ax5.YLim = [0 1];
    ax5.XTick = []; ax5.YTick = []; 
    axis off
    titleTXT = [num2str(myCellMetrics.dayInfo(idx(idxTag_none(1))))];
    sgtitle(titleTXT, 'Interpreter','none', 'Fontsize', 10);
    set(gca,'FontSize',8.5);
end
end