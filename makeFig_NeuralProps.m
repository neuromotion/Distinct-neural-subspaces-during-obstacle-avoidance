function subfuncs = makeFig_NeuralProps(varargin)

% Code to make single neural properties figure for paper
% 
% Layout:
% A--------------------------------------------------------  B-------------
% |  Boomer example neuron       Starbuck example neuron     | Boomer pop
% |                                                          | phase plot
% |  Raster (basic walk)         Raster (basic walk)         |
% |  Raster (walking obs step)   Raster (walking obs step)   | Starbuck pop
% |  PSTH (walking obs step)     PSTH (walking obs step)     | phase plot
% 
% C-------------------   D--------------------   E---------------------
% |  Boomer phase        |  Boomer depth of      |  Boomer mean FR
% |  shifts              |  mod shifts           |  shifts
% |                      |                       |
% |                      |                       |
% |  Starbuck phase      |  Starbuck depth of    |  Starbuck mean FR
% |  shifts              |  mod shifts           |  shifts
% |                      |                       |

% expose subfunctions in case any other figures want to use them
subfuncs.makeRasterAndPETH = @makeRasterAndPETH;
subfuncs.testUnimodal = @testUnimodal;
subfuncs.getNeuralProperties = @getNeuralProperties;

if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 11.5, 11.2])

figMargins = 0.05;
subplotGap = 0.1;

%% Subplot A -- Example neuron firing activity

% load in data
load('./Data/TrialsDataBoomer.mat')

% only using leg array
trialsDataBoomer = trialsLegM1;

% set example neuron
exampleNeuron = 43;

% get spike times
% find the indices for each task type
badTrials = filterTrials(trialsDataBoomer,90,5);
walkTrialInds = find(strcmpi('Walk',string({trialsDataBoomer.Task})));
walkTrialInds = setdiff(walkTrialInds, badTrials);

walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataBoomer.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataBoomer(walkTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataBoomer(walkObsTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);

% for later, when calculating spks/s also get the distribution of gait
% cycle durations for each step
% first basic walking
walkTrialDurations_Boomer = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkTrialInds).SpikeCounts})*10;

% then walking obstacle
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkObsTrialInds).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataBoomer(walkObsTrialInds).TrialEvents) walkObsTotalDurations];

steps_Boomer = -3:3;
walkObsTrialDurations_Boomer = diff(walkObsEvents(:,1:2:end),1,2)*10;

% get neural stats
[walkPETH{1}, walkObsPETH{1}, neuralStats{1}] = getNeuralProperties(trialsDataBoomer,...
    walkTrialInds, walkObsTrialInds, walkTrialDurations_Boomer, walkObsTrialDurations_Boomer, steps_Boomer);


% now do the actual plotting
colorMap = parula(8);
makeRasterAndPETH([figMargins, 5+subplotGap+figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), {'Basic Walking', 'Obstacle Stride'}, mean(walkTrialDurations_Boomer), mean(walkObsTrialDurations_Boomer(:,4)),...
    exampleNeuron, neuralStats{1});

% Now do Starbuck 

% load in data
load('./Data/TrialsDataStarbuck.mat')

% only using leg array
trialsDataStarbuck = trialsLegM1;

% set example neuron
exampleNeuron = 26;

% get spike trains
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkTrialInds = find(strcmpi('Walk',string({trialsDataStarbuck.Task})));
walkTrialInds = setdiff(walkTrialInds, badTrials);

walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataStarbuck.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataStarbuck(walkTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataStarbuck(walkObsTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);

% for later, when calculating spks/s also get the distribution of gait
% cycle durations for each step
% first basic walking
walkTrialDurations_Starbuck = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkTrialInds).SpikeCounts})*10;

% then walking obstacle
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkObsTrialInds).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataStarbuck(walkObsTrialInds).TrialEvents) walkObsTotalDurations];

steps_Starbuck = -3:2;
walkObsTrialDurations_Starbuck = diff(walkObsEvents(:,1:2:end),1,2)*10;

% get neural stats
[walkPETH{2}, walkObsPETH{2}, neuralStats{2}] = getNeuralProperties(trialsDataStarbuck,...
    walkTrialInds, walkObsTrialInds, walkTrialDurations_Starbuck, walkObsTrialDurations_Starbuck, steps_Starbuck);

makeRasterAndPETH([figMargins + 3.2 + subplotGap, 5+subplotGap+figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Starbuck),...
    mean(walkObsTrialDurations_Starbuck(:,4)), exampleNeuron, neuralStats{2});


%% Subplot B -- Population phase distribution

% Now, sort neurons by preferred phase in basic walking
% first Boomer
% only use unimodal neurons and neurons with dispersion above 0.15 and
% neurons who pass the rayleigh test
dispersionCutoff = 0.15;
goodNeurons = intersect(find(neuralStats{1}.isUnimodal(:,1)), find(neuralStats{1}.dispersion(:,1) > dispersionCutoff));
goodNeurons = intersect(goodNeurons, find(neuralStats{1}.p_circTestR(:,1)<(0.05/50)));

% now sort by preferred phase (ascending)
[~, indsPerm] = sort(neuralStats{1}.preferedPhase(goodNeurons,1));
sortedPETH = walkPETH{1}(goodNeurons(indsPerm), :);

% smooth out PETH by convolving with gaussian
countsSmooth = convGauss(repmat(sortedPETH, 1, 3), 10, 20);
countsSmooth = countsSmooth(:,length(sortedPETH)+1:length(sortedPETH)*2);

% finally normalize to max FR
walkPhaseBoomer = countsSmooth./repmat(max(countsSmooth,[],2), 1, size(countsSmooth,2));

% make plot
boomerWalkPhaseH = axes('Units','inches','OuterPosition',...
    [3.4*2 + subplotGap*2 - 0.2, 5 + subplotGap*2 + figMargins + 2.3, 2.2, 2.2]);

imagesc(walkPhaseBoomer)
ylabel('Neuron')
box off
set(boomerWalkPhaseH, 'YTick', [], 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XColor', 'none', 'XLim', [0 100])

% next do walking obstacle step
sortedPETH = squeeze(walkObsPETH{1}(goodNeurons(indsPerm), :, 4));
countsSmooth = convGauss(repmat(sortedPETH, 1, 3), 10, 20);
countsSmooth = countsSmooth(:,length(sortedPETH)+1:length(sortedPETH)*2);
walkObsPhaseBoomer = countsSmooth./repmat(max(countsSmooth,[],2), 1, size(countsSmooth,2));

% make plot
boomerWalkObsPhaseH = axes('Units','inches','OuterPosition',...
    [3.4*2 + subplotGap*2 + 1.8, 5 + subplotGap*2 + figMargins + 2.3, 2.2, 2.2]);

imagesc(walkObsPhaseBoomer)
ylabel('Neuron')
box off
set(boomerWalkObsPhaseH, 'YTick', [], 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02], ...
    'LineWidth', 2, 'YColor','none', 'XLim', [0 100], 'XColor', 'none')


% finally do starbuck
goodNeurons = intersect(find(neuralStats{2}.isUnimodal(:,1)), find(neuralStats{2}.dispersion(:,1) > dispersionCutoff));
goodNeurons = intersect(goodNeurons, find(neuralStats{2}.p_circTestR(:,1)<(0.05/42)));

[~, indsPerm] = sort(neuralStats{2}.preferedPhase(goodNeurons,1));

sortedPETH = walkPETH{2}(goodNeurons(indsPerm), :);
countsSmooth = convGauss(repmat(sortedPETH, 1, 3), 10, 20);
countsSmooth = countsSmooth(:,length(sortedPETH)+1:length(sortedPETH)*2);

walkPhaseStarbuck = countsSmooth./repmat(max(countsSmooth,[],2), 1, size(countsSmooth,2));

% make plot
starbuckWalkPhaseH = axes('Units','inches','OuterPosition',...
    [3.4*2 + subplotGap*2 - 0.2, 5 + subplotGap + figMargins, 2.2, 2.3]);

imagesc(walkPhaseStarbuck)
ylabel('Neuron')
xlabel('Gait Cycle %')
box off
set(starbuckWalkPhaseH, 'YTick', [], 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', [0 100])

% next do walking obstacle step
sortedPETH = squeeze(walkObsPETH{2}(goodNeurons(indsPerm), :, 4));
countsSmooth = convGauss(repmat(sortedPETH, 1, 3), 10, 20);
countsSmooth = countsSmooth(:,length(sortedPETH)+1:length(sortedPETH)*2);
walkObsPhaseStarbuck = countsSmooth./repmat(max(countsSmooth,[],2), 1, size(countsSmooth,2));

% make plot
starbuckWalkObsPhaseH = axes('Units','inches','OuterPosition',...
    [3.4*2 + subplotGap*2 + 1.8, 5 + subplotGap + figMargins, 2.2, 2.3]);

imagesc(walkObsPhaseStarbuck)
ylabel('Neuron')
xlabel('Gait Cycle %')
box off
set(starbuckWalkObsPhaseH, 'YTick', [], 'FontSize',11, 'YColor','none', 'XLim', [0 100],...
    'TickDir','out','TickLength', [0.02 0.02], 'LineWidth', 2)

% Labels
text(12, -31, 'Obstacle Stride', 'FontSize', 13)
text(-91, -31, 'Basic Walking', 'FontSize', 13)


%% Subplot C -- Preferred phase changes

% first do boomer
goodNeurons = intersect(find(neuralStats{1}.isUnimodal(:,1)), find(neuralStats{1}.dispersion(:,1) > dispersionCutoff));
goodNeurons = intersect(goodNeurons, find(neuralStats{1}.p_circTestR(:,1)<(0.05/50)));

% make plot
boomerPhaseShiftH = axes('Units','inches','OuterPosition',...
    [figMargins - 0.3, figMargins + 1.8, 5, 3.5]);

walkPhase = neuralStats{1}.preferedPhase(goodNeurons,1);
walkObsPhase = neuralStats{1}.preferedPhase(goodNeurons,5);

% since phase is circular, change the walk obs phase to be greater than 100
% or less than 0 if it will bring it closer to the walk phase
walkObsPhaseMod = [walkObsPhase walkObsPhase+100 walkObsPhase-100];
[~, bestMod] = min(abs(walkObsPhaseMod - walkPhase),[],2);
for iNeuron = 1:length(bestMod)
    walkObsPhase(iNeuron) = walkObsPhaseMod(iNeuron,bestMod(iNeuron));
end

% plot by significance
pValues = neuralStats{1}.preferedPhasePValue(goodNeurons);
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(walkPhase(~sigNeurs), walkObsPhase(~sigNeurs), 300, 'marker', '.');
hold on
scatter(walkPhase(sigNeurs), walkObsPhase(sigNeurs), 300, 'marker', '.', 'MarkerFaceColor', 'r');
line([0 1000], [0 1000], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
title('Preferred Phase (%)')
box off
set(boomerPhaseShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'xlim', [0 115], 'ylim', [0 115])

% next, starbuck
goodNeurons = intersect(find(neuralStats{2}.isUnimodal(:,1)), find(neuralStats{2}.dispersion(:,1) > dispersionCutoff));
goodNeurons = intersect(goodNeurons, find(neuralStats{2}.p_circTestR(:,1)<(0.05/42)));

% make plot
starbuckPhaseShiftH = axes('Units','inches','OuterPosition',...
    [figMargins - 0.3, -0.5 , 5, 3.5]);

walkPhase = neuralStats{2}.preferedPhase(goodNeurons,1);
walkObsPhase = neuralStats{2}.preferedPhase(goodNeurons,5);

% since phase is circular, change the walk obs phase to be greater than 100
% or less than 0 if it will bring it closer to the walk phase
walkPhaseMod = [walkPhase walkPhase+100];
[~, bestMod] = min(abs(walkPhaseMod - walkObsPhase),[],2);
for iNeuron = 1:length(bestMod)
    walkPhase(iNeuron) = walkPhaseMod(iNeuron,bestMod(iNeuron));
end

pValues = neuralStats{2}.preferedPhasePValue(goodNeurons);
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(walkPhase(~sigNeurs), walkObsPhase(~sigNeurs), 300, 'marker', '.');
hold on
scatter(walkPhase(sigNeurs), walkObsPhase(sigNeurs), 300, 'marker', '.', 'MarkerFaceColor', 'r');
line([0 1000], [0 1000], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
box off
set(starbuckPhaseShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'xlim', [0 115], 'ylim', [0 115])


%% Subplot D -- depth of modulation changes

% first do boomer
goodNeurons = intersect(find(neuralStats{1}.isUnimodal(:,1)), find(neuralStats{1}.dispersion(:,1) > dispersionCutoff));

% make plot
boomerDoMShiftH = axes('Units','inches','OuterPosition',...
    [figMargins+3.1, figMargins + 1.8, 5, 3.5]);

% plot by significance
pValues = neuralStats{1}.depthModsPValue;
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(neuralStats{1}.depthOfMod(~sigNeurs,1), neuralStats{1}.depthOfMod(~sigNeurs, 5), 300, 'marker', '.');
hold on
scatter(neuralStats{1}.depthOfMod(sigNeurs,1), neuralStats{1}.depthOfMod(sigNeurs, 5), 300, 'marker', '.', 'MarkerFaceColor','r');
autoXLim = get(gca, 'XLim');
autoYLim = get(gca, 'YLim');
line([0 1000], [0 1000], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
title('Depth of Modulation (spks/s)')
box off
set(boomerDoMShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', autoXLim, 'YLim', autoYLim)

% next, starbuck
goodNeurons = intersect(find(neuralStats{2}.isUnimodal(:,1)), find(neuralStats{2}.dispersion(:,1) > dispersionCutoff));

% make plot
starbuckDoMShiftH = axes('Units','inches','OuterPosition',...
    [figMargins + 3.1, -0.5 , 5, 3.5]);

pValues = neuralStats{2}.depthModsPValue;
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(neuralStats{2}.depthOfMod(~sigNeurs,1), neuralStats{2}.depthOfMod(~sigNeurs, 5), 300, 'marker', '.');
hold on
scatter(neuralStats{2}.depthOfMod(sigNeurs,1), neuralStats{2}.depthOfMod(sigNeurs, 5), 300, 'marker', '.', 'MarkerFaceColor','r');
autoXLim = get(gca, 'XLim');
autoYLim = get(gca, 'YLim');
line([0 1000], [0 1000], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
box off
set(starbuckDoMShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', autoXLim, 'YLim', autoYLim)


%% Subplot E -- Mean discharge rate changes

% first do boomer
goodNeurons = intersect(find(neuralStats{1}.isUnimodal(:,1)), find(neuralStats{1}.dispersion(:,1) > dispersionCutoff));

% make plot
boomerMeanRateShiftH = axes('Units','inches','OuterPosition',...
    [figMargins+6.5, figMargins + 1.8, 5, 3.5]);

% plot by significance
pValues = neuralStats{1}.meanFRsPValue;
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(neuralStats{1}.meanRate(~sigNeurs,1), neuralStats{1}.meanRate(~sigNeurs, 5), 300, 'marker', '.');
hold on
scatter(neuralStats{1}.meanRate(sigNeurs,1), neuralStats{1}.meanRate(sigNeurs, 5), 300, 'marker', '.','MarkerFaceColor','r');
autoXLim = get(gca, 'XLim');
autoYLim = get(gca, 'YLim');
line([0 1000], [0 1000], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
title('Mean Discharge Rate (spks/s)')
box off
set(boomerMeanRateShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', autoXLim, 'YLim', autoYLim)

% next, starbuck
goodNeurons = intersect(find(neuralStats{2}.isUnimodal(:,1)), find(neuralStats{2}.dispersion(:,1) > dispersionCutoff));

% make plot
starbuckMeanRateShiftH = axes('Units','inches','OuterPosition',...
    [figMargins + 6.5, -0.5 , 5, 3.5]);

pValues = neuralStats{2}.meanFRsPValue;
cutoffFDR = FDRcutoff(pValues,0.05);
sigNeurs = pValues < cutoffFDR;

scatter(neuralStats{2}.meanRate(~sigNeurs,1), neuralStats{2}.meanRate(~sigNeurs, 5), 300, 'marker', '.');
hold on
scatter(neuralStats{2}.meanRate(sigNeurs,1), neuralStats{2}.meanRate(sigNeurs, 5), 300, 'marker', '.','MarkerFaceColor','r');
autoXLim = get(gca, 'XLim');
autoYLim = get(gca, 'YLim');
line([0 100], [0 100], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([0 1000], [10 1010], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])
line([0 1000], [-10 990], 'linestyle', '--', 'linewidth', 1, 'color', [0.5 0.5 0.5 0.5])

xlabel('Basic Walking')
ylabel('Obstacle Stride')
box off
set(starbuckMeanRateShiftH, 'FontSize',11, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', autoXLim, 'YLim', autoYLim)


%% Subplot F -- Change timings for Boomer (taken from previous time course figure)
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

% defined from walking data
dutyPercent = 67;

subFuncsTimeCourse = makeFig_TimeCourse(false);

[kinDataWalkObs, frDataWalkObs, spikeCountsWalkObs, obsMoveWalkObs] = subFuncsTimeCourse.getData(trialsDataBoomer, 'WalkingObstacle', dutyPercent);

%do PCA on the neural data
catFRWalkObs = frDataWalkObs(:,:);

[~, catPCWalkObs, neurVarWalkObs] = pca(catFRWalkObs');

% divide back into trials and steps
tmp = mat2cell(catPCWalkObs', size(catPCWalkObs,2), repmat(100*size(frDataWalkObs,3), 1, size(frDataWalkObs,4)));
tmp = cellfun(@(x) mat2cell(x, size(catPCWalkObs,2), repmat(100,1,size(frDataWalkObs,3))), tmp, 'un', 0);
tmp = cellfun(@(x) cat(3, x{:}), tmp, 'un', 0);
pcWalkObs = cat(4, tmp{:});

nNeurDims = 15;
% get mahalanobis distance on neural PCs, using first step as reference
neurMahalRef = repmat(subFuncsTimeCourse.getMahal(squeeze(pcWalkObs(1:nNeurDims,:,1,:)), squeeze(pcWalkObs(1:nNeurDims,:,1,:)), size(pcWalkObs,4)), 1, 1, size(pcWalkObs,3));
for iStep = 1:size(pcWalkObs,3)
    neurMahalWalkObs_Boomer(:,:,iStep) = subFuncsTimeCourse.getMahal(squeeze(pcWalkObs(1:nNeurDims,:,1,:)), squeeze(pcWalkObs(1:nNeurDims,:,iStep,:)), 0);
end
neurMahalWalkObs_Boomer(:,:,1) = neurMahalRef(:,:,1);

% also get for kinematics
kinMahalRef = repmat(subFuncsTimeCourse.getMahal(squeeze(kinDataWalkObs(:,:,1,:)), squeeze(kinDataWalkObs(:,:,1,:)), size(kinDataWalkObs,4)), 1, 1, size(kinDataWalkObs,3));
for iStep = 1:size(kinDataWalkObs,3)
    kinMahalWalkObs_Boomer(:,:,iStep) = subFuncsTimeCourse.getMahal(squeeze(kinDataWalkObs(:,:,1,:)), squeeze(kinDataWalkObs(:,:,iStep,:)), 0);
end
kinMahalWalkObs_Boomer(:,:,1) = kinMahalRef(:,:,1);

% find delay
[kinNeurDiffxCorr_Boomer, kinNeurDiffLags_Boomer] = crosscorr(mean(neurMahalWalkObs_Boomer(:,:,4)),mean(kinMahalWalkObs_Boomer(:,:,4)));
[~, maxCorrInd] = max(kinNeurDiffxCorr_Boomer);
changeDelay_Boomer = kinNeurDiffLags_Boomer(maxCorrInd);

% plotting
boomerMahalsH = axes('Units','inches','OuterPosition',...
    [subplotGap-0.9, figMargins-0.1, 5.5, 3]);

yyaxis left
shadedErrorBar(1:100,mean(squeeze(neurMahalWalkObs_Boomer(:,:,4))),std(squeeze(neurMahalWalkObs_Boomer(:,:,4))),...
    {'color',[0 0.4470 0.7410],'linewidth',2,'linestyle','-'},1)
hold on
mahalLegendH(1) = plot(1:100, mean(squeeze(neurMahalWalkObs_Boomer(:,:,4))), '.-', 'color',[0 0.4470 0.7410], 'LineWidth',2);
ylabel('Mahalanobis Dist')

yyaxis right
shadedErrorBar(1:100,mean(squeeze(kinMahalWalkObs_Boomer(:,:,4))),std(squeeze(kinMahalWalkObs_Boomer(:,:,4))),...
    {'color',[0.8500, 0.3250, 0.0980],'linewidth',2,'linestyle','-'},1)
hold on
mahalLegendH(2) = plot(1:100, mean(squeeze(kinMahalWalkObs_Boomer(:,:,4))), '.-', 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',2);
xlabel('Obstacle Step Gait %')

box off
set(boomerMahalsH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'ylim', [0, 6e4])
legend(mahalLegendH, 'Neural', 'Kinematics', 'box', 'off')

boomerXCorrH = axes('Units','inches','OuterPosition',...
    [subplotGap+2.4, figMargins-0.1, 4, 3]);

plot(kinNeurDiffLags_Boomer, kinNeurDiffxCorr_Boomer, 'k', 'linewidth', 2)
hold on
line(repmat(changeDelay_Boomer,1,2), get(gca,'ylim'), 'linewidth', 2, 'color','r','linestyle','--')

box off
set(boomerXCorrH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [-20 0 changeDelay_Boomer, 20])



%% Subplot  -- Change timings for Starbuck
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

% defined from walking data
dutyPercent = 69;

[kinDataWalkObs, frDataWalkObs, obsMoveWalkObs] = subFuncsTimeCourse.getData(trialsDataStarbuck, 'WalkingObstacle', dutyPercent);

%do PCA on the neural data
catFRWalkObs = frDataWalkObs(:,:);

[~, catPCWalkObs, neurVarWalkObs] = pca(catFRWalkObs');

% divide back into trials and steps
tmp = mat2cell(catPCWalkObs', size(catPCWalkObs,2), repmat(100*size(frDataWalkObs,3), 1, size(frDataWalkObs,4)));
tmp = cellfun(@(x) mat2cell(x, size(catPCWalkObs,2), repmat(100,1,size(frDataWalkObs,3))), tmp, 'un', 0);
tmp = cellfun(@(x) cat(3, x{:}), tmp, 'un', 0);
pcWalkObs = cat(4, tmp{:});

nNeurDims = 15;
% get mahalanobis distance on neural PCs, using first step as reference
neurMahalRef = repmat(subFuncsTimeCourse.getMahal(squeeze(pcWalkObs(1:nNeurDims,:,1,:)), squeeze(pcWalkObs(1:nNeurDims,:,1,:)), size(pcWalkObs,4)), 1, 1, size(pcWalkObs,3));
for iStep = 1:size(pcWalkObs,3)
    neurMahalWalkObs_Starbuck(:,:,iStep) = subFuncsTimeCourse.getMahal(squeeze(pcWalkObs(1:nNeurDims,:,1,:)), squeeze(pcWalkObs(1:nNeurDims,:,iStep,:)), 0);
end
neurMahalWalkObs_Starbuck(:,:,1) = neurMahalRef(:,:,1);

% also get for kinematics
kinMahalRef = repmat(subFuncsTimeCourse.getMahal(squeeze(kinDataWalkObs(:,:,1,:)), squeeze(kinDataWalkObs(:,:,1,:)), size(kinDataWalkObs,4)), 1, 1, size(kinDataWalkObs,3));
for iStep = 1:size(kinDataWalkObs,3)
    kinMahalWalkObs_Starbuck(:,:,iStep) = subFuncsTimeCourse.getMahal(squeeze(kinDataWalkObs(:,:,1,:)), squeeze(kinDataWalkObs(:,:,iStep,:)), 0);
end
kinMahalWalkObs_Starbuck(:,:,1) = kinMahalRef(:,:,1);

% find delay
[kinNeurDiffxCorr_Starbuck, kinNeurDiffLags_Starbuck] = crosscorr(mean(neurMahalWalkObs_Starbuck(:,:,4)),mean(kinMahalWalkObs_Starbuck(:,:,4)));
[~, maxCorrInd] = max(kinNeurDiffxCorr_Starbuck);
changeDelay_Starbuck = kinNeurDiffLags_Starbuck(maxCorrInd);

% plot
starbuckMahalsH = axes('Units','inches','OuterPosition',...
    [subplotGap + 4.6, figMargins-0.1, 5.5, 3]);

yyaxis left
shadedErrorBar(1:100,mean(squeeze(neurMahalWalkObs_Starbuck(:,:,4))),std(squeeze(neurMahalWalkObs_Starbuck(:,:,4))),...
    {'color',[0 0.4470 0.7410],'linewidth',2,'linestyle','-'},1)

yyaxis right
shadedErrorBar(1:100,mean(squeeze(kinMahalWalkObs_Starbuck(:,:,4))),std(squeeze(kinMahalWalkObs_Starbuck(:,:,4))),...
    {'color',[0.8500, 0.3250, 0.0980],'linewidth',2,'linestyle','-'},1)
xlabel('Obstacle Step Gait %')

box off
set(starbuckMahalsH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2)

starbuckXCorrH = axes('Units','inches','OuterPosition',...
    [subplotGap+8, figMargins-0.1, 4, 3]);

plot(kinNeurDiffLags_Starbuck, kinNeurDiffxCorr_Starbuck, 'k', 'linewidth', 2)
hold on
line(repmat(changeDelay_Starbuck,1,2), get(gca,'ylim'), 'linewidth', 2, 'color','r','linestyle','--')

box off
set(starbuckXCorrH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [-20 0 changeDelay_Starbuck, 20])


%% Subfunctions for creating plots
function [raster1H, raster2H, pethH] = makeRasterAndPETH(botLeftPos, spiketimes1, spiketimes2, plotColor1, plotColor2, legendLabels,...
    stepLength1, stepLength2, neuralInd, neuralStats)
% make raster plot of spike trains and right under it plot the PSTH

% Calc axes positions. PSTH plot will be 1.5 inches high, raster plot will
% be 2 inches high, and 0.15 inch gap in between. Both will be 3 inches long.

raster1Height = 1.5;
raster2Height = 1.5;
psthHeight = 1.9;
gapHeight = 0.0;
plotLength = 3.2;

% make first raster plot
raster1H = axes('Units','inches','OuterPosition',...
    [botLeftPos(1), botLeftPos(2) + psthHeight + gapHeight + raster2Height, plotLength, raster1Height],...
    'FontSize',12, 'LineWidth',2, 'TickDir','out', 'XColor', 'none');

rasterplot(spiketimes1,'times','|',[], [], 'Color', plotColor1);
box off
ylabel('Trial')
xlim([0 100])

% make second raster plot
raster2H = axes('Units','inches','OuterPosition',...
    [botLeftPos(1), botLeftPos(2) + psthHeight + gapHeight, plotLength, raster2Height],...
    'FontSize',12, 'LineWidth',2, 'TickDir','out', 'XColor', 'none');

rasterplot(spiketimes2,'times','|', [], [], 'Color', plotColor2);
box off
ylabel('Trial')
xlim([0 100])

% make peth
pethH = axes('Units','inches','OuterPosition',...
    [botLeftPos(1), botLeftPos(2), plotLength, psthHeight],...
    'FontSize', 12, 'LineWidth', 2, 'TickDir', 'out');

% calculate bin counts
counts1 = histcounts(cat(1,spiketimes1{:}), 0:1:100)/length(spiketimes1);
counts2 = histcounts(cat(1,spiketimes2{:}), 0:1:100)/length(spiketimes2);

%convert to spks/sec based on average step length corresponding to 100%
counts1 = counts1*(100/(stepLength1/1000));
counts2 = counts2*(100/(stepLength2/1000));

% smooth a little using gaussian (repeat at ends to avoid edge artifacts)
counts1Smooth = convGauss(repmat(counts1, 1, 3), 10, 20);
counts1Smooth = counts1Smooth(length(counts1)+1:length(counts1)*2);
counts2Smooth = convGauss(repmat(counts2, 1, 3), 10, 20);
counts2Smooth = counts2Smooth(length(counts2)+1:length(counts2)*2);

pethLine(1) = plot(0:1:99,counts1Smooth, 'LineWidth', 2, 'Color', plotColor1);
hold on
pethLine(2) = plot(0:1:99,counts2Smooth, 'LineWidth', 2, 'Color', plotColor2);
ylabel('spikes/s')
box off
xlabel('Gait cycle %')

set(pethH, 'FontSize', 12, 'LineWidth', 2, 'TickDir', 'out');

% show dispersion for the two peths
title(['r_w_a_l_k = ' num2str(round(neuralStats.dispersion(neuralInd, 1)*100)/100)...
    ',   r_o_b_s_t_a_c_l_e = ' num2str(round(neuralStats.dispersion(neuralInd, 5)*100)/100)], 'FontSize', 10)

% also put line at the neuron's preferrred phase
walkPhase = neuralStats.preferedPhase(neuralInd, 1);
walkObsPhase = neuralStats.preferedPhase(neuralInd, 5);

autoYLim = get(gca, 'YLim');
line([walkPhase walkPhase], autoYLim, 'LineStyle', '--', 'linewidth', 1, 'color', plotColor1)
line([walkObsPhase walkObsPhase], autoYLim, 'LineStyle', '--', 'linewidth', 1, 'color', plotColor2)

if ~isempty(legendLabels)
    legend(pethLine, legendLabels, 'Box', 'off', 'LineWidth', 2, 'Location', 'best', 'FontSize', 10)
end



function [walkPETH, walkObsPETH, neuralStats] = getNeuralProperties(trialsData,...
    walkTrialInds, walkObsTrialInds, walkTrialDurations, walkObsTrialDurations, obstacleSteps)
% get neuron PETHs and calculate single neuron properties (unimodality,
% test for circular uniformity, dispersion, preferred phase, mean firing
% rate, and depth of modulation)

% for circular statistics, the bin centers 0.5 - 99.5 are need to be
% converted to radians
binRads = (0.5:1:99.5)*2*pi/100;
binSpacingRads = 1/100*2*pi;

% now go through each neuron and get stats
for iNeuron = 1:length(trialsData(1).SpikeTimes)
    
    %calculate PETH for basic walking first
    spikeTimes = cellfun(@(x) x{iNeuron}(x{iNeuron} > 0 & x{iNeuron} < 100),...
        {trialsData(walkTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);
    spikeCounts = histcounts(cat(1,spikeTimes{:}), 0:1:100)/length(spikeTimes);
    %convert to spks/sec based on average step length corresponding to 100%
    walkPETH(iNeuron,:) = spikeCounts*(100/mean(walkTrialDurations/1000));
    thisNeuronPETH = walkPETH(iNeuron,:);
    
    %save spike times for permutation test later
    walkSpikeTimes = spikeTimes;
    
    %if no spikes at all, set everything to nan
    if isempty(cat(1,spikeTimes{:}))
        neuralStats.isUnimodal(iNeuron,:) = 1;
        neuralStats.p_circTestR(iNeuron,:) = NaN;
        neuralStats.p_circTestO(iNeuron,:) = NaN;
        neuralStats.dispersion(iNeuron,:) = NaN;
        neuralStats.preferedPhase(iNeuron,:) = NaN;
        neuralStats.meanRate(iNeuron,:) = NaN;
        neuralStats.depthOfMod(iNeuron,:) = NaN;
        continue
    end
    
    %test bimodality
    countsSmooth = convGauss(repmat(thisNeuronPETH, 1, 3), 10, 20);
    countsSmooth = countsSmooth(length(thisNeuronPETH)+1:length(thisNeuronPETH)*2);
    neuralStats.isUnimodal(iNeuron,1) = testUnimodal(countsSmooth);
    
    %do non-uniformity tests
    neuralStats.p_circTestR(iNeuron, 1) = circ_rtest(cat(1,spikeTimes{:})/100*2*pi);
    neuralStats.p_circTestO(iNeuron, 1) = circ_otest(cat(1,spikeTimes{:})/100*2*pi);
    neuralStats.dispersion(iNeuron, 1) = circ_r(binRads,...
        spikeCounts*length(spikeTimes), binSpacingRads, 2);
    
    %get the preferred phase
    meanRad = circ_mean(binRads', thisNeuronPETH')*100/(2*pi);
    if meanRad < 0
        meanRad = meanRad + 100;
    end
    neuralStats.preferedPhase(iNeuron, 1) = meanRad;
    
    %get the mean firing rate
    neuralStats.meanRate(iNeuron, 1) = mean(thisNeuronPETH);
    
    %get depth of modulation, which is max-min firing rate
    neuralStats.depthOfMod(iNeuron, 1) = prctile(thisNeuronPETH,95) - prctile(thisNeuronPETH,5);
    
%     plotting to check
%     figure;
%     plot(thisNeuronPETH)
%     hold on
%     plot(countsSmooth, 'LineWidth', 2)
%     line(repmat(neuralStats..preferedPhase(iNeuron, 1), 1, 2), get(gca,'ylim'),'linestyle', '--', 'color','r','linewidth',2)
%     title(['p_r = ' num2str(neuralStats.p_circTestR(iNeuron, 1)) ', p_o = ' num2str(neuralStats.p_circTestO(iNeuron,1)), ...
%         ', r = ' num2str(neuralStats.dispersion(iNeuron,1)) ', Uni = ' num2str(neuralStats.isUnimodal(iNeuron, 1))])
    
    %now do obstacle walking
    if isempty(walkObsTrialInds)
        %don't want walk obs, just do walk
        walkObsPETH = [];
        continue
    end
    
    for iStep = 1:length(obstacleSteps)
        
        spikeTimes = cellfun(@(x) x{iNeuron}(x{iNeuron} > (iStep-1)*100 & x{iNeuron} < iStep*100) - (iStep-1)*100,...
            {trialsData(walkObsTrialInds).GaitNormalizedSpikeTimes}, 'un', 0);
        spikeCounts = histcounts(cat(1,spikeTimes{:}), 0:1:100)/length(spikeTimes);
        walkObsPETH(iNeuron,:,iStep) = spikeCounts*(100/mean(walkObsTrialDurations(:,iStep)/1000));
        stepPETH = squeeze(walkObsPETH(iNeuron,:,iStep));
        
        %save spike times of the obstacle step (assuming it's step 5) for permutation test
        if obstacleSteps(iStep) == 0
            walkObsSpikeTimes = spikeTimes;
        end
        
        %test bimodality
        countsSmooth = convGauss(repmat(stepPETH, 1, 3), 10, 20);
        countsSmooth = countsSmooth(length(stepPETH)+1:length(stepPETH)*2);
        neuralStats.isUnimodal(iNeuron, iStep+1) = testUnimodal(countsSmooth);
        
        %do non-uniformity tests
        neuralStats.p_circTestR(iNeuron, iStep+1) = circ_rtest(cat(1,spikeTimes{:})/100*2*pi);
        neuralStats.p_circTestO(iNeuron, iStep+1) = circ_otest(cat(1,spikeTimes{:})/100*2*pi);
        neuralStats.dispersion(iNeuron, iStep+1) = circ_r(binRads,...
            spikeCounts*length(spikeTimes), binSpacingRads, 2);
        
        %get the preferred phase
        meanRad = circ_mean(binRads', stepPETH')*100/(2*pi);
        if meanRad < 0
            meanRad = meanRad + 100;
        end
        neuralStats.preferedPhase(iNeuron, iStep+1) = meanRad;
        
        %get the mean firing rate
        neuralStats.meanRate(iNeuron, iStep+1) = mean(stepPETH);
        
        %get depth of modulation, which is max-min firing rate
        neuralStats.depthOfMod(iNeuron, iStep+1) = prctile(stepPETH,95) - prctile(stepPETH,5);
        
    end
    
    %do permutation test to test for significant differences between
    %obstacle and basic strides
    nPerms = 1000;
    for iPerm = 1:nPerms
        
        nWalkTrials = length(walkSpikeTimes);
        nWalkObsTrials = length(walkObsSpikeTimes);
        
        permInds = randperm(nWalkTrials+nWalkObsTrials);
        allSpikeTimes = [walkSpikeTimes walkObsSpikeTimes];
        
        permWalkSpikeTimes = allSpikeTimes(permInds(1:nWalkTrials));
        permWalkObsSpikeTimes = allSpikeTimes(permInds(nWalkTrials+1:end));
        
        %calculate PETHs
        meanBinSize = (mean(walkTrialDurations)+mean(walkObsTrialDurations(:,obstacleSteps==0)))/2;
        spikeCounts = histcounts(cat(1,permWalkSpikeTimes{:}), 0:1:100)/length(permWalkSpikeTimes);
        walkPermPETH = spikeCounts*(100/(meanBinSize/1000));
        
        %and for obstacle step
        spikeCounts = histcounts(cat(1,permWalkObsSpikeTimes{:}), 0:1:100)/length(permWalkObsSpikeTimes);
        walkObsPermPETH = spikeCounts*(100/(meanBinSize/1000));
        
        permPETHs = {walkPermPETH walkObsPermPETH};
        
        for iCond = 1:2
            %find the difference in:
            %preferred phase
            meanRad = circ_mean(binRads', permPETHs{iCond}')*100/(2*pi);
            if meanRad < 0
                meanRad = meanRad + 100;
            end
            preferedPhases(iCond) = meanRad;
            
            %mean FR
            meanFRs(iCond) = mean(permPETHs{iCond});
            
            %depth of modulation
            depthMods(iCond) = prctile(permPETHs{iCond},95) - prctile(permPETHs{iCond},5);
        end
        
        %add to the distributions
        meanFRsDiffPerm(iPerm) = diff(meanFRs);
        depthModsDiffPerm(iPerm) = diff(depthMods);
        
        
        % since phase is circular, change the walk obs phase to be greater than 100
        % or less than 0 if it will bring it closer to the walk phase
        walkObsPhaseMod = [preferedPhases(2) preferedPhases(2)+100 preferedPhases(2)-100];
        [~, bestMod] = min(abs(walkObsPhaseMod - preferedPhases(1)),[],2);
        preferedPhases(2) = walkObsPhaseMod(bestMod);
        preferedPhaseDiffPerm(iPerm) = diff(preferedPhases);
        
    end
    
    %finally, calculate the p-values using the null distribution from above
    meanFRsDiff = neuralStats.meanRate(iNeuron, find(obstacleSteps==0)+1) - neuralStats.meanRate(iNeuron,1);
    neuralStats.meanFRsPValue(iNeuron) = (sum(abs(meanFRsDiffPerm) > abs(meanFRsDiff))+1)/(nPerms+1);
    
    depthModsDiff = neuralStats.depthOfMod(iNeuron, find(obstacleSteps==0)+1) - neuralStats.depthOfMod(iNeuron,1);
    neuralStats.depthModsPValue(iNeuron) = (sum(abs(depthModsDiffPerm) > abs(depthModsDiff))+1)/(nPerms+1);
    
    walkObsPhase = neuralStats.preferedPhase(iNeuron, find(obstacleSteps==0)+1);
    walkPhase = neuralStats.preferedPhase(iNeuron,1);
    walkObsPhaseMod = [walkObsPhase walkObsPhase+100 walkObsPhase-100];
    [~, bestMod] = min(abs(walkObsPhaseMod - walkPhase),[],2);
    walkObsPhase = walkObsPhaseMod(bestMod);
    preferredPhaseDiff = walkObsPhase - walkPhase;
    neuralStats.preferedPhasePValue(iNeuron) = (sum(abs(preferedPhaseDiffPerm) > abs(preferredPhaseDiff))+1)/(nPerms+1);
    
end



function isUnimodal = testUnimodal(smoothedPETH)
% function to determine whether the PETH of neurons can be considered
% unimodal or not

% get the max and min values
[maxVal maxInd] = max(smoothedPETH);
[minVal, minInd] = min(smoothedPETH);

if maxVal == minVal
    isUnimodal = true;
    return
end

% segment timecourse into rising segment (min to max) and dropping segment
% (max down to min)
tmp = circshift(smoothedPETH, 1-maxInd);
[~, minShifted] = min(tmp);
dropSeg = tmp(1:minShifted);

tmp = circshift(smoothedPETH, 1-minInd);
[~, maxShifted] = max(tmp);
riseSeg = tmp(1:maxShifted);

% find any peaks in the dropping seg and any troughs in the rising seg
[~, pkInds, pkWidth, pkProm] = findpeaks(dropSeg);
[~, trInds, trWidth, trProm] = findpeaks(-riseSeg);

% any peaks that are greater than 50% of the main amplitude and more than
% 10% of the gait cycle will be counted as another mode
mainAmp = maxVal-minVal;
ampThresh = mainAmp * 0.5;
widthThresh = 10;

if any([pkProm > ampThresh & pkWidth > widthThresh, ...
        trProm > ampThresh & trWidth > widthThresh])
    isUnimodal = false;
else
    isUnimodal = true;
end
    


% 
