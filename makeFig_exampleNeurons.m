function makeFig_exampleNeurons

% Code to make rasters and PETHs of various example neurons
% 
% Layout:
% A---------------------------  B-------------------------  C--------------------------
% |  Boomer example neurons     |                           |
% |                             |                           |
% |    No Change(30,9,16,17)    |       Complex (19,12)        |      Large Change (32)                     
% |  Raster (basic walk)        | Raster (basic walk)       |  Raster (basic walk)      
% |  Raster (walking obs step)  | Raster (walking obs step) |  Raster (walking obs step)   
% |  PSTH (walking obs step)    | PSTH (walking obs step)   |  PSTH (walking obs step)
% 
% D---------------------------  E-------------------------  F--------------------------
% |  Starbuck example neurons   |                           |
% |                             |                           |
% |       Bimodal (11)          |  Weakly modulated (31)    |      Weird Change (26)                     
% |  Raster (basic walk)        | Raster (basic walk)       |  Raster (basic walk)      
% |  Raster (walking obs step)  | Raster (walking obs step) |  Raster (walking obs step)   
% |  PSTH (walking obs step)    | PSTH (walking obs step)   |  PSTH (walking obs step)

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 11.3, 11.5])

figMargins = 0.05;
subplotGap = 0.1;

load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

% use the subfunction from makeFig_NeuralProps to get neural properties
% get the trial indices
subFuncs = makeFig_NeuralProps(1);

badTrials = filterTrials(trialsDataBoomer,90,5);
walkTrialInds{1} = find(strcmpi('Walk',string({trialsDataBoomer.Task})));
walkTrialInds{1} = setdiff(walkTrialInds{1}, badTrials);

walkObsTrialInds{1} = find(strcmpi('WalkingObstacle',string({trialsDataBoomer.Task})));
walkObsTrialInds{1} = setdiff(walkObsTrialInds{1}, badTrials);

steps_Boomer = -3:3;
walkTrialDurations_Boomer = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkTrialInds{1}).SpikeCounts})*10;
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkObsTrialInds{1}).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataBoomer(walkObsTrialInds{1}).TrialEvents) walkObsTotalDurations];
walkObsTrialDurations_Boomer = diff(walkObsEvents(:,1:2:end),1,2)*10;

[walkPETH{1}, walkObsPETH{1}, neuralStats{1}] = subFuncs.getNeuralProperties(trialsDataBoomer,...
    walkTrialInds{1}, walkObsTrialInds{1}, walkTrialDurations_Boomer, walkObsTrialDurations_Boomer, steps_Boomer);

% next Starbuck
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkTrialInds{2} = find(strcmpi('Walk',string({trialsDataStarbuck.Task})));
walkTrialInds{2} = setdiff(walkTrialInds{2}, badTrials);

walkObsTrialInds{2} = find(strcmpi('WalkingObstacle',string({trialsDataStarbuck.Task})));
walkObsTrialInds{2} = setdiff(walkObsTrialInds{2}, badTrials);

steps_Starbuck = -3:2;
walkTrialDurations_Starbuck = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkTrialInds{2}).SpikeCounts})*10;
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkObsTrialInds{2}).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataStarbuck(walkObsTrialInds{2}).TrialEvents) walkObsTotalDurations];
walkObsTrialDurations_Starbuck = diff(walkObsEvents(:,1:2:end),1,2)*10;

[walkPETH{2}, walkObsPETH{2}, neuralStats{2}] = subFuncs.getNeuralProperties(trialsDataStarbuck,...
    walkTrialInds{2}, walkObsTrialInds{2}, walkTrialDurations_Starbuck, walkObsTrialDurations_Starbuck, steps_Starbuck);



%% Subplot A -- Boomer example neuron - little difference between obstacle step

exampleNeuron = 30;

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataBoomer(walkTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataBoomer(walkObsTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
colorMap = parula(8);
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins, 5.2+subplotGap+figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), {'Basic Walking', 'Walking Obstacle'}, mean(walkTrialDurations_Boomer),...
    mean(walkObsTrialDurations_Boomer(:,4)), exampleNeuron, neuralStats{1});


%% Subplot B -- Boomer example neuron - complex waveform

exampleNeuron = 19; 

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataBoomer(walkTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataBoomer(walkObsTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins+3.75, 5.2+subplotGap+figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Boomer),...
    mean(walkObsTrialDurations_Boomer(:,4)), exampleNeuron, neuralStats{1});


%% Subplot C -- Boomer example neuron - large difference between obstacle step

exampleNeuron = 32;

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataBoomer(walkTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataBoomer(walkObsTrialInds{1}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins+7.5, 5.2+subplotGap+figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Boomer),...
    mean(walkObsTrialDurations_Boomer(:,4)), exampleNeuron, neuralStats{1});

%% Subplot D -- Starbuck example neuron - bimodal

exampleNeuron = 11;

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataStarbuck(walkTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataStarbuck(walkObsTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins, figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Starbuck),...
    mean(walkObsTrialDurations_Starbuck(:,4)), exampleNeuron, neuralStats{2});


%% Subplot E -- Starbuck example neuron - weakly modulated neuron

exampleNeuron = 31;

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataStarbuck(walkTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataStarbuck(walkObsTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins+3.75, figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Starbuck),...
    mean(walkObsTrialDurations_Starbuck(:,4)), exampleNeuron, neuralStats{2});


%% Subplot F -- Starbuck example neuron - complex difference in obstacle step

exampleNeuron = 13;

% get the spike times, only during the gait cycle
walkSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 0 & x{exampleNeuron} < 100),...
    {trialsDataStarbuck(walkTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);
walkObsSpikeTimes = cellfun(@(x) x{exampleNeuron}(x{exampleNeuron} > 300 & x{exampleNeuron} < 400)-300,...
    {trialsDataStarbuck(walkObsTrialInds{2}).GaitNormalizedSpikeTimes}, 'un', 0);

% do the actual plotting
[rast1H, rast2H, pethH] = subFuncs.makeRasterAndPETH([figMargins+7.5, figMargins], walkSpikeTimes, walkObsSpikeTimes,...
    colorMap(1,:), colorMap(5,:), [], mean(walkTrialDurations_Starbuck),...
    mean(walkObsTrialDurations_Starbuck(:,4)), exampleNeuron, neuralStats{2});


% 
