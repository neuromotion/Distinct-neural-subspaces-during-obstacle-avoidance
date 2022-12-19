function subfuncs = makefig_NeuralDataSupp(varargin)

% Code to make neural data processing supplementary figure
% 
% Layout:
% A--------------------------  B--------------------------
% |  Boomer implant location   | Starbuck implant location
% |                                                          
% |  
% 
% C--------------------------   D-------------------------
% |  Example neural             |  Example kinematics
% |  normalization              |  normalization
% |
% 
% E--------------------------   F-------------------------
% |  Boomer neuron r-metric     |  Starbuck neuron r-metric
% |  distribution               |  distribution
% |                             |

% expose subfunctions in case any other figures want to use them
subfuncs.makeRasterAndPETH = @makeRasterAndPETH;
subfuncs.testUnimodal = @testUnimodal;
subfuncs.getNeuralProperties = @getNeuralProperties;

if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 11.5, 11.2])


%% Subplot A -- Example raster normalization
% use example from Boomer
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

% set example neuron
exampleNeuron = 4;

% only use basic walking trials
badTrials = filterTrials(trialsDataBoomer,90,5);
walkTrialIndsBoomer = find(strcmpi('Walk',string({trialsDataBoomer.Task})));
walkTrialIndsBoomer = setdiff(walkTrialIndsBoomer, badTrials);

% get normalized spike data
spikeTimesNorm = cellfun(@(x) x{exampleNeuron}, {trialsDataBoomer(walkTrialIndsBoomer).GaitNormalizedSpikeTimes}, 'un', 0);
% don't use pre and post-trial
spikeTimesNorm = cellfun(@(x) x(x>0 & x<100), spikeTimesNorm, 'un', 0);

% get unnormalized spike data
spikeTimesUnnorm = cellfun(@(x) x{exampleNeuron}, {trialsDataBoomer(walkTrialIndsBoomer).SpikeTimes}, 'un', 0);
% don't use pre and post-trial
spikeTimesStop = {trialsDataBoomer(walkTrialIndsBoomer).SpikeTimesStop};
spikeTimesUnnorm = cellfun(@(x,y) x(x>0 & x<y), spikeTimesUnnorm, spikeTimesStop, 'un', 0);

% get toe off event times
toeOffTimes = cellfun(@(x) x(2), {trialsDataBoomer(walkTrialIndsBoomer).TrialEvents})*10;

% sort by ascending spike times stop
spikeTimesStop = [spikeTimesStop{:}];
[spikeTimesStopOrdered, orderInds] = sort(spikeTimesStop);

spikeTimesUnnormOrdered = spikeTimesUnnorm(orderInds);
toeOffTimesOrdered = toeOffTimes(orderInds);
spikeTimesNormOrdered = spikeTimesNorm(orderInds);

% make raster plot
unnormalizedRasterH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [-0.8, 5.7, 6.4, 3.2]);

rasterplot(spikeTimesUnnormOrdered,'Times','|')

colorMap = lines(2);
% add lines for toe off and step stop
for iTrial = 1:length(toeOffTimesOrdered)
    hold on;
    line(repmat(toeOffTimesOrdered(iTrial),1, 2), [iTrial-0.5, iTrial+0.5], 'color', colorMap(1,:), 'linewidth', 3)
    line(repmat(spikeTimesStopOrdered(iTrial),1, 2), [iTrial-0.5, iTrial+0.5], 'color', colorMap(2,:), 'linewidth', 3)
end

set(unnormalizedRasterH, 'FontSize',12, 'LineWidth',2, 'TickDir','out', 'YColor', 'none', 'XLim', [0 1220], 'XAxisLocation', 'top');
xlabel('Time (ms)')
box off

% now make raster plot for normalized
% make raster plot
normalizedRasterH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [-0.8, 4.1, 6.4, 3.2]);

rasterplot(spikeTimesNormOrdered,'Times','|')
hold on
line([67 67], [0 51], 'color', colorMap(1,:), 'linewidth', 3);
line([100 100], [0 51], 'color', colorMap(2,:), 'linewidth', 3);

set(normalizedRasterH, 'FontSize',12, 'LineWidth',2, 'TickDir','out', 'YColor', 'none', 'XLim', [0 105]);
xlabel('Gait Cycle Percentage')
box off


%% Subplot B -- Example kinematics normalization

% only plot 5 trials
nTrialsPlot = 5;

% get normalized kinematics (use toe height)
toeHeightNormalized = cellfun(@(x) x.PinkyToe(2,:), {trialsDataBoomer(walkTrialIndsBoomer).GaitNormalizedKinematics}, 'un', 0);
toeHeightNormalizedOrdered = toeHeightNormalized(orderInds(1:nTrialsPlot));

% get unnormalized kinematics
toeHeightUnnormalized = cellfun(@(x) x.PinkyToe(2,:), {trialsDataBoomer(walkTrialIndsBoomer).Kinematics}, 'un', 0);
toeHeightUnormalizedOrdered = toeHeightUnnormalized(orderInds(1:nTrialsPlot));

% plot unnormalized kinematics
unnormalizedKinH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [3.9, 5.5, 6.4, 3.4]);

colorMap = lines(nTrialsPlot+2);
for iTrial = 1:length(toeHeightUnormalizedOrdered)
    hold on
    plot(toeHeightUnormalizedOrdered{iTrial}+0.02*(iTrial-1), 'color', colorMap(iTrial+2,:))
    plot(toeOffTimesOrdered(iTrial)/10, toeHeightUnormalizedOrdered{iTrial}(toeOffTimesOrdered(iTrial)/10)+0.02*(iTrial-1),...
        '.','MarkerSize', 15, 'color', colorMap(1,:))
    plot(length(toeHeightUnormalizedOrdered{iTrial}), toeHeightUnormalizedOrdered{iTrial}(end)+0.02*(iTrial-1),...
        '.','MarkerSize', 15, 'color', colorMap(2,:))
end

set(unnormalizedKinH, 'FontSize',12, 'LineWidth',2, 'TickDir','out', 'YColor', 'none', 'XAxisLocation', 'top', ...
    'xLim', [0 110], 'XTickLabel', string(num2cell(0:200:1000)));
xlabel('Time (ms)')

% plot normalized kinematics
normalizedKinH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [3.9, 4.1, 6.4, 3.4]);

for iTrial = 1:length(toeHeightUnormalizedOrdered)
    hold on
    plot(toeHeightNormalizedOrdered{iTrial}+0.02*(iTrial-1), 'color', colorMap(iTrial+2,:))
    plot(67, toeHeightNormalizedOrdered{iTrial}(67)+0.02*(iTrial-1), '.', 'MarkerSize', 15,'color', colorMap(1,:))
    plot(101, toeHeightNormalizedOrdered{iTrial}(101)+0.02*(iTrial-1), '.','MarkerSize', 15, 'color', colorMap(2,:))
end

set(normalizedKinH, 'FontSize',12, 'LineWidth',2, 'TickDir','out', 'YColor', 'none', 'XLim', [1 105]);
xlabel('Gait Cycle Percentage')


%% Subplot C -- Raster to PETH to circular dispersion

% get PETH
walkTrialDurations_Boomer = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkTrialIndsBoomer).SpikeCounts})*10;
counts = histcounts(cat(1,spikeTimesNorm{:}), 0:1:100)/length(spikeTimesNorm) * ...
    (100/(mean(walkTrialDurations_Boomer)/1000));

countsSmooth = convGauss(repmat(counts, 1, 3), 10, 20);
countsSmooth = countsSmooth(length(counts)+1:length(counts)*2);

% get circular statistics
% convert from percentage to radians
binRads = (0.5:1:99.5)*2*pi/100;
binSpacingRads = 1/100*2*pi;
meanRad = circ_mean(binRads', countsSmooth');
dispersion = circ_r(binRads, counts*length(spikeTimesNorm), binSpacingRads, 2);

% plot PETH
pethH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [-0.7, 1.6, 6.4 3.6]);

plot(1:100, countsSmooth, 'linewidth', 3)
hold on
line([0 100], repmat(mean(countsSmooth),1,2), 'linestyle', '--', 'linewidth', 2, 'color', 'k')
set(pethH, 'FontSize',12, 'LineWidth',2, 'TickDir','out', 'xlim', [0 105]);
box off
xlabel('Gait Cycle Percentage')
ylabel('Firing Rate (spks/s)')


% next plot in polar coordinates
polarH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [4, 1.2, 6.4 3.8]);

polarPlotH = polarplot((1:100)/100*2*pi, countsSmooth, 'LineWidth', 3);
hold on
polarplot(repmat(meanRad,1,2), [0, dispersion*max(countsSmooth)], 'LineWidth', 3)
text(100,100,{['r = ' num2str(dispersion) ' ,'], ['ang = ' num2str(meanRad/pi*180)]})

set(polarPlotH.Parent, 'FontSize',12, 'LineWidth',2, 'RTickLabel', {});


%% Subplot D -- Neural dispersion distribution for Boomer

% get neural stats
neuralPropFuncs = makeFig_NeuralProps(false);
[~, ~, neuralStatsBoomer] = neuralPropFuncs.getNeuralProperties(trialsDataBoomer,...
    walkTrialIndsBoomer, [], walkTrialDurations_Boomer);

boomerDispersionH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [-0.7, -0.6, 6.4, 3.4]);
    
histogram(neuralStatsBoomer.dispersion,0:0.075:1, 'EdgeAlpha', 0)
hold on
line([0.15 0.15], get(gca,'YLim'), 'linestyle', '--', 'color', 'r', 'linewidth', 2)

set(boomerDispersionH, 'FontSize',12, 'LineWidth',2, 'TickDir','out');
box off
xlabel('Dispersion (r)')
ylabel('Count')



%% Subplot E -- Neural dispersion distribution for Starbuck

load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

% only use basic walking trials
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkTrialIndsStarbuck = find(strcmpi('Walk',string({trialsDataStarbuck.Task})));
walkTrialIndsStarbuck = setdiff(walkTrialIndsStarbuck, badTrials);

walkTrialDurations_Starbuck = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkTrialIndsStarbuck).SpikeCounts})*10;

[~, ~, neuralStatsStarbuck] = neuralPropFuncs.getNeuralProperties(trialsDataStarbuck,...
    walkTrialIndsStarbuck, [], walkTrialDurations_Starbuck);

% make histogram
starbuckDispersionH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [3.9, -0.6, 6.4, 3.4]);

histogram(neuralStatsStarbuck.dispersion,0:0.075:1, 'EdgeAlpha', 0)
hold on
line([0.15 0.15], get(gca,'YLim'), 'linestyle', '--', 'color', 'r', 'linewidth', 2)

set(starbuckDispersionH, 'FontSize',12, 'LineWidth',2, 'TickDir','out');
box off
xlabel('Dispersion (r)')
ylabel('Count')



% 
