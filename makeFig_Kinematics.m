function subfuncs = makeFig_Kinematics(varargin)

% Code to make figure showing kinematics of obstacle stepping for paper
% 
% Layout:
% A---------------------     B---------------------    B-------------------   C--------------------   
% |  Boomer                  |                          |  Step height         |  Step duration
% |  Stick diagram basic     | Stick diagram obstacle   |  Boomer              |  Boomer        
% |  walking                 | step                     |                      |                       
% |                          |                          |                      |
% |  Starbuck                |                          |  Starbuck            |  Starbuck
% |  Stick diagram basic     | Stick diagram obstacle   |                      |
% |  walking                 | step                     |                      |
% 

% expose subfunctions in case any other figures want to use them
subfuncs.makeStickFig = @makeStickFig;

% if something only wants to get the subfunctions, don't actually run the
% script to make the figure
if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 15.5, 8.6])

figMargins = 0.05;
subplotGap = 0.1;

jointNames = string({'Crest','Hip','Knee','Ankle','FootKnuckle','PinkyToe'});
legendNames = string({'Crest', 'Hip', 'Knee', 'Ankle', '5th MT Head', '5th Digit Tip'});

%% Subplot A -- Stick diagram for a basic walking step

% load in data
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

% get basic walking trials
badTrials = filterTrials(trialsDataBoomer,90,5);
walkTrialInds = find(strcmpi('Walk',string({trialsDataBoomer.Task})));
walkTrialInds = setdiff(walkTrialInds, badTrials);

% example trial to plot stick diagram:
boomerExampleTrial = 1;

% make stick figure
plotColors = parula(6);
timeSpacing = 6;

boomerWalkStickH = axes('Units','inches','OuterPosition',...
    [-1.7, 2, 7, 7]);
makeStickFig(trialsDataBoomer(walkTrialInds(boomerExampleTrial)).Kinematics, jointNames, plotColors, timeSpacing, true, legendNames)

% also plot the toe trajectory
toeX = trialsDataBoomer(walkTrialInds(boomerExampleTrial)).Kinematics.PinkyToe(1,:);
toeY = trialsDataBoomer(walkTrialInds(boomerExampleTrial)).Kinematics.PinkyToe(2,:);
hold on
plot(toeX,toeY, 'linewidth', 3, 'color', [0.8500, 0.3250, 0.0980])
box off
axis off

% now do Starbuck
starbuckExampleTrial = 1;
boomerWalkStickH = axes('Units','inches','OuterPosition',...
    [-1.7, -1.5, 7, 7]);
makeStickFig(trialsDataStarbuck(walkTrialInds(starbuckExampleTrial)).Kinematics, jointNames, plotColors, timeSpacing, false)
box off
axis off

toeX = trialsDataStarbuck(walkTrialInds(boomerExampleTrial)).Kinematics.PinkyToe(1,:);
toeY = trialsDataStarbuck(walkTrialInds(boomerExampleTrial)).Kinematics.PinkyToe(2,:);
hold on
plot(toeX,toeY, 'linewidth', 3, 'color', [0.8500, 0.3250, 0.0980])
box off
axis off


%% Subplot B -- Stick diagram for a Obstacle step

% first Boomer
walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataBoomer.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

% example trial to plot stick diagram:
boomerExampleTrial = 10;

% make stick figure
timeSpacing = 6;

boomerWalkObsStickH = axes('Units','inches','OuterPosition',...
    [2, 2, 7, 7]);

% get the data
kinData = trialsDataBoomer(walkObsTrialInds(boomerExampleTrial)).Kinematics;

% get the obstacle step
obsStepStartInd = trialsDataBoomer(walkObsTrialInds(boomerExampleTrial)).TrialEvents(...
    strcmpi(trialsDataBoomer(walkObsTrialInds(boomerExampleTrial)).TrialEventsLabels, 'Obstacle Step Limb Strike'));

obsStepStopInd = trialsDataBoomer(walkObsTrialInds(boomerExampleTrial)).TrialEvents(...
    strcmpi(trialsDataBoomer(walkObsTrialInds(boomerExampleTrial)).TrialEventsLabels, '1 Obstacle Step Limb Strike'));

for iJoint = 1:length(jointNames)
    
    kinData.(jointNames{iJoint}) = kinData.(jointNames{iJoint})(:, obsStepStartInd:obsStepStopInd);
    
end

makeStickFig(kinData, jointNames, plotColors, timeSpacing, false)
box off
axis off

% also plot the toe trajectory
toeX = kinData.PinkyToe(1,:);
toeY = kinData.PinkyToe(2,:);
hold on
plot(toeX,toeY, 'linewidth', 3, 'color', [0.8500, 0.3250, 0.0980])
box off
axis off

% then starbuck
walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataStarbuck.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

% example trial to plot stick diagram:
starbuckExampleTrial = 10;

% get the data
kinData = trialsDataStarbuck(walkObsTrialInds(starbuckExampleTrial)).Kinematics;

% get the obstacle step
obsStepStartInd = trialsDataStarbuck(walkObsTrialInds(starbuckExampleTrial)).TrialEvents(...
    strcmpi(trialsDataStarbuck(walkObsTrialInds(starbuckExampleTrial)).TrialEventsLabels, 'Obstacle Step Limb Strike'));

obsStepStopInd = trialsDataStarbuck(walkObsTrialInds(starbuckExampleTrial)).TrialEvents(...
    strcmpi(trialsDataStarbuck(walkObsTrialInds(starbuckExampleTrial)).TrialEventsLabels, '1 Obstacle Step Limb Strike'));

for iJoint = 1:length(jointNames)
    
    kinData.(jointNames{iJoint}) = kinData.(jointNames{iJoint})(:, obsStepStartInd:obsStepStopInd);
    
end

% make stick figure
timeSpacing = 6;

starbuckWalkStickH = axes('Units','inches','OuterPosition',...
    [2, -1.5, 7, 7]);
makeStickFig(kinData, jointNames, plotColors, timeSpacing, false)
box off
axis off
toeX = kinData.PinkyToe(1,:);
toeY = kinData.PinkyToe(2,:);
hold on
plot(toeX,toeY, 'linewidth', 3, 'color', [0.8500, 0.3250, 0.0980])


%% Subplot C -- Step height across steps

% get heights for Boomer
walkTrialInds = find(strcmpi('Walk',string({trialsDataBoomer.Task})));
walkTrialInds = setdiff(walkTrialInds, badTrials);
walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataBoomer.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
nSteps = length(trialsDataBoomer(walkObsTrialInds(1)).TrialEvents)/2 - 1;

[walkDurationBoomer, walkHeightBoomer, walkObsDurationBoomer, walkObsHeightBoomer] = ...
    getStepProps(trialsDataBoomer, walkTrialInds, walkObsTrialInds, nSteps);

% make bar plot
colorMap = parula(8);
boomerStepHeightsH = axes('Units','inches','OuterPosition',...
    [6.3, 3.8, 6.3, 4]);

barHeights = [mean(walkHeightBoomer) mean(walkObsHeightBoomer)];
barSEMs = [std(walkHeightBoomer) std(walkObsHeightBoomer)];
boomerHeightBarH = bar([1 3:3+nSteps-1], barHeights, 'EdgeAlpha', 0, 'FaceColor', 'flat');

% colors
for iBar = 1:size(boomerHeightBarH.CData, 1)
    boomerHeightBarH.CData(iBar,:) = colorMap(iBar,:);
end

% plot error bars (std)
hold on
errorbar([1 3:3+nSteps-1], barHeights, barSEMs, '.k', 'MarkerSize', 1, 'LineWidth', 2)

% and line to indicate obstacle height
line(get(gca, 'xlim'), [7.62 7.62], 'lineStyle', '--', 'linewidth', 1, 'color', 'r')

% make pretty
box off
set(boomerStepHeightsH, 'FontSize',13, 'TickDir','out','TickLength', [0.02 0.02], 'ylim', [0 21], ...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'XTickLabel', {'Basic Walk', 'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2', 'Obs Stride + 3'})
ylabel('Stride height (cm)')


% now do the same thing for Starbuck
walkTrialInds = find(strcmpi('Walk',string({trialsDataStarbuck.Task})));
walkTrialInds = setdiff(walkTrialInds, badTrials);
walkObsTrialInds = find(strcmpi('WalkingObstacle',string({trialsDataStarbuck.Task})));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
nSteps = length(trialsDataStarbuck(walkObsTrialInds(1)).TrialEvents)/2;

[walkDurationStarbuck, walkHeightStarbuck, walkObsDurationStarbuck, walkObsHeightStarbuck] = ...
    getStepProps(trialsDataStarbuck, walkTrialInds, walkObsTrialInds, nSteps);

barHeights = [mean(walkHeightStarbuck) mean(walkObsHeightStarbuck)];
barSEMs = [std(walkHeightStarbuck) std(walkObsHeightStarbuck)];

% make plot
starbuckStepHeightsH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [6.3, 0.2, 6.3, 4]);
starbuckHeightBarH = bar([1 3:3+nSteps-1], barHeights, 'EdgeAlpha', 0, 'FaceColor', 'flat');

% colors
for iBar = 1:size(starbuckHeightBarH.CData, 1)
    starbuckHeightBarH.CData(iBar,:) = colorMap(iBar,:);
end

% plot error bars (std)
hold on
errorbar([1 3:3+nSteps-1], barHeights, barSEMs, '.k', 'MarkerSize', 1, 'LineWidth', 2)

% and line to indicate obstacle height
line(get(gca, 'xlim'), [7.62 7.62], 'lineStyle', '--', 'linewidth', 1, 'color', 'r')

% make pretty
box off
set(starbuckStepHeightsH, 'FontSize',13, 'TickDir','out','TickLength', [0.02 0.02], 'ylim', [0 12.5], 'ytick', 0:2:12, ...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'XTickLabel', {'Basic Walk', 'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2'})
ylabel('Stride height (cm)')


%% Subplot D -- Step durations across steps

% first Boomer
boomerStepDurationH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [10.3, 3.8, 6.3, 4]);

nSteps = length(trialsDataBoomer(walkObsTrialInds(1)).TrialEvents)/2 - 1;

barHeights = [mean(walkDurationBoomer) mean(walkObsDurationBoomer)];
barSEMs = [std(walkDurationBoomer) std(walkObsDurationBoomer)];
boomerDurationBarH = bar([1 3:3+nSteps-1], barHeights, 'EdgeAlpha', 0, 'FaceColor', 'flat');

% colors
for iBar = 1:size(boomerDurationBarH.CData, 1)
    boomerDurationBarH.CData(iBar,:) = colorMap(iBar,:);
end

% plot error bars (std)
hold on
errorbar([1 3:3+nSteps-1], barHeights, barSEMs, '.k', 'MarkerSize', 1, 'LineWidth', 2)

% make pretty
box off
set(boomerStepDurationH, 'FontSize',13, 'TickDir','out','TickLength', [0.02 0.02], 'ylim', [0 1300], 'ytick', 0:200:1200, ...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'XTickLabel', {'Basic Walk', 'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2', 'Obs Stride + 3'})
ylabel('Stride duration (ms)')


% next starbuck
starbuckStepDurationH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [10.3, 0.2, 6.3, 4]);

nSteps = length(trialsDataStarbuck(walkObsTrialInds(1)).TrialEvents)/2;

barHeights = [mean(walkDurationStarbuck) mean(walkObsDurationStarbuck)];
barSEMs = [std(walkDurationStarbuck) std(walkObsDurationStarbuck)];
starbuckDurationBarH = bar([1 3:3+nSteps-1], barHeights, 'EdgeAlpha', 0, 'FaceColor', 'flat');

% colors
for iBar = 1:size(starbuckDurationBarH.CData, 1)
    starbuckDurationBarH.CData(iBar,:) = colorMap(iBar,:);
end

% plot error bars (std)
hold on
errorbar([1 3:3+nSteps-1], barHeights, barSEMs, '.k', 'MarkerSize', 1, 'LineWidth', 2)

% make pretty
box off
set(starbuckStepDurationH, 'FontSize',13, 'TickDir','out','TickLength', [0.02 0.02], 'ylim', [0 1300], 'ytick', 0:200:1200, ...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'XTickLabel', {'Basic Walk', 'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2', 'Obs Stride + 3'})
ylabel('Stride duration (ms)')




function makeStickFig(kinematics, jointNames, plotColors, timeSpacing, makeLegend, legendNames, timesToPlot)
% Plot stick diagram of a limb

nTimePoints = size(kinematics.(jointNames{1}),2);

% always plot the time of maximum toe height
[~, maxToeTime] = max(kinematics.PinkyToe(2,:));
if ~isempty(timeSpacing)
    timesToPlot = sort(unique([maxToeTime:-timeSpacing:1 maxToeTime:timeSpacing:nTimePoints]));
end

% set max toe height stick figure to non-transparent
maxToeInd = find(timesToPlot == maxToeTime);

hold on;
for iTime = 1:length(timesToPlot)
    
    if iTime == maxToeInd
        trans = 1;
    else
        trans = 0.2;
    end
    
    for iJoint = 1:length(jointNames)-1
        
        %plot limb segments as lines
        line([kinematics.(jointNames{iJoint})(1, timesToPlot(iTime)) kinematics.(jointNames{iJoint+1})(1, timesToPlot(iTime))],...
            [kinematics.(jointNames{iJoint})(2, timesToPlot(iTime)) kinematics.(jointNames{iJoint+1})(2, timesToPlot(iTime))],...
            'Color', [0 0 0 trans], 'LineWidth', 2.5);
        
        %also plot joints as dots
        h = scatter(kinematics.(jointNames{iJoint})(1, timesToPlot(iTime)), kinematics.(jointNames{iJoint})(2, timesToPlot(iTime)), 40,....
            'Marker','o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [plotColors(iJoint, :)], ...
            'MarkerEdgeAlpha', trans, 'MarkerFaceAlpha', trans);
        
        %for legend
        if iTime == 1
            jointLegendH(iJoint) = h;
        end
        
    end
    
    %plot last joint
    h = scatter(kinematics.(jointNames{end})(1, timesToPlot(iTime)), kinematics.(jointNames{end})(2, timesToPlot(iTime)), 40, ....
        'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', plotColors(length(jointNames), :), ...
        'MarkerEdgeAlpha', trans, 'MarkerFaceAlpha', trans);
    
    if iTime == 1
        jointLegendH(length(jointNames)) = h;
    end
    
end

axis square

% also put a scale bar
xlims = get(gca,'xLim');
ylims = get(gca,'yLim');
line([xlims(1) xlims(1)+0.1], [ylims(1) ylims(1)], 'color', 'k', 'linewidth', 4)
line([xlims(1) xlims(1)], [ylims(1) ylims(1)+0.1], 'color', 'k', 'linewidth', 4)
text(xlims(1), ylims(1)-0.02, '10 cm', 'FontSize', 14)

% legend
if makeLegend
    legendH = legend(jointLegendH, legendNames, 'Box', 'off', 'FontSize', 14);
end



function [walkDuration, walkHeight, walkObsDuration, walkObsHeight] = getStepProps(trialsData, walkTrialInds, walkObsTrialInds, nSteps)
% funtion to calculat the step heights and durations for basic walking and
% all of the walking obstacle steps

% first get basic walking
for iTrial = 1:length(walkObsTrialInds)
    %stride duration (in ms)
    walkDuration(iTrial) = size(trialsData(walkTrialInds(iTrial)).SpikeCounts, 2)*10;
    
    %step height (toe) (in cm)
    toeKin = trialsData(walkTrialInds(iTrial)).Kinematics.PinkyToe;
    walkHeight(iTrial) = (max(toeKin(2,:)) - min(toeKin(2,:)))*100;
end

% then each of the walking obstacle steps
for iTrial = 1:length(walkObsTrialInds)
    
    for iStep = 1:nSteps
        
        [trialStartInd, trialStartEvent, trialEndInd, trialEndEvent] = extractGaitInds(trialsData, 'WalkingObstacle', '', false, iStep-4);
        if isempty(trialEndEvent)
            stepInds(1) = trialsData(walkObsTrialInds(iTrial)).TrialEvents(trialStartEvent);
            stepInds(2) = size(trialsData(walkObsTrialInds(iTrial)).SpikeCounts,2);
        else
            stepInds = trialsData(walkObsTrialInds(iTrial)).TrialEvents([trialStartEvent trialEndEvent]);
            %don't include the next foot strike as part of the current step
            stepInds(2) = stepInds(2) - 1;
        end
        
        %stride duration (in ms)
        walkObsDuration(iTrial, iStep) = diff(stepInds)*10;
       
        %step height (toe) (in cm)
        toeKin = trialsData(walkObsTrialInds(iTrial)).Kinematics.PinkyToe(:,stepInds(1):stepInds(2));
        walkObsHeight(iTrial, iStep) = (max(toeKin(2,:)) - min(toeKin(2,:)))*100;
        
    end    
    
end


% 
