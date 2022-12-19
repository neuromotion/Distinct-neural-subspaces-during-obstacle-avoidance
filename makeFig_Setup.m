function subfuncs = makeFig_Setup(varargin)
% Make experimental setup and example data figure
% 
% Layout:
% A------------------------  B---------------------
% |  Recording setup and     | Behavioral paradigm
% |  obstacle rig            | 
% |                          | 
% 
% C------------------------------------------------
% |  Example data from a trial from Boomer
% |
% |  population raster
% |  toe height
% |  toe horizontal position
% |  stepping diagram
% |  obstacle position

if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 10.1, 10])

figMargins = 0.05;
subplotGap = 0.1;


%% Subplot A -- Rig and setup
% don't need to do anything here, will just import the graphic



%% Subplot B -- Paradigm 
% don't need to do anything here, will just import the graphic



%% Subplot C -- Example recording dataset
% load in Boomer data
load('./Data/TrialsDataBoomer.mat')

% example trial to use
exTrial = 134;

% only use up to 2nd step after obstacle step
trialLength = trialsLegM1(exTrial).TrialEvents(13);
trialTimes = 0.01:0.01:trialLength/100;

% make raster plot
exampleRasterH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [figMargins-0.4, 1.8 + subplotGap*4 + figMargins, 11, 2.5]);

% remove last step
spikeData = trialsLegM1(exTrial).SpikeTimes;
for iNeuron = 1:length(spikeData)
    spikeData{iNeuron}(spikeData{iNeuron} > trialLength*10) = [];
end
rasterplot(fliplr(spikeData),'Times','|')

xlim([0 trialLength*10])
box off
axis off

% make toe height plot
exampleToeHeightH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [figMargins-0.4, 1.2 + subplotGap*4 + figMargins, 11, 2.3]);

toeHeight = trialsLegM1(exTrial).Kinematics.PinkyToe(2,1:trialLength);

plot(trialTimes, toeHeight, 'LineWidth', 2, 'color', [0 0.4470 0.7410])

box off
set(exampleToeHeightH, 'TickDir', 'out', 'XColor', 'none', 'LineWidth', 2, ...
    'FontSize', 12, 'YTick', [0 0.2], 'YLim', [0 0.22], 'xlim', [0 trialLength/100])
ylabel('Toe Y (m)', 'FontSize', 11)

% make toe x position
exampleToePositionH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [figMargins-0.4, 0.5 + subplotGap*4 + figMargins, 11, 2.3]);

toePosition = trialsLegM1(exTrial).Kinematics.PinkyToe(1,1:trialLength);

plot(trialTimes, toePosition, 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980])

box off
set(exampleToePositionH, 'TickDir', 'out', 'XColor', 'none', 'LineWidth', 2, ...
    'FontSize', 12, 'YTick', [-0.2 0.3], 'YLim', [-0.3 0.3], 'xlim', [0 trialLength/100])
ylabel('Toe X (m)', 'FontSize', 11)


% stepping diagram
exampleStepDiagramH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [figMargins-0.4, 0.3 + figMargins, 11, 2.3]);

% get limb gait events
startFrame = trialsLegM1(exTrial).startFrame;

right_footStrike = [64464 64563 64665 64747 64851 64940 65033 65129]-startFrame;
right_footOff = [64528 64637 64725 64809 64910 65002 65100 65195]-startFrame;

for iStep = 1:length(right_footStrike)
    patch([repmat(right_footStrike(iStep)/100,1,2) repmat(right_footOff(iStep)/100,1,2)],...
        [0 1 1 0], [0.4940, 0.1840, 0.5560], 'EdgeColor','none')
end

right_handStrike = [startFrame 64516 64616 64756 64855 64956 65070 65178]-startFrame;
right_handOff = [64487 64588 64698 64828 64929 65040 65151 65256]-startFrame;

for iStep = 1:length(right_handStrike)
    patch([repmat(right_handStrike(iStep)/100,1,2) repmat(right_handOff(iStep)/100,1,2)],...
        [1 2 2 1], [0.4660, 0.6740, 0.1880], 'EdgeColor','none')
end

left_footStrike = [startFrame 64512 64616 64709 64804 64893 64982 65078 65175]-startFrame;
left_footOff = [64483 64584 64680 64772 64865 64956 65052 65146 65249]-startFrame;

for iStep = 1:length(left_footStrike)
    patch([repmat(left_footStrike(iStep)/100,1,2) repmat(left_footOff(iStep)/100,1,2)],...
        [3.2 4.2 4.12 3.2], [0.4940, 0.1840, 0.5560], 'EdgeColor','none')
end

left_handStrike = [64466 64564 64671 64812 64905 65010 65124 65231]-startFrame;
left_handOff = [64536 64639 64757 64876 64983 65096 65201 65312]-startFrame;

for iStep = 1:length(left_handStrike)
    patch([repmat(left_handStrike(iStep)/100,1,2) repmat(left_handOff(iStep)/100,1,2)],...
        [2.2 3.2 3.2 2.2], [0.4660, 0.6740, 0.1880], 'EdgeColor','none')
end

xlim([0 trialLength/100])
box off
axis off


% obstacle position
exampleObstaclePosH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [figMargins-0.4, figMargins-0.4, 11, 2.1]);

obsPos = trialsLegM1(exTrial).ObstaclePosition(1:trialLength);
plot(trialTimes, obsPos, 'linewidth',3, 'color', [0.3 0.3 0.3])
box off

set(exampleObstaclePosH, 'TickDir', 'out', 'LineWidth', 2, ...
    'FontSize', 12, 'YTick', [0 1.3], 'YLim', [-0.1 1.3], 'xlim', [0 trialLength/100],...
    'XTick', [0:5 trialLength/100])
ylabel({'Obstacle', 'Position (m)'}, 'FontSize', 11)
xlabel('Time (s)')

% also add dotted lines separating gait cycles
stepEnds = trialsLegM1(exTrial).TrialEvents(3:2:12)/100;
for iStep = 1:length(stepEnds)
    line(repmat(stepEnds(iStep),1,2), [0 2], 'linestyle','--','linewidth',2,'color','r')
end

% add line indicating when each limb crosses the obstacle
% manually found for this particular trial
leadingLimbCrossTime = (64788 - trialsLegM1(exTrial).startFrame)/100;
laggingLimbCrossTime = (64823 - trialsLegM1(exTrial).startFrame)/100;

line(repmat(leadingLimbCrossTime,1,2), [0 2], 'linestyle','--','linewidth',2,'color','k')
line(repmat(laggingLimbCrossTime,1,2), [0 2], 'linestyle','--','linewidth',2,'color','k')


% 
