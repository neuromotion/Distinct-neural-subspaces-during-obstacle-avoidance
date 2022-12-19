function subfuncs = makeFig_PLDS(varargin)
% Code to plot PLDS trajectories for paper
% 
% Layout:
% A------------------------  E--------------------- 
% |  Boomer PLDS             | Starbuck PLDS        
% |  PLDS plot with          | PLDS plot with view   
% |  view showing similar    | showing similar      
% |  trajectories            | trajectories          
% 
% B------------------------  F--------------------- 
% |  Boomer PLDS             | Starbuck              
% |  PLDS plot with          | PLDS plot with view  
% |  view showing diverging  | showing divering     
% |  trajectories            | trajectories          
% 
% C------------------------  G---------------------
% | Boomer kinematics        | Starbuck kinematics
% | Projection               | Projection
% |                          |
% 
% D------------------------  H--------------------- 
% | Boomer Spread            | Starbuck Spread 
% | metrics                  | metrics
% |                          |


% expose subfunctions in case any other figures want to use them
subfuncs.splitSteps = @splitSteps;
subfuncs.plotTrajs = @plotTrajs;
subfuncs.findBestOverlapView = @findBestOverlapView;
subfuncs.samplePLDS = @samplePLDS;
subfuncs.calcVAF = @calcVAF;
subfuncs.calcSpreadMetric = @calcSpreadMetric;

if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 9, 12.8])

figMargins = 0.05;
subplotGap = 0.1;


%% Subplot A -- Boomer PLDS trajectories step invarient proj

% first boomer, split into steps
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

minDim_Boomer = 9;
load(['./Data/PLDS/Boomer_PLDS_Dim' num2str(minDim_Boomer)])

trajBoomer = sepPLDSTraj{1};

% defined from walking data
dutyPercent = 67;

% get the duty cycle percentage
badTrials = filterTrials(trialsDataBoomer,90,5);
walkTrialInds = find(cellfun(@(x) strcmpi(x,'walk'), {trialsDataBoomer.Task}));
walkTrialInds = setdiff(walkTrialInds, badTrials);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsDataBoomer.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
    
walkTrialDurations_Boomer = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkTrialInds).SpikeCounts});
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataBoomer(walkObsTrialInds).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataBoomer(walkObsTrialInds).TrialEvents) walkObsTotalDurations];

steps_Boomer = -3:3;
walkObsTrialDurations_Boomer = diff(walkObsEvents(:,1:2:end),1,2)*10;

splitTrajs = splitSteps(trajBoomer, trialsDataBoomer, 7, dutyPercent);
% don't use 3rd step after obstacle
splitTrajs = splitTrajs([1 3 2 4:minDim_Boomer],:,1:6,:);

% find good overlap camera angle
meanTrajs = mean(splitTrajs, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuth, elevation] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadNeur_Boomer = calcSpreadMetric(splitTrajs, azimuth, elevation);

% make plot
boomerPLDSIndH = axes('Units','inches','OuterPosition',...
    [-0.4, 7.8, 5, 5.3]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;

meanTracesH = plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent);

xlabel('LD1')
ylabel('LD3')
zlabel('LD2')
view([azimuth elevation])

set(boomerPLDSIndH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)

% legend
stepLabels = join([repmat("Step",7,1),string(num2cell(-3:3))']);
legend(meanTracesH, stepLabels, 'box', 'off', 'location', 'best', 'NumColumns', 1)



%% Subplot B -- Boomer PLDS step dependent proj
boomerPLDSDepH = axes('Units','inches','OuterPosition',...
    [-0.4, 5, 5, 4.7]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;

% spread metric
spreadDep_Boomer = calcSpreadMetric(splitTrajs, 70, 12);

plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent)

xlabel('LD1')
ylabel('LD3')
zlabel('LD2')
view([70 12])

set(boomerPLDSDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)

%% Subplot C -- Boomer kinematic PCA trajectories

% get Kinematics
badTrials = filterTrials(trialsDataBoomer,90,5);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsDataBoomer.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

allKins = {trialsDataBoomer(walkObsTrialInds).Kinematics};

% concatenate joint posisitons into rows of a matrix
allKins = cellfun(@struct2cell, allKins, 'un', 0);
allKins = cellfun(@(x) cat(1,x{:}), allKins, 'un', 0);

% concatenate across trials
nTrialPoints = cellfun(@(x) size(x, 2), allKins);
allKinsCat = cat(2, allKins{:});

% run PCA
[kinematicProjs, kinematicPCs, kinematicsVars] = pca(allKinsCat');

% split back into trials
kinPCsTrials = mat2cell(kinematicPCs', size(kinematicPCs, 2), nTrialPoints);

% split into steps and normalize
dutyPercent = 67;
splitTrajsKin = splitSteps(kinPCsTrials, trialsDataBoomer, 7, dutyPercent);
splitTrajsKin = splitTrajsKin([1 3 2 4:end],:,1:6,:);

% get good overlap angle
meanTrajs = mean(splitTrajsKin, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuthKin, elevationKin] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadKin_Boomer = calcSpreadMetric(splitTrajsKin, azimuthKin, elevationKin);

% also, for normalizing spread later on, get the amplitude of the top 3
% components
normPCA_Boomer = norm([max(meanTrajs(1,:))-min(meanTrajs(1,:)), max(meanTrajs(2,:))-min(meanTrajs(2,:)), ...
    max(meanTrajs(3,:))-min(meanTrajs(3,:))]);

% make plot
boomerKinPCAH = axes('Units','inches','OuterPosition',...
    [-0.4, 2, 5, 5]);

meanTracesH = plotTrajs(splitTrajsKin, 1, 2, 3, plotColors(2:end,:), dutyPercent);
xlabel('PC1')
ylabel('PC3')
zlabel('PC2')

view(-20, 15);
set(boomerKinPCAH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)



%% Subplot D -- Comparing overlap metrics for Boomer

overlapPlotH = axes('Units','inches','OuterPosition',...
    [-0.4, -0.4, 5, 4]);

% get average and std for bar plots
% have two plots, neural inv to neural dep and neural inv to kin invar
meanSpreads = [mean(spreadNeur_Boomer), mean(spreadDep_Boomer), mean(spreadNeur_Boomer), mean(spreadKin_Boomer)];
stdSpreads = [std(spreadNeur_Boomer), std(spreadDep_Boomer), std(spreadNeur_Boomer), std(spreadKin_Boomer)];

% plot each of the trials change in spread metric separately
line(repmat([1;2],1,length(spreadNeur_Boomer)), [spreadNeur_Boomer; spreadDep_Boomer],'color',[0.3 0.3 0.3 0.5],'linewidth',2)
line(repmat([4;5],1,length(spreadNeur_Boomer)), [spreadNeur_Boomer; spreadKin_Boomer],'color',[0.3 0.3 0.3 0.5],'linewidth',2)

% plot the average spreads of neural plds vs kin pca
line([0.7 1.3], repmat(meanSpreads(1),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([1.7 2.3], repmat(meanSpreads(2),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([1 2], meanSpreads(1:2),'color',[0, 0.4470, 0.7410],'linewidth',2)

% same for starbuck
line([3.7 4.3], repmat(meanSpreads(3),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([4.7 5.3], repmat(meanSpreads(4),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([4 5], meanSpreads(3:4),'color',[0, 0.4470, 0.7410],'linewidth',2)

% add signifiance stars for significance tests (5% alpha lvl)
p_BoomerNeur = signrank(spreadDep_Boomer,spreadNeur_Boomer,'tail','right');
if p_BoomerNeur < 0.05
    text(1.4, 0.46,'*','FontSize',18)
    line([1 2], [0.45 0.45],'color','k','linewidth',2)
end

p_BoomerKin = signrank(spreadKin_Boomer,spreadNeur_Boomer,'tail','right');
if p_BoomerKin < 0.05
    text(4.4, 0.46,'*','FontSize',18)
    line([4 5], [0.45 0.45],'color','k','linewidth',3)
end

box off
ylabel('Spread Index')
xlim([0 6])
set(overlapPlotH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2,...
    'XTick',[1 2 4 5], 'XTickLabel',{'PLDS Inv', 'PLDS Dep', 'PLDS Inv', 'Kin Inv'}, 'XTickLabelRotation', 45)



%% Subplot E -- Starbuck PLDS trajectories step invarient proj

% first starbuck, split into steps
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

minDim_Starbuck = 9;
load(['./Data/PLDS/Starbuck_PLDS_Dim' num2str(minDim_Starbuck)])

trajStarbuck = sepPLDSTraj{1};

% defined from walking data
dutyPercent = 69;

% get the duty cycle percentage
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkTrialInds = find(cellfun(@(x) strcmpi(x,'walk'), {trialsDataStarbuck.Task}));
walkTrialInds = setdiff(walkTrialInds, badTrials);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsDataStarbuck.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
    
walkTrialDurations_Starbuck = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkTrialInds).SpikeCounts});
walkObsTotalDurations = cellfun(@(x) size(x, 2), {trialsDataStarbuck(walkObsTrialInds).SpikeCounts})';
walkObsEvents = [cat(1,trialsDataStarbuck(walkObsTrialInds).TrialEvents) walkObsTotalDurations];

steps_Starbuck = -3:3;
walkObsTrialDurations_Starbuck = diff(walkObsEvents(:,1:2:end),1,2)*10;

splitTrajs = splitSteps(trajStarbuck, trialsDataStarbuck, 6, dutyPercent);

% get the camera angle with the best overlap
meanTrajs = mean(splitTrajs, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuth, elevation] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadNeur_Starbuck = calcSpreadMetric(splitTrajs, azimuth, elevation);

% make plot
starbuckPLDSIndH = axes('Units','inches','OuterPosition',...
    [4.2, 7.8, 5, 5.3]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;

meanTracesH = plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent);

xlabel('LD1')
ylabel('LD2')
zlabel('LD3')
view([azimuth elevation])

set(starbuckPLDSIndH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)


%% Subplot F -- Starbuck PLDS step dependent proj
starbuckPLDSDepH = axes('Units','inches','OuterPosition',...
    [4.2, 5, 5, 4.7]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;


% spread metric
spreadDep_Starbuck = calcSpreadMetric(splitTrajs, 78, 2);

plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent)

xlabel('LD1')
ylabel('LD2')
zlabel('LD3')
view([78 2])

set(starbuckPLDSDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)

%% Subplot G -- Starbuck kinematic PCA trajectories

% get Kinematics
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsDataStarbuck.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);

allKins = {trialsDataStarbuck(walkObsTrialInds).Kinematics};

% concatenate joint posisitons into rows of a matrix
allKins = cellfun(@struct2cell, allKins, 'un', 0);
allKins = cellfun(@(x) cat(1,x{:}), allKins, 'un', 0);

% concatenate across trials
nTrialPoints = cellfun(@(x) size(x, 2), allKins);
allKinsCat = cat(2, allKins{:});

% run PCA
[kinematicProjs, kinematicPCs, kinematicsVars] = pca(allKinsCat');

% split back into trials
kinPCsTrials = mat2cell(kinematicPCs', size(kinematicPCs, 2), nTrialPoints);

% split into steps and normalize
dutyPercent = 67;
splitTrajsKin = splitSteps(kinPCsTrials, trialsDataStarbuck, 6, dutyPercent);
splitTrajsKin = splitTrajsKin(:,:,1:6,:);

% get good overlap angle
meanTrajs = mean(splitTrajsKin, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuthKin, elevationKin] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadKin_Starbuck = calcSpreadMetric(splitTrajsKin, azimuthKin, elevationKin);

% also, for normalizing spread later on, get the amplitude of the top 3
% components
normPCA_Starbuck = norm([max(meanTrajs(1,:))-min(meanTrajs(1,:)), max(meanTrajs(2,:))-min(meanTrajs(2,:)), ...
    max(meanTrajs(3,:))-min(meanTrajs(3,:))]);

% make plot
starbuckKinPCAH = axes('Units','inches','OuterPosition',...
    [4.2, 2, 5, 5]);

meanTracesH = plotTrajs(splitTrajsKin, 1, 2, 3, plotColors(2:end,:), dutyPercent);
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

view(5, -16);
set(starbuckKinPCAH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)



%% Subplot H -- Comparing overlap metrics for Starbuck

overlapPlotH = axes('Units','inches','OuterPosition',...
    [4.2, -0.4, 5, 4]);

% get average and std for bar plots
% have two plots, neural inv to neural dep and neural inv to kin invar
meanSpreads = [mean(spreadNeur_Starbuck), mean(spreadDep_Starbuck), mean(spreadNeur_Starbuck), mean(spreadKin_Starbuck)];
stdSpreads = [std(spreadNeur_Starbuck), std(spreadDep_Starbuck), std(spreadNeur_Starbuck), std(spreadKin_Starbuck)];

% plot each of the trials change in spread metric separately
line(repmat([1;2],1,length(spreadNeur_Starbuck)), [spreadNeur_Starbuck; spreadDep_Starbuck],'color',[0.3 0.3 0.3 0.5],'linewidth',2)
line(repmat([4;5],1,length(spreadNeur_Starbuck)), [spreadNeur_Starbuck; spreadKin_Starbuck],'color',[0.3 0.3 0.3 0.5],'linewidth',2)

% plot the average spreads of neural plds vs kin pca
line([0.7 1.3], repmat(meanSpreads(1),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([1.7 2.3], repmat(meanSpreads(2),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([1 2], meanSpreads(1:2),'color',[0, 0.4470, 0.7410],'linewidth',2)

% same for starbuck
line([3.7 4.3], repmat(meanSpreads(3),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([4.7 5.3], repmat(meanSpreads(4),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([4 5], meanSpreads(3:4),'color',[0, 0.4470, 0.7410],'linewidth',2)

% add signifiance stars for significance tests (5% alpha lvl)
p_StarbuckNeur = signrank(spreadDep_Starbuck,spreadNeur_Starbuck,'tail','right');
if p_StarbuckNeur < 0.05
    text(1.4, 0.46,'*','FontSize',18)
    line([1 2], [0.45 0.45],'color','k','linewidth',2)
end

p_StarbuckKin = signrank(spreadKin_Starbuck,spreadNeur_Starbuck,'tail','right');
if p_StarbuckKin < 0.05
    text(4.4, 0.46,'*','FontSize',18)
    line([4 5], [0.45 0.45],'color','k','linewidth',3)
end

box off
ylabel('Spread Index')
xlim([0 6])
set(overlapPlotH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2,...
    'XTick',[1 2 4 5], 'XTickLabel',{'PLDS Inv', 'PLDS Dep', 'PLDS Inv', 'Kin Inv'}, 'XTickLabelRotation', 45)




function splitTrajs = splitSteps(trajs, trialsData, nSteps, dutyPercent)
% functiton for normalizing trajectories to gait cycle, and splitting the
% gait cycle into another dimension in the data

if ~isempty(dutyPercent)
    
    % do normalization
    % get the trial indices
    badTrials = filterTrials(trialsData,90,5);
    walkObsInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsData.Task}));
    walkObsInds = setdiff(walkObsInds, badTrials);
    
    assert(length(walkObsInds) == length(trajs));
    
    for iTrial = 1:length(trajs)
        
        trajs{iTrial} = timeNormalize(trajs{iTrial}, ...
            [trialsData(walkObsInds(iTrial)).TrialEvents(1:nSteps*2) size(trialsData(walkObsInds(iTrial)).SpikeCounts,2)+1],...
            sort([dutyPercent:100:nSteps*100 0:100:nSteps*100]), 0:nSteps*100);
        
    end
    
end

% now split into step cycles for each trial
trajs = cellfun(@(x) mat2cell(x(:, 1:(100*nSteps)), size(x, 1), repmat(100,1,nSteps)), trajs, 'un', 0);
% then concatenate along 3rd dimension for the steps
trajs = cellfun(@(x) cat(3, x{:}), trajs,'un',0);

% finally concatenate along 4th dimenion for trials
splitTrajs = cat(4,trajs{:});



function meanTracesH = plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors, dutyPercent)
% function for plotting multi-trial trajectories

hold on;
for iStep = 1:size(splitTrajs, 3)
    
%     % plot each individual trial thinly
%     for iTrial = 1:size(splitTrajs, 4)
%         
%         plot3(squeeze(splitTrajs(dim1, :, iStep, iTrial)), squeeze(splitTrajs(dim2, :, iStep, iTrial)), ...
%             squeeze(splitTrajs(dim3, :, iStep, iTrial)), 'Color', [plotColors(iStep,:) 0.2], 'LineWidth', 0.5)
%         
%     end
    
    %now plot the average
    trialMeans = mean(splitTrajs, 4);
    
    meanTracesH(iStep) = plot3(squeeze(trialMeans(dim1, :, iStep)), squeeze(trialMeans(dim2, :, iStep)), ...
        squeeze(trialMeans(dim3, :, iStep)), 'Color', plotColors(iStep,:), 'LineWidth', 3);
    
    %plot dots at the stance-swing transitions
    plot3(squeeze(trialMeans(dim1, dutyPercent, iStep)), squeeze(trialMeans(dim2, dutyPercent, iStep)), ...
        squeeze(trialMeans(dim3, dutyPercent, iStep)), 'o', 'MarkerFaceColor', plotColors(iStep,:),...
        'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerSize', 10);
    
    if iStep ~= size(splitTrajs, 3)
    %plot dots at the step transitions
    plot3(squeeze(trialMeans(dim1, 100, iStep)), squeeze(trialMeans(dim2, 100, iStep)), ...
        squeeze(trialMeans(dim3, 100, iStep)), 'd', 'MarkerFaceColor', plotColors(iStep,:),...
        'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerSize', 10)
    end
    
end



function [azimuth, elevation] = findBestOverlapView(meanTrajs)
% function to find the best rotational angle that has the most overlap
% in the trajectories between all the steps

% exhaustive search
for iAz = 1:180
    for iEl = 1:180
        
        %get projection matrix
        projV = viewmtx(iAz, iEl);
        
        %project 
        for iStep = 1:size(meanTrajs, 3)
            
            projTraj(:,:,iStep) = projV(1:2,1:3) * squeeze(meanTrajs(:,:,iStep));
            
        end
        
        %calculate euclidean distance
        spreadMetric(iAz,iEl,:) = calcSpreadMetric(projTraj, [], []);
%         meanProjTraj = mean(projTraj, 3);
%         diffs = projTraj - repmat(meanProjTraj, 1, 1, size(projTraj,3));
%         euclidDists = squeeze(sqrt(sum(diffs.^2)));
%         spreadMetric(iAz,iEl,:) = max(euclidDists,[],2);
        
    end
end

% find the angles with the smallest max spreadMetric across the gait
% cycle
meanSpreads = max(spreadMetric,[],3);

minSpread = min(meanSpreads(:));
[azimuth, elevation] = find(meanSpreads==minSpread);


function spreadMetric = calcSpreadMetric(allTrajs, azimuth, elevation)
% function to calculate the spread metric of each trial given a camera view

for iTrial = 1:size(allTrajs,4)
    
%     if size(allTrajs,4)>1
%         figure('OuterPosition',[100 100 1000 500])
%     end
            
    %project
    for iStep = 1:size(allTrajs, 3)
        
        if isempty(azimuth)
            projTraj(:,:,iStep) = squeeze(allTrajs(:,:,iStep,iTrial));
        else
            projV = viewmtx(azimuth, elevation);
            projTraj(:,:,iStep) = projV(1:2,1:3) * squeeze(allTrajs(1:3,:,iStep,iTrial));
        end
            
%         %make plots
%         if size(allTrajs,4)>1
%             
%             subplot(1,2,1)
%             hold on;
%             if iStep == 1
%                 lineWidth = 3;
%             else
%                 lineWidth = 1;
%             end
%             plot(squeeze(projTraj(1,:,iStep)),squeeze(projTraj(2,:,iStep)),'LineWidth', lineWidth)
%             
%         end
        
    end
    
    %calculate euclidean distance 
    %go through each time point of the obstacle step and find best alignment of
    %other steps
    for iTime = 1:100
        
        basePoint = squeeze(projTraj(:,iTime,4));
        
        %only consider points within +/- 10 time points around the current
        %point
        dists = sqrt(sum( (projTraj(:,:,[1:3 5:end]) - repmat(basePoint,1,100, size(projTraj,3)-1)).^2 , 1));
        dists = dists(:,max(iTime-10,1):min(iTime+10,100),:);
        [minDists, minDistTimes] = min(dists,[],2);
%         for iStep = 1:5
%             minDistPoints(iStep,:) = squeeze(projTraj(:,minDistTimes(1,1,iStep),iStep+1));
%         end
        maxDist(iTime) = max(minDists);
        
    end
    
    %normalize by size of bounding box of the projection
    normFactor = sqrt(sum((max(max(projTraj,[],2),[],3)-min(min(projTraj,[],2),[],3)).^2));
    
    maxDistSort = sort(maxDist(20:80),'descend');
    spreadMetric(iTrial) = maxDistSort(6)/normFactor;
%     meanProjTraj = mean(projTraj(:,:,[1 3 4]), 3);
%     diffs = projTraj(:,:,[1 3 4]) - repmat(meanProjTraj, 1, 1, 3);
%     euclidDists = max(squeeze(sqrt(sum(diffs.^2,1))),[],2);
%     spreadMetric(iTrial) = mean(euclidDists); 

%     if size(allTrajs,4)>1
%         subplot(1,2,2)
%         plot(maxDist)
%     end
end



function [sampleSpikeTrain sampleSpikeTimes] = samplePLDS(traj, model, dims, nSamples)
% function to resample population spike trains given a PLDS model

for iTrial = 1:length(traj)
    
    %get intensity function (lambda)
    lambdas = model.C(:,dims) * traj{iTrial}(dims,:) + model.d;

    %sample from poisson
    for iSample = 1:nSamples
        
        sampleSpikeTrain{iSample, iTrial} = poissrnd(exp(lambdas));
        
    end
    
end



function vaf = calcVAF(realSpikes, modelSpikes, trialsData, dutyPercent)
% get PETHs

% normalize to gait cycle
% get the trial indices
badTrials = filterTrials(trialsData,90,5);
walkObsInds = find(cellfun(@(x) strcmpi(x,'walkingObstacle'), {trialsData.Task}));
walkObsInds = setdiff(walkObsInds, badTrials);

nSteps = length(trialsData(walkObsInds(1)).TrialEvents)/2;

assert(length(walkObsInds) == length(realSpikes));

% do normalization for each trial
for iTrial = 1:length(realSpikes)
    
    realSpikes{iTrial} = timeNormalize(realSpikes{iTrial}, ...
        [trialsData(walkObsInds(iTrial)).TrialEvents size(trialsData(walkObsInds(iTrial)).SpikeCounts,2)+1],...
        sort([dutyPercent:100:nSteps*100 0:100:nSteps*100]), 0:nSteps*100);
    
    modelSpikes{iTrial} = timeNormalize(modelSpikes{iTrial}, ...
        [trialsData(walkObsInds(iTrial)).TrialEvents size(trialsData(walkObsInds(iTrial)).SpikeCounts,2)+1],...
        sort([dutyPercent:100:nSteps*100 0:100:nSteps*100]), 0:nSteps*100);
    
end
    
% sum bins across trials
countsReal = sum(cat(3,realSpikes{:}),3);
countsModel = sum(cat(3,modelSpikes{:}),3);

% smooth a little using gaussian (repeat at ends to avoid edge artifacts)
smoothedCountsReal = convGauss(repmat(countsReal, 1, 3), 10, 30);
smoothedCountsReal = smoothedCountsReal(:,size(countsReal,2)+1:size(countsReal,2)*2);

smoothedCountsModel = convGauss(repmat(countsModel, 1, 3), 10, 30);
smoothedCountsModel = smoothedCountsModel(:,size(countsModel,2)+1:size(countsModel,2)*2);

% also mean center
smoothedCountsReal = smoothedCountsReal' - mean(smoothedCountsReal');
smoothedCountsModel = smoothedCountsModel' - mean(smoothedCountsModel');

% calculate VAF
vaf = (norm(smoothedCountsReal,'fro') - norm(smoothedCountsReal - smoothedCountsModel, 'fro')) /...
    norm(smoothedCountsReal,'fro');



% 

