function makeFig_ModelPLDS
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


figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 9, 9])

figMargins = 0.05;
subplotGap = 0.1;

load('PLDS_Spreads.mat')

%% Subplot A -- Boomer PLDS trajectories step invarient proj

% first boomer, split into steps
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

minDim_Boomer = 9;
load(['./Data/PLDS/BoomerModel_PLDS.mat'])

trajBoomer = sepPLDSTraj{1};

% defined from walking data
dutyPercent = 67;

steps_Boomer = -3:3;

splitTrajs = splitSteps(trajBoomer, trialsDataBoomer, 7, dutyPercent);
% don't use 3rd step after obstacle
splitTrajs = splitTrajs([1 3 2 4:minDim_Boomer],:,1:6,:);

% get good overlap angle
meanTrajs = mean(splitTrajs, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuth, elevation] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadModel_Boomer = calcSpreadMetric(splitTrajs, azimuth, elevation);

% make plot
boomerPLDSIndH = axes('Units','inches','OuterPosition',...
    [-0.4, 4.3, 5, 5]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;

meanTracesH = plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent);

xlabel('LD1')
ylabel('LD3')
zlabel('LD2')
view([44 74])

set(boomerPLDSIndH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)

% legend
stepLabels = join([repmat("Step",7,1),string(num2cell(-3:3))']);
legend(meanTracesH, stepLabels, 'box', 'off', 'location', 'best', 'NumColumns', 1)


%% Subplot B -- Comparing overlap metrics for Boomer

overlapPlotH = axes('Units','inches','OuterPosition',...
    [4.2, 4.8, 4, 4]);

% get average and std for bar plots
% have two plots, neural inv to neural dep and neural inv to kin invar
meanSpreads = [mean(spreadNeur_Boomer), mean(spreadModel_Boomer)];
stdSpreads = [std(spreadNeur_Boomer), std(spreadModel_Boomer)];

% plot each of the trials change in spread metric separately
line(repmat([1;2],1,length(spreadNeur_Boomer)), [spreadNeur_Boomer; spreadModel_Boomer],'color',[0.3 0.3 0.3 0.5],'linewidth',2)

% plot the average spreads of neural plds vs kin pca
line([0.7 1.3], repmat(meanSpreads(1),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([1.7 2.3], repmat(meanSpreads(2),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([1 2], meanSpreads(1:2),'color',[0, 0.4470, 0.7410],'linewidth',2)

% add signifiance stars for significance tests (5% alpha lvl)
p_BoomerNeur = signrank(spreadModel_Boomer,spreadNeur_Boomer,'tail','right');
if p_BoomerNeur < 0.05
    text(1.4, 0.46,'*','FontSize',18)
    line([1 2], [0.45 0.45],'color','k','linewidth',2)
end

box off
ylabel('Spread Index')
xlim([0 3])
set(overlapPlotH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2,...
    'XTick',[1 2 4 5], 'XTickLabel',{'PLDS Inv', 'PLDS Dep', 'PLDS Inv', 'Kin Inv'}, 'XTickLabelRotation', 45)



%% Subplot C -- Starbuck PLDS trajectories step invarient proj

% first starbuck, split into steps
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

minDim_Starbuck = 9;
load(['./Data/PLDS/StarbuckModel_PLDS.mat'])

trajStarbuck = sepPLDSTraj{1};

% defined from walking data
dutyPercent = 69;

% get the duty cycle percentage
splitTrajs = splitSteps(trajStarbuck, trialsDataStarbuck, 6, dutyPercent);

% get the camera angle with the best overlap
meanTrajs = mean(splitTrajs, 4);
meanTrajs = meanTrajs(1:3, :, :);
[azimuth, elevation] = findBestOverlapView(meanTrajs);

% get spread metric across trials
spreadModel_Starbuck = calcSpreadMetric(splitTrajs, azimuth, elevation);

% make plot
starbuckPLDSIndH = axes('Units','inches','OuterPosition',...
    [-0.4, -0.4, 5, 5]);

plotColors = parula(8);
plotColors(5,:) = [0.25 0.25 0.25];
dim1 = 1;
dim2 = 2;
dim3 = 3;

meanTracesH = plotTrajs(splitTrajs, dim1, dim2, dim3, plotColors(2:end,:), dutyPercent);

xlabel('LD1')
ylabel('LD2')
zlabel('LD3')
view([12 21])

set(starbuckPLDSIndH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2)


%% Subplot D -- Comparing overlap metrics for Starbuck

overlapPlotH = axes('Units','inches','OuterPosition',...
    [4.2, 0, 4, 4]);

% get average and std for bar plots
% have two plots, neural inv to neural dep and neural inv to kin invar
meanSpreads = [mean(spreadNeur_Starbuck), mean(spreadModel_Starbuck)];
stdSpreads = [std(spreadNeur_Starbuck), std(spreadModel_Starbuck)];

% plot each of the trials change in spread metric separately
line(repmat([1;2],1,length(spreadNeur_Starbuck)), [spreadNeur_Starbuck; spreadModel_Starbuck],'color',[0.3 0.3 0.3 0.5],'linewidth',2)

% plot the average spreads of neural plds vs kin pca
line([0.7 1.3], repmat(meanSpreads(1),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
hold on
line([1.7 2.3], repmat(meanSpreads(2),1,2),'color',[0, 0.4470, 0.7410],'linewidth',3)
line([1 2], meanSpreads(1:2),'color',[0, 0.4470, 0.7410],'linewidth',2)


% add signifiance stars for significance tests (5% alpha lvl)
p_StarbuckNeur = signrank(spreadModel_Starbuck,spreadNeur_Starbuck,'tail','right');
if p_StarbuckNeur < 0.05
    text(1.4, 0.36,'*','FontSize',18)
    line([1 2], [0.35 0.35],'color','k','linewidth',2)
end

box off
ylabel('Spread Index')
xlim([0 3])
set(overlapPlotH, 'FontSize',12, 'TickDir','out','TickLength', [0.03 0.03], 'LineWidth', 2,...
    'XTick',[1 2 4 5], 'XTickLabel',{'PLDS Inv', 'PLDS Dep', 'PLDS Inv', 'Kin Inv'}, 'XTickLabelRotation', 45)




function splitTrajs = splitSteps(trajs, trialsData, nSteps, dutyPercent)
% functiton for normalizing trajectories to gait cycle, and splitting the
% gait cycle into another dimension in the data

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

