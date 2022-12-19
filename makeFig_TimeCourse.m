function subfuncs = makeFig_TimeCourse(varargin)

% Code to make figure showing kinematics of obstacle stepping for paper
% 
% Layout:
% A-----------------------------------------------------
% |  Boomer Decoding
% |  Toe height timeseries error            MSE by step
% |
% 
% B-----------------------------------------------------
% |  Starbuckk Decoding
% |  Toe height timeseries error            MSE by step
% |  
% 
% C------------------------ D---------------------------
% | Boomer                  | Starbuck mahal dists
% | mahal dists  cross-corr | mahal dists    cross-corr
% |

% expose subfunctions in case any other figures want to use them
subfuncs.getData = @getData;
subfuncs.doDecoding = @doDecoding;
subfuncs.getMahal = @getMahal;

% if something only wants to get the subfunctions, don't actually run the
% script to make the figure
if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 12, 8.2])

figMargins = 0.05;
subplotGap = 0.1;

% first do boomer
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

% defined from walking data
dutyPercent = 67;

[kinDataWalkObs, frDataWalkObs, spikeCountsWalkObs, obsMoveWalkObs] = getData(trialsDataBoomer, 'WalkingObstacle', dutyPercent);
% don't use 3rd step after obstacle
kinDataWalkObs = kinDataWalkObs(:,:,1:6,:);
frDataWalkObs = frDataWalkObs(:,:,1:6,:);

% decoding
neurStepCatWalkObs = permute(frDataWalkObs,[1 4 2 3]);
neurStepCatWalkObs = permute(neurStepCatWalkObs(:,:,:), [1 3 2]);
kinStepCatWalkObs = permute(kinDataWalkObs,[1 4 2 3]);
kinStepCatWalkObs = permute(kinStepCatWalkObs(:,:,:), [1 3 2]);
nTrials = size(neurStepCatWalkObs,3);
nSteps = size(frDataWalkObs,3);

% do all trained CV
[performance, trainWeights, estKins] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'CV', 0, nTrials);
allTrainErr_Boomer = cat(3,performance.absErr{:});
% get metrics on individual steps
[performance, trainWeights, estKins] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
allTrainMSE_Boomer = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
allTrainMSE_Boomer = cat(3,allTrainMSE_Boomer{:});

% train on obstacle step, decode rest
[~, obsWeights] = doDecoding(neurStepCatWalkObs(:,291:400,:), kinStepCatWalkObs(:,291:400,:), 10, 'Train', 0);
[performance, ~, ~] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'Test', 0, [], obsWeights);
obsTrainErr_Boomer =  mat2cell(performance.absErr, size(performance.absErr,1), repmat(size(performance.absErr,2)/nTrials, 1, nTrials));
obsTrainErr_Boomer = cat(3,obsTrainErr_Boomer{:});

% CV for the obstacle step
[performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,291:400,:), kinStepCatWalkObs(:,291:400,:), 10, 'CV', 0, nTrials);
obsTrainErrCV = cat(3,performance.absErr{:});
obsTrainErr_Boomer(:,301:400,:) = obsTrainErrCV;
obsTrainCVMSE_Boomer = [performance.MSE{:}];

% also get MSE on each step individually
for iStep = 1:size(frDataWalkObs,3)
    for iTrial = 1:nTrials
        [performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial),...
            kinStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial), 10, 'Test', 0, [], obsWeights);
        obsTrainMSE_Boomer(:,iTrial,iStep) = performance.MSE;
    end
end
obsTrainMSE_Boomer(:,:,4) = obsTrainCVMSE_Boomer;

% train on unobstructed walk step, decode rest
[~, walkWeights] = doDecoding(neurStepCatWalkObs(:,1:100,:), kinStepCatWalkObs(:,1:100,:), 10, 'Train', 0);
[performance, ~, ~] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'Test', 0, [], walkWeights);
walkTrainErr_Boomer =  mat2cell(performance.absErr,size(performance.absErr,1), repmat(size(performance.absErr,2)/nTrials, 1, nTrials));
walkTrainErr_Boomer = cat(3,walkTrainErr_Boomer{:});

% CV for the walk step
[performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1:100,:), kinStepCatWalkObs(:,1:100,:), 10, 'CV', 0, nTrials);
walkTrainErrCV = cat(3,performance.absErr{:});
walkTrainErr_Boomer(:,1:90,:) = walkTrainErrCV;
walkTrainCVMSE_Boomer = [performance.MSE{:}];

% also get MSE on each step individually
for iStep = 1:size(frDataWalkObs,3)
    for iTrial = 1:nTrials
        [performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial),...
            kinStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial), 10, 'Test', 0, [], walkWeights);
        walkTrainMSE_Boomer(:,iTrial,iStep) = performance.MSE;
    end
end
walkTrainMSE_Boomer(:,:,1) = walkTrainCVMSE_Boomer;


% next, look at decoding with dPCA components
subfuncs = makeFig_dPCA(false);

% run dPCA and analysis on neural data
[spikeDataAverage_Boomer, WNeur_Boomer, VNeur_Boomer, whichMargNeur_Boomer, explVarNeur_Boomer, trajInvNeur_Boomer,...
    trajDepNeur_Boomer, prinAnglesNeur_Boomer, allAnglesNeur_Boomer] = subfuncs.dpcaAnalysis(frDataWalkObs, [5 5], size(frDataWalkObs,4), false);
trajInvNeurCat_Boomer = squeeze(mat2cell(trajInvNeur_Boomer, 5, 100, ones(1,6), size(trajInvNeur_Boomer,4)));
trajInvNeurCat_Boomer = squeeze(cat(2,trajInvNeurCat_Boomer{:}));
trajDepNeurCat_Boomer = squeeze(mat2cell(trajDepNeur_Boomer, 5, 100, ones(1,6), size(trajInvNeur_Boomer,4)));
trajDepNeurCat_Boomer = squeeze(cat(2,trajDepNeurCat_Boomer{:}));

% train an test (CV) with just step invarient
[performance, trainWeights, estKins] = doDecoding(trajInvNeurCat_Boomer, kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
invCVMSE_Boomer = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
invCVMSE_Boomer = cat(3,invCVMSE_Boomer{:});


% train an test (CV) with just step dependent
[performance, trainWeights, estKins] = doDecoding(cat(1, trajDepNeurCat_Boomer, trajInvNeurCat_Boomer), kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
depCVMSE_Boomer = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
depCVMSE_Boomer = cat(3,depCVMSE_Boomer{:});



% Now do Starbuck
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

% defined from walking data
dutyPercent = 69;

[kinDataWalkObs, frDataWalkObs, spikeCountsWalkObs, obsMoveWalkObs] = getData(trialsDataStarbuck, 'WalkingObstacle', dutyPercent);

% decoding
neurStepCatWalkObs = permute(frDataWalkObs,[1 4 2 3]);
neurStepCatWalkObs = permute(neurStepCatWalkObs(:,:,:), [1 3 2]);
kinStepCatWalkObs = permute(kinDataWalkObs,[1 4 2 3]);
kinStepCatWalkObs = permute(kinStepCatWalkObs(:,:,:), [1 3 2]);
nTrials = size(neurStepCatWalkObs,3);
nSteps = size(frDataWalkObs,3);

% do all trained CV
[performance, trainWeights, estKins] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'CV', 0, nTrials);
allTrainErr_Starbuck = cat(3,performance.absErr{:});
% get metrics on individual steps
[performance, trainWeights, estKins] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
allTrainMSE_Starbuck = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
allTrainMSE_Starbuck = cat(3,allTrainMSE_Starbuck{:});

% train on obstacle step, decode rest
[~, obsWeights] = doDecoding(neurStepCatWalkObs(:,291:400,:), kinStepCatWalkObs(:,291:400,:), 10, 'Train', 0);
[performance, ~, ~] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'Test', 0, [], obsWeights);
obsTrainErr_Starbuck =  mat2cell(performance.absErr, size(performance.absErr,1), repmat(size(performance.absErr,2)/nTrials, 1, nTrials));
obsTrainErr_Starbuck = cat(3,obsTrainErr_Starbuck{:});

% CV for the obstacle step
[performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,291:400,:), kinStepCatWalkObs(:,291:400,:), 10, 'CV', 0, nTrials);
obsTrainErrCV = cat(3,performance.absErr{:});
obsTrainErr_Starbuck(:,301:400,:) = obsTrainErrCV;
obsTrainCVMSE_Starbuck = [performance.MSE{:}];

% also get MSE on each step individually
for iStep = 1:size(frDataWalkObs,3)
    for iTrial = 1:nTrials
        [performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial),...
            kinStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial), 10, 'Test', 0, [], obsWeights);
        obsTrainMSE_Starbuck(:,iTrial,iStep) = performance.MSE;
    end
end
obsTrainMSE_Starbuck(:,:,4) = obsTrainCVMSE_Starbuck;

% train on unobstructed walk step, decode rest
[~, walkWeights] = doDecoding(neurStepCatWalkObs(:,1:100,:), kinStepCatWalkObs(:,1:100,:), 10, 'Train', 0);
[performance, ~, ~] = doDecoding(neurStepCatWalkObs, kinStepCatWalkObs, 10, 'Test', 0, [], walkWeights);
walkTrainErr_Starbuck =  mat2cell(performance.absErr,size(performance.absErr,1), repmat(size(performance.absErr,2)/nTrials, 1, nTrials));
walkTrainErr_Starbuck = cat(3,walkTrainErr_Starbuck{:});

% CV for the walk step
[performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1:100,:), kinStepCatWalkObs(:,1:100,:), 10, 'CV', 0, nTrials);
walkTrainErrCV = cat(3,performance.absErr{:});
walkTrainErr_Starbuck(:,1:90,:) = walkTrainErrCV;
walkTrainCVMSE_Starbuck = [performance.MSE{:}];

% also get MSE on each step individually
for iStep = 1:size(frDataWalkObs,3)
    for iTrial = 1:nTrials
        [performance, ~, ~] = doDecoding(neurStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial),...
            kinStepCatWalkObs(:,1+(iStep-1)*100:iStep*100,iTrial), 10, 'Test', 0, [], walkWeights);
        walkTrainMSE_Starbuck(:,iTrial,iStep) = performance.MSE;
    end
end
walkTrainMSE_Starbuck(:,:,1) = walkTrainCVMSE_Starbuck;


% run dPCA and analysis on neural data
[spikeDataAverage_Starbuck, WNeur_Starbuck, VNeur_Starbuck, whichMargNeur_Starbuck, explVarNeur_Starbuck, trajInvNeur_Starbuck,...
    trajDepNeur_Starbuck, prinAnglesNeur_Starbuck, allAnglesNeur_Starbuck] = subfuncs.dpcaAnalysis(frDataWalkObs, [5 5], size(frDataWalkObs,4), false);
trajInvNeurCat_Starbuck = squeeze(mat2cell(trajInvNeur_Starbuck, 5, 100, ones(1,6), size(trajInvNeur_Starbuck,4)));
trajInvNeurCat_Starbuck = squeeze(cat(2,trajInvNeurCat_Starbuck{:}));
trajDepNeurCat_Starbuck = squeeze(mat2cell(trajDepNeur_Starbuck, 5, 100, ones(1,6), size(trajInvNeur_Starbuck,4)));
trajDepNeurCat_Starbuck = squeeze(cat(2,trajDepNeurCat_Starbuck{:}));

% train an test (CV) with just step invarient
[performance, trainWeights, estKins] = doDecoding(trajInvNeurCat_Starbuck, kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
invCVMSE_Starbuck = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
invCVMSE_Starbuck = cat(3,invCVMSE_Starbuck{:});


% train an test (CV) with just step dependent
[performance, trainWeights, estKins] = doDecoding(cat(1, trajDepNeurCat_Starbuck, trajInvNeurCat_Starbuck), kinStepCatWalkObs, 10, 'CV', nSteps, nTrials);
depCVMSE_Starbuck = cellfun(@(x) [x.MSE{:}], performance, 'un', 0);
depCVMSE_Starbuck = cat(3,depCVMSE_Starbuck{:});



%% Subplot A -- boomer absolute toe error
boomerDecodingErrH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [subplotGap-1, 4.5 + subplotGap*2 + figMargins, 12.5, 3]);

colorMap = parula(8);

errlegendH(1) = plot(11:600, mean(squeeze(walkTrainErr_Boomer(14,:,:))')*100, 'color', colorMap(1,:), 'linewidth', 2);
hold on;
shadedErrorBar(11:600, mean(squeeze(walkTrainErr_Boomer(14,:,:))')*100, std(squeeze(walkTrainErr_Boomer(14,:,:))')*100,...
    {'color',colorMap(1,:), 'linewidth', 2},1)

errlegendH(2) = plot(11:600, mean(squeeze(obsTrainErr_Boomer(14,:,:))')*100, 'color', colorMap(4,:), 'linewidth', 2);
shadedErrorBar(11:600, mean(squeeze(obsTrainErr_Boomer(14,:,:))')*100, std(squeeze(obsTrainErr_Boomer(14,:,:))')*100,...
    {'color', colorMap(4,:), 'linewidth', 2},1)

legend(errlegendH, 'Pre-Obs step trained', 'Obs step trained', 'Box', 'off')

box off
set(boomerDecodingErrH, 'FontSize',12, 'TickDir','out','TickLength', [0.01 0.01],...
    'LineWidth', 2)

ylabel('Decoding error (cm)')
xlabel('Gait %')


%% Subplot B -- Starbuck absolute toe error
starbuckDecodingErrH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [subplotGap-1, 2.3 + subplotGap*2 + figMargins, 12.5, 3]);

shadedErrorBar(11:600, mean(squeeze(walkTrainErr_Starbuck(14,:,:))')*100, std(squeeze(walkTrainErr_Starbuck(14,:,:))')*100,...
    {'color',colorMap(1,:), 'linewidth', 2},1)
hold on
shadedErrorBar(11:600, mean(squeeze(obsTrainErr_Starbuck(14,:,:))')*100, std(squeeze(obsTrainErr_Starbuck(14,:,:))')*100,...
    {'color',colorMap(4,:), 'linewidth', 2},1)

box off
set(starbuckDecodingErrH, 'FontSize',12, 'TickDir','out','TickLength', [0.01 0.01],...
    'LineWidth', 2)

ylabel('Decoding error (cm)')
xlabel('Gait %')


%% Subplot C -- Boomer 
boomerDecodingStepsH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [subplotGap-1, 0.1, 7.2, 3.4]);

% for each of the steps, show pre-obs trained, obs trained, all trained,
% inv dPCA trained, and all dPCA trained
nBarTypes = 5;
xVals = setdiff(1:(nBarTypes+1)*nSteps, nBarTypes+1:nBarTypes+1:(nBarTypes+1)*nSteps);

% concatentate all values across kin variables and trials
walkTrainCat = reshape(walkTrainMSE_Boomer,size(walkTrainMSE_Boomer,1)*size(walkTrainMSE_Boomer,2),size(walkTrainMSE_Boomer,3));
obsTrainCat = reshape(obsTrainMSE_Boomer,size(obsTrainMSE_Boomer,1)*size(obsTrainMSE_Boomer,2),size(obsTrainMSE_Boomer,3));
allTrainCat = reshape(allTrainMSE_Boomer,size(allTrainMSE_Boomer,1)*size(allTrainMSE_Boomer,2),size(allTrainMSE_Boomer,3));
invTrainCat = reshape(invCVMSE_Boomer,size(invCVMSE_Boomer,1)*size(invCVMSE_Boomer,2),size(invCVMSE_Boomer,3));
depTrainCat = reshape(depCVMSE_Boomer,size(depCVMSE_Boomer,1)*size(depCVMSE_Boomer,2),size(depCVMSE_Boomer,3));

barMeans = cat(1, squeeze(mean(walkTrainMSE_Boomer(14,:,:)))', squeeze(mean(obsTrainMSE_Boomer(14,:,:)))', squeeze(mean(allTrainMSE_Boomer(14,:,:)))',...
    squeeze(mean(invCVMSE_Boomer(14,:,:)))', squeeze(mean(depCVMSE_Boomer(14,:,:)))');

barStds = cat(1, squeeze(std(walkTrainMSE_Boomer(14,:,:)))', squeeze(std(obsTrainMSE_Boomer(14,:,:)))', squeeze(std(allTrainMSE_Boomer(14,:,:)))',...
    squeeze(std(invCVMSE_Boomer(14,:,:)))', squeeze(std(depCVMSE_Boomer(14,:,:)))');

barMeans = barMeans(:);
barStds = barStds(:);

boomerBarH = bar(xVals, barMeans, 'EdgeAlpha', 0, 'FaceColor', 'flat');

% colors
barColorMap = lines(nBarTypes);
barColorMap(4,:) = barColorMap(2,:);
barColorMap(1,:) = colorMap(1,:);
barColorMap(2,:) = colorMap(4,:);

iColor = 1;
for iBar = 1:size(boomerBarH.CData, 1)
    boomerBarH.CData(iBar,:) = barColorMap(iColor,:);
    if iColor == nBarTypes
        iColor = 1;
    else
        iColor = iColor+1;
    end
end

% legend
hold on
for iLeg = 1:nBarTypes
    legH(iLeg) = plot(0,0,'color', barColorMap(iLeg,:));
end

legend(legH, {'Walk Trained', 'Obstacle Trained', 'All Trained', 'Invarient dPCA trained', 'All dPCA trained'},'FontSize',12,'Box','off')

% error bars
errorbar(xVals, barMeans, barStds, '.k', 'MarkerSize', 1, 'LineWidth', 2)

set(boomerDecodingStepsH, 'FontSize', 13, 'XTick', [3:6:33], 'TickDir', 'out', 'XTickLabel', {'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2', 'Obs Stride + 3'}, 'LineWidth', 2, 'XTickLabelRotation', 25)
box off
ylabel('Decoding MSE')


%% Subplot D -- Starbuck MSE plots
starbuckDecodingStepsH = axes('Units','inches','PositionConstraint','innerposition','OuterPosition',...
    [4.6, 0.1, 7.2, 3.4]);

% for each of the steps, show pre-obs trained, obs trained, all trained,
% inv dPCA trained, and all dPCA trained
nBarTypes = 5;
xVals = setdiff(1:(nBarTypes+1)*nSteps, nBarTypes+1:nBarTypes+1:(nBarTypes+1)*nSteps);

% concatentate all values across kin variables and trials
walkTrainCat = reshape(walkTrainMSE_Starbuck,size(walkTrainMSE_Starbuck,1)*size(walkTrainMSE_Starbuck,2),size(walkTrainMSE_Starbuck,3));
obsTrainCat = reshape(obsTrainMSE_Starbuck,size(obsTrainMSE_Starbuck,1)*size(obsTrainMSE_Starbuck,2),size(obsTrainMSE_Starbuck,3));
allTrainCat = reshape(allTrainMSE_Starbuck,size(allTrainMSE_Starbuck,1)*size(allTrainMSE_Starbuck,2),size(allTrainMSE_Starbuck,3));
invTrainCat = reshape(invCVMSE_Starbuck,size(invCVMSE_Starbuck,1)*size(invCVMSE_Starbuck,2),size(invCVMSE_Starbuck,3));
depTrainCat = reshape(depCVMSE_Starbuck,size(depCVMSE_Starbuck,1)*size(depCVMSE_Starbuck,2),size(depCVMSE_Starbuck,3));

barMeans = cat(1, squeeze(mean(walkTrainMSE_Starbuck(14,:,:)))', squeeze(mean(obsTrainMSE_Starbuck(14,:,:)))', squeeze(mean(allTrainMSE_Starbuck(14,:,:)))',...
    squeeze(mean(invCVMSE_Starbuck(14,:,:)))', squeeze(mean(depCVMSE_Starbuck(14,:,:)))');

barStds = cat(1, squeeze(std(walkTrainMSE_Starbuck(14,:,:)))', squeeze(std(obsTrainMSE_Starbuck(14,:,:)))', squeeze(std(allTrainMSE_Starbuck(14,:,:)))',...
    squeeze(std(invCVMSE_Starbuck(14,:,:)))', squeeze(std(depCVMSE_Starbuck(14,:,:)))');

barMeans = barMeans(:);
barStds = barStds(:);

starbuckBarH = bar(xVals, barMeans, 'EdgeAlpha', 0, 'FaceColor', 'flat');
hold on

% colors
barColorMap = lines(nBarTypes);
barColorMap(4,:) = barColorMap(2,:);
barColorMap(1,:) = colorMap(1,:);
barColorMap(2,:) = colorMap(4,:);

iColor = 1;
for iBar = 1:size(starbuckBarH.CData, 1)
    starbuckBarH.CData(iBar,:) = barColorMap(iColor,:);
    if iColor == nBarTypes
        iColor = 1;
    else
        iColor = iColor+1;
    end
end

% error bars
errorbar(xVals, barMeans, barStds, '.k', 'MarkerSize', 1, 'LineWidth', 2)

set(starbuckDecodingStepsH, 'FontSize', 13, 'XTick', [3:6:33], 'TickDir', 'out', 'XTickLabel', {'Obs Stride - 3', 'Obs Stride - 2', ...
    'Obs Stride - 1', 'Obs Stride', 'Obs Stride + 1', 'Obs Stride + 2', 'Obs Stride + 3'}, 'LineWidth', 2, 'XTickLabelRotation', 25)
box off
ylabel('Decoding MSE')




function [kinData, frData, spikeCounts, obsMovePhase] = getData(trialsData, task, dutyPercent)

badTrials = filterTrials(trialsData,90,5);
taskTrialInds = find(strcmpi(task,string({trialsData.Task})));
taskTrialInds = setdiff(taskTrialInds, badTrials);

jointNames = fieldnames(trialsData(1).Kinematics);

% get spike data and smooth
% for gait normalization, there could be different number of pre and
% post-trial data, get the min of all trials
nPrePoints = min(cellfun(@(x) size(x,2), {trialsData(taskTrialInds).GaitNormalizedPreTrialSpikeCounts}));
nPostPoints = min(cellfun(@(x) size(x,2), {trialsData(taskTrialInds).GaitNormalizedPostTrialSpikeCounts}));
nSteps = size(trialsData(taskTrialInds(1)).GaitNormalizedSpikeCounts,2)/100;

for iTrial = 1:length(taskTrialInds)
    
    %get neural data
    allSpikes = [trialsData(taskTrialInds(iTrial)).GaitNormalizedPreTrialSpikeCounts(:,1:nPrePoints), ...
        trialsData(taskTrialInds(iTrial)).GaitNormalizedSpikeCounts, ...
        trialsData(taskTrialInds(iTrial)).GaitNormalizedPostTrialSpikeCounts(:, end-nPostPoints+1:end)];
    
    %smooth by convolving with gaussian
    smoothedFR = convGauss(allSpikes, 10, 4*10)*100;
    smoothedFR = smoothedFR(:, nPrePoints+1:end-nPostPoints);
    
    %divide into steps
    for iStep = 1:nSteps
        spikeCounts(:,:,iStep,iTrial) = allSpikes(:, (1+(iStep-1)*100:iStep*100) + nPrePoints);
        frData(:,:,iStep,iTrial) = smoothedFR(:,1+(iStep-1)*100:iStep*100);
    end
    
    %next get kinematics (don't use crest)
    kinCell = struct2cell(trialsData(taskTrialInds(iTrial)).GaitNormalizedKinematics);
    allKin = cat(1, kinCell{2:end});
    allKin = allKin(:,1:100*nSteps);
    
    %split into steps
    for iStep = 1:nSteps
        kinData(:,:,iStep,iTrial) = allKin(:,1+(iStep-1)*100:iStep*100);
    end
    
    %also get obstacle
    if strcmpi(task,'Walk')
        %no obstacle moving during basic walking
        continue
    end
    obsPos = trialsData(taskTrialInds(iTrial)).ObstaclePosition;
    obsMoveInd(iTrial) = find(diff(obsPos) > 0, 1);
    %get it in terms of normalized gait phase
    limbStrikes = [trialsData(taskTrialInds(iTrial)).TrialEvents(...
        contains(string(trialsData(taskTrialInds(iTrial)).TrialEventsLabels),'Limb Strike')), ...
        size(trialsData(taskTrialInds(iTrial)).SpikeCounts,2)];
    
    limbOffs = trialsData(taskTrialInds(iTrial)).TrialEvents(...
        contains(string(trialsData(taskTrialInds(iTrial)).TrialEventsLabels),'Limb Off'));
    
    %make sure we got the correct number of gait events
    assert(length(limbStrikes) == nSteps+1 && length(limbOffs) == nSteps, ...
        'Incorrect number of limb strike/off events found!');
    
    %determine which step the obstacle move was in
    obsMoveStep = find(limbStrikes < obsMoveInd(iTrial), 1, 'last');
    
    %convert to gait phase (starting index (1) is 0 %, ending index (length of the step) is
    %100)
    obsMoveStepPhase = interp1([limbStrikes(obsMoveStep)+1 limbStrikes(obsMoveStep+1) limbOffs(obsMoveStep)+1]-limbStrikes(obsMoveStep), ...
        [0 100 60], obsMoveInd(iTrial)-limbStrikes(obsMoveStep)+1);
    obsMovePhase(iTrial) = obsMoveStepPhase + 100*(obsMoveStep-1);
    
end

kinData = squeeze(kinData);
frData = squeeze(frData);


function mahalDist = getMahal(reference, data, nFolds)

if nFolds > 0
    cvFoldsInds = divideBlocks(randperm(size(reference,3)), nFolds);
end

for iPhase = 1:size(data,2)
    
    if nFolds > 0
        
        %do cross-validation
        nTrials = size(reference,3);
        for iFold = 1:nFolds
            
            testTrials = cvFoldsInds{iFold};
            trainTrials = setdiff(1:nTrials, testTrials);
            
            mahalDist(iFold,iPhase) = mahal(squeeze(reference(:,iPhase,testTrials))', squeeze(reference(:,iPhase,trainTrials))');
        end
        
    else
    
        %no cross-validation, just compare data to reference
        mahalDist(:,iPhase) = mahal(squeeze(data(:,iPhase,:))', squeeze(reference(:,iPhase,:))');
        
    end
    
end


function [decodingPerformance, trainWeights, estKin] = doDecoding(neur, kin, nHist, type, nSteps, nFolds, weights)

% add history
for iTrial = 1:size(neur,3)
    neurHist(:,:,iTrial) = addHistory(squeeze(neur(:,:,iTrial)), nHist, 2);
end

% remove the first nHist points (0's seem to affect the decoding)
neurHist(:,1:nHist,:) = [];
kin(:,1:nHist,:) = [];

if strcmpi(type,'Train')
    
    %no cross validation, just get training weights
    trainNeur = neurHist(:,:)-mean(neurHist(:,:),2);
    trainKin = kin(:,:)-mean(kin(:,:),2);
    
    trainWeights = trainKin/trainNeur;
    decodingPerformance = [];
    estKin = [];
    
elseif strcmpi(type, 'Test')
    
    %no cross validation, just use given weights to decode and get
    %performance metrics
    testNeur = neurHist(:,:)-mean(neurHist(:,:),2);
    realKin = kin(:,:)-mean(kin(:,:),2);
    
    %decode
    estKin = weights*testNeur;
    
    %get performance metrics
    trainWeights = [];
    decodingPerformance = calcPerformanceMetrics(estKin, realKin);
    
elseif strcmpi(type, 'CV')

    % do cross vaidation
    % divide into test and training sets
    nTrials = size(neur,3);
    cvFoldsInds = divideBlocks(randperm(nTrials), nFolds);
    
    for iFold = 1:nFolds
        
        testTrials = cvFoldsInds{iFold};
        trainTrials = setdiff(1:nTrials, testTrials);
        
        %get data and mean center
        trainNeur = neurHist(:,:,trainTrials);
        trainNeur = trainNeur(:,:)-mean(trainNeur(:,:),2);
        trainKin = kin(:,:,trainTrials);
        trainKinMeans = mean(trainKin(:,:),2);
        trainKin = trainKin(:,:)-trainKinMeans;
        
        testNeur = neurHist(:,:,testTrials);
        testNeur = testNeur(:,:)-mean(testNeur(:,:),2);
        testKin = kin(:,:,testTrials);
        realKin{iFold} = testKin(:,:)-mean(testKin(:,:),2);
        
        %get weights
        trainWeights{iFold} = trainKin/trainNeur;
        
        %decode
        estKin{iFold} = trainWeights{iFold}*testNeur;
        
    end
    
    %get performance metrics
    if nSteps <= 0
        %get metrics for all steps concatenated
        decodingPerformance = calcPerformanceMetrics(estKin, realKin);
    else
        %get metrics for individual steps
        for iStep = 1:nSteps
            if iStep == 1 
                %lost first 10 time points due to history
                estKinStep = cellfun(@(x) x(:,1:90), estKin, 'un', 0);
                realKinStep = cellfun(@(x) x(:,1:90), realKin, 'un', 0);
            else
                estKinStep = cellfun(@(x) x(:,(1+(iStep-1)*100:iStep*100)-10), estKin, 'un', 0);
                realKinStep = cellfun(@(x) x(:,(1+(iStep-1)*100:iStep*100)-10), realKin, 'un', 0);
            end
            decodingPerformance{iStep} = calcPerformanceMetrics(estKinStep, realKinStep);
        end
    end
    
end

    
% 
