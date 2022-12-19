function subfuncs = makeFig_dPCA(varargin)
% Code to plot dPCA analysis
% 
% Layout:
% A-------------------------------------------------   B-------------------------------------------------
% |  Boomer                                            |  Starbuck
% |  Top two components for step invariant subspace    |  Top two components for step invariant subspace
% |                                                    |
% |  Top two components for step depedent subspace     |  Top two components for step depedent subspace 
%                   
% C---------------------  D-------------------------   E------------------------   F--------------------
% |  Principle angles     |                            | Rotation fit              |
% |  Starbuck             | Boomer                     | Starbuck                  | Boomer      
% |                       |                            | Inv rotations             | Inv rotations
% |                       |                            | Dep rotations  R2 ratio   | Dep rotations  R2 ratio
% |                       |                            |                           |

% expose subfunctions in case any other figures want to use them
subfuncs.splitSteps = @splitSteps;
subfuncs.dPCAProj = @dPCAProj;
subfuncs.getRotationMetric = @getRotationMetric;
subfuncs.dpcaAnalysis = @dpcaAnalysis;
subfuncs.plotTrajVecField = @plotTrajVecField;

if nargin == 1
    return
end

figure('Color','w', 'Units', 'inches',...
    'OuterPosition',[2, 0.5, 14.2, 8.5])

figMargins = 0.05;
subplotGap = 0.1;

% first run dPCA analysis on neural as well as kinematic data
subFuncs = makeFig_PLDS(false);

% first Boomer
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

load('./Data/ModelDataBoomer.mat')
modelDataBoomer = modeledData;

[spikeData, kinData] = splitSteps(trialsDataBoomer, 'WalkingObstacle');
[kinModelSpikeData] = splitSteps(modelDataBoomer, 'WalkingObstacle',0);

% don't use 3rd step after obstacle
spikeData = spikeData(:,:,1:6,:);
kinData = kinData(:,:,1:6,:);
kinModelSpikeData = kinModelSpikeData(:,:,1:6,:);

% run dPCA and analysis on neural data
[spikeDataAverage_Boomer, WNeur_Boomer, VNeur_Boomer, whichMargNeur_Boomer, explVarNeur_Boomer, trajInvNeur_Boomer,...
    trajDepNeur_Boomer, prinAnglesNeur_Boomer, allAnglesNeur_Boomer, canonCorrNeur_Boomer] = dpcaAnalysis(spikeData, [5 5], size(spikeData,4), false);

% run with kin-modeled spike data
[modelDataAverage_Boomer, WModel_Boomer, VModel_Boomer, whichMargModel_Boomer, explVarModel_Boomer, trajInvModel_Boomer,...
    trajDepModel_Boomer, prinAnglesModel_Boomer, allAnglesModel_Boomer, canonCorrModel_Boomer] = dpcaAnalysis(kinModelSpikeData, [5 5], size(kinModelSpikeData,4), false);

%run with just basic steps for control noise level measurement
[noiseDataAverage_Boomer, WNoise_Boomer, VNoise_Boomer, whichMargNoise_Boomer, explVarNoise_Boomer, trajInvNoise_Boomer,...
    trajDepNoise_Boomer, prinAnglesNoise_Boomer, allAnglesNoise_Boomer, canonCorrNoise_Boomer] = dpcaAnalysis(spikeData(:,:,[1 6],:), [5 5], size(spikeData,4), false);

allNoDepNoise = squeeze(cat(4,trajDepNoise_Boomer(:,:,1,:), trajDepNoise_Boomer(:,:,2,:)));
upperNoDepNoise_Boomer = mean(allNoDepNoise,3) + 2*std(allNoDepNoise,[],3);
lowerNoDepNoise_Boomer = mean(allNoDepNoise,3) - 2*std(allNoDepNoise,[],3);

% get significance cutoff (99% and 1%)
for iDim = 1:5
    
    allDimNoise = allNoDepNoise(iDim,:,:);
    noiseCutoff(iDim) = prctile(allNoDepNoise(:),99.99);
    meanDeps = mean(trajDepNeur_Boomer,4);
    h_Boomer(iDim,:,:) = meanDeps(iDim,:,:)>noiseCutoff(iDim) | meanDeps(iDim,:,:)<-1*noiseCutoff(iDim);
    
    %get significant signals for noise analysis
    allSteps = zeros(1,600);
    allSigInds = {};
    allSigVals = {};
    for iStep = 1:6
        if any(h_Boomer(iDim,:,iStep))
            
            sig = squeeze(meanDeps(iDim,:,iStep));
            sigInds = squeeze(h_Boomer(iDim,:,iStep));
            allSigInds{iStep} = find(sigInds)+100*(iStep-1);
            allSigVals{iStep} = sig(sigInds);
            
        end
    end
%     xq = setdiff(1:600,[0 600 allSigInds{:}]);
%     yq = pchip([1:100 501:600 allSigInds{:}],[zeros(1,200) allSigVals{:}],xq);
%     allSteps([1:100 501:600 allSigInds{:}]) = [zeros(1,200) allSigVals{:}];
%     allsteps(xq) = yq;

    allSteps([allSigInds{:}]) = [allSigVals{:}];
    allSteps = reshape(allSteps,100,6);
    
    %get noise and signal components
    trajDepExtNoise_Boomer(iDim,:,:,:) = squeeze(trajDepNeur_Boomer(iDim,:,:,:)) - repmat(allSteps,1,1,size(trajDepNeur_Boomer,4));
    trajDepExtSig_Boomer(iDim,:,:,:) = repmat(allSteps,1,1,size(trajDepNeur_Boomer,4));
    
end

% get inv noise and signal components
meanInvs = mean(mean(trajInvNeur_Boomer,4),3);
trajInvExtSig_Boomer = repmat(meanInvs,1,1,6,size(trajDepNeur_Boomer,4));
trajInvExtNoise_Boomer = trajInvNeur_Boomer - trajInvExtSig_Boomer;

% % do t-test
% for iDim = 1:5
%     for iPhase = 1:100
%         for iStep = 1:6
%             [hRaw_Boomer(iDim,iPhase,iStep), p_Boomer(iDim,iPhase,iStep)] = ttest2(squeeze(trajDepNeur_Boomer(iDim,iPhase,iStep,:)),squeeze(allNoDepNoise(iDim,iPhase,:)));
%         end
%     end
% end
% 
% benHochCutoff = FDRcutoff(p_Boomer(:),0.01,false);
% benHochCutoff = 0.0001/3000;
% h_Boomer = p_Boomer<benHochCutoff;

% do it with kinematics too
[kinDataAverage_Boomer, WKin_Boomer, VKin_Boomer, whichMargKin_Boomer, explVarKin_Boomer, trajInvKin_Boomer,...
    trajDepKin_Boomer, prinAnglesKin_Boomer, allAnglesKin_Boomer, canonCorrKin_Boomer] = dpcaAnalysis(kinData, [5 5], size(spikeData,4),false);

% get obstacle stepping
[obsSpikeData, obsKinData] = splitSteps(trialsDataBoomer, 'Obstacle');

[trajInv, trajDep, prinAngles, allAngles] = dPCAProj(permute(obsSpikeData, [1 3 2 4]), WNeur_Boomer, VNeur_Boomer, whichMargNeur_Boomer, false);

% get the rotational fit
[anglesInvNeur_Boomer, MskewInvNeur_Boomer, skewRatioInvNeur_Boomer] = getRotationMetric(trajInvNeur_Boomer);
[anglesDepNeur_Boomer, MskewDepNeur_Boomer, skewRatioDepNeur_Boomer] = getRotationMetric(trajDepNeur_Boomer);

[anglesInvKin_Boomer, MskewInvKin_Boomer, skewRatioInvKin_Boomer] = getRotationMetric(trajInvKin_Boomer);
[anglesDepKin_Boomer, MskewDepKin_Boomer, skewRatioDepKin_Boomer] = getRotationMetric(trajDepKin_Boomer);

% compare timing with kinematics
subFuncsTimeCourse = makeFig_TimeCourse(false);
depMahal_Boomer = subFuncsTimeCourse.getMahal(squeeze(trajDepNeur_Boomer(:,:,1,:)), squeeze(trajDepNeur_Boomer(:,:,4,:)), 0);
invMahal_Boomer = subFuncsTimeCourse.getMahal(squeeze(trajInvNeur_Boomer(:,:,1,:)), squeeze(trajInvNeur_Boomer(:,:,4,:)), 0);
kinMahal_Boomer = subFuncsTimeCourse.getMahal(squeeze(kinData(:,:,1,:)), squeeze(kinData(:,:,4,:)), 0);
% find delay
[kinDepDiffxCorr_Boomer, kinDepDiffLags_Boomer] = crosscorr(mean(depMahal_Boomer),mean(kinMahal_Boomer));
[~, maxCorrInd] = max(kinDepDiffxCorr_Boomer);
depDelay_Boomer = kinDepDiffLags_Boomer(maxCorrInd);

% to get confidence intervals, do 500 bootstraps of trials
nBootstraps = 500;
nTrials = size(spikeData, 4);
for iShuff = 1:nBootstraps
    
    bootInds_Boomer{iShuff} = bootstrp(1, @(x) x, 1:nTrials);
    
    spikeDataShuff{iShuff} = spikeData(:,:,:,bootInds_Boomer{iShuff});
    kinDataShuff{iShuff} = kinData(:,:,:,bootInds_Boomer{iShuff});
    modelDataShuff{iShuff} = kinModelSpikeData(:,:,:,bootInds_Boomer{iShuff});
    
    % run dPCA and analysis on bootstrapped data
    [spikeDataAverageShuff_Boomer{iShuff}, WNeurShuff_Boomer{iShuff}, VNeurShuff_Boomer{iShuff}, whichMargNeurShuff_Boomer{iShuff},...
        explVarNeurShuff_Boomer{iShuff}, trajInvNeurShuff_Boomer{iShuff}, trajDepNeurShuff_Boomer{iShuff},...
        prinAnglesNeurShuff_Boomer{iShuff}, allAnglesNeurShuff_Boomer{iShuff}, canonCorrNeurShuff_Boomer{iShuff}] = ...
        dpcaAnalysis(spikeDataShuff{iShuff}, [5 5], size(spikeData,4), false);

    % kinematics too
    [kinDataAverageShuff_Boomer{iShuff}, WKinShuff_Boomer{iShuff}, VKinShuff_Boomer{iShuff}, whichMargKinShuff_Boomer{iShuff},...
        explVarKinShuff_Boomer{iShuff}, trajInvKinShuff_Boomer{iShuff}, trajDepKinShuff_Boomer{iShuff},...
        prinAnglesKinShuff_Boomer{iShuff}, allAnglesKinShuff_Boomer{iShuff}, canonCorrKKinShuff_Boomer{iShuff}] = ...
        dpcaAnalysis(kinDataShuff{iShuff}, [5 5], size(spikeData,4), false);

    % and the kinematic-modeled spike data
    [modelDataAverageShuff_Boomer{iShuff}, WModelShuff_Boomer{iShuff}, VModelShuff_Boomer{iShuff}, whichMargModelShuff_Boomer{iShuff},...
        explVarModelShuff_Boomer{iShuff}, trajInvModelShuff_Boomer{iShuff}, trajDepModelShuff_Boomer{iShuff},...
        prinAnglesModelShuff_Boomer{iShuff}, allAnglesModelShuff_Boomer{iShuff}, canonCorrKModelShuff_Boomer{iShuff}] = ...
        dpcaAnalysis(modelDataShuff{iShuff}, [5 5], size(spikeData,4), false);
    
    % get the rotational fit
    [anglesInvNeurShuff_Boomer{iShuff}, MskewInvNeurShuff_Boomer{iShuff}, skewRatioInvNeurShuff_Boomer(iShuff)] = ...
        getRotationMetric(trajInvNeurShuff_Boomer{iShuff});
    [anglesDepNeurShuff_Boomer{iShuff}, MskewDepNeurShuff_Boomer{iShuff}, skewRatioDepNeurShuff_Boomer(iShuff)] = ...
        getRotationMetric(trajDepNeurShuff_Boomer{iShuff});
    
    [anglesInvKinShuff_Boomer{iShuff}, MskewInvKinShuff_Boomer{iShuff}, skewRatioInvKinShuff_Boomer(iShuff)] = ...
        getRotationMetric(trajInvKinShuff_Boomer{iShuff});
    [anglesDepKinShuff_Boomer{iShuff}, MskewDepKinShuff_Boomer{iShuff}, skewRatioDepKinShuff_Boomer(iShuff)] = ...
        getRotationMetric(trajDepKinShuff_Boomer{iShuff});
    
end


% for pairs of bootstraps, calculate the principle angles between thier
% subspaces for the null distribution of similar subspaces
for iShuff = 1:nBootstraps/2
    
    %first inv subspaces
    invProj1 = VNeurShuff_Boomer{iShuff}(:,whichMargNeurShuff_Boomer{iShuff}==2);
    invProj2 = VNeurShuff_Boomer{iShuff+nBootstraps/2}(:,whichMargNeurShuff_Boomer{iShuff+nBootstraps/2}==2);
    
    %calc angle
    [~, S, ~] = svd(invProj1'*invProj2);
    tmp = diag(S);
    prinAnglesNull_Boomer(:,iShuff) = acosd(tmp);
    
    %also get dep subspaces
    depProj1 = VNeurShuff_Boomer{iShuff}(:,whichMargNeurShuff_Boomer{iShuff}==1);
    depProj2 = VNeurShuff_Boomer{iShuff+nBootstraps/2}(:,whichMargNeurShuff_Boomer{iShuff+nBootstraps/2}==1);
    
    %calc angle
    [~, S, ~] = svd(depProj1'*depProj2);
    tmp = diag(S);
    prinAnglesNull_Boomer(:,iShuff+nBootstraps/2) = acosd(tmp);
end

% permutation test for differences in R2 ratio, do shuffle of labels
nPermShuffs = 500;
for iShuff = 1:nPermShuffs
    
    mergedTrajsNeur = cat(4, trajInvNeur_Boomer, trajDepNeur_Boomer);
    mergedTrajsKin = cat(4, trajInvKin_Boomer, trajDepKin_Boomer);
    permIndShuff{iShuff} = randperm(size(mergedTrajsNeur,4));
    
    %first do neural components
    mergedTrajsNeur = mergedTrajsNeur(:,:,:,permIndShuff{iShuff});
    trajInvNeurPerm_Boomer = mergedTrajsNeur(:,:,:,1:size(mergedTrajsNeur,4)/2);
    trajDepNeurPerm_Boomer = mergedTrajsNeur(:,:,:,size(mergedTrajsNeur,4)/2+1:end);
    
    [~, ~, skewRatioInvNeurPerm] = getRotationMetric(trajInvNeurPerm_Boomer);
    [~, ~, skewRatioDepNeurPerm] = getRotationMetric(trajDepNeurPerm_Boomer);
    
    skewRatioDiffNeur_Boomer(iShuff) = skewRatioInvNeurPerm - skewRatioDepNeurPerm;
    
    %next to kinematics
    mergedTrajsKin = mergedTrajsKin(:,:,:,permIndShuff{iShuff});
    trajInvKinPerm_Boomer = mergedTrajsKin(:,:,:,1:size(mergedTrajsKin,4)/2);
    trajDepKinPerm_Boomer = mergedTrajsKin(:,:,:,size(mergedTrajsKin,4)/2+1:end);
    
    [~, ~, skewRatioInvKinPerm] = getRotationMetric(trajInvNeurPerm_Boomer);
    [~, ~, skewRatioDepKinPerm] = getRotationMetric(trajDepNeurPerm_Boomer);
    
    skewRatioDiffKin_Boomer(iShuff) = skewRatioInvKinPerm - skewRatioDepKinPerm;
    
end

% next do starbuck
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

load('./Data/ModelDataStarbuck.mat')
modelDataStarbuck = modeledData;

[spikeData, kinData] = splitSteps(trialsDataStarbuck, 'WalkingObstacle');
[kinModelSpikeData] = splitSteps(modelDataStarbuck, 'WalkingObstacle',0);

% neural
[spikeDataAverage_Starbuck, WNeur_Starbuck, VNeur_Starbuck, whichMargNeur_Starbuck, explVarNeur_Starbuck, trajInvNeur_Starbuck,...
    trajDepNeur_Starbuck, prinAnglesNeur_Starbuck, allAnglesNeur_Starbuck, canonCorrNeur_Starbuck] = dpcaAnalysis(spikeData, [5 5], size(spikeData,4), false);

% kinematics-modeled neural
[modelDataAverage_Starbuck, WModel_Starbuck, VModel_Starbuck, whichMargModel_Starbuck, explVarModel_Starbuck, trajInvModel_Starbuck,...
    trajDepModel_Starbuck, prinAnglesModel_Starbuck, allAnglesModel_Starbuck, canonCorrModel_Starbuck] = dpcaAnalysis(kinModelSpikeData, [5 5], size(spikeData,4), false);

% get no-dep noise level
[noiseDataAverage_Starbuck, WNoise_Starbuck, VNoise_Starbuck, whichMargNoise_Starbuck, explVarNoise_Starbuck, trajInvNoise_Starbuck,...
    trajDepNoise_Starbuck, prinAnglesNoise_Starbuck, allAnglesKin_Starbuck, canonCorrNoise_Starbuck] = dpcaAnalysis(spikeData(:,:,[1 6],:), [5 5], size(spikeData,4), false);

allNoDepNoise = squeeze(cat(4,trajDepNoise_Starbuck(:,:,1,:), trajDepNoise_Starbuck(:,:,2,:)));
upperNoDepNoise_Starbuck = mean(allNoDepNoise,3) + 2*std(allNoDepNoise,[],3);
lowerNoDepNoise_Starbuck = mean(allNoDepNoise,3) - 2*std(allNoDepNoise,[],3);

% get significance cutoff (99% and 1%)
for iDim = 1:5
    
    allDimNoise = allNoDepNoise(iDim,:,:);
    noiseCutoff(iDim) = prctile(allNoDepNoise(:),99.99);
    meanDeps = mean(trajDepNeur_Starbuck,4);
    h_Starbuck(iDim,:,:) = meanDeps(iDim,:,:)>noiseCutoff(iDim) | meanDeps(iDim,:,:)<-1*noiseCutoff(iDim);
    
    %get significant signals for noise analysis
    allSteps = zeros(1,600);
    allSigInds = {};
    allSigVals = {};
    for iStep = 1:6
        if any(h_Starbuck(iDim,:,iStep))
            
            sig = squeeze(meanDeps(iDim,:,iStep));
            sigInds = squeeze(h_Starbuck(iDim,:,iStep));
            allSigInds{iStep} = find(sigInds)+100*(iStep-1);
            allSigVals{iStep} = sig(sigInds);
            
        end
    end
%     xq = setdiff(1:600,[0 600 allSigInds{:}]);
%     yq = pchip([1:100 501:600 allSigInds{:}],[zeros(1,200) allSigVals{:}],xq);
%     allSteps([1:100 501:600 allSigInds{:}]) = [zeros(1,200) allSigVals{:}];
%     allsteps(xq) = yq;

    allSteps([allSigInds{:}]) = [allSigVals{:}];
    allSteps = reshape(allSteps,100,6);
    
    %get noise and signal components
    trajDepExtNoise_Starbuck(iDim,:,:,:) = squeeze(trajDepNeur_Starbuck(iDim,:,:,:)) - repmat(allSteps,1,1,size(trajDepNeur_Starbuck,4));
    trajDepExtSig_Starbuck(iDim,:,:,:) = repmat(allSteps,1,1,size(trajDepNeur_Starbuck,4));
    
end

% get inv noise and signal components
meanInvs = mean(mean(trajInvNeur_Starbuck,4),3);
trajInvExtSig_Starbuck = repmat(meanInvs,1,1,6,size(trajDepNeur_Starbuck,4));
trajInvExtNoise_Starbuck = trajInvNeur_Starbuck - trajInvExtSig_Starbuck;


% % do t-test
% for iDim = 1:5
%     for iPhase = 1:100
%         for iStep = 1:6
%             [hRaw_Starbuck(iDim,iPhase,iStep), p_Starbuck(iDim,iPhase,iStep)] = ttest2(squeeze(trajDepNeur_Starbuck(iDim,iPhase,iStep,:)),squeeze(allNoDepNoise(iDim,iPhase,:)));
%         end
%     end
% end
% 
% benHochCutoff = FDRcutoff(p_Starbuck(:),0.01,false);
% benHochCutoff = 0.05/3000;
% h_Starbuck = p_Starbuck<benHochCutoff;

% kinematics 
[kinDataAverage_Starbuck, WKin_Starbuck, VKin_Starbuck, whichMargKin_Starbuck, explVarKin_Starbuck, trajInvKin_Starbuck,...
    trajDepKin_Starbuck, prinAnglesKin_Starbuck, allAnglesKin_Starbuck, canonCorrKin_Starbuck] = dpcaAnalysis(kinData, [5 5], size(spikeData,4), false);

% get rotational fit
[anglesInvNeur_Starbuck, MskewInvNeur_Starbuck, skewRatioInvNeur_Starbuck] = getRotationMetric(trajInvNeur_Starbuck);
[anglesDepNeur_Starbuck, MskewDepNeur_Starbuck, skewRatioDepNeur_Starbuck] = getRotationMetric(trajDepNeur_Starbuck);

[anglesInvKin_Starbuck, MskewInvKin_Starbuck, skewRatioInvKin_Starbuck] = getRotationMetric(trajInvKin_Starbuck);
[anglesDepKin_Starbuck, MskewDepKin_Starbuck, skewRatioDepKin_Starbuck] = getRotationMetric(trajDepKin_Starbuck);

% compare timing with kinematics
depMahal_Starbuck = subFuncsTimeCourse.getMahal(squeeze(trajDepNeur_Starbuck(:,:,1,:)), squeeze(trajDepNeur_Starbuck(:,:,4,:)), 0);
invMahal_Starbuck = subFuncsTimeCourse.getMahal(squeeze(trajInvNeur_Starbuck(:,:,1,:)), squeeze(trajInvNeur_Starbuck(:,:,4,:)), 0);
kinMahal_Starbuck = subFuncsTimeCourse.getMahal(squeeze(kinData(:,:,1,:)), squeeze(kinData(:,:,4,:)), 0);
% find delay
[kinDepDiffxCorr_Starbuck, kinDepDiffLags_Starbuck] = crosscorr(mean(depMahal_Starbuck),mean(kinMahal_Starbuck));
[~, maxCorrInd] = max(kinDepDiffxCorr_Starbuck);
depDelay_Starbuck = kinDepDiffLags_Starbuck(maxCorrInd);

% bootstrap
nTrials = size(spikeData, 4);
for iShuff = 1:nBootstraps
    
    bootInds_Starbuck{iShuff} = bootstrp(1, @(x) x, 1:nTrials);
    
    spikeDataShuff{iShuff} = spikeData(:,:,:,bootInds_Starbuck{iShuff});
    kinDataShuff{iShuff} = kinData(:,:,:,bootInds_Starbuck{iShuff});
    modelDataShuff{iShuff} = kinModelSpikeData(:,:,:,bootInds_Starbuck{iShuff});
    
    % run dPCA and analysis on shuffled data
    [spikeDataAverageShuff_Starbuck{iShuff}, WNeurShuff_Starbuck{iShuff}, VNeurShuff_Starbuck{iShuff}, whichMargNeurShuff_Starbuck{iShuff},...
        explVarNeurShuff_Starbuck{iShuff}, trajInvNeurShuff_Starbuck{iShuff}, trajDepNeurShuff_Starbuck{iShuff},...
        prinAnglesNeurShuff_Starbuck{iShuff}, allAnglesNeurShuff_Starbuck{iShuff}, canonCorrNeurShuff_Starbuck{iShuff}] = ...
        dpcaAnalysis(spikeDataShuff{iShuff}, [5 5], size(spikeData,4), false);
    
    % run dPCA and analysis on kin-modeled data
    [modelDataAverageShuff_Starbuck{iShuff}, WModelShuff_Starbuck{iShuff}, VModelShuff_Starbuck{iShuff}, whichMargModelShuff_Starbuck{iShuff},...
        explVarModelShuff_Starbuck{iShuff}, trajInvModelShuff_Starbuck{iShuff}, trajDepModelShuff_Starbuck{iShuff},...
        prinAnglesModelShuff_Starbuck{iShuff}, allAnglesModelShuff_Starbuck{iShuff}, canonCorrModelShuff_Starbuck{iShuff}] = ...
        dpcaAnalysis(modelDataShuff{iShuff}, [5 5], size(spikeData,4), false);
    
    % kinematics too
    [kinDataAverageShuff_Starbuck{iShuff}, WKinShuff_Starbuck{iShuff}, VKinShuff_Starbuck{iShuff}, whichMargKinShuff_Starbuck{iShuff},...
        explVarKinShuff_Starbuck{iShuff}, trajInvKinShuff_Starbuck{iShuff}, trajDepKinShuff_Starbuck{iShuff},...
        prinAnglesKinShuff_Starbuck{iShuff}, allAnglesKinShuff_Starbuck{iShuff}, canonCorrKinShuff_Starbuck{iShuff}] = ...
        dpcaAnalysis(kinDataShuff{iShuff}, [5 5], size(spikeData,4), false);

    % get the rotational fit
    [anglesInvNeurShuff_Starbuck{iShuff}, MskewInvNeurShuff_Starbuck{iShuff}, skewRatioInvNeurShuff_Starbuck(iShuff)] = ...
        getRotationMetric(trajInvNeurShuff_Starbuck{iShuff});
    [anglesDepNeurShuff_Starbuck{iShuff}, MskewDepNeurShuff_Starbuck{iShuff}, skewRatioDepNeurShuff_Starbuck(iShuff)] = ...
        getRotationMetric(trajDepNeurShuff_Starbuck{iShuff});
    
    [anglesInvKinShuff_Starbuck{iShuff}, MskewInvKinShuff_Starbuck{iShuff}, skewRatioInvKinShuff_Starbuck(iShuff)] = ...
        getRotationMetric(trajInvKinShuff_Starbuck{iShuff});
    [anglesDepKinShuff_Starbuck{iShuff}, MskewDepKinShuff_Starbuck{iShuff}, skewRatioDepKinShuff_Starbuck(iShuff)] = ...
        getRotationMetric(trajDepKinShuff_Starbuck{iShuff});
    
end

% for pairs of bootstraps, calculate the principle angles between thier
% subspaces for the null distribution of similar subspaces
for iShuff = 1:nBootstraps/2
    
    %first inv subspaces
    invProj1 = VNeurShuff_Starbuck{iShuff}(:,whichMargNeurShuff_Starbuck{iShuff}==2);
    invProj2 = VNeurShuff_Starbuck{iShuff+nBootstraps/2}(:,whichMargNeurShuff_Starbuck{iShuff+nBootstraps/2}==2);
    
    %calc angle
    [~, S, ~] = svd(invProj1'*invProj2);
    tmp = diag(S);
    prinAnglesNull_Starbuck(:,iShuff) = acosd(tmp);
    
    %also get dep subspaces
    depProj1 = VNeurShuff_Starbuck{iShuff}(:,whichMargNeurShuff_Starbuck{iShuff}==1);
    depProj2 = VNeurShuff_Starbuck{iShuff+nBootstraps/2}(:,whichMargNeurShuff_Starbuck{iShuff+nBootstraps/2}==1);
    
    %calc angle
    [~, S, ~] = svd(depProj1'*depProj2);
    tmp = diag(S);
    prinAnglesNull_Starbuck(:,iShuff+nBootstraps/2) = acosd(tmp);
end

% permutation test for differences in R2 ratio, do shuffle of labels
for iShuff = 1:nPermShuffs
    
    mergedTrajsNeur = cat(4, trajInvNeur_Starbuck, trajDepNeur_Starbuck);
    mergedTrajsKin = cat(4, trajInvKin_Starbuck, trajDepKin_Starbuck);
    permIndShuff{iShuff} = randperm(size(mergedTrajsNeur,4));
    
    %first do neural components
    mergedTrajsNeur = mergedTrajsNeur(:,:,:,permIndShuff{iShuff});
    trajInvNeurPerm_Starbuck = mergedTrajsNeur(:,:,:,1:size(mergedTrajsNeur,4)/2);
    trajDepNeurPerm_Starbuck = mergedTrajsNeur(:,:,:,size(mergedTrajsNeur,4)/2+1:end);
    
    [~, ~, skewRatioInvNeurPerm] = getRotationMetric(trajInvNeurPerm_Starbuck);
    [~, ~, skewRatioDepNeurPerm] = getRotationMetric(trajDepNeurPerm_Starbuck);
    
    skewRatioDiffNeur_Starbuck(iShuff) = skewRatioInvNeurPerm - skewRatioDepNeurPerm;
    
    %next to kinematics
    mergedTrajsKin = mergedTrajsKin(:,:,:,permIndShuff{iShuff});
    trajInvKinPerm_Starbuck = mergedTrajsKin(:,:,:,1:size(mergedTrajsKin,4)/2);
    trajDepKinPerm_Starbuck = mergedTrajsKin(:,:,:,size(mergedTrajsKin,4)/2+1:end);
    
    [~, ~, skewRatioInvKinPerm] = getRotationMetric(trajInvNeurPerm_Starbuck);
    [~, ~, skewRatioDepKinPerm] = getRotationMetric(trajDepNeurPerm_Starbuck);
    
    skewRatioDiffKin_Starbuck(iShuff) = skewRatioInvKinPerm - skewRatioDepKinPerm;
    
end



%% Subplot A -- dPCA components for Boomer
% first plot step invarient trajectories

invInds = find(whichMargNeur_Boomer==2);
colorMap = parula(8);

% dim 1
dim1BoomerInvH = axes('Units','inches','OuterPosition',...
    [figMargins-1.2, 7.3 + subplotGap*5 + figMargins, 6, 3]);

trajInvMeans_Boomer = mean(trajInvNeur_Boomer,4);
for iStep = 1:size(trajInvMeans_Boomer,3)
    
    legendH(iStep) = plot(trajInvMeans_Boomer(1,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
end

line([0 0], [-0.9 -0.8], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Boomer.componentVar(invInds(1))*10)/10;
set(dim1BoomerInvH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% legend
stepLabels = join([repmat("Step",6,1),string(num2cell(-3:2))']);
legend(legendH, stepLabels, 'box', 'off', 'location', 'best', 'NumColumns', 1, 'FontSize', 10, 'box', 'off')

% next, dim2
dim2BoomerInvH = axes('Units','inches','OuterPosition',...
    [figMargins+2.1, 7.3 + subplotGap*5 + figMargins, 6, 3]);

for iStep = 1:size(trajInvMeans_Boomer,3)
    
    plot(trajInvMeans_Boomer(2,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
end

line([0 0], [-0.5 -0.4], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Boomer.componentVar(invInds(2))*10)/10;
set(dim2BoomerInvH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% next, do step dependent components
depInds = find(whichMargNeur_Boomer==1);

% dim 1
dim1BoomerDepH = axes('Units','inches','OuterPosition',...
    [figMargins-1.2, 5 + subplotGap*5 + figMargins, 6, 3]);

trajDepMeans_Boomer = mean(trajDepNeur_Boomer,4);
for iStep = 1:size(trajDepMeans_Boomer,3)
    
    plot(trajDepMeans_Boomer(1,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
    plot(find(h_Boomer(1,:,iStep)),ones(1,sum(h_Boomer(1,:,iStep)))*0.45+0.04*iStep,'.','color',colorMap(iStep+1,:));
end

line([2 2], [-0.3 -0.2], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Boomer.componentVar(depInds(1))*10)/10;
set(dim1BoomerDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% next, dim2
dim2BoomerDepH = axes('Units','inches','OuterPosition',...
    [figMargins+2.1, 5 + subplotGap*5 + figMargins, 6, 3]);

for iStep = 1:size(trajDepMeans_Boomer,3)
    
    plot(trajDepMeans_Boomer(2,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
    plot(find(h_Boomer(2,:,iStep)),ones(1,sum(h_Boomer(2,:,iStep)))*0.20+0.01*iStep,'.','color',colorMap(iStep+1,:));
end

line([2 2], [-0.2 -0.1], 'linewidth', 3, 'color', 'k')


box off
varPerct = round(explVarNeur_Boomer.componentVar(depInds(2))*10)/10;
set(dim2BoomerDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100], 'ylim', [-0.5 0.5])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)


%% Subplot B -- Starbuck dPCA components
% first plot step invarient trajectories

invInds = find(whichMargNeur_Starbuck==2);

% dim 1
dim1StarbuckInvH = axes('Units','inches','OuterPosition',...
    [figMargins+5.7, 7.3 + subplotGap*5 + figMargins, 6, 3]);

trajInvMeans_Starbuck = mean(trajInvNeur_Starbuck,4);
for iStep = 1:size(trajInvMeans_Starbuck,3)
    
    plot(trajInvMeans_Starbuck(1,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
end

line([0 0], [-0.9 -0.8], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Starbuck.componentVar(invInds(1))*10)/10;
set(dim1StarbuckInvH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% next, dim2
dim2StarbuckInvH = axes('Units','inches','OuterPosition',...
    [figMargins+9, 7.3 + subplotGap*5 + figMargins, 6, 3]);

for iStep = 1:size(trajInvMeans_Starbuck,3)
    
    plot(trajInvMeans_Starbuck(2,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
end

line([0 0], [-0.5 -0.4], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Starbuck.componentVar(invInds(2))*10)/10;
set(dim2StarbuckInvH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% next, do step dependent components
depInds = find(whichMargNeur_Starbuck==1);

% dim 1
dim1StarbuckDepH = axes('Units','inches','OuterPosition',...
    [figMargins+5.7, 5 + subplotGap*5 + figMargins, 6, 3]);

trajDepMeans_Starbuck = mean(trajDepNeur_Starbuck,4);
for iStep = 1:size(trajDepMeans_Starbuck,3)
    
    plot(trajDepMeans_Starbuck(1,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
    plot(find(h_Starbuck(1,:,iStep)),ones(1,sum(h_Starbuck(1,:,iStep)))*0.26+0.04*iStep,'.','color',colorMap(iStep+1,:));
end

line([2 2], [-0.3 -0.2], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Starbuck.componentVar(depInds(1))*10)/10;
set(dim1StarbuckDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)

% next, dim2
dim2StarbuckDepH = axes('Units','inches','OuterPosition',...
    [figMargins+9, 5 + subplotGap*5 + figMargins, 6, 3]);

for iStep = 1:size(trajDepMeans_Starbuck,3)
    
    plot(trajDepMeans_Starbuck(2,:,iStep), 'color', colorMap(iStep+1,:), 'LineWidth', 2);
    hold on
    
    plot(find(h_Starbuck(2,:,iStep)),ones(1,sum(h_Starbuck(2,:,iStep)))*0.1+0.02*iStep,'.','color',colorMap(iStep+1,:));
end

line([2 2], [-0.2 -0.1], 'linewidth', 3, 'color', 'k')

box off
varPerct = round(explVarNeur_Starbuck.componentVar(depInds(2))*10)/10;
set(dim2StarbuckDepH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [0 50 100], 'YLim', [-0.3 0.46])
title(['Var = ' num2str(varPerct) '%'], 'FontSize', 12)


%% Subplot C -- Cross-corr timing between kinematic changes and Stride-dependent mode changes

% plot Boomer
depKinTimingBoomerH = axes('Units','inches','OuterPosition',...
    [figMargins-1, figMargins+2.7, 5.2, 3.5]);
plot(kinDepDiffLags_Boomer, kinDepDiffxCorr_Boomer, 'k', 'linewidth', 2)
hold on
line(repmat(depDelay_Boomer,1,2), get(gca,'ylim'), 'linewidth', 2, 'color','r','linestyle','--')

box off
set(depKinTimingBoomerH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [-20 0 depDelay_Boomer, 20])

% then Starbuck
depKinTimingStarbuckH = axes('Units','inches','OuterPosition',...
    [figMargins+5.8, figMargins+2.7, 5.2, 3.5]);
plot(kinDepDiffLags_Starbuck, kinDepDiffxCorr_Starbuck, 'k', 'linewidth', 2)
hold on
line(repmat(depDelay_Starbuck,1,2), get(gca,'ylim'), 'linewidth', 2, 'color','r','linestyle','--')

box off
set(depKinTimingStarbuckH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'YColor', 'none', 'XTick', [-20 0 depDelay_Starbuck, 20])


%% Subplot D -- Principle angles between step invarient vs step dependent subspaces

% first boomer
prinAngleBoomerH = axes('Units','inches','OuterPosition',...
    [figMargins+1.8, figMargins+2.7, 5.8, 3.5]);
hold on;

% get confidence interval from bootstrap
allPrinAnglesNeurBoomer = cat(2,prinAnglesNeurShuff_Boomer{:});
upperLim = prctile(allPrinAnglesNeurBoomer', 97.5);
lowerLim = prctile(allPrinAnglesNeurBoomer', 2.5);

% plot neural
prinAngleH(1) = plot(prinAnglesNeur_Boomer, 'color', [0 0.4470 0.7410], 'linewidth', 2);
shadedErrorBar(1:5, prinAnglesNeur_Boomer, [upperLim - prinAnglesNeur_Boomer'; prinAnglesNeur_Boomer' - lowerLim],...
    {'color', [0 0.4470 0.7410]}, true)

% get confidence interval for kinematics
allPrinAnglesKinBoomer = cat(2,prinAnglesKinShuff_Boomer{:});
upperLim = prctile(allPrinAnglesKinBoomer', 97.5);
lowerLim = prctile(allPrinAnglesKinBoomer', 2.5);

% plot kinematics
prinAngleH(2) = plot(prinAnglesKin_Boomer, 'color', [0.8500, 0.3250, 0.0980], 'linewidth', 2);
shadedErrorBar(1:5, prinAnglesKin_Boomer, [upperLim - prinAnglesKin_Boomer'; prinAnglesKin_Boomer' - lowerLim],...
    {'color', [0.8500, 0.3250, 0.0980]}, true)

% and kinematics-modeled neural
allPrinAnglesModelBoomer = cat(2,prinAnglesModelShuff_Boomer{:});
upperLim = prctile(allPrinAnglesModelBoomer', 97.5);
lowerLim = prctile(allPrinAnglesModelBoomer', 2.5);

% plot kinematic-modeled spikes
prinAngleH(3) = plot(prinAnglesModel_Boomer, 'color', [0.5, 0.5, 0.5], 'linewidth', 2);
shadedErrorBar(1:5, prinAnglesModel_Boomer, [upperLim - prinAnglesModel_Boomer'; prinAnglesModel_Boomer' - lowerLim],...
    {'color', [0.5, 0.5, 0.5]}, true)

% plot null distirbution 95% percentile
plot(prctile(prinAnglesNull_Boomer',97.5),'--', 'linewidth', 2, 'color', [0.8 0.8 0.8])

box off
ylabel('Degrees')
xlabel('Principle Angle')

set(prinAngleBoomerH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', [0.5 5.5])

legend(prinAngleH, 'Neural', 'Kinematics', 'Kin-model', 'box', 'off')

% now do starbuck
prinAngleStarbuckH = axes('Units','inches','OuterPosition',...
    [figMargins+8.5, figMargins+2.7, 5.8, 3.5]);
hold on;

% get confidence interval from bootstrap
allPrinAnglesNeurStarbuck = cat(2,prinAnglesNeurShuff_Starbuck{:});
upperLim = prctile(allPrinAnglesNeurStarbuck', 97.5);
lowerLim = prctile(allPrinAnglesNeurStarbuck', 2.5);

% plot neural
plot(prinAnglesNeur_Starbuck, 'color', [0 0.4470 0.7410], 'linewidth', 2)
shadedErrorBar(1:5, prinAnglesNeur_Starbuck, [upperLim - prinAnglesNeur_Starbuck'; prinAnglesNeur_Starbuck' - lowerLim],...
    {'color', [0 0.4470 0.7410]}, true)

% get confidence interval for kinematics
allPrinAnglesKiinStarbuck = cat(2,prinAnglesKinShuff_Starbuck{:});
upperLim = prctile(allPrinAnglesKiinStarbuck', 97.5);
lowerLim = prctile(allPrinAnglesKiinStarbuck', 2.5);

% plot kinematics
plot(prinAnglesKin_Starbuck, 'color', [0.8500, 0.3250, 0.0980], 'linewidth', 2)
shadedErrorBar(1:5, prinAnglesKin_Starbuck, [upperLim - prinAnglesKin_Starbuck'; prinAnglesKin_Starbuck' - lowerLim],...
    {'color', [0.8500, 0.3250, 0.0980]}, true)

% and kinematics-modeled neural
allPrinAnglesModelStarbuck = cat(2,prinAnglesModelShuff_Starbuck{:});
upperLim = prctile(allPrinAnglesModelStarbuck', 97.5);
lowerLim = prctile(allPrinAnglesModelStarbuck', 2.5);

% plot kinematic-modeled spikes
prinAngleH(3) = plot(prinAnglesModel_Starbuck, 'color', [0.5, 0.5, 0.5], 'linewidth', 2);
shadedErrorBar(1:5, prinAnglesModel_Starbuck, [upperLim - prinAnglesModel_Starbuck'; prinAnglesModel_Starbuck' - lowerLim],...
    {'color', [0.5, 0.5, 0.5]}, true)

% plot null distirbution 95% percentile
plot(prctile(prinAnglesNull_Starbuck',97.5),'--', 'linewidth', 2, 'color', [0.8 0.8 0.8])

box off
ylabel('Degrees')
xlabel('Principle Angle')

set(prinAngleStarbuckH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XLim', [0.5 5.5])

%% Subplot E -- Rotational strength of components

% first boomer
invSlopeFieldBoomerH = axes('Units','inches','OuterPosition',...
    [figMargins+0, figMargins+0.2, 2.2, 2.2]);

plotTrajVecField(trajInvNeur_Boomer, MskewInvNeur_Boomer, colorMap)
box off
axis off

depSlopeFieldBoomerH = axes('Units','inches','OuterPosition',...
    [figMargins+2.2, figMargins+0.2, 2.2, 2.2]);

plotTrajVecField(trajDepNeur_Boomer, MskewDepNeur_Boomer, colorMap)
box off
axis off

skewFitBoomerH = axes('Units','inches','OuterPosition',...
    [figMargins+4.5, figMargins+0.1, 2.4, 2.2]);

ratios = [skewRatioInvNeur_Boomer skewRatioDepNeur_Boomer  skewRatioInvKin_Boomer skewRatioDepKin_Boomer];
boomerBar = bar([1 2 4 5], ratios, 'EdgeAlpha', 0, 'FaceColor', 'flat');
hold on;

% colors
for iBar = 1:size(boomerBar.CData, 1)
    if iBar < 3
        boomerBar.CData(iBar,:) = [0 0.4470 0.7410];
    else
        boomerBar.CData(iBar,:) = [0.8500, 0.3250, 0.0980];
    end
end

% get confidence interval from bootstrap
topErrBar = [prctile(skewRatioInvNeurShuff_Boomer, 97.5), ...
    prctile(skewRatioDepNeurShuff_Boomer, 97.5), prctile(skewRatioInvKinShuff_Boomer, 97.5), prctile(skewRatioDepKinShuff_Boomer, 97.5)];
botErrBar = [prctile(skewRatioInvNeurShuff_Boomer, 2.5), ...
    prctile(skewRatioDepNeurShuff_Boomer, 2.5), prctile(skewRatioInvKinShuff_Boomer, 2.5), prctile(skewRatioDepKinShuff_Boomer, 2.5)];

errorbar([1 2 4 5], ratios, topErrBar-ratios, botErrBar-ratios, '.k', 'LineWidth', 2)

box off
ylabel('R^2_S_k_e_w / R^2_F_u_l_l')

set(skewFitBoomerH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'xTickLabel', {'Inv', 'Dep', 'Inv', 'Dep'}, 'ylim', [0 1])

% add significance stars for t-test
% get pvalue for neural inv-dep
observedDiff = skewRatioInvNeur_Boomer - skewRatioDepNeur_Boomer;
p_neurBoomer = (sum(abs(skewRatioDiffNeur_Boomer) > abs(observedDiff))+1) / (nPermShuffs+1);
% p_neurBoomer = ranksum(skewRatioInvNeurShuff_Boomer, skewRatioDepNeurShuff_Boomer);

if p_neurBoomer < 0.05
    text(1.3, 0.95, '*','FontSize',20)
    line([1 2], [0.93 0.93],'color','k','linewidth',2)
end
    
% finally, get p-value for kin inv-dep
observedDiff = skewRatioInvKin_Boomer - skewRatioDepKin_Boomer;
p_kinBoomer = (sum(abs(skewRatioDiffKin_Boomer) > abs(observedDiff))+1) / (nPermShuffs+1);
% p_kinBoomer = ranksum(skewRatioInvKinShuff_Boomer, skewRatioDepKinShuff_Boomer);

if p_kinBoomer < 0.05
    text(4.3, 0.55, '*','FontSize',20)
    line([4 5], [0.53 0.53],'color','k','linewidth',2)
end
    
% next, do starbuck
invSlopeFieldStarbuckH = axes('Units','inches','OuterPosition',...
    [figMargins+7.1, figMargins+0.2, 2.2, 2.2]);

plotTrajVecField(trajInvNeur_Starbuck, MskewInvNeur_Starbuck, colorMap)
box off
axis off

depSlopeFieldStarbuckH = axes('Units','inches','OuterPosition',...
    [figMargins+9.3, figMargins+0.2, 2.2, 2.2]);

plotTrajVecField(trajDepNeur_Starbuck, MskewDepNeur_Starbuck, colorMap)
box off
axis off

skewFitStarbuckH = axes('Units','inches','OuterPosition',...
    [figMargins+11.6, figMargins+0.1, 2.4, 2.2]);

ratios = [skewRatioInvNeur_Starbuck skewRatioDepNeur_Starbuck skewRatioInvKin_Starbuck skewRatioDepKin_Starbuck];
starbuckBar = bar([1 2 4 5], ratios, 'EdgeAlpha', 0, 'FaceColor', 'flat');
hold on;

% colors
for iBar = 1:size(starbuckBar.CData, 1)
    if iBar < 3
        starbuckBar.CData(iBar,:) = [0 0.4470 0.7410];
    else
        starbuckBar.CData(iBar,:) = [0.8500, 0.3250, 0.0980];
    end
end

% get confidence interval from bootstrap
topErrBar = [prctile(skewRatioInvNeurShuff_Starbuck, 97.5), ...
    prctile(skewRatioDepNeurShuff_Starbuck, 97.5), prctile(skewRatioInvKinShuff_Starbuck, 97.5), prctile(skewRatioDepKinShuff_Starbuck, 97.5)];
botErrBar = [prctile(skewRatioInvNeurShuff_Starbuck, 2.5), ...
    prctile(skewRatioDepNeurShuff_Starbuck, 2.5), prctile(skewRatioInvKinShuff_Starbuck, 2.5), prctile(skewRatioDepKinShuff_Starbuck, 2.5)];

errorbar([1 2 4 5], ratios, topErrBar-ratios, botErrBar-ratios, '.k', 'LineWidth', 2)

% add significance stars from bootstrap
% get p-value for neural inv-dep
observedDiff = skewRatioInvNeur_Starbuck - skewRatioDepNeur_Starbuck;
p_neurStarbuck = (sum(abs(skewRatioDiffNeur_Starbuck) > abs(observedDiff))+1) / (nPermShuffs+1);
% p_neurStarbuck = ranksum(skewRatioInvNeurShuff_Starbuck, skewRatioDepNeurShuff_Starbuck);

if p_neurStarbuck < 0.05
    text(1.3, 1.02, '*','FontSize',20)
    line([1 2], [1 1],'color','k','linewidth',2)
end
    
observedDiff = skewRatioInvKin_Starbuck - skewRatioDepKin_Starbuck;
p_kinStarbuck = (sum(abs(skewRatioDiffKin_Starbuck) > abs(observedDiff))+1) / (nPermShuffs+1);
% p_kinStarbuck = ranksum(skewRatioInvKinShuff_Starbuck, skewRatioDepKinShuff_Starbuck);

if p_kinStarbuck < 0.05
    text(4.3, 0.45, '*','FontSize',20)
    line([4 5], [0.43 0.43],'color','k','linewidth',2)
end
    
box off
ylabel('R^2_S_k_e_w / R^2_F_u_l_l')

set(skewFitStarbuckH, 'FontSize',12, 'TickDir','out','TickLength', [0.02 0.02],...
    'LineWidth', 2, 'XTickLabelRotation', 45, 'xTickLabel', {'Inv', 'Dep', 'Inv', 'Dep'}, 'ylim', [0 1])



function [spikeData, kinData, spikeDataUnsmoothed] = splitSteps(trialsData, task, getKin, badTrials)

if nargin == 2
    getKin = 1;
end

% get trials
if nargin == 2
    badTrials = filterTrials(trialsData,90,5);
else
    badTrials = [];
end

walkObsTrialInds = find(cellfun(@(x) strcmpi(x, task), {trialsData.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
    
% go through each trial, get normalized data
for iTrial = 1:length(walkObsTrialInds)
    
    %get spike counts
    preTrialCounts = trialsData(walkObsTrialInds(iTrial)).GaitNormalizedPreTrialSpikeCounts;
    counts = trialsData(walkObsTrialInds(iTrial)).GaitNormalizedSpikeCounts;
    postTrialCounts = trialsData(walkObsTrialInds(iTrial)).GaitNormalizedPostTrialSpikeCounts;
    
    %smooth with gaussian
    smoothedCounts = convGauss([preTrialCounts counts postTrialCounts], 10, 40);
    smoothedCounts = smoothedCounts(:, size(preTrialCounts,2)+1 : end-size(postTrialCounts,2));
        
    nSteps = size(counts,2)/100;

    %concatenate all joints
    if getKin
        kin = struct2cell(trialsData(walkObsTrialInds(iTrial)).GaitNormalizedKinematics);
        kin = cat(1,kin{2:end});
        kin = kin(:,1:100*nSteps);
    end
    
    %divide into steps
    for iStep = 1:nSteps
        
        spikeData(:, :, iStep, iTrial) = smoothedCounts(:, 1 + 100*(iStep-1) : 100*(iStep));
        spikeDataUnsmoothed(:, :, iStep, iTrial) = counts(:, 1 + 100*(iStep-1) : 100*(iStep));
        if getKin
            kinData(:, :, iStep, iTrial) = kin(:, 1 + 100*(iStep-1) : 100*(iStep));
        end
        
    end
   
end

% spikeData = squeeze(spikeData);
% spikeDataUnsmoothed = squeeze(spikeDataUnsmoothed);
% kinData = squeeze(kinData);


function [trajInv, trajDep, prinAngles, allAngles, canonR] = dPCAProj(inputData, W, V, marges, useSurr)    
    
% first, project
compMeans = mean(inputData(:,:),2);
inputDataCen = inputData - repmat(compMeans,1,1,size(inputData,3));

if useSurr

    % for each step
        for iStep = 1:size(inputData,2)
            
            % get the step invarient part
            invInds = find(marges == 2);
            invW = W(:, invInds);
            
            trajInv(:,:,iStep) = invW' * squeeze(inputDataCen(:,iStep,:));
            
            %step depedent part
            depInds = find(marges == 1);
            depW = W(:, depInds);
            
            trajDep(:,:,iStep) = depW' * squeeze(inputDataCen(:,iStep,:));
            
        end
        
else
    
    % split into trials
    nTrials = size(inputData,4);
    for iTrial = 1:nTrials
        
        % for each step
        for iStep = 1:size(inputData,2)
            
            % get the step invarient part
            invInds = find(marges == 2);
            invW = W(:, invInds);
            
            trajInv(:,:,iStep,iTrial) = invW' * squeeze(inputDataCen(:, iStep, 1+100*(iTrial-1):100*iTrial));
            
            %step depedent part
            depInds = find(marges == 1);
            depW = W(:, depInds);
            
            trajDep(:,:,iStep,iTrial) = depW' * squeeze(inputDataCen(:,iStep, 1+100*(iTrial-1):100*iTrial));
            
        end
        
    end
    
end

% get principle angle between subspaces
[~, S, ~] = svd(V(:,invInds)'*V(:,depInds));
tmp = diag(S);
prinAngles = acosd(tmp);

% alternatively, get angles between each of the dimensions
for iInv = 1:length(invInds)
    for iDep = 1:length(depInds)
        
        allAngles(iInv,iDep) = acosd(V(:,invInds(iInv))' * V(:,depInds(iDep)));
        
    end
end
allAngles = allAngles(:,:);

% or, use the actual data and look at cannonical correlation analysis
[~,~,canonR] = canoncorr(trajInv(:,:)',trajDep(:,:)');

% permute axes
% for i = 1:500
%     permuteV = V(:,randperm(size(V,2)));
%     [~, permuteS, ~] = svd(permuteV(:,1:length(invInds))'*V(:,length(invInds)+1:end));
%     canonR = acosd(diag(permuteS));
% end

function [angles, Mskew, R2Ratio] = getRotationMetric(trajs)
% rotation metric, compare the angle of the position vs velocity of the
% top two average components

% get averages
meanTrajs = mean(trajs,4);
meanTrajs = meanTrajs(:,:);

for iTime = 1:size(meanTrajs,2)-1
    
    pos = meanTrajs(1:2,iTime);
    vel = diff(meanTrajs(1:2,iTime:iTime+1),[],2);
    angles(iTime) = acosd(pos'*vel/(norm(pos)*norm(vel)));
    
end

% fit jPCA
% get all the trajs concatenated across trials and steps
neuronMeans = mean(trajs(:,:),2);

tmp = trajs(:,1:99,:,:);
catTrajs = tmp(:,:);
% also mean center it
% neuronMeans = mean(catTrajs,2);
catTrajs = catTrajs - repmat(neuronMeans,1,size(catTrajs,2));

% get change in trajs
dTrajs = diff(trajs-repmat(neuronMeans,1,size(trajs,2),size(trajs,3),size(trajs,4)),[],2);
catdTrajs = dTrajs(:,:);

Mskew = skewSymRegress(catdTrajs',catTrajs')';
Mfull = catdTrajs/catTrajs;

% calculate R2s
catdTrajsVar = sum(reshape(catdTrajs-repmat(mean(catdTrajs')',1,size(catdTrajs,2)),1,numel(catdTrajs)).^2);
MskewVar = sum(reshape(catdTrajs - Mskew*catTrajs,1,numel(catdTrajs)).^2);
MfullVar = sum(reshape(catdTrajs - Mfull*catTrajs,1,numel(catdTrajs)).^2);
R2Skew = (catdTrajsVar - MskewVar)/catdTrajsVar;
R2Full = (catdTrajsVar - MfullVar)/catdTrajsVar;
R2Ratio = R2Skew/R2Full;


function [spikeDataAverage, W, V, whichMarg, explVar, trajInv, trajDep, prinAngles, allAngles, canonR] = dpcaAnalysis(spikeData, nDims, nTrials, useSurr)

% avrage across trials
if (useSurr)
    
    spikeDataAverage = spikeData;
    % get Cnoise
    neurCnoise = dpca_getNoiseCovariance(spikeDataAverage, ...
        repmat(spikeData, 1, 1, 1, nTrials), size(spikeData,4), 'simultaneous', false, 'type', 'averaged');
    
else
    spikeDataAverage = permute(mean(spikeData,4), [1 3 2]);
    
    % get Cnoise
    neurCnoise = dpca_getNoiseCovariance(spikeDataAverage, ...
        permute(spikeData, [1 3 2 4]), size(spikeData,4), 'simultaneous', true);
end

% run dPCA
% marginalizations will be phase and phase-step
combinedParams = {{1, [1,2]}, {[2]}};
[W, V, whichMarg] = dpca(spikeDataAverage, nDims,'combinedParams', combinedParams, 'Cnoise', neurCnoise);

% get explained variance
explVar = dpca_explainedVariance(spikeDataAverage, W, V, 'combinedParams', combinedParams);

% get projections and principle angles
if useSurr
    [trajInv, trajDep, prinAngles, allAngles, canonR] = dPCAProj(spikeData, W, V, whichMarg, useSurr);
else
    [trajInv, trajDep, prinAngles, allAngles, canonR] = dPCAProj(permute(spikeData, [1 3 2 4]), W, V, whichMarg, useSurr);
end


function plotTrajVecField(traj, Mskew, colorMap)
% make verification plots
meanTrajs = mean(traj,4);

hold on;
for iStep = 1:size(meanTrajs,3)
    plot(meanTrajs(1,:,iStep), meanTrajs(2,:,iStep), 'color', colorMap(iStep,:), 'LineWidth', 1.5);
end

meanTrajs = meanTrajs(:,:);
xRange = get(gca, 'xlim');
yRange = get(gca, 'ylim');

[vecfieldX vecfieldY] = meshgrid(xRange(1)*2:(diff(xRange)/10):xRange(2)*2, yRange(1)*2:(diff(yRange)/10):yRange(2)*2);
pos = cat(3,vecfieldX,vecfieldY);

for iX = 1:size(pos,1)
    for iY = 1:size(pos,2)

        vel = Mskew(1:2,1:2)*squeeze(pos(iX,iY,:));
        vecfieldU(iX, iY) = vel(1);
        vecfieldV(iX, iY) = vel(2);
        
    end
end

hold on
quiver(vecfieldX, vecfieldY, vecfieldU, vecfieldV, 1.5, 'color', [0.8500, 0.3250, 0.0980])
set(gca,'XLim',xRange,'YLim',yRange)



function prinAngles = randSubspaceAngles(nPerms, nMainDims, nSubDims)

for iPerm = 1:nPerms

    basis = {};
    for iSubspace = 1:2
        %random inital vector
        basis{iSubspace} = rand(nMainDims,1);
        basis{iSubspace} = basis{iSubspace}/norm(basis{iSubspace});
        
        for iSubDim = 1:nSubDims-1
            
            %random vector again
            newInit = rand(nMainDims,1);
            
            %project on to all previous basis vectors
            newInitOrth = newInit - sum(repmat(basis{iSubspace}' * newInit, 1, nMainDims)' .* basis{iSubspace},2);
            basis{iSubspace}(:,end+1) = newInitOrth/norm(newInitOrth);
            
        end
    end
    
    %get principle angles
    [~, S, ~] = svd(basis{1}'*basis{2});
    Svalues = diag(S);
    prinAngles(iPerm,:) = acosd(Svalues);
    
end


% 
