% script to add simulated gait-noramlized spike counts (pre and post) and
% unnormalized spike counts to the Boomer and Starbuck leg arrays from
% kineamtics to use as control for dPCA and PLDS analysis

% neural simluation based off of Gallego et al 2018 and Perich et al 2018:
% y = a+sum(b*k)+u
% y is the mean poisson rate
% a is randomly sampled baseline firing rate between 0 and 0.1
% b are the weights sampled from standard normal distribution
% k is the kinematics
% u is randomly noise (zero mean, 0.05 std)

% first do Boomer
clear

% load in data
load('./Data/TrialsDataBoomer.mat')
trialsDataBoomer = trialsLegM1;

nNeurons = size(trialsDataBoomer(1).SpikeCounts,1);
nKinematics = 15; %joint angles
dutyPercent = 67;
nSteps = 7;

% sample weights
kinWeights = randn(nNeurons,nKinematics);

% sample baseline
baselines = rand(nNeurons,1)*0.1;

% get mean firing rates across all walking and walking obstacle trials
badTrials = filterTrials(trialsDataBoomer,90,5);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x, 'WalkingObstacle'), {trialsDataBoomer.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
walkTrialInds = find(cellfun(@(x) strcmpi(x, 'Walk'), {trialsDataBoomer.Task}));
walkTrialInds = setdiff(walkTrialInds, badTrials);
allTrialSpikeCounts = {trialsDataBoomer([walkObsTrialInds walkTrialInds]).SpikeCounts};
allTrialSpikeCounts = cat(2,allTrialSpikeCounts{:});
meanSpikeCounts = mean(allTrialSpikeCounts,2);

stdSpikeCounts = std(allTrialSpikeCounts-meanSpikeCounts,[],2);

% also get all kinematics to make them positive and scaled down
allTrialKin = cellfun(@(x) struct2cell(x),{trialsDataBoomer([walkObsTrialInds walkTrialInds]).Kinematics}, 'un', 0);
% allTrialKin = cellfun(@(x) getJointAngle(x), allTrialKin,'un',0);
allTrialKin = cellfun(@(x) cat(1,x{:}), allTrialKin,'un',0);
allTrialKin = cat(2,allTrialKin{:});
allTrialKin(1:3,:) = [];
allTrialKinMeans = mean(allTrialKin,2);
allTrialKinAmps = max(allTrialKin,[],2) - min(allTrialKin,[],2);
alltrialLambdas = kinWeights*((allTrialKin - allTrialKinMeans)./repmat(allTrialKinAmps,1,size(allTrialKin,2)));
lambdaMins = min(alltrialLambdas');

% go through each trial and simulate
iModelTrial = 1;
for iTrial = 1:length(trialsDataBoomer)
    
    %only do walk and walk obstacle
    if ~any(iTrial==walkTrialInds) & ~any(iTrial==walkObsTrialInds)
        continue
    end
    
    %get trial lengths
    preTrialLength = size(trialsDataBoomer(iTrial).PreTrialSpikeCounts,2);
    trialLength = size(trialsDataBoomer(iTrial).SpikeCounts,2);
    postTrialLength = size(trialsDataBoomer(iTrial).PostTrialSpikeCounts,2);
    
    %events for normalization
    if length(trialsDataBoomer(iTrial).TrialEvents)==2
        events = [trialsDataBoomer(iTrial).TrialEvents size(trialsDataBoomer(iTrial).SpikeCounts,2)+1];
        gaitPercentages = [0 dutyPercent 100];
    else
        events = [trialsDataBoomer(iTrial).TrialEvents size(trialsDataBoomer(iTrial).SpikeCounts,2)+1];
        gaitPercentages = sort([dutyPercent:100:nSteps*100 0:100:nSteps*100]);
    end
    
    %get kinematics
    kins = {struct2cell(trialsDataBoomer(iTrial).PreTrialKinematics),...
        struct2cell(trialsDataBoomer(iTrial).Kinematics),...
        struct2cell(trialsDataBoomer(iTrial).PostTrialKinematics)};
    
%     kins = cellfun(@(x) cellfun(@(y) y(1:2,:), x, 'un',0), kins, 'un', 0);
    kins = cellfun(@(x) cat(1,x{2:end}), kins, 'un', 0);
%     kins = cellfun(@(x) getJointAngle(x), kins,'un',0);
    
    sampleSpikeTrain = {};
    for iSeg = 1:length(kins)
        
        nTimePoints = size(kins{iSeg},2);
        
        %generate noise
        segNoise = randn(nNeurons, nTimePoints)*0.01;
        
        %calc poisson rate
        timeVarRates = kinWeights*((kins{iSeg} - allTrialKinMeans)./repmat(allTrialKinAmps,1,size(kins{iSeg},2)));
        lambda = repmat(baselines,1,nTimePoints) + timeVarRates - lambdaMins' + segNoise;
        lambda(lambda<0) = 0;
        lambda(isnan(lambda)) = 0;
        
        %sample spikes
        %sample from poisson
        sampleSpikeTrain{iSeg} = poissrnd(lambda);
        
        if iSeg == 2
            allLambdas{iModelTrial} = lambda;
        end
        
    end
    
    %now time normalize
    if strcmpi(trialsDataBoomer(iTrial).Task,'Walk')
        normalizedData = timeNormalize([sampleSpikeTrain{:}],events+size(sampleSpikeTrain{1},2),gaitPercentages,-7:107);
        modeledData(iModelTrial).Task = 'Walk';
        modeledData(iModelTrial).GaitNormalizedPreTrialSpikeCounts = normalizedData(:,1:8);
        modeledData(iModelTrial).GaitNormalizedSpikeCounts = normalizedData(:,9:108);
        modeledData(iModelTrial).GaitNormalizedPostTrialSpikeCounts = normalizedData(:,109:end);
    else
        normalizedData = timeNormalize([sampleSpikeTrain{:}],events+size(sampleSpikeTrain{1},2),gaitPercentages,-7:100*nSteps+7);
        modeledData(iModelTrial).Task = 'WalkingObstacle';
        modeledData(iModelTrial).GaitNormalizedPreTrialSpikeCounts = normalizedData(:,1:8);
        modeledData(iModelTrial).GaitNormalizedSpikeCounts = normalizedData(:,9:100*nSteps+8);
        modeledData(iModelTrial).GaitNormalizedPostTrialSpikeCounts = normalizedData(:,100*nSteps+9:end);
    end
    
    iModelTrial = iModelTrial+1;
    
end

save('./Data/ModelDataBoomer','modeledData');

% next do starbuck
clear

% load in data
load('./Data/TrialsDataStarbuck.mat')
trialsDataStarbuck = trialsLegM1;

nNeurons = size(trialsDataStarbuck(1).SpikeCounts,1);
nKinematics = 15; %hip-toe 3 dims
dutyPercent = 69;
nSteps = 6;

% sample weights
kinWeights = randn(nNeurons,nKinematics);

% sample baseline
baselines = rand(nNeurons,1)*0.1;

% get mean firing rates across all walking and walking obstacle trials
badTrials = filterTrials(trialsDataStarbuck,90,5);
walkObsTrialInds = find(cellfun(@(x) strcmpi(x, 'WalkingObstacle'), {trialsDataStarbuck.Task}));
walkObsTrialInds = setdiff(walkObsTrialInds, badTrials);
walkTrialInds = find(cellfun(@(x) strcmpi(x, 'Walk'), {trialsDataStarbuck.Task}));
walkTrialInds = setdiff(walkTrialInds, badTrials);
allTrialSpikeCounts = {trialsDataStarbuck([walkObsTrialInds walkTrialInds]).SpikeCounts};
allTrialSpikeCounts = cat(2,allTrialSpikeCounts{:});
meanSpikeCounts = mean(allTrialSpikeCounts,2);

stdSpikeCounts = std(allTrialSpikeCounts-meanSpikeCounts,[],2);

% also get all kinematics to make them positive and scaled down
allTrialKin = cellfun(@(x) struct2cell(x),{trialsDataStarbuck([walkObsTrialInds walkTrialInds]).Kinematics}, 'un', 0);
% allTrialKin = cellfun(@(x) getJointAngle(x), allTrialKin,'un',0);
allTrialKin = cellfun(@(x) cat(1,x{:}), allTrialKin,'un',0);
allTrialKin = cat(2,allTrialKin{:});
allTrialKin(1:3,:) = [];
allTrialKinMeans = mean(allTrialKin,2);
allTrialKinAmps = max(allTrialKin,[],2) - min(allTrialKin,[],2);
alltrialLambdas = kinWeights*((allTrialKin - allTrialKinMeans)./repmat(allTrialKinAmps,1,size(allTrialKin,2)));
lambdaMins = min(alltrialLambdas');

allLambdas = {};
modeledData = {};
% go through each trial and simulate
iModelTrial = 1;
for iTrial = 1:length(trialsDataStarbuck)
    
    %only do walk and walk obstacle
    if ~any(iTrial==walkTrialInds) & ~any(iTrial==walkObsTrialInds)
        continue
    end
    
    %get trial lengths
    preTrialLength = size(trialsDataStarbuck(iTrial).PreTrialSpikeCounts,2);
    trialLength = size(trialsDataStarbuck(iTrial).SpikeCounts,2);
    postTrialLength = size(trialsDataStarbuck(iTrial).PostTrialSpikeCounts,2);
    
    %events for normalization
    if length(trialsDataStarbuck(iTrial).TrialEvents)==2
        events = [trialsDataStarbuck(iTrial).TrialEvents size(trialsDataStarbuck(iTrial).SpikeCounts,2)+1];
        gaitPercentages = [0 dutyPercent 100];
    else
        events = [trialsDataStarbuck(iTrial).TrialEvents size(trialsDataStarbuck(iTrial).SpikeCounts,2)+1];
        gaitPercentages = sort([dutyPercent:100:nSteps*100 0:100:nSteps*100]);
    end
    
    %get kinematics
    kins = {struct2cell(trialsDataStarbuck(iTrial).PreTrialKinematics),...
        struct2cell(trialsDataStarbuck(iTrial).Kinematics),...
        struct2cell(trialsDataStarbuck(iTrial).PostTrialKinematics)};
    
%     kins = cellfun(@(x) cellfun(@(y) y(1:2,:), x, 'un',0), kins, 'un', 0);
    kins = cellfun(@(x) cat(1,x{2:end}), kins, 'un', 0);
%     kins = cellfun(@(x) getJointAngle(x), kins,'un',0);
    
    sampleSpikeTrain = {};
    for iSeg = 1:length(kins)
        
        nTimePoints = size(kins{iSeg},2);
        
        %generate noise
        segNoise = randn(nNeurons, nTimePoints)*0.01;
        
        %calc poisson rate
        timeVarRates = kinWeights*((kins{iSeg} - allTrialKinMeans)./repmat(allTrialKinAmps,1,size(kins{iSeg},2)));
        lambda = repmat(baselines,1,nTimePoints) + timeVarRates - lambdaMins' + segNoise;
        lambda(lambda<0) = 0;
        lambda(isnan(lambda)) = 0;
        
        %sample spikes
        %sample from poisson
        sampleSpikeTrain{iSeg} = poissrnd(lambda);
        
        if iSeg == 2
            allLambdas{iModelTrial} = lambda;
        end
        
    end
    
    %now time normalize
    if strcmpi(trialsDataStarbuck(iTrial).Task,'Walk')
        normalizedData = timeNormalize([sampleSpikeTrain{:}],events+size(sampleSpikeTrain{1},2),gaitPercentages,-7:107);
        modeledData(iModelTrial).Task = 'Walk';
        modeledData(iModelTrial).GaitNormalizedPreTrialSpikeCounts = normalizedData(:,1:8);
        modeledData(iModelTrial).GaitNormalizedSpikeCounts = normalizedData(:,9:108);
        modeledData(iModelTrial).GaitNormalizedPostTrialSpikeCounts = normalizedData(:,109:end);
    else
        normalizedData = timeNormalize([sampleSpikeTrain{:}],events+size(sampleSpikeTrain{1},2),gaitPercentages,-7:100*nSteps+7);
        modeledData(iModelTrial).Task = 'WalkingObstacle';
        modeledData(iModelTrial).GaitNormalizedPreTrialSpikeCounts = normalizedData(:,1:8);
        modeledData(iModelTrial).GaitNormalizedSpikeCounts = normalizedData(:,9:100*nSteps+8);
        modeledData(iModelTrial).GaitNormalizedPostTrialSpikeCounts = normalizedData(:,100*nSteps+9:end);
    end
    
    iModelTrial = iModelTrial+1;
    
end

save('./Data/ModelDataStarbuck','modeledData');



function jointAngles = getJointAngle(kins)

for iKin = 1:length(kins)-2
    
    vect1 = kins{iKin+1} - kins{iKin};
    vect2 = kins{iKin+2} - kins{iKin+1};
    jointAngles(iKin,:) = acos(diag(vect1'*vect2)' ./ (sqrt(sum(vect1.^2,1)) .* sqrt(sum(vect2.^2,1))));
    
end

end

% 
