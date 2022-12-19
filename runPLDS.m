% fit PLDS model on spiking data. Treadmill data must be in './Data/'
% directory. Can run on both real collected spike counts or modeled spike
% counts. Will save the fitted low-dimensional latent states

clear; close all;

runOnModeledData = false;

if ~runOnModeledData
    %run on actual recorded data
    
    %set parameters
    dataset = 'Starbuck';
    limb = 'Leg';
    tasks = {'Walk','Obstacle','WalkingObstacle'};
    normalization = '';
    onlySwing = false;
    gaussStd = 5;
    dim=15;
    
    if strcmpi(dataset,'Boomer')
        
        if strcmpi(limb, 'Arm')
            load('./Data/TrialsDataBoomer.mat')
            trialsData = trialsArmM1;
            jointNames = {'Shoulder','Elbow','Wrist','HandKnuckle','PinkyFinger'};
        elseif strcmpi(limb, 'Leg')
            load('./Data/TrialsDataBoomer.mat')
            trialsData = trialsLegM1;
            walkingObstacleStep = -3:3;
            jointNames = {'Hip','Knee','Ankle','FootKnuckle','PinkyToe'};
        else
            error('limb has to be arm or leg')
        end
        hasKin = true;
        
    elseif strcmpi(dataset, 'Starbuck')
        
        if strcmpi(limb, 'Arm')
            load('./Data/TrialsDataStarbuck.mat')
            trialsData = trialsArmM1;
            walkingObstacleStep = -2:3;
            jointNames = {'Shoulder','Elbow','Wrist','HandKnuckle','PinkyFinger'};
        elseif strcmpi(limb, 'Leg')
            load('./Data/TrialsDataStarbuck.mat')
            trialsData = trialsLegM1;
            walkingObstacleStep = -3:2;
            jointNames = {'Hip','Knee','Ankle','FootKnuckle','PinkyToe'};
        else
            error('limb has to be arm or leg')
        end
        
        hasKin = false;
        
    else
        
        error('dataset has to be Starbuck or Boomer')
        
    end
    
else
    %run on modeled data
    
    %set parameters
    load('./Data/ModelDataBoomer.mat')
    normalization = 'GaitNormalized';
    tasks = {'WalkingObstacle'};
    trialsData = modeledData;
    dim = 9;
    
end

nNeurons = size(trialsData(1).([normalization 'SpikeCounts']),1);

% remove trials with lots of dropped signal, and also neurons with
% extremely low FRs
% badTrials = filterTrials(trialsData,90,5);
% excludedNeurons = [3 14 21 33];

%go through each task (basic walking, walking with obstacle, ect) and run
%the PLDS model
for iTask = 1:length(tasks)
    
    allTaskInds = find(cellfun(@(x) strcmpi(x,tasks{iTask}), {trialsData.Task}));
%     taskTrialInds{iTask} = setdiff(allTaskInds, badTrials);
    taskTrialInds{iTask} = allTaskInds;
    
    for iTrial=1:length(taskTrialInds{iTask})
        
        trialInd = taskTrialInds{iTask}(iTrial);
        spikeData{iTask}{iTrial} = trialsData(trialInd).([normalization 'SpikeCounts']);
%         spikeData{iTask}{iTrial}(excludedNeurons,:) = [];
        
        gaussStd = 5;
        preTrialSpikes = trialsData(trialInd).([normalization 'PreTrialSpikeCounts'])(:, end-gaussStd+1:end);
        postTrialSpikes = trialsData(trialInd).([normalization 'PostTrialSpikeCounts'])(:, 1:gaussStd);
        smoothedFR = convGauss([preTrialSpikes spikeData{iTask}{iTrial} postTrialSpikes], 10, gaussStd*10);
        firingRates{iTask}{iTrial} = smoothedFR(:,gaussStd+1:end-gaussStd);
        
    end

    [PLDSTrajTrain{iTask},PLDSTrajTest{iTask},LL{iTask}] = RunPLDSModel(spikeData{iTask},spikeData{iTask},dim,60,3600*20);

end
    
% can also run across all tasks combined
% [combPLDSTraj{iTask}, combPLDSModel{iTask}, combLL{iTask}] = RunPLDSModel(cat(2,spikeData{:}),cat(2,spikeData{:}),dim,60,3600);

% save
if ~runOnModeledData
    save([dataset '_PLDS'],'PLDSTrajTrain','sepLL');%,'combPLDSTraj','combPLDSModel','combLL')
else
    save([ './Data//PLDS/BoomerModel_PLDS_Unnormalized'],'PLDSTrajTrain','sepLL');
end


% 

