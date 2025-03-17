function resultEye = eyeRegressionFunc(subnoStr,saveText,dirs)

%basic pupil regression
% set paths
basePath = dirs.basePath;
eyeDir =    dirs.eyeDir;
behaveDir = dirs.behaveDir;
% smuPath =   dirs.smuPath;

cd(basePath);
addpath(genpath(basePath))
% addpath(genpath(smuPath)) 

set(groot,'defaultfigureposition',[400 250 900 750])

% setting certain parameters
blinkWindow = 150;
leftArea = 4;
rightArea = 7;
nTrials = 300;
nBlockTrials = 120;
nPracticeTrials = 60;
firstTrial = [61,62,181,182];
notFirst = true(nTrials,1);
notFirst(firstTrial) = false;

baselineTime = 1000;
timeBefore = 1000; %must be >= baselineTime
timeAfter = 4000;

badSubs = 30;
thresh = 0.2;
doBaseline = 1;

subs = [2079,2035,2033,2045,2099,2062,2044,2075,2057,2034,2048,2066,2104,2041,2086,2039,2074,2071,2084,2073,2056,2098,2082,2095,2042,20280,2037,2105,2063,2094 ... 
        2096,2103,2070,2097,2090,2100,2055,2081,2083,2087,2092,20380,2069,2050,2053,2093,2058,2080,2046,2030,2036,2060,2040,2101,2065,2068,2106,2043,2088,2047]; 


if doBaseline == 1
    titles = ["Intercept","MeanSTP","MeanSTPCPOB","MeanEntropy","Condition","SumSinColor","SumCosColor","Baseline 1000 ms"];
else
    titles = ["Intercept","MeanSTP","MeanSTPCPOB","MeanEntropy","Condition","SumSinColor","SumCosColor"];
end

    %load eyedata
    subno = str2num(subnoStr);
    eyeMat = sprintf('%s.mat',subnoStr);
    file1 = fullfile(eyeDir,eyeMat);
    load(file1)

    eyeData = data;

    %interpolate single eye data
    if size(eyeData,2)==5
       eyeData(:,6:8)=0;
       eyeData(:,8)=eyeData(:,5);
       eyeData(:,5:7)=eyeData(:,2:4);
    end
    %Remove bad eye from data with one bad eye
    if ismember(subno,[2030,2063])
        eyeData(:,2:4)=eyeData(:,5:7);
    end
    if ismember(subno,[2035 2057 2058 2062 2071 2083 2087 2098 2101 2103 2106])
        eyeData(:,5:7)=eyeData(:,2:4);
    end
    %load behavioral data
    behaveMat = sprintf('allSubCombined/%s_allBlockData.mat',subnoStr);
    file2 = fullfile(behaveDir,behaveMat);
    load(file2)

    % identify blinks change from 0 to NaN
    % eyeData(:,rightArea) = eyeData(:,leftArea);
    eyeData((eyeData(:,leftArea)==0),leftArea)=NaN;
    eyeData((eyeData(:,rightArea)==0),rightArea)=NaN;
    
    % make sure other eye is also NaN when one eye is closed
    eyeData(isnan(eyeData(:,leftArea)),rightArea)=NaN;
    eyeData(isnan(eyeData(:,rightArea)),leftArea)=NaN;
    
    %NaN blink window
    blinks = find(isnan(eyeData(:,leftArea)));
    isBlink = isnan(eyeData(:,leftArea));

    for t=1:length(blinks)
        eyeData(max(blinks(t)-blinkWindow,1):min(blinks(t)+blinkWindow,length(eyeData)),leftArea)= NaN;
        eyeData(max(blinks(t)-blinkWindow,1):min(blinks(t)+blinkWindow,length(eyeData)),rightArea)= NaN;
    end
     
    %interpolate blink window 
    [value,indices] = fillmissing(eyeData(:,leftArea),'linear', 'EndValues', 'none');
    eyeData(indices,leftArea)=value(indices);
    [value,indices] = fillmissing(eyeData(:,rightArea),'linear', 'EndValues', 'none');
    eyeData(indices,rightArea)=value(indices);
    
    %calculate mean pupil area
    meanArea=(eyeData(:,leftArea)+eyeData(:,rightArea))/2;
    
    %add meanArea to eyeData
    eyeData(:,9) = meanArea;
    eyeData(:,10) = isBlink;

    %cut out data before start and after end of session
    startTime = find(eyeData(:,8) == 15);
    endTime = find(eyeData(:,8) == 14);
    eyeData = eyeData(startTime(end):endTime(1),:);
    
    % Find start of each trial (stimOn) in terms of eyeData timesteps
    
    trialStartAll = [0; find(eyeData(:,8) == 4)];
    trialStart = [];
    for i = 2:(length(trialStartAll))
        diffStart = trialStartAll(i,1) - trialStartAll(i-1,1);
        if diffStart > 5
            trialStart = [trialStart trialStartAll(i)];
        end
    end
    trialStart = trialStart';

    %establish time window in terms of eyeData timesteps
    timeWindow = length([-timeBefore:timeAfter]);
    allTrialAreas = zeros([length(trialStart),timeWindow]);
    allTrialBlinks = zeros([length(trialStart),timeWindow]);


    %fill trial areas based on trial start times
    endData = max(find(~isnan(eyeData(:,4))));
    for i = 1:height(allTrialAreas)
        if trialStart(i) + timeAfter <= endData
            allTrialAreas(i,:) = eyeData(trialStart(i)-timeBefore:trialStart(i)+timeAfter,9);
            allTrialBlinks(i,:) = eyeData(trialStart(i)-timeBefore:trialStart(i)+timeAfter,10);
        else
            allTrialAreas(i,:) = nan;
            allTrialBlinks(i,:) = 1;
        end
    end
    badBlinks = mean(allTrialBlinks,2)>thresh;

    if ismember(subno,[20280,20380])
        allTrialAreas(1:4,:) = [];
        badBlinks = [true(56,1);badBlinks];
    end

    %normalize trial area
    allTrialAreasNorm = reshape(nanzscore(allTrialAreas(:)),size(allTrialAreas));

    %calculate mean baseline for each trial
    if ismember(subno,[20280,20380])
        allBaselineMeansNorm = [nan(nPracticeTrials,1);mean(allTrialAreasNorm(:,(timeBefore-baselineTime+1:timeBefore+1)),2)];
        allBaselineMeans = [nan(nPracticeTrials,1);mean(allTrialAreas(:,(timeBefore-baselineTime+1:timeBefore+1)),2)];
    else
        allBaselineMeansNorm = mean(allTrialAreasNorm(:,(timeBefore-baselineTime+1:timeBefore+1)),2);
        allBaselineMeans = mean(allTrialAreas(:,(timeBefore-baselineTime+1:timeBefore+1)),2);
    end
    
    if ismember(subno,[20280,20380])
        allTrialAreasNorm = [zeros(nPracticeTrials,length(allTrialAreasNorm));allTrialAreasNorm];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    PREPROCESSING END    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%     %create logical list of surprise trials
%     surpriseTrials = [0,0;((alldata.surprise(2:end,:) == 1) & abs(circ_dist(deg2rad(alldata.colorArray(2:end,:)),deg2rad(alldata.colorArray(1:end-1,:)))) > 1)];
    
    surpriseTrials = alldata.surprise;
    surpriseTrials = sum(surpriseTrials(:,:),2);
    surpriseTrials = surpriseTrials >= 1;
    surpriseTrials(61,181) = 1;

    %define block 3 and 4 areas in terms of oddball/changepoint depending
    %on condition
    context = zeros(nTrials,1);
    if alldata.condition(3)==1
        changepointBlockTrials = 61:180;
        oddballBlockTrials = 181:nTrials;
        context(61:180) = 1;
        context(181:nTrials) = -1;
    elseif alldata.condition(3)==2
        changepointBlockTrials = 181:nTrials;
        oddballBlockTrials= 61:180;
        context(61:180) = -1;
        context(181:nTrials) = 1;
    else
        disp('something is very wrong with the alldata.condition variable')
        keyboard
    end
 
    %define colors and find sin/cos of left/right colors
    colors = alldata.colorArray;

    sinLeftColor = sin(deg2rad(colors(:,1)));
    sinRightColor = sin(deg2rad(colors(:,2)));
    cosLeftColor = cos(deg2rad(colors(:,1)));
    cosRightColor = cos(deg2rad(colors(:,2)));

    %add trial number
    trialNum = [1:nTrials]';
    trialNumInBlock = [zeros(1,60),1:120,1:120]';
    
    %CALCULATING LEARNING RATE
    
    behaveMat2 = sprintf('subCombined/%s_3and4BlockData.mat',subnoStr);
    file3 = fullfile(behaveDir,behaveMat2);
    alldataShort = load(file3);

    predictions = alldataShort.alldata.pred;
    outcomes = alldataShort.alldata.est;
    newBlock = nBlockTrials + 1;
    
    %run CLR function (done twice, one for left and right stimulus
    [LR1,UP1,PE1] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polar');
    [LR2,UP2,PE2] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polar');
    
    %concatenate LR UP and PE
    LR = [LR1,LR2;nan,nan];
    UP = [UP1,UP2;nan,nan];
    PE = [PE1,PE2];

    %truncate LR range from 0 to 1
    newUpdate = UP;
    overupdate = abs(UP(:,:)) > abs(PE(:,:)) & sign(UP) == sign(PE);
    newUpdate(overupdate) = UP(overupdate) - 2.*(UP(overupdate)-PE(overupdate));

    LR = UP./PE;

    LR(LR > 1) = 1;
    LR(LR < 0) = 0;
    
    %define LR in terms of blocks
    block3LR = LR(1:nBlockTrials,:);
    block4LR = [LR(nBlockTrials+1:239,:);nan,nan];
    
    %find LR for surprise trials in each block
    surpriseNew = [0,0;(abs(circ_dist(outcomes(2:end,:),outcomes(1:end-1,:))) > 1 & alldataShort.alldata.surprise(2:end,:) == 1)];

    block3surpriseLR = block3LR(find(surpriseNew(1:nBlockTrials,:) == 1));
    block4surpriseLR = block4LR(find(surpriseNew(nBlockTrials+1:end,:) == 1));
    
    %divide LR into changepoint and oddball groups based on condition
    if alldataShort.alldata.condition(1) == 1
        oddballLR = rmmissing(block4surpriseLR);
        changepointLR = rmmissing(block3surpriseLR);
    else
        oddballLR = rmmissing(block3surpriseLR);
        changepointLR = rmmissing(block4surpriseLR);
    end

    %Compute avg/max LR for a given trial
    meanLR = [nan(nPracticeTrials,1);nanmean(LR,2)];
    maxLR = [nan(nPracticeTrials,1);max(LR,[],2)];

    %add mean oddball, changepoint, and total learning rate to group
    %variable

    %add perceptual bias as a regressor
    % bias = alldata.estErr(:,:)./alldata.predictErr(:,:);
    % bias(bias > 1) = 1;
    % bias(bias < 0) = 0;
    % meanBias = mean(bias,2);
    % minBias = min(bias,[],2);

    %load surprise previously calculated using model
    allModelData = load(fullfile(behaveDir,['allModelData',saveText,'/', subnoStr, '_allBlockData.mat']));
    allModelData = allModelData.allDataStruct;

    meanSurpriseCP = mean(reshape(allModelData.surpriseCP',[nBlockTrials,2]),2);
    maxSurpriseCP = max(reshape(allModelData.surpriseCP',[nBlockTrials,2]),[],2);

    meanSurpriseOB = mean(reshape(allModelData.surpriseOB',[nBlockTrials,2]),2);
    maxSurpriseOB = max(reshape(allModelData.surpriseOB',[nBlockTrials,2]),[],2);

    if alldata.condition(1) == 1
    meanSurprise = [zeros(nPracticeTrials,1);(meanSurpriseCP);(meanSurpriseOB)];
    maxSurprise = [zeros(nPracticeTrials,1);zscore(maxSurpriseCP);zscore(maxSurpriseOB)];
    else
    meanSurprise = [zeros(nPracticeTrials,1);(meanSurpriseOB);(meanSurpriseCP)];
    maxSurprise = [zeros(nPracticeTrials,1);zscore(maxSurpriseOB);zscore(maxSurpriseCP)];
    end

    maxSTPCPOB = zscore(maxSurprise);
    maxSTPCPOB(oddballBlockTrials) = -maxSTPCPOB(oddballBlockTrials);
    meanSTPCPOB = (meanSurprise);
    meanSTPCPOB(oddballBlockTrials) = -meanSTPCPOB(oddballBlockTrials);

    % meanLRSTP = zscore(meanSurprise).*nanzscore(meanLR);
    %load entropy previously calculated using model

    meanEntropyCP = mean(reshape(allModelData.entropyCP',[nBlockTrials,2]),2);
    maxEntropyCP = max(reshape(allModelData.entropyCP',[nBlockTrials,2]),[],2);

    meanEntropyOB = mean(reshape(allModelData.entropyOB',[nBlockTrials,2]),2);
    maxEntropyOB = max(reshape(allModelData.entropyOB',[nBlockTrials,2]),[],2);

    if alldata.condition(1) == 1
    meanEntropy = [zeros(nPracticeTrials,1);zscore(meanEntropyCP);zscore(meanEntropyOB)];
    maxEntropy = [zeros(nPracticeTrials,1);maxEntropyCP;maxEntropyOB];
    else
    meanEntropy = [zeros(nPracticeTrials,1);zscore(meanEntropyOB);zscore(meanEntropyCP)];
    maxEntropy = [zeros(nPracticeTrials,1);maxEntropyOB;maxEntropyCP];
    end

    % Find trials since surprising event
    trialsSinceSurprise = zeros(nTrials,1);
    for i = 1:nTrials
        if i == 1
            trialsSinceSurprise(i) = 0;
        elseif surpriseTrials(i) == 1
            trialsSinceSurprise(i) = 0;
        else
            trialsSinceSurprise(i) = trialsSinceSurprise(i-1)+1;
        end
    end

    meanSTPMinus1 = [0;meanSurprise(1:end-1)];
    meanSTPCPOBMinus1 = [0;meanSTPCPOB(1:end-1)];

    meanSTPTrialNum = meanSurprise.*(trialNum-nPracticeTrials);
    meanSTPCPOBTrialNum = meanSTPCPOB.*(trialNumInBlock);
    
    meanEntropyCPOB = zscore(meanEntropy);
    meanEntropyCPOB(oddballBlockTrials) = -meanEntropyCPOB(oddballBlockTrials);

    meanPE = nanmean(abs(alldata.predictErr),2);
    meanRE = nanmean(abs(alldata.estErr),2);

    %Calculate regressed bias and add to allData
    % biasXes = zeros(2*nBlockTrials,2*nBlockTrials+1);
    % biasXes(:,1) = 1;
    % for i = 1:2*nBlockTrials
    %     biasXes(2*i-1,i+1) = alldataShort.alldata.predictErr(i,1);
    %     biasXes(2*i,i+1) = alldataShort.alldata.predictErr(i,2);
    % end
    % biasYcol = reshape(alldataShort.alldata.estErr',[],1);
    % regBias = regress(biasYcol,biasXes);
    % regBias(regBias>1) = 1;
    % regBias(regBias<0) = 0;
    % 
    % regBias = [zeros(60,1);regBias(2:end)];

    %new regressed way of calculating learning
    % updateXes = zeros(2*nBlockTrials,2*nBlockTrials+1);
    % updateXes(:,1) = 1;
    % for i = 1:2*nBlockTrials
    %     updateXes(2*i-1,i+1) = PE(i,1);
    %     updateXes(2*i,i+1) = PE(i,2);
    % end
    % updateYcol = reshape(alldataShort.alldata.predUpdate',[],1);
    % regLRs = regress(updateYcol,updateXes);
    % regLRs(regLRs>1) = 1;
    % regLRs(regLRs<0) = 0;
    % regLRs = [zeros(60,1);regLRs(2:end)];

    %ACTUAL REGRESSION PART
    %define xes for regression
    if doBaseline==1
        xes = [ones(nTrials,1),meanSurprise,meanSTPCPOB,meanEntropy,context,sinLeftColor+sinRightColor,cosLeftColor+cosRightColor,allBaselineMeans];
    else
        xes = [ones(nTrials,1),meanSurprise,meanSTPCPOB,meanEntropy,context,sinLeftColor+sinRightColor,cosLeftColor+cosRightColor];        
    end

    isPredTrial=isfinite(meanLR);

    toi = isPredTrial & ~badBlinks & notFirst;

    selX= xes(toi,:);

    %zscore xes
    zSelX=[selX(:,1), zscore(selX(:,2:end))];

    B = [];
    %run actual regression
    y=allTrialAreasNorm;
    for i =1:timeWindow
        [B(:,i),BINT,R,RINT,STATS] = regress(y(toi,i), zSelX);
    end
        [BBaseline] = regress(zSelX(:,end),zSelX(:,1:end-1));

resultEye.B = B;
resultEye.BBaseline = BBaseline;

