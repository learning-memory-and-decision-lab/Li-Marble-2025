function result = eyeRegressionPredResp(subnoStr,saveText,dirs)

% Prediction Phase Pupil Regression (locked to PredResp1)

% set paths
basePath = dirs.basePath;
eyeDir =    dirs.eyeDir;
behaveDir = dirs.behaveDir;
% smuPath =   dirs.smuPath;

% addpath(genpath(smuPath))
addpath(genpath(basePath))
cd(basePath);

set(groot,'defaultfigureposition',[400 250 900 750])

% setting certain parameters
blinkWindow = 150;
timeBefore = 2000;
timeAfter = 4000;
baselineTime = 2000; %use this for doBaseline == 1 or 2
baselineTimeStart = 2000; %use these for doBaseline == 3 time is ms before predResp
baselineTimeEnd   = 500;
leftArea = 4;
rightArea = 7;
nTrials = 476;
nBlockTrials = 238;
nAllTrials = 480;
nPracticeTrials = 60;
thresh = 0.2;
firstResp = 1;
secondResp = 1;
doBaseline = 3; % 1 = preResp baseline, 2 = preStim baseline, 3 = prePredResponse baseline
noPredTrials = [239,240,479,480];

gaussMu = 0; 
gaussSd = 1; 
sdWindow = 150; 

if doBaseline == 0
    titles = ["Intercept","STP","Entropy"];
else
    titles = ["Intercept","STP","STPCPOB","Entropy","Condition","Baseline"];
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
    
    if ~ismember(subno,[20280 2030 2079 2103 7777])
        trialStartAll = [0; find(eyeData(:,8) == 11)];
    else
        trial1StartAll = [0; find(eyeData(:,8) == 11)];
        trial2StartAll = [0; find(eyeData(:,8) == 13)];
        trialStartAll = sort([trial1StartAll;trial2StartAll]);
    end

    trialStart = [];
    for i = 2:(length(trialStartAll))
        diffStart = trialStartAll(i,1) - trialStartAll(i-1,1);
        if diffStart > 5
            trialStart = [trialStart trialStartAll(i)];
        end
    end
    trialStart = trialStart';
 
    %find start of pred phase
    if ~ismember(subno,[20280 2030 2079 2103 7777])
        predStartAll = [0; find(eyeData(:,8) == 7)];
    else
        pred1StartAll = [0; find(eyeData(:,8) == 7)];
        pred2StartAll = [0; find(eyeData(:,8) == 9)];
        predStartAll = sort([pred1StartAll;pred2StartAll]);
    end

    predStart = [];
    for i = 2:(length(predStartAll))
        diffStart = predStartAll(i,1) - predStartAll(i-1,1);
        if diffStart > 5
            predStart = [predStart predStartAll(i)];
        end
    end
    predStart = predStart';
    
    if ismember(subno,[20280,20380])
        predStart(1:8) = [];
        predStart(noPredTrials) = [];
    else 
        predStart(1:2*nPracticeTrials) = [];
        predStart(noPredTrials) = [];
    end

    %establish time window in terms of eyeData timesteps
    timeWindow = length([-timeBefore:timeAfter]);
    allTrialAreas = zeros([nTrials,timeWindow]);
    allPredAreas  = zeros([nTrials,baselineTime+1]);

    %fill trial areas based on trial start times
    endData = max(find(~isnan(eyeData(:,4))));
    for i = 1:length(trialStart)
        if trialStart(i) + timeAfter <= endData
            allTrialAreas(i,:) = eyeData(trialStart(i)-timeBefore:trialStart(i)+timeAfter,9);
            allPredAreas(i,:) = eyeData(predStart(i)-baselineTime:predStart(i),9);
            allTrialBlinks(i,:) = eyeData(trialStart(i)-timeBefore:trialStart(i)+timeAfter,10);
        else
            allTrialAreas(i,:) = nan;
            allPredAreas(i,:) = nan;
            allTrialBlinks(i,:) = 1;
        end
    end
        badBlinks = mean(allTrialBlinks,2)>thresh;
    allPredAreasNorm = reshape(nanzscore(allPredAreas(:)),size(allPredAreas));
    allTrialAreasNorm = reshape(nanzscore(allTrialAreas(:)),size(allTrialAreas));
    % allTrialAreasNorm = [reshape(nanzscore(allTrialAreas(1:end/2,:)),nBlockTrials,timeBefore+timeAfter+1);reshape(nanzscore(allTrialAreas(end/2+1:end,:)),nBlockTrials,timeBefore+timeAfter+1)];
%     imagesc(allTrialBlinks)
%     title(subnoStr)
%         fig = gcf;
%         scName = append(subnoStr,'.png');
%         scLoc = append(scDir,scName);
%         saveas(fig,scLoc)
    

    %get baseline areas
    stimStartAll = [0; find(eyeData(:,8) == 4)];
    stimStart = [];
    for i = 2:(length(stimStartAll))
        diffStim = stimStartAll(i,1) - stimStartAll(i-1,1);
        if diffStim > 5 
            stimStart = [stimStart stimStartAll(i)];
        end
    end

    stimStart = stimStart';



    baselineWindow = length([0:baselineTime]);
    if doBaseline == 2
        allBaselineAreas = zeros([length(stimStart),baselineWindow]);
        for i=1:length(stimStart)
            allBaselineAreas(i,:) = eyeData(stimStart(i)-baselineTime:stimStart(i),9);
        end

        if ~ismember(subno,[20280,20380])
            allBaselineAreas([1:60,180,300],:) = [];
        else
            allBaselineAreas([1:4,124,244],:) = [];
        end
        allBaselineAreasNorm = reshape(nanzscore(allBaselineAreas(:)),size(allBaselineAreas));
        allBaselineAreasDoubled = nan(2*size(allBaselineAreas,1),size(allBaselineAreas,2));
        allBaselineAreasDoubled(1:2:end,:) = allBaselineAreasNorm;
        allBaselineAreasDoubled(2:2:end,:) = allBaselineAreasNorm;
        allBaselineAreasNorm = allBaselineAreasDoubled;
    end
    
    
    if doBaseline == 1
        allBaselineAreas = zeros([length(stimStart),baselineWindow]);
        for i = 1:length(predStart)
            allBaselineAreas(i,:) = allPredAreasNorm(i,:);
        end
        allBaselineAreasNorm = allBaselineAreas;
    end

    if doBaseline == 3
        allBaselineAreas = zeros([length(stimStart),length(baselineTimeEnd:baselineTimeStart)]);
        for i = 1:length(predStart)
            allBaselineAreas(i,:) = allTrialAreasNorm(i,timeBefore-baselineTimeStart+1:timeBefore-baselineTimeEnd+1);
        end
        allBaselineAreasNorm = allBaselineAreas;
    end

    %create logical list of surprise trials
    surpriseTrials = [0,0;((alldata.surprise(2:end,:) == 1) > 1)];

    %define block 3 and 4 areas in terms of oddball/changepoint depending
    %on condition
    context = zeros(nTrials,1);
    isChangepoint = zeros(nTrials,1);
    isOddball = zeros(nTrials,1);
    if alldata.condition(3)==1
        changepointBlockTrials = 1:nTrials/2-2;
        oddballBlockTrials = (nTrials/2-1):nTrials-4;
        context(1:nTrials/2) = 1;
        context((nTrials/2+1):nTrials) = -1;
        isChangepoint(1:nBlockTrials) = 1;
        isOddball((nBlockTrials+1):end) = 1;
    elseif alldata.condition(3)==2
        changepointBlockTrials = (nTrials/2-1):nTrials-4;
        oddballBlockTrials = 1:nTrials/2-2;
        context(1:nTrials/2) = -1;
        context((nTrials/2+1):nTrials) = 1;
        isOddball(1:nBlockTrials) = 1;
        isChangepoint((nBlockTrials+1):end) = 1;
    else
        disp('something is very wrong with the alldata.condition variable')
        keyboard
    end

    %define colors and find sin/cos of left/right colors
    colors = alldata.colorArray(61:end,:);

    sinColor = reshape(sin(deg2rad(colors(:,:))),[nAllTrials,1]);
    cosColor = reshape(cos(deg2rad(colors(:,:))),[nAllTrials,1]);
    
    sinColor(noPredTrials) = [];
    cosColor(noPredTrials) = [];

    %add trial number
    trialNum = [1:nTrials]';
    
    %CALCULATING LEARNING RATE
    
    behaveMat2 = sprintf('subCombined/%s_3and4BlockData.mat',subnoStr);
    file3 = fullfile(behaveDir,behaveMat2);
    alldataShort = load(file3);

    predictions = alldataShort.alldata.pred;
    outcomes = deg2rad(alldataShort.alldata.colorArray);
    newBlock = nBlockTrials + 1;

    %run CLR function (done twice, one for left and right stimulus
    [~,UP1,PE1] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polar');
    [~,UP2,PE2] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polar');
    
    %concatenate UP and PE
    UP = [UP1,UP2;nan,nan];
    UP = reshape(UP',nAllTrials,1);
    PE = [PE1,PE2];
    PE = reshape(PE',nAllTrials,1);
    %truncate LR range from 0 to 1
    newUpdate = UP;
    overupdate = abs(UP(:,:)) > abs(PE(:,:)) & sign(UP) == sign(PE);
    newUpdate(overupdate) = UP(overupdate) - 2.*(UP(overupdate)-PE(overupdate));

    LR = UP./PE;

    LR(LR > 1) = 1;
    LR(LR < 0) = 0;
    
    %define LR in terms of blocks
    block3LR = LR(1:nTrials/2,:);
    block4LR = LR(nTrials/2+1:nTrials,:);
    
    %find LR for surprise trials in each block
%     surpriseNew = [0,0;(abs(circ_dist(outcomes(2:end,:),outcomes(1:end-1,:))) > 1 & alldataShort.alldata.surprise(2:end,:) == 1)];
% 
%     block3surpriseLR = block3LR(find(surpriseNew(1:nBlockTrials,:) == 1));
%     block4surpriseLR = block4LR(find(surpriseNew(nBlockTrials+1:end,:) == 1)); 
    %divide LR into changepoint and oddball groups based on condition
%     if alldataShort.alldata.condition(1) == 1
%         oddballLR = rmmissing(block4surpriseLR);
%         changepointLR = rmmissing(block3surpriseLR);
%     else
%         oddballLR = rmmissing(block3surpriseLR);
%         changepointLR = rmmissing(block4surpriseLR);
%     end
    %Compute avg/max LR for a given trial
%     meanLR = nanmean(LR,2);
%     maxLR = max(LR,[],2);

    %add mean oddball, changepoint, and total learning rate to group
    %variable

    %add perceptual bias as a regressor
    % bias = alldata.estErr(:,:)./alldata.predictErr(:,:);
    % bias(1:60,:) = [];
    % bias = reshape(bias',nAllTrials,1);
    % bias(noPredTrials) = [];
    % bias(bias > 1) = 1;
    % bias(bias < 0) = 0;
%     meanBias = mean(bias,2);
%     minBias = min(bias,[],2);
    
    %STP * CP/OB
%     CPOB = double(surpriseTrials);
%     CPOB(oddballBlockTrials) = -CPOB(oddballBlockTrials);
     
    %load surprise previously calculated using model
    allModelData = load(fullfile(behaveDir,['allModelData',saveText,'/', subnoStr, '_allBlockData.mat']));
    allModelData = allModelData.allDataStruct;

    surpriseCP = allModelData.surpriseCP;
    surpriseOB = allModelData.surpriseOB;
    
    surpriseCP = reshape(surpriseCP,[],2)';
    surpriseCP(:,end) = [];
    surpriseCP = reshape(surpriseCP,[],1);

    surpriseOB = reshape(surpriseOB,[],2)';
    surpriseOB(:,end) = [];
    surpriseOB = reshape(surpriseOB,[],1);

    if alldata.condition(1) == 1
    STP = [(surpriseCP);(surpriseOB)];
    else
    STP = [(surpriseOB);(surpriseCP)];
    end

    STPCPOB = zscore(STP);
    STPCPOB(find(context==-1)) = -STPCPOB(find(context==-1));

    %load entropy previously calculated using model

    entropyCP = allModelData.entropyCP;
    entropyOB = allModelData.entropyOB;

    if alldata.condition(1) == 1
    entropy = [nanzscore(entropyCP);nanzscore(entropyOB)];
    else
    entropy = [nanzscore(entropyOB);nanzscore(entropyCP)];
    end
    
    LR(noPredTrials) = [];
    entropy(noPredTrials) = [];

    surpriseTrials = reshape(surpriseTrials(61:end,:)',[nAllTrials,1]);
    surpriseTrials(noPredTrials) = [];
    %mean 2 xes across each trial

    for i = 1:nTrials
        if mod(i,2) == 1
            meanSTP(i,:) = mean([STP(i),STP(i+1)]);
            meanEntropy(i,:) = mean([entropy(i),entropy(i+1)]);
            % meanLR(i,:) = mean([LR(i),LR(i+1)]);
            meanSTPCPOB(i,:) = mean([STPCPOB(i),STPCPOB(i+1)]);
            % meanBias(i,:) = mean([bias(i),bias(i+1)]);
        elseif mod(i,2) == 0
            meanSTP(i,:) = mean([STP(i),STP(i-1)]);
            meanEntropy(i,:) = mean([entropy(i),entropy(i-1)]);
            % meanLR(i,:) = mean([LR(i),LR(i-1)]);
            meanSTPCPOB(i,:) = mean([STPCPOB(i),STPCPOB(i-1)]);
            % meanBias(i,:) = mean([bias(i),bias(i-1)]);
        end
    end

    meanSTPMinus1 = [0;0;meanSTP(1:end-2)];
    meanSTPCPOBMinus1 = [0;0;meanSTPCPOB(1:end-2)];

    %ACTUAL REGRESSION PART
    %define xes for regression
    
    if doBaseline == 0
        %xes = [ones(nTrials,1),meanSTP,meanEntropy];
        xes = [ones(nTrials,1),STP,entropy];
    else
        %xes = [ones(nTrials,1),meanSTP,meanEntropy,nanmean(allBaselineAreasNorm,2)];
        xes = [ones(nTrials,1),STP,STPCPOB,entropy,context,nanmean(allBaselineAreasNorm,2)];
    end

    %remove non-prediction trials and bad eye Trials
    isPredTrial=isfinite(LR);

    %specify which trials were predicted first (because trials in other arrays are left right left right)
    isFirstEye = logical(mod(1:nAllTrials,2))';
    isFirstPred = isFirstEye;
    alldataShort.alldata.chosenTargPredict = [alldataShort.alldata.chosenTargPredict(2:end,:);nan,nan];
    badBlinksBehav = badBlinks;
    for i = 1:length(alldataShort.alldata.chosenTargPredict)
       if alldataShort.alldata.chosenTargPredict(i,1) ~= 1
           isFirstPred(2*i) = 1;
           isFirstPred(2*i-1) = 0;
       end
    end
    
    chosenTargPredict = alldataShort.alldata.chosenTargPredict;
    chosenTargPredict([120,240],:) = [];
    xesSwitched =xes;

    for i = 1:size(chosenTargPredict,1)
       if chosenTargPredict(i,1) ~= 1
           xesSwitched(2*i,:) = xes(2*i-1,:);
           xesSwitched(2*i-1,:) = xes(2*i,:);          
       end
    end

    isFirstPred(noPredTrials) = [];
    isFirstEye(noPredTrials) = [];    
    isSecondEye = logical(~isFirstEye);
    isSecondPred = logical(~isFirstPred);
    noY = [zeros(nBlockTrials-2,1);1;1;zeros(nBlockTrials-2,1);1;1];

    if firstResp == 1
        goodAndPred = isFirstPred;
    elseif secondResp == 1
        goodAndPred = isSecondPred;
    end

    %zscore xes
    zX=[xesSwitched(:,1), nanzscore(xesSwitched(:,2:end-1)), xesSwitched(:,end)];
    
    zSelX= zX;
    
    B = [];
    %select good trial areas 
    if firstResp == 1 && secondResp == 0
        isY = isFirstEye;
    elseif secondResp == 1 && firstResp == 0
        isY = isSecondEye;
    elseif firstResp == 1 && secondResp == 1
        isY = isFirstEye | isSecondEye;
    end
        
    y=allTrialAreasNorm; 
    toi = isY & ~badBlinks & max(abs(allTrialAreasNorm),[],2)<8;
    %run actual regression
    for i =1:timeWindow
        [B(:,i)] = regress(y(toi,i), zSelX(toi,:));
    end

    %run regression on derivative
    gaussX = linspace(-4*gaussSd,4*gaussSd,8*sdWindow);
    gaussY = 1/(2*pi*gaussSd)*exp(-(gaussX-gaussMu).^2/(2*gaussSd^2));
    for i = 1:nTrials
        convY(i,:) = conv(y(i,:),gaussY);
    end
    y = convY;
    y = diff(y,1,2);
    y = y(:,sdWindow*4+1:timeWindow+sdWindow*4);
    for i =1:timeWindow
        [diffB(:,i)] = regress(y(isY & ~badBlinks,i), zSelX(isY & ~badBlinks,1:end-1));
    end
    %save subject coefficients
    result.B = B;
    result.diffB = diffB;
%%