%eye_analysis('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/Meera022620.mat',...
%200, 4000, 150);
function result = eyeAnalysis(path, before, after, blinkWindow,subNum,baseDir,DP,sharedFuncDir)
%% Loading variables and data

% whichComp=4;
showPlots=false;

% if whichComp==1 % iMac in lab manager office
% 
%     baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
%     modelDataDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined';
%     DP=sprintf('/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/%s_allBlockData.mat',subNum);
%     addpath(genpath(baseDir))
%     % Add entire set of subdirectories to the path
%     addpath(genpath("~/Dropbox (Brown)/sharedMatlabUtilities"));
% 
% elseif whichComp==2 
% 
%     baseDir='/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry';
%     modelDataDir='/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry/behaveData/subCombined';
%     DP=sprintf('/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry/behaveData/subCombined/%s_3and4BlockData.mat',subNum);
% 
%     addpath(genpath(baseDir))
%     addpath(genpath("/Users/ttli/Dropbox/ALP_projectFolder/sharedMatlabUtilities"));
% 
% 
% elseif whichComp==3 
% 
%     baseDir='/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry';
%     modelDataDir='/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry/behaveData/subCombined';
%     DP=sprintf('/Users/ttli/Dropbox/ALP_projectFolder/vwm_task_Harry/behaveData/subCombined/%s_allBlockData.mat',subNum);
%     addpath(genpath(baseDir))
% 
% 
%     addpath(genpath("'/Users/ttli/Dropbox/ALP_projectFolder/sharedMatlabUtilities"));
% 
% elseif whichComp==4 
% 
%     baseDir='/Users/gulli/Dropbox (Brown)/ALP_Summer2022/vwm_task_Harry';
%     modelDataDir='/Users/gulli/Dropbox (Brown)/ALP_Summer2022/vwm_task_Harry/behaveData/subCombined';
%     DP=sprintf('/Users/gulli/Dropbox (Brown)/ALP_Summer2022/vwm_task_Harry/behaveData/subCombined/%s_3and4BlockData.mat',subNum);
%     addpath(genpath(baseDir))
% 
% 
%     addpath(genpath("'/Users/gulli/Dropbox (Brown)/ALP_Summer2022/sharedMatlabUtilities"));
% 
% end
%keyboard
% if whichComp>2
%     datapath=fullfile(baseDir,'EEG_data'); %eeg data folder
% else
%     datapath=fullfile(baseDir,'eeg_data');
% end

addpath(genpath(baseDir))
% addpath(genpath(sharedFuncDir))

% Types of events
start           =15;
endBlock        =14;

block           =1;
instructionsOn  =2;
instructionsOff =3;
stimOn          =4;
stimOff         =5;
cue1on          =6;
responseMade1   =7;
cue2on          =8;
responseMade2   =9;
predCue1on      =10;
predResp1       =11;
predCue2on      =12;
predResp2       =13;

%eye data structure columns
timeStamp=1;
gazeLx=2;
gazeLy=3;
left_area = 4;
gazeRx=5;
gazeRy=6;
right_area = 7;
event = 8;

behaveDataDir=fullfile(baseDir, 'behaveData');
cd(baseDir)
% Loading data
% vars = load('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/martin022520_10.mat');
% vars = load('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/Meera022620.mat');

vars = load(path);

dataEye = vars.data;

dataEyeCopy=dataEye;

%DP="/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/2001_allBlockData.mat";
%DP=['/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/Meera1.mat';
%vars = shared_variables(DP);
vars = shared_variables(DP);


nTrials = 120; %400;%vars.nTrials;

nTrialsRand=40; %400;

% keyboard
%load behave data for random block
% fn=fullfile(behaveDataDir,'subCombined',[subNum,'_randBlockData.mat']);
% randData=load(fn);
% randData=randData.alldata;

%load behave data for condition blocks
fn=fullfile(behaveDataDir,'subCombined',[subNum,'_3and4BlockData.mat']);
condData=load(fn);
condData=condData.alldata;

%randData=condData;


% stimOn=0;
% stimOff=500; % ish (might be 548?)ALSO maskOn
% maskOff=1500; % mask on for 1s
% cue1on=2400; % retention interval for 900ms

%imagesc(dataEye(:,4)')

preprocessBefore=0;

runPred=0;

makeVid=0;

screenX=1920;
screenY=1080;

targ1Xpos=800;
targ2Xpos=1120;

targ1Ypos=540;
targ2Ypos=540;

rectSize= 120;
plotGaze=1;
checkColorPupil=0;
oneSubOnly=0;

power=0;
%%

% behaveFiles=dir(fullfile(behaveDataDir, ['vwm_multiStim_color_', subNum, '*']));
% behavFNs={behaveFiles.name};
%
% for block=3:length(behavFNs)
%
%     load(fullfile(behaveDataDir, behavFNs{block}))
%     %         fileName=sprintf('vwm_multiStim_color_%s%s.mat',sn_str,num2str(block));
%     %         fileDir=fullfile(dataDir,fileName);
%     %         load(fileDir)
%
%     nTrial=size(data.reportedValue,1);
%
%     estErr=nan(nTrial,nStim);
%     est=nan(nTrial,nStim);
%
%     predictErr=nan(nTrial,nStim);
%     pred=nan(nTrial,nStim);
%     predDeg=nan(nTrial,nStim);
%
%
%     trueEst=nan(nTrial,nStim);
%     truePred=nan(nTrial,nStim);
%
%     predEstDiff=nan(nTrial,nStim);
%     predUpdate=nan(nTrial,nStim);
%
%     surpriseTrial=nan(nTrial,nStim);
%
%
%     %        if ~data.prefs.doFeedback
%     for r=1:nTrial
%         for c=1:nStim
%             %if data.chosenTarg(r,c)==1
%             if data.chosenTargEst(r,c)==1
%                 estDeg(r,c)=data.reportedValue(r,1);
%                 est(r,c)=deg2rad(data.reportedValue(r,1));
%                 estErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,1)),deg2rad(data.presentedColor(r,1)));
%                 trueEst(r,c)=deg2rad(data.presentedColor(r,1));
%             else
%                 estDeg(r,c)=data.reportedValue(r,2);
%                 est(r,c)=deg2rad(data.reportedValue(r,2));
%                 estErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,2)),deg2rad(data.presentedColor(r,2)));
%                 trueEst(r,c)=deg2rad(data.presentedColor(r,2));
%             end
%         end
%     end
%     %         else
%     %             for r=1:nTrial
%     %                 for c=1:nStim
%     %
%     %                 end
%     %             end
%     %end
%
%     if data.prefs.doPredict
%         predV=data.predictedValue;
%         predVUpdated=cat(1,predV(2:end,:),[0 0]);
%         data.chosenTargPredict=cat(1,[nan nan],data.chosenTargPredict);
%         for r=2:nTrial
%             for c=1:nStim
%                 %if data.chosenTarg(r,c)==1
%                 if data.chosenTargPredict(r,c)==1
%                     predDeg(r,c)=data.predictedValue(r-1,1);
%                     pred(r,c)=deg2rad(data.predictedValue(r-1,1));
%                     predictErr(r,c)=circ_dist(deg2rad(data.predictedValue(r-1,1)),deg2rad(data.toPredict(r-1,1)));
%                     truePred(r,c)=deg2rad(data.toPredict(r-1,1));
%                     predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,1)),deg2rad(data.predictedValue(r-1,1)));
%                     predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,1)),deg2rad(predV(r-1,1)));
%                     surpriseTrial(r,c)=data.surprise(r-1,1);
%                 else
%                     predDeg(r,c)=data.predictedValue(r-1,2);
%                     pred(r,c)=deg2rad(data.predictedValue(r-1,2));
%                     predictErr(r,c)=circ_dist(deg2rad(data.predictedValue(r-1,2)),deg2rad(data.toPredict(r-1,2)));
%                     truePred(r,c)=deg2rad(data.toPredict(r-1,2));
%                     predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,2)),deg2rad(data.predictedValue(r-1,2)));
%                     predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,2)),deg2rad(predV(r-1,2)));
%                     surpriseTrial(r,c)=data.surprise(r-1,2);
%                 end
%             end
%         end
%     else
%         predVUpdated=nan(nTrial,2);
%         data.chosenTargPredict=nan(nTrial,2);
%         predDeg=nan(nTrial,2);
%         pred=nan(nTrial,2);
%         predictErr=nan(nTrial,2);
%         truePred=nan(nTrial,2);
%         predEstDiff=nan(nTrial,2);
%         predUpdate=nan(nTrial,2);
%         surpriseTrial=nan(nTrial,2);
%     end
%
%     blockData=data;
%
%     blockData = rmfield(blockData, 'prefs');
%     blockData = rmfield(blockData, 'totScore');
%     blockData = rmfield(blockData, 'initRT');
%     if data.prefs.doPredict
%         blockData = rmfield(blockData, 'predictedValue');
%         blockData = rmfield(blockData, 'toPredict');
%     end
%
%     blockData.block=ones(length(blockData.colorArray),1).*block;
%     blockData.trueEst=trueEst;
%     blockData.truePred=truePred;
%     blockData.est=est;
%     blockData.estErr=estErr;
%     blockData.pred=pred;
%     blockData.predictErr=predictErr; %objective
%     blockData.predEstDiff=predEstDiff; %subjective
%     blockData.predUpdate=predUpdate;
%     blockData.surpriseTrial=surpriseTrial;
%     blockData.predDeg=predDeg;
%     blockData.estDeg=estDeg;
%
%     blockData=straightStruct(blockData);
%
%     if ~exist('alldata')
%         alldata=blockData;
%     else
%         alldata=catBehav(blockData, alldata, true);
%     end
%
%     fn=fullfile(behaveDataDir,'subCombined',[subNum,'_allBlockData.mat']);
%     save(fn,'alldata');
% end
%
%% Copying one eye's data and using it for the other if only one eye

if size(dataEye,2) == 5
    dataEye(:,6:8) = 0;
    dataEye(:,8) = dataEye(:,5);
    dataEye(:,5:7) = dataEye(:,2:4);
    
end
%% Removing bad eye from eye data with bad eye
subno = str2double(subNum);
if ismember(subno,[2030,2063,3014])
    dataEye(:,2:4)=dataEye(:,5:7);
end
if ismember(subno,[2035,2057,2058,2062,2071,2083,2087,2098,2101,2103,2106,3018])
    dataEye(:,5:7)=dataEye(:,2:4);
end
%% Finding start of each block
%keyboard
% Figurign out at that timestamp each block begins
blocks = dataEye(dataEye(:,event) == block, 1);
blockChange = nan(4, 1); %The starts will be stored here, we have 4 blocks.but only saving for rand,cp,ob blocks
blockNumber = 1;
previous = blocks(1);
place = 1;
for i = 1:size(blocks, 1)
    if blocks(i) - previous ~= 1 || i == 1
        blockChange(place) = blocks(i);
        place = place + 1;
    end
    
    previous = blocks(i);
    
end
%% Getting data for rand and conditions only stim on locked

% Getting the rows that have stimulus information for all 3 blocks
stimulusData = dataEye((dataEye(:,event) == stimOn), :);
%random_indices = boolean((stimulusData(:,1) > blockChange(1)) & (stimulusData(:,1) < blockChange(2)));


rand_indices = (stimulusData(:,1)>=blockChange(1))&(stimulusData(:,1) < blockChange(2));
randEyeData=stimulusData(rand_indices, :);
%randEyeData=stimulusData(:, :);

if condData.condition(1)
    CP_indices = (stimulusData(:,1) >= blockChange(2)) & (stimulusData(:,1) < blockChange(3));
    OB_indices = stimulusData(:,1) >= blockChange(3) ;
else
    OB_indices = (stimulusData(:,1) >= blockChange(2)) & (stimulusData(:,1) < blockChange(3));
    CP_indices = stimulusData(:,1) >= blockChange(3) ;
end
% Getting the rows that have the stimulus information for changepoints and
% odball condition only

stimulusData = stimulusData(CP_indices | OB_indices, :);


%% Getting data for conditions only pred1 on locked

% % Getting the rows that have stimulus information for all 3 blocks
% pred1Data = dataEye((dataEye(:,event) == predCue1on), :);
% %random_indices = boolean((pred1Data(:,1) > blockChange(1)) & (pred1Data(:,1) < blockChange(2)));
% % NEED TO COUNTERBALANCE HERE
% CP_indices_pred = (pred1Data(:,1) >= blockChange(2)) & (pred1Data(:,1) < blockChange(3));
% OB_indices_pred = pred1Data(:,1) >= blockChange(3) ;
% 
% % Getting the rows that have the stimulus information for changepoints and
% % odball condition only
% 
% pred1Data = pred1Data(CP_indices_pred | OB_indices_pred, :);
% CP_indices_pred = (pred1Data(:,1) >= blockChange(2)) & (pred1Data(:,1) < blockChange(3));
% OB_indices_pred = pred1Data(:,1) >= blockChange(3) ;
%% Finding start of each stimulus

% We identify the start of each stimulus in stimStart
previous = stimulusData(1,1);
stimStart = nan(nTrials * 2, 1);
stimStart(1) = stimulusData(1,1);
stimCount = 2;

for i = 1:size(stimulusData, 1)
    %     Sometimes there is a 0 within the stimulus, we need to skip it
    if stimulusData(i,1) - previous > 10 && i ~= 1
        stimStart(stimCount) = stimulusData(i,1);
        stimCount = stimCount +  1;
    end
    previous = stimulusData(i,1);
end

%now do the same for rand block
previousRand = randEyeData(1,1);
stimStartRand = nan(nTrialsRand, 1);
stimStartRand(1) = randEyeData(1,1);
stimCount = 2;

for i = 1:size(randEyeData, 1)
    %     Sometimes there is a 0 within the stimulus, we need to skip it
    if randEyeData(i,1) - previousRand > 10 && i ~= 1
        stimStartRand(stimCount) = randEyeData(i,1);
        stimCount = stimCount +  1;
    end
    previousRand = randEyeData(i,1);
end

%% Finding end of each stimulus
% stimOff = 5;
% stimOffEvents  = dataEye(find(dataEye(:,event) == stimOff & dataEye(:,1) > stimStart(1)), 1);
% endStim = [stimOffEvents(1,1)];
% previous = stimOffEvents(1,1);
% for i = 1:size(stimOffEvents, 1)
%     %     Sometimes there is a 0 within the stimulus, we need to skip it
%     if stimOffEvents(i,1) - previous > 10 && i ~= 1
%         endStim = [endStim; stimOffEvents(i, 1)];
%     end
%     previous = stimOffEvents(i,1);
% end
% meanEndStim = mean(endStim - stimStart);

% %now do the same for rand block
% stimOff = 5;
% stimOffEventsRand  = dataEye(dataEye(:,event) == stimOff & dataEye(:,1) > stimStartRand(1)& dataEye(:,1) < stimStart(1));
% %stimOffEventsRand  = dataEye(dataEye(:,event) == stimOff );
% endStimRand=nan(nTrialsRand,1);
% endStimRand(1) = stimOffEventsRand(1,1);
% previousRand = stimOffEventsRand(1,1);
% stimCount=2;
% for i = 1:size(stimOffEventsRand, 1)
%     %     Sometimes there is a 0 within the stimulus, we need to skip it
%     if stimOffEventsRand(i,1) - previousRand > 10 && i ~= 1
%         endStimRand(stimCount) = stimOffEventsRand(i, 1);
%         stimCount=stimCount+1;
%     end
%     previousRand = stimOffEventsRand(i,1);
% end
% meanEndStimRand = mean(endStimRand - stimStartRand);




%% Finding start of each pred1

% We identify the start of each stimulus in pred1
% previous_pred = pred1Data(1,1);
% predStart = nan(nTrials * 2, 1);
% predStart(1) = pred1Data(1,1);
% predCount = 2;
% 
% for i = 1:size(pred1Data, 1)
%     %     Sometimes there is a 0 within the stimulus, we need to skip it
%     if pred1Data(i,1) - previous_pred > 1 && i ~= 1
%         predStart(predCount) = pred1Data(i,1);
%         predCount = predCount +  1;
%         %disp(i)
%     end
%     previous_pred = pred1Data(i,1);
% end

%% Finding start of each cue1
% cue1  = dataEye(find(dataEye(:,event) == predCue1on & dataEye(:,1) > stimStart(1)), 1);
% cue1Start = [cue1(1,1)];
% previous = cue1(1,1);
% for i = 1:size(cue1, 1)
%     %     Sometimes there is a 0 within the stimulus, we need to skip it
%     if cue1(i,1) - previous > 10 && i ~= 1
%         cue1Start = [cue1Start; cue1(i, 1)];
%     end
%     previous = cue1(i,1);
% end
% %cue1mean = mean(cue1Start - stimStart);
% predStart=cue1Start;

%% interpolate left and right seperately using raw continuous data


dataEye((dataEye(:,4)==0),4)=NaN;
dataEye((dataEye(:,7)==0),7)=NaN;

dataEye(isnan(dataEye(:,4)),7)=NaN;
dataEye(isnan(dataEye(:,7)),4)=NaN;

blinks = find(isnan(dataEye(:,4)));
isBlink=isnan(dataEye(:,4));

%create an array of logical blinks and put it into 5000x180 form to nan out
%the data before regression

for t=1:length(blinks)
    dataEye(max(blinks(t)-blinkWindow,1):min(blinks(t)+blinkWindow,length(dataEye)),4)= NaN;
    dataEye(max(blinks(t)-blinkWindow,1):min(blinks(t)+blinkWindow,length(dataEye)),7)= NaN;
    
end

[value,indices] = fillmissing(dataEye(:,4),'linear');
dataEye(indices,4)=value(indices);
[value,indices] = fillmissing(dataEye(:,7),'linear');
dataEye(indices,7)=value(indices);

meanArea=(dataEye(:,4)+dataEye(:,7))/2;

%keyboard
dataEyeOriginal=dataEye;

%% Taking a chunk of data for each stimulus


% Getting an interval of size stimSize corresponding to each of the
% stimuli. We need to get info even beyond the end of the stimuli because
% pupil needs time to change size.
%
% stimData = same columns as data + one column for stimNumber (9) + one column
% for unitsAfterSimulus (10)

% stimBefore = before;
% stimAfter = after;
% stimSize = stimBefore + stimAfter;
% stimData = nan(2*nTrials*stimSize, size(dataEye, 2) + 2);
% stimNumber = [1:nTrials, 1:nTrials];
% for i = 0: (2*nTrials-1)
%     keyboard
%     stimData(i*stimSize + 1:(i+1)*stimSize, 1:8) = dataEye(dataEye(:,1) >= stimStart(i+1) - stimBefore & ...
%                                     dataEye(:,1) < stimStart(i+1) + stimAfter, :);
%
%
%     stimData((i*stimSize + 1):(i+1)*stimSize, 9) = ones(stimSize, 1) .* stimNumber(i+1);
%     stimData((i*stimSize + 1):(i+1)*stimSize, 10) =  -stimBefore+1:stimAfter;  %1:stimSize;
% end
%keyboard

%using the interpolated raw data, select out averaged eye data and isBlink
stimBefore = before;
stimAfter = after;
stimSize = stimBefore + stimAfter;
dataEye(:,size(dataEye,2)+1)=meanArea;
dataEye(:,size(dataEye,2)+1)=isBlink;
stimData = nan(2*nTrials*stimSize, size(dataEye, 2) + 2);
stimNumber = [1:nTrials, 1:nTrials];
meanAreaMat=nan(stimSize,2*nTrials);
isBlinkMat=nan(stimSize,2*nTrials);
%keyboard
for i = 0: (2*nTrials-1)
    
    meanAreaMat(:,i+1)=dataEye(dataEye(:,1) >= stimStart(i+1) - stimBefore & ...
        dataEye(:,1) < stimStart(i+1) + stimAfter, 9);
    
    isBlinkMat(:,i+1)=dataEye(dataEye(:,1) >= stimStart(i+1) - stimBefore & ...
        dataEye(:,1) < stimStart(i+1) + stimAfter, 10);
    
    stimData(i*stimSize + 1:(i+1)*stimSize, 1:10) = dataEye(dataEye(:,1) >= stimStart(i+1) - stimBefore & ...
        dataEye(:,1) < stimStart(i+1) + stimAfter, :);
    
    
    stimData((i*stimSize + 1):(i+1)*stimSize, 11) = ones(stimSize, 1) .* stimNumber(i+1);
    stimData((i*stimSize + 1):(i+1)*stimSize, 12) =  -stimBefore+1:stimAfter;  %1:stimSize;
end
% figure()
% imagesc(meanAreaMat')
% figure()
% imagesc(isBlinkMat')

%take the mean of isBlinkMat on the time axis, if greater than 50%, is bad
%pupil trial,turn pupil (Y) data to nan before regression
trialQuality=mean(isBlinkMat);
isGoodTrial=(trialQuality<0.5);
isBadTrial=(trialQuality>0.5);
isGoodSub=((mean(isGoodTrial))>0.5);

meanAreaMat(:,isBadTrial)=nan;

%set a threshold for both participants and trials (less than 50% good
%trials)
if ~isGoodSub
    p=sprintf('subject %s has too many blink trials',subNum);
    disp(p);
    result=[];
    return
end







avgBaseline = nanmean(meanAreaMat(1:before,:));

zArea=nanzscore(meanAreaMat);

%calculate diameter from area

% A=pi*r^2
% r=sqrt(A/pi)
% d=2*sqrt(A/pi)

%check paper for usage of area vs diameter

diameterMat=2*sqrt(((meanAreaMat)./pi));
zDiam=nanzscore(diameterMat);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hello
%using the interpolated raw data, select out averaged eye data and isBlink
%for rand block
% stimBefore = before;
% stimAfter = after;
% stimSize = stimBefore + stimAfter;
% dataEyeRand=dataEyeOriginal;
% dataEyeRand(:,size(dataEyeOriginal,2)+1)=meanArea;
% dataEyeRand(:,size(dataEyeRand,2)+1)=isBlink;
% stimDataRand = nan(nTrialsRand*stimSize, size(dataEyeOriginal, 2) + 2);
% stimNumberRand = [1:nTrialsRand];
% meanAreaMatRand=nan(stimSize,nTrialsRand);
% isBlinkMatRand=nan(stimSize,nTrialsRand);
% %keyboard
% for i = 0: (nTrialsRand-1)
%     
%     
%     meanAreaMatRand(:,i+1)=dataEyeRand(dataEyeRand(:,1) >= stimStartRand(i+1) - stimBefore & ...
%         dataEyeRand(:,1) < stimStartRand(i+1) + stimAfter, 9);
%     
%     isBlinkMatRand(:,i+1)=dataEyeRand(dataEye(:,1) >= stimStartRand(i+1) - stimBefore & ...
%         dataEyeRand(:,1) < stimStartRand(i+1) + stimAfter, 10);
%     
%     stimDataRand(i*stimSize + 1:(i+1)*stimSize, 1:10) = dataEyeRand(dataEye(:,1) >= stimStartRand(i+1) - stimBefore & ...
%         dataEyeRand(:,1) < stimStartRand(i+1) + stimAfter, :);
%     
%     
%     stimDataRand((i*stimSize + 1):(i+1)*stimSize, 11) = ones(stimSize, 1) .* stimNumberRand(i+1);
%     stimDataRand((i*stimSize + 1):(i+1)*stimSize, 12) =  -stimBefore+1:stimAfter;  %1:stimSize;
%     
%     
% end
% 
% 
% %pupil trial,turn pupil (Y) data to nan before regression
% trialQualityRand=mean(isBlinkMatRand);
% isGoodTrialRand=(trialQualityRand<0.5);
% isBadTrialRand=(trialQualityRand>0.5);
% isGoodSubRand=((mean(isGoodTrialRand))>0.5);
% 
% %store this and look later in regression to see if isbadtrial is true 
% 
% meanAreaMatRand(:,isBadTrialRand)=nan;
% 
% %set a threshold for both participants and trials (less than 50% good
% %trials)
% % if ~isGoodSubRand
% %     p=sprintf('subject %s has too many blink trials',subNum);
% %     disp(p);
% %     result=[];
% %     return
% % end
% 
% zAreaRand=nanzscore(meanAreaMatRand);



%% Taking a chunk of data for each pred1
% if runPred
%     predBefore = before;
%     predAfter = after;
%     predSize =predBefore + predAfter;
%     predData = nan(2*(nTrials-1)*predSize, size(dataEye, 2) + 2);
%     predNumber = [1:nTrials-1, 1:nTrials-1];
%     for j = 0: (2*nTrials-3)
%         disp(j)
%         predData(j*predSize + 1:(j+1)*predSize, 1:8) = dataEye(dataEye(:,1) >= predStart(j+1) - predBefore & ...
%             dataEye(:,1) < predStart(j+1) + predAfter, :);
% 
% 
%         predData((j*predSize + 1):(j+1)*predSize, 9) = ones(predSize, 1) .* predNumber(j+1);
%         predData((j*predSize + 1):(j+1)*predSize, 10) =  -predBefore+1:predAfter;  %1:stimSize;
%     end
% end


%% check saccades during fixation period (-500-0ms)

    %from stim locked data,get the left right average of gaze position first
    %     gazeLRx = (stimData(:, gazeLx)+stimData(:, gazeRx))/2;
    %     gazeLRy = (stimData(:, gazeLy)+stimData(:, gazeRy))/2;
    %
    %     %interpolate gaze at zero position
    %     gazeLRx(gazeLRx==0)=NaN;
    %     gazeLRy(gazeLRy==0)=NaN;
    %
    %     [value,indices] = fillmissing(gazeLRx(:,1),'linear');
    %     gazeLRx(indices,1)=value(indices);
    %     [value,indices] = fillmissing(gazeLRy(:,1),'linear');
    %     gazeLRy(indices,1)=value(indices);
    %
    %
    %     %insert back to stimData matrix
    %     stimData(:,end+1)=gazeLRx;
    %     stimData(:,end+1)=gazeLRy;
    %
    %     %then reshape data so we have time x trial
    %     gazeLRxMat = reshape(gazeLRx, [stimSize, 2*nTrials]);
    %     gazeLRyMat = reshape(gazeLRy, [stimSize, 2*nTrials]);
    %
    %     %in these matrices, row stimBefore is stim start
    %     %fixation period starts at row: stimBefore-500 till stimBefore-1
    %     gazeLRxFix = gazeLRxMat(stimBefore-500:stimBefore-1,:);
    %     gazeLRyFix = gazeLRyMat(stimBefore-500:stimBefore-1,:);
    %
    %     %save the x y coordinates in the final results for plotting across all
    %     %subs later
    %
    %
    %
    %     %426c monitor resolution: screenX=1920 screenY=1080
    %
    %
    %
    % %     figure()
    % %     plot(gazeLRxFix(400:end,100),gazeLRyFix(400:end,100),'x')
    % %     xlim([0 1920])
    % %     ylim([0 1080])
    % %     set(gca,'Ydir','reverse')
    % %
    %     for tr=1:size(gazeLRxFix,2)
    %         trialNum=num2str(tr);
    %         %make a vid for left right average saccade over fixation period
    %         mkdir([baseDir,'/ET_data/',subNum])
    %         vidName=[baseDir,'/ET_data/',subNum,'/',sprintf('fixate%strial%s.mp4',subNum,trialNum)]
    %
    %         vidObj=VideoWriter(vidName,'MPEG-4'); %use the VideoWriter object to create a video file
    %         % input is file name,can be .avi or .mp4
    %
    %
    %         open(vidObj); %open the file to write you video data in
    %         if makeVid
    %         for ti=1:size(gazeLRxFix,1) %loop through the number of frames you want to get in the video,
    %             %here I want each frame to be the topoplot of each time point, so I'm
    %             %looping through time
    %
    %
    %             %now I'm just plotting whatever I want in this frame (iteration of the
    %             %loop, naming this plot whatever
    %             stimLockTime=-(501-ti);
    %
    %             currPlot=plot(gazeLRxFix(ti,tr),gazeLRyFix(ti,tr),'x');
    %             set(gca, 'XLim',[0 screenX])
    %             set(gca, 'YLim',[0 screenY])
    %             set(gca,'Ydir','reverse')
    %             set(get(gca,'title'),'string',sprintf('trial%s time%s',trialNum,num2str(stimLockTime)))
    %
    %             %then I capture this figure as my movieframe using the getframe
    %             %function, naming the frame whatever
    %             currFrame=getframe(gcf);
    %
    %
    %             %then I write video data (frame) to my movie file using writeVideo
    %             writeVideo(vidObj,currFrame);
    %
    %
    %             %then you need to clear your figure and delete before you move on (since you're gonna be using the
    %             %same plot variable name in the next iteration of loop
    %             clf
    %             delete(currPlot)
    %
    %
    %
    %         end
    %         close(vidObj)
    %         end
    %
    %
    %
    %
    %
    %     end
    %
    %optional for rand block
    %
    %     gazeLRx = (stimDataRand(:, gazeLx)+stimDataRand(:, gazeRx))/2;
    %     gazeLRy = (stimDataRand(:, gazeLy)+stimDataRand(:, gazeRy))/2;
%
 %   
%     gazeRxDat = stimDataRand(:, gazeRx);
%     gazeRyDat = stimDataRand(:, gazeRy);
%     
%     %interpolate gaze at zero position
%     gazeRxDat(gazeRxDat==0)=NaN;
%     gazeRyDat(gazeRyDat==0)=NaN;
%     
%     [value,indices] = fillmissing(gazeRxDat(:,1),'linear');
%     gazeRxDat(indices,1)=value(indices);
%     [value,indices] = fillmissing(gazeRyDat(:,1),'linear');
%     gazeRyDat(indices,1)=value(indices);
%     
%     gazeRxMat = reshape(gazeRxDat, [stimSize, nTrials]);
%     gazeRyMat = reshape(gazeRyDat, [stimSize, nTrials]);
%     
%     
%     
%     
%     gazeLxDat = stimDataRand(:, gazeLx);
%     gazeLyDat = stimDataRand(:, gazeLy);
%     
%     gazeLxDat(gazeLxDat==0)=NaN;
%     gazeLyDat(gazeLyDat==0)=NaN;
%     
%     [value,indices] = fillmissing(gazeLxDat(:,1),'linear');
%     gazeLxDat(indices,1)=value(indices);
%     [value,indices] = fillmissing(gazeLyDat(:,1),'linear');
%     gazeLyDat(indices,1)=value(indices);
%     
%     gazeLxMat = reshape(gazeLxDat, [stimSize, nTrials]);
%     gazeLyMat = reshape(gazeLyDat, [stimSize, nTrials]);
%     
%     gazeLRxMat=(gazeLxMat+gazeRxMat)/2;
%     gazeLRyMat=(gazeLyMat+gazeRyMat)/2;
%     
    % for non rand blocks
    %keyboard
    
    gazeRxDat = stimData(:, gazeRx);
    gazeRyDat = stimData(:, gazeRy);
    
    %interpolate gaze at zero position
    gazeRxDat(gazeRxDat==0)=NaN;
    gazeRyDat(gazeRyDat==0)=NaN;
    
    [value,indices] = fillmissing(gazeRxDat(:,1),'linear');
    gazeRxDat(indices,1)=value(indices);
    [value,indices] = fillmissing(gazeRyDat(:,1),'linear');
    gazeRyDat(indices,1)=value(indices);
    
    gazeRxMat = reshape(gazeRxDat, [stimSize, 2*nTrials]);
    gazeRyMat = reshape(gazeRyDat, [stimSize, 2*nTrials]);
    
    
    
    
    gazeLxDat = stimData(:, gazeLx);
    gazeLyDat = stimData(:, gazeLy);
    
    gazeLxDat(gazeLxDat==0)=NaN;
    gazeLyDat(gazeLyDat==0)=NaN;
    
    [value,indices] = fillmissing(gazeLxDat(:,1),'linear');
    gazeLxDat(indices,1)=value(indices);
    [value,indices] = fillmissing(gazeLyDat(:,1),'linear');
    gazeLyDat(indices,1)=value(indices);
    
    gazeLxMat = reshape(gazeLxDat, [stimSize, 2*nTrials]);
    gazeLyMat = reshape(gazeLyDat, [stimSize, 2*nTrials]);
    
    gazeLRxMat=(gazeLxMat+gazeRxMat)/2;
    gazeLRyMat=(gazeLyMat+gazeRyMat)/2;
    
    
    gazeLRxDat=[gazeLxDat,stimData(:, gazeLx),gazeRxDat,stimData(:, gazeRx)];
    
    
    
    %% check attentional bias during stim presentation (0-500ms)
    
    %in these matrices, row stimBefore is stim start
    %stim end is at stimBefore+500
    gazeLRxStimOn = gazeLRxMat(stimBefore:stimBefore+500-1,:);
    gazeLRyStimOn = gazeLRyMat(stimBefore:stimBefore+500-1,:);
    
    
    % gazeLxStimOn = mean(gazeLxMat(stimBefore:stimBefore+500-1,:),1);
    % gazeRxStimOn = mean(gazeRxMat(stimBefore:stimBefore+500-1,:),1);
   
    
    
    %save the x y coordinates in the final results for plotting across all
    %subs later
    
    
    
    %for two stims, if gaze towards left,score is negative |center-gaze|,then right side get the positive of the same score
    gazeAttention1=nan(size(gazeLRxMat,2),2);
    
    % awayEstErr=[];
    % awayEst=[];
    % towardsEstErr=[];
    % towardsEst=[];
    %keyboard
    for tr=1:size(gazeLRxMat,2)
        trialNum=num2str(tr);
        %make a vid for left right average saccade over stim on period
        %         mkdir([baseDir,'/ET_data/',subNum])
        %         vidName=[baseDir,'/ET_data/',subNum,'/',sprintf('stimOn%strial%s.mp4',subNum,trialNum)]
        %
        %         vidObj=VideoWriter(vidName,'MPEG-4'); %use the VideoWriter object to create a video file
        %         % input is file name,can be .avi or .mp4
        %
        %         if makeVid
        %         open(vidObj); %open the file to write you video data in
        %
        %         for ti=1:size(gazeLRxFix,1) %loop through the number of frames you want to get in the video,
        %             %here I want each frame to be the topoplot of each time point, so I'm
        %             %looping through time
        %
        %
        %             %now I'm just plotting whatever I want in this frame (iteration of the
        %             %loop, naming this plot whatever
        %             stimLockTime=ti-1;
        %
        %             currPlot=plot(gazeLRxStimOn(ti,tr),gazeLRyStimOn(ti,tr),'x');
        %             set(gca, 'XLim',[0 screenX])
        %             set(gca, 'YLim',[0 screenY])
        %             set(gca,'Ydir','reverse')
        %             set(get(gca,'title'),'string',sprintf('trial%s time%s',trialNum,num2str(stimLockTime)))
        %
        %             %then I capture this figure as my movieframe using the getframe
        %             %function, naming the frame whatever
        %             currFrame=getframe(gcf);
        %
        %
        %             %then I write video data (frame) to my movie file using writeVideo
        %             writeVideo(vidObj,currFrame);
        %
        %
        %             %then you need to clear your figure and delete before you move on (since you're gonna be using the
        %             %same plot variable name in the next iteration of loop
        %             clf
        %             delete(currPlot)
        %
        %         end
        %
        %         close(vidObj)
        %         end
        %here we need to come up with a way to get the gazePunishScore
        %per trial, average over the 500ms stimOn period?
        
        centerX=screenX*0.5;
        meanGazeX(tr,1)=mean(gazeLRxStimOn(:,tr));
        meanGazeY(tr,1)=mean(gazeLRyStimOn(:,tr));
        
        score=abs(centerX-meanGazeX(tr,1));
        if centerX>meanGazeX(tr,1)
            gazeAttention1(tr,1)=-score;
            gazeAttention1(tr,2)=score;
            % towardsEst=[towardsEst;condData.est(tr,1)];
            % towardsEstErr=[towardsEstErr;condData.estErr(tr,1)];
            % awayEst=[awayEst;condData.est(tr,2)];
            %  awayEstErr=[awayEstErr;condData.estErr(tr,2)];
        else
            gazeAttention1(tr,1)=score;
            gazeAttention1(tr,2)=-score;
            % towardsEst=[towardsEst;condData.est(tr,2)];
            % towardsEstErr=[towardsEstErr;condData.estErr(tr,2)];
            % awayEst=[awayEst;condData.est(tr,1)];
            %  awayEstErr=[awayEstErr;condData.estErr(tr,1)];
       
        end
        
    end
    
    %save this for behavioral regression
    gazeAttention=zscore(gazeAttention1,0,'all');
    
    % gazeLRxTrialAve=mean(gazeLRxStimOn,1);
    % 
    % gazeAttentionL=mean(zscore(gazeLRxStimOn,0,1),1);
    % gazeAttentionR=-(gazeAttentionL);
    % 
    % gazeAttention1=[gazeAttentionL;gazeAttentionR]';
    
    
    
    if plotGaze
%         figure()
%         hold on
%         plot(meanGazeX,meanGazeY,'x','markerSize',15)
%         xlim([0 screenX])
%         ylim([0 screenY])
%         plot(targ2Xpos+100,targ2Ypos,'s','markerSize',150,'lineWidth',5)
%         plot(targ1Xpos-100,targ2Ypos,'s','markerSize',150,'lineWidth',5)
%         xlabel('horizontal pixel','fontSize',50)
%         set(gca,'FontSize',40)
%         ylabel('vertical pixel','fontSize',50)
%         set(gca,'FontSize',40)
%         
%         
%         
%         plot(targ1Xpos,targ1Ypos,'s','markerSize',70,'lineWidth',5)
%         plot(targ2Xpos,targ2Ypos,'s','markerSize',70,'lineWidth',10)
%         
%         figure()
%         title('percpetual error distribution based on attention')
%         subplot(2,1,1)
%         title('percpetual error distribution based on attention')
%         hold on
%         histogram(towardsEstErr,15)
%         xlim([-1 1])
%         ylim([0 40])
%         legend({'attention towards target histogram'})
%         subplot(2,1,2)
%         hold on
%         histogram(awayEstErr,15)
%         xlim([-1 1])
%         ylim([0 40])
%         legend({'attention away from target histogram'})
        
        
        
    end

%% Preparing data for regressions

% b_kl = nan(90, stimSize);
% b_s = nan(90, stimSize);
% % for j = 1:30
% b_kl = nan(1,stimSize);
% b_s = nan(1,stimSize);
% Analysing the behavioral data to get kl and surprise



% sharedVar = shared_variables(DP);
% sharedVar.replace = true;
% dataCP = learningAndPerception(sharedVariables);
% dataOB = oddBall_Inference(sharedVariables);
% close all

% dataCP = modelCP(sharedVar);
% dataOB = modelOB(sharedVar);
%close all
% 
% dataCP2 = learningAndPerception(sharedVar,2);
% dataOB2 = oddBall_Inference(sharedVar,2);
% close all
% 
% klCP=[dataCP.klCP(1:nTrials);dataCP.klCP(nTrials+1:end)];
% klOB=[dataOB.klOB(1:nTrials);dataOB.klOB(nTrials+1:end)];
% 
% surpriseCP=[dataCP.surpriseCP(1:nTrials);dataCP.surpriseCP(nTrials+1:end)];
% surpriseOB=[dataOB.surpriseOB(1:nTrials);dataOB.surpriseOB(nTrials+1:end)];
% 
% condNumCP=ones(length(klCP),1);
% condNumOB=ones(length(klOB),1).*-1;
% 
% %
% if condData.condition(1,1)==1
%     condNum=[(condNumCP);(condNumOB)];
%     kl = [klCP;klOB];
%     surprise = [surpriseCP;surpriseOB];
% %     kl = [mean(klCP,1),mean(klOB,1)]';
% %     surprise = [mean(surpriseCP,1),mean(surpriseOB,1)]';
% 
% else
%     condNum=[(condNumOB);(condNumCP)];
%     kl = [klOB;klCP];
%     surprise = [surpriseOB;surpriseCP];
% %     kl = [mean(klOB,1),mean(klCP,1)]';
% %     surprise = [mean(surpriseOB,1),mean(surpriseCP,1)]';
% 
% end


%condNum=[(condNumCP);(condNumOB)];

% kl(1)=[];
% kl(50)=[];
%
% surprise(1)=[];
% surprise(50)=[];

% if ~preprocessBefore %we changed to new method to preprocess blink data before we get the average eye data
%     % Getting average area
%     % from stim locked data
%     rArea = (stimData(:, left_area)+stimData(:, right_area))/2;
%     %aveArea = stimData(:, 9);
% 
%     %from pred locked data
%     %rArea = (predData(:, left_area)+predData(:, right_area))/2;
%     %rArea = stimData(:, left_area);
% 
% 
%     rArea( rArea == 0) = NaN;
%     %rArea = reshape(rArea, [stimSize, 2*(nTrials-1)]);
%     rArea = reshape(rArea, [stimSize, 2*nTrials]);
%     rArea(rArea<1000)=NaN;
% 
%     %raw area value
%     rawArea = rArea;
%     % surpriseCond=zscore(surprise.*condNum);
%     % 
%     % [rho,pval]=corr(surpriseCond,kl);
%     % [rhoCP,pvalCP]=corr(surpriseCond(1:50),kl(1:50));
%     % [rhoOB,pvalOB]=corr(surpriseCond(51:end),kl(51:end));
%     % 
% 
% 
% 
%     %%%%%%%%%%% AREA %%%%%%%%%%%%
%     for trial = 1:2*(nTrials)
%         blinks = find(isnan(rArea(:,trial)));
%         if isempty(blinks)
%         else
%             for i = 1:size(blinks,1)
%                 rArea(max(blinks(i)-blinkWindow, 1):min(blinks(i)+blinkWindow, stimSize), trial) = NaN;
%             end
%         end
%     end
% 
%     %rArea=reshape(rArea,
% 
%     for i = 1:2*(nTrials)
%         [value,indices] = fillmissing(rArea(:,i),'linear','SamplePoints',1:stimSize);
%         rArea(indices,i) = value(indices);
%     end
% 
% 
%     zArea = reshape(nanzscore(reshape(rArea, [2*(nTrials)*stimSize, 1])), [stimSize, 2*(nTrials)]);
% 
% 
%     avgBaseline = mean(zArea(1:before,:));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %calculate diameter from raw area
% 
%     % A=pi*r^2
%     % r=sqrt(A/pi)
%     % d=2*sqrt(A/pi)
% 
%     %check paper for usage of area vs diameter
% 
%     diameter=2*sqrt(((rawArea)./pi));
% 
%     %diameter derivative
%     % matD=diameter(2:end,:)-diameter(1:end-1,:);
%     % zMat=zeros(1,2*nTrials);
%     % diffDiam=cat(1,zMat,matD);
% 
% 
% 
% 
%     %%%%%%%%%% DIAMETER %%%%%%%%%%%
%     for trial = 1:2*(nTrials)
%         blinks = find(isnan(diameter(:,trial)));
%         if isempty(blinks)
%         else
%             for i = 1:size(blinks,1)
%                 diameter(max(blinks(i)-blinkWindow, 1):min(blinks(i)+blinkWindow, stimSize), trial) = NaN;
%             end
%         end
%     end
% 
%     for i = 1:2*(nTrials)
%         [value,indices] = fillmissing(diameter(:,i),'linear','SamplePoints',1:stimSize);
%         diameter(indices,i) = value(indices);
%     end
%     nonzDiameter = reshape(reshape(diameter, [2*(nTrials)*stimSize, 1]), [stimSize, 2*(nTrials)]);
% 
%     zDiameter = reshape(nanzscore(reshape(diameter, [2*(nTrials)*stimSize, 1])), [stimSize, 2*(nTrials)]);
%     avgBaseline_d = mean(zDiameter(1:before,:));
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 
%     % gaussKernal=normpdf(-30:30,0,8);
%     % gaussKernal=gaussKernal./sum(gaussKernal);
%     % %plot(gaussKernal)
%     % 
%     % for i=1:size(diameter,2)
%     %     smoothData(:,i)=conv(diameter(:,i),gaussKernal,'same');
%     % end
%     % 
%     % 
%     % diffDiam=diff(smoothData);
%     % 
%     % zMat=zeros(1,2*(nTrials));
%     % diffDiam=cat(1,zMat,diffDiam);
%     % 
%     % 
%     % %%%%%%%%%%%%%% DIAMETER DERIVATIVE %%%%%%%%%%%%%%
%     % for trial = 1:2*(nTrials)
%     %     blinks = find(isnan(diffDiam(:,trial)));
%     %     if isempty(blinks)
%     %     else
%     %         for i = 1:size(blinks)
%     %             diffDiam(max(blinks(i)-blinkWindow, 1):min(blinks(i)+blinkWindow, stimSize), trial) = NaN;
%     %         end
%     %     end
%     % end
%     % 
%     % for i = 1:2*(nTrials)
%     %     [value,indices] = fillmissing(diffDiam(:,i),'linear','SamplePoints',1:stimSize);
%     %     diffDiam(indices,i) = value(indices);
%     % end
%     % 
%     % zDiffDiam=reshape(nanzscore(reshape(diffDiam, [2*(nTrials)*stimSize, 1])), [stimSize, 2*(nTrials)]);
%     % avgBaseline_diff = mean(zDiffDiam(1:before,:));
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% end
% 
% 
% % kl = [mean([dataCP.kl(1:nTrials); dataCP.kl(nTrials + 1:end)]),...
% %         mean([dataOB.kl(1:nTrials); dataOB.kl(nTrials + 1:end)])]';
% % surprise = [mean([dataCP.surprise(1:nTrials); dataCP.surprise(nTrials + 1:end)]),...
% %         mean([dataOB.surprise(1:nTrials); dataOB.surprise(nTrials + 1:end)])]';
% 
% %keyboard

%% show pupil has an equal response to the different colors used in this task

%divide colors into n colors sampled every x degrees
%keyboard

% if oneSubOnly
% for n=1:10
% 
%     x=360./n;
% 
%     colorBoundaries=0:x:360;
% 
%     isColor=nan(nTrialsRand,n);
% 
%     for trial=1:nTrialsRand
% 
%         trialC1=randData.presentedColor(trial,1);
%         trialC2=randData.presentedColor(trial,2);
% 
%         for c=1:n
% 
%             if (colorBoundaries(c+1)>trialC1 && trialC1>colorBoundaries(c)) || (colorBoundaries(c+1)>trialC2 && trialC2>colorBoundaries(c))
%                 isColor(trial,c)=1;
%             else
%                 isColor(trial,c)=0;
% 
%             end
%         end
% 
%     end
% 
%     %
%     colorRegressXes=[ones(nTrialsRand,1),isColor];
% 
%     %stimOnRandEye=meanAreaMatRand((before+1):(before+after),:);
%     stimOnRandEye=meanAreaMat((before+1):(before+after),:);
% 
%     eyeColorY=mean(stimOnRandEye,1)';
% 
% 
% 
%     %[h,p]=vartest2(colorRegressXes);
% 
% 
%     [bEyeColor, bintEyeColor,R,RINT,STATSeyeColor] = regress(zscore(eyeColorY),colorRegressXes);
% 
%     %second value in stats give f statistic
% 
%     fValues(n,1,1)=x;
%     fValues(n,2,1)=STATSeyeColor(2);
% 
%     sigV=icdf('f',0.95,n,size(eyeColorY,1));
% 
% end
% end
% 

%keyboard
% if checkColorPupil
% 
% colorsSin=sin(randData.presentedColor);
% colorsCos=cos(randData.presentedColor);
% 
% %colorPupilXes=[ones(length(colorsSin),1),colorsSin,colorsCos];
% 
% % for i = -stimBefore+1:stimAfter
% %     
% %     colorPupilY=zAreaRand(i+stimBefore,:)';
% %     %colorPupilY=zArea(i+stimBefore,:)';
% %     [bEyeColor, bintEyeColor,R,RINT,STATSeyeColor] = regress(colorPupilY,colorPupilXes);
% %     
% %     fValueSinCos(i+stimBefore)=STATSeyeColor(2);
% %     bEyeColorAll(:,i+stimBefore)=bEyeColor(:,1);
% %     %     bEyeColor(i+stimBefore)=bEyeColor;
% %     %     bintEyeColor(i+stimBefore,1)=bintEyeColor(;
% %     %
% %     
% % end
% 
% %we are adding sine and cosine into 
% %small effect and inconsistent across subjects 
% 
% %05/06/2022
% % test 3 subs with luminance 50, show binned color plot SEM across subs, in
% % RR say we're piloting new data with actual task with new luminance 
% % continue collectikng data only if we can reproduce previous results with
% % new luminance 
% %if we something totally different we wont collect more....
% 
% % we did see an effect in one subs, reviewer is right
% % we increased background luminance and 
% % can't exclude some particpant will have this effect, so we include the
% % sine and cosine in the model
% 
% %color difference just use beta bint in regress with one term of circdist(hypothesis is difference affect pupil
% %but we show it doens't)
% 
% %if oneSubOnly
% % figure()
% % plot(-stimBefore+1:stimAfter,fValueSinCos,'-')
% % hold on
% % ylabel('f values color effect on pupil (using sin and cos)')
% % xlabel('time')
% % xlim([-200 2000])
% % ylim([0 6])
% % sigV=icdf('f',0.95,5,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% 
% 
% 
% 
% windowOfInterest=zAreaRand(750:1300,:);
% meanWindowPupil=mean(windowOfInterest,1)';
% 
% 
% colorBins=[0:36:360]';
% colorsPerTrial=randData.colorArray;
% 
% pupilByColorBins=struct();
% 
% 
% for t=1:length(meanWindowPupil)
%     for b=1:length(colorBins)-1
% 
%         if colorsPerTrial(t,1)>colorBins(b)&&colorsPerTrial(t,1)<colorBins(b+1)
%             binName=num2str(colorBins(b));
%             fieldName=sprintf('binStart%s',binName);
%             if ~isfield(pupilByColorBins,fieldName)
%                 pupilByColorBins.(fieldName)=meanWindowPupil(t)
%             else
%                 pupilByColorBins.(fieldName)=[pupilByColorBins.(fieldName);meanWindowPupil(t)]
%             end
% 
%         end
%         if colorsPerTrial(t,2)>colorBins(b)&&colorsPerTrial(t,2)<colorBins(b+1)
%             binName=num2str(colorBins(b))
%             fieldName=sprintf('binStart%s',binName);
%             if ~isfield(pupilByColorBins,fieldName)
%                 pupilByColorBins.(fieldName)=meanWindowPupil(t)
%             else
%                 pupilByColorBins.(fieldName)=[pupilByColorBins.(fieldName);meanWindowPupil(t)]
%             end
% 
%         end
%     end
% end
% 
% 
% %figure 9 of rebuttal, binned color pupil effect 
% 
% for b=1:length(colorBins)-1
%     binName=num2str(colorBins(b));
%     fieldName=sprintf('binStart%s',binName);
% 
%     binValues=getfield(pupilByColorBins,fieldName);
%     binValuesAll(b,1:length(binValues))=binValues;
% 
% end
% 
% binValuesAll(binValuesAll==0)=nan;
% 
% binValuesAllMean=nanmean(binValuesAll,2);
% 
% binValuesAllStd=nanstd(binValuesAll,0,2)
% binValuesAllN=sum(~isnan(binValuesAll),2)
% binValuesAllSem=binValuesAllStd./sqrt(binValuesAllN)
% 
% bin95ConfInterval=binValuesAllSem.*norminv([.975], 0, 1)
% CI=[binValuesAllMean+bin95ConfInterval, binValuesAllMean-bin95ConfInterval];
% 
% 
% % figure()
% % bar(colorBins(1:end-1),binValuesAllMean)
% % hold on
% % plot([colorBins(1:end-1), colorBins(1:end-1)]', CI', '-k')
% % %errorbar(colorBins(1:end-1),binValuesAllMean,CI(:,1),CI(:,2),'k','linestyle','none');
% % xlabel('color angle','fontsize',40)
% % ylabel('average pupil size','fontsize',40)
% % set(gca,'fontsize',30)
% % 
% % %figure for difference between color space
% 
% absoluteColorDiff=abs(circ_dist(randData.colorArray(:,1),randData.colorArray(:,2)));
% 
% %run nested f tests reduced model with just baseline, rich has more terms
% 
% % Input
% % ~~~~~
% % y   column vector of dependent variable, call its length m.
% % Arich   m by p  matrix of independent variables. There are p regressors
% % Ared    m by q  matrix of independent variables. q < p
% % alpha   serves 3 purposes -- criterion for hypothesis testing for our
% %         model comparison and also to test for lack of fit of the richer
% %         model. Also serves as basis for CI on the regressors.
% % alpha_inc  the nested regression is only performed if the rich model is 
% %           achieves at least this level of significance.  default = 0.2
% % Output
% % ~~~~~~
% % h    1 if richer model offers significantly better fit. NaN if neither
% %      model achieves significance at p<=alpha
% % p    p value to reject H0: reduced model is just as good
% % F    F ratio used to compute probability
% % df   2-vector: [p-q, m-p]
% % brich   p-vector: lsq estimates for richer model
% % bred    q-vector: lsq estimates for reduced model
% % bintrich  p by 2 matrix containing CI for regressors, rich model
% % bintred  p by 2 matrix containing CI for regressors, reduced model
% % statrich  rsquared value, F and p for the rich model (alone)  (3-vector)
% % statred   rsquared value, F and p for the reduced model (alone)
% 
% 
% %find baseline -200ms-0ms
% baselinePupil=nanmean(zAreaRand(stimBefore-200:stimBefore-1,:),1)';
% 
% Ared=[ones(length(randData.colorArray),1),baselinePupil];
% ArichSinCos=[ones(length(randData.colorArray),1),baselinePupil,colorsSin,colorsCos];
% ArichDiff=[ones(length(randData.colorArray),1),baselinePupil,absoluteColorDiff];
% 
% alpha=0.05;
% alpha_inc=1;
% 
% 
% 
% 
% 
% % x=normrnd(0,1,100,1);
% % y=x+normrnd(0,1,100,1);
% % z=normrnd(0,1,100,1);
% % 
% % [h,p,F,df,brich,bred,bintrich,bintred,statrich,statred] = nestedF(y,[ones(100,1), x, z],[ones(100,1), x],.05,1);
% 
% 
% 
% clear h2_stored p2_stored F2_stored fRich2_stored fRed2_stored
% for i = -stimBefore+1:stimAfter
% 
% 
%     colorPupilY=zAreaRand(i+stimBefore,:)';
%     colorPupilY=colorPupilY(isGoodTrialRand);
% 
%     ArichSinCosReg=ArichSinCos;
%     ArichDiffReg=ArichDiff;
%     AredReg=Ared;
% 
%     ArichSinCosReg=ArichSinCosReg(isGoodTrialRand,:);
%     ArichDiffReg=ArichDiffReg(isGoodTrialRand,:);
%     AredReg=AredReg(isGoodTrialRand,:);
% 
%     %color difference
%     [bEyeColorD, bintEyeColorD,RD,RINTD,STATSeyeColorD] = regress(colorPupilY,ArichDiffReg);
% 
% 
%     fValueDiff(i+stimBefore)=STATSeyeColorD(2);
%     bEyeColorDiff(:,i+stimBefore)=bEyeColorD(3,1);
%     bintEyeColorDiff1(:,i+stimBefore)=bintEyeColorD(3,1);
%     bintEyeColorDiff2(:,i+stimBefore)=bintEyeColorD(3,2);
% 
%     %sin cosine
%     [h2,p2,F2,df2,brich2,bred2,bintrich2,bintred2,statrich2,statred2] = nestedF(colorPupilY,ArichSinCosReg,AredReg,alpha,alpha_inc);
%     h2_stored(i+stimBefore)=h2;
%     p2_stored(i+stimBefore)=p2;
%     F2_stored(i+stimBefore)=F2;
%     fRich2_stored(i+stimBefore)=statrich2(2);
%     fRed2_stored(i+stimBefore)=statred2(2);
% 
% end
% 
% 
% % figure()
% % %subplot(3,1,1)
% % plot(-stimBefore+1:stimAfter,F2_stored,'-','linewidth',2)
% % hold on
% % ylabel('f values color effect on pupil (using sin and cos)','fontsize',20)
% % xlabel('time','fontsize',20)
% % xlim([-400 2000])
% % title('luminance 50')
% % ylim([0 6])
% % sigV=icdf('f',0.95,6,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% % legend('sub3')
% % set(gca,'fontsize',25)
% % 
% % 
% % figure()
% % plot(-stimBefore+1:stimAfter,bEyeColorDiff(3,:),'-')
% % hold on
% % ylabel('beta values for color difference effect on pupil')
% % title('circular difference between two stim colors')
% % xlabel('time')
% % %ylim([0 6])
% % %sigV=icdf('f',0.95,3,sum(isGoodTrialRand));
% % yline(0,'--');
% % 
% % figure()
% % shadedErrorBar(-stimBefore+1:stimAfter,bEyeColorDiff,[bintEyeColorDiff1;bintEyeColorDiff2])
% % xlim([0 2000])
% % 
% 
% 
% % %xline(750)
% % %xline(1300)
% % 
% % 
% save('sub902F','F2_stored', 'fValueDiff','bEyeColorDiff','bintEyeColorDiff1','bintEyeColorDiff2') 
% % 
% % 
% end
% % 
% % figure()
% % shadedErrorBar(-stimBefore+1:stimAfter,bEyeColorDiff,[bintEyeColorDiff1;bintEyeColorDiff2])
% %
% % figure()
% % %subplot(3,1,1)
% % plot(-stimBefore+1:stimAfter,fRich2_stored,'-','linewidth',2)
% % hold on
% % ylabel('f values color effect on pupil (using sin and cos)','fontsize',20)
% % xlabel('time','fontsize',20)
% % xlim([-400 2000])
% % ylim([0 6])
% % sigV=icdf('f',0.95,6,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% % figure()
% % %subplot(3,1,2)
% % plot(-stimBefore+1:stimAfter,fRed2_stored,'-','linewidth',2)
% % hold on
% % ylabel('f values absolute color difference effect on pupil','fontsize',20)
% % xlabel('time','fontsize',20)
% % xlim([-400 2000])
% % %ylim([0 6])
% % sigV=icdf('f',0.95,3,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% % 
% % subplot(3,1,3)
% % plot(-stimBefore+1:stimAfter,fRed,'-','linewidth',2)
% % hold on
% % ylabel('f values baseline only','fontsize',20)
% % xlabel('time','fontsize',20)
% % xlim([-400 2000])
% % %ylim([0 6])
% % sigV=icdf('f',0.95,2,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% % 
% % 
% % figure()
% % plot(-stimBefore+1:stimAfter,bEyeColorAll(1,:))
% % hold on
% % xlim([-200 2000])
% % xline(750)
% % xline(1300)
% % xlabel('time')
% % ylabel('intercept')
% % 
% % figure()
% % plot(-stimBefore+1:stimAfter,mean(zAreaRand,2))
% % hold on
% % %xlim([-200 2000])
% % xline(750)
% % xline(1300)
% % xlabel('time')
% % ylabel('normalized pupil area averaged across trials')
% % 
% % figure()
% % %subplot(2,1,1)
% % plot(-stimBefore+1:stimAfter,fValueSinCos,'-','linewidth',2)
% % hold on
% % ylabel('f values color effect on pupil (using sin and cos)','fontsize',20)
% % xlabel('time','fontsize',20)
% % xlim([-400 2000])
% % ylim([0 6])
% % sigV=icdf('f',0.95,4,400);
% % yline(sigV,'--');
% % xline(750)
% % xline(1300)
% % set(gca,'fontsize',20)
% % 
% % 
% % %subplot(2,1,2)
% % plot(-stimBefore+1:stimAfter,mean(zAreaRand,2),'linewidth',2)
% % hold on
% % xlim([-400 2000])
% % xline(750)
% % xline(1300)
% % xlabel('time','fontsize',20)
% % ylabel('normalized pupil area averaged across trials','fontsize',20)
% % 
% % 
% % 
% % 
% % 
% % histogram(meanWindowPupil,10)
% % end
% %take average 
% %end
% %
% 
% %% Regression for both conditions together
% bes_intercept=nan(stimSize, 1);
% bes_kl = nan(stimSize, 1);
% bes_surprise = nan(stimSize, 1);
% bes_cond= nan(stimSize, 1);
% 
% bes_intercept_d=nan(stimSize, 1);
% bes_kl_d=nan(stimSize, 1);
% bes_surprise_d=nan(stimSize, 1);
% bes_cond_d= nan(stimSize, 1);
% 
% intercept_params=nan(stimSize, 2);
% kl_params = nan(stimSize, 2);
% surprise_params =nan(stimSize, 2);
% cond_params=nan(stimSize, 2);
% 
% intercept_params_d=nan(stimSize, 2);
% kl_params_d=nan(stimSize, 2);
% surprise_params_d=nan(stimSize, 2);
% cond_params_d=nan(stimSize, 2);
% 
% bes_intercept_diff=nan(stimSize, 1);
% bes_kl_diff=nan(stimSize, 1);
% bes_surprise_diff=nan(stimSize, 1);
% bes_cond_diff= nan(stimSize, 1);
% 
% intercept_params_diff=nan(stimSize, 2);
% kl_params_diff = nan(stimSize, 2);
% surprise_params_diff =nan(stimSize, 2);
% cond_params_diff=nan(stimSize, 2);
% 
% %xes = [ones(size(kl)), zscore(kl), zscore(surprise), avgBaseline'];
% 

% xes = [ones(size(kl)), zscore(surprise.*condNum), zscore(surprise), condNum,nanzscore(avgBaseline')];
% % xes_d = [ones(size(kl)), zscore(surprise.*condNum), zscore(surprise), condNum,nanzscore(avgBaseline_d')];
% % xes_diff = [ones(size(kl)), zscore(surprise.*condNum), zscore(surprise), condNum,nanzscore(avgBaseline_diff')];
% 
% 
% 
% 
% for i = -stimBefore+1:stimAfter
%     
%     [b, bint] = regress(zArea(i+stimBefore,:)', xes);
%     [b_d, bint_d] = regress(zDiam(i+stimBefore,:)', xes);
%     %     [b, bint] = regress(zArea(i+stimBefore,:)', xes);
%     %     [b_d, bint_d] = regress(zDiameter(i+stimBefore,:)', xes_d);
%     %     [b_diff, bint_diff] = regress(zDiffDiam(i+stimBefore,:)', xes_diff);
%     
%     intercept_params(i+stimBefore, 1) = bint(1,1);
%     intercept_params(i+stimBefore, 2) = bint(1,2);
%     learning_params(i+stimBefore, 1) = bint(2,1);
%     learning_params(i+stimBefore, 2) = bint(2,2);
%     st_params(i+stimBefore, 1) = bint(3, 1);
%     st_params(i+stimBefore, 2) = bint(3, 2);
%     cond_params(i+stimBefore, 1) = bint(4, 1);
%     cond_params(i+stimBefore, 2) = bint(4, 2);
%     bes_intercept(i+stimBefore) =b(1);
%     bes_learning(i+stimBefore) = b(2);
%     bes_st(i+stimBefore) = b(3);
%     bes_cond(i+stimBefore) = b(4);
%     
%     intercept_params_d(i+stimBefore, 1) = bint_d(1,1);
%     intercept_params_d(i+stimBefore, 2) = bint_d(1,2);
%     kl_params_d(i+stimBefore, 1) = bint_d(2,1);
%     kl_params_d(i+stimBefore, 2) = bint_d(2,2);
%     surprise_params_d(i+stimBefore, 1) = bint_d(3, 1);
%     surprise_params_d(i+stimBefore, 2) = bint_d(3, 2);
%     bes_intercept_d(i+stimBefore) =b_d(1);
%     bes_kl_d(i+stimBefore) =b_d(2);
%     bes_surprise_d(i+stimBefore) =b_d(3);
%     bes_cond_d(i+stimBefore) = b_d(4);
%     cond_params_d(i+stimBefore, 1) = bint_d(4, 1);
%     cond_params_d(i+stimBefore, 2) = bint_d(4, 2);
%     %
%     %     intercept_params_diff(i+stimBefore, 1) = bint_diff(1,1);
%     %     intercept_params_diff(i+stimBefore, 2) = bint_diff(1,2);
%     %     kl_params_diff(i+stimBefore, 1) = bint_diff(2,1);
%     %     kl_params_diff(i+stimBefore, 2) = bint_diff(2,2);
%     %     surprise_params_diff(i+stimBefore, 1) = bint_diff(3, 1);
%     %     surprise_params_diff(i+stimBefore, 2) = bint_diff(3, 2);
%     %     bes_intercept_diff(i+stimBefore) =b_diff(1);
%     %     bes_kl_diff(i+stimBefore) =b_diff(2);
%     %     bes_surprise_diff(i+stimBefore) =b_diff(3);
%     %     bes_cond_diff(i+stimBefore) = b_diff(4);
%     %     cond_params_diff(i+stimBefore, 1) = bint_diff(4, 1);
%     %     cond_params_diff(i+stimBefore, 2) = bint_diff(4, 2);
%     
% end
% % xes = [ones(size(kl)), zscore(kl), zscore(surprise)];
% %
% % for i = -stimBefore+1:stimAfter
% %     [b, bint] = regress(gradient(rArea(i+stimBefore,:))', xes);
% %     kl_params(i+stimBefore, 1) = bint(2,1);
% %     kl_params(i+stimBefore, 2) = bint(2,2);
% %     surprise_params(i+stimBefore, 1) = bint(3, 1);
% %     surprise_params(i+stimBefore, 2) = bint(3, 2);
% %     bes_kl(i+stimBefore) = b(2);
% %     bes_surprise(i+stimBefore) = b(3);
% % end
% 
% %% Regression seperately for each condition
% % kl = zscore(kl);
% % surprise = zscore(surprise);
% % cp_bes_kl = nan(stimSize, 1);
% % cp_bes_surprise = nan(stimSize, 1);
% % cp_kl_params = nan(stimSize, 2);
% % cp_surprise_params = nan(stimSize, 2);
% % cpArea = zArea(:,1:nTrials);
% % cpxes = [ones(nTrials, 1), zscore(kl(1:nTrials)), zscore(surprise(1:nTrials))];
% %
% % ob_bes_kl = nan(stimSize, 1);
% % ob_bes_surprise = nan(stimSize, 1);
% % ob_kl_params = nan(stimSize, 2);
% % ob_surprise_params = nan(stimSize, 2);
% % obArea = zArea(:, nTrials + 1:end);
% % obxes = [ones(nTrials, 1), zscore(kl(nTrials + 1:end)), zscore(surprise(nTrials + 1:end))];
% % for i = 1:stimSize
% %     [b, bint] = regress(cpArea(i,:)', cpxes);
% %     cp_bes_kl(i) = b(2);
% %     cp_bes_surprise(i) = b(3);
% %     cp_kl_params(i,1) = bint(2,1);
% %     cp_kl_params(i,2) = bint(2,2);
% %     cp_surprise_params(i,1) = bint(3,1);
% %     cp_surprise_params(i,2) = bint(3,2);
% %
% %     [b, bint] = regress(obArea(i,:)', obxes);
% %     ob_bes_kl(i) = b(2);
% %     ob_bes_surprise(i) = b(3);
% %     ob_kl_params(i,1) = bint(2,1);
% %     ob_kl_params(i,2) = bint(2,2);
% %     ob_surprise_params(i,1) = bint(3,1);
% %     ob_surprise_params(i,2) = bint(3,2);
% %
% % end
% 
% %% Average
% 
% 
% means = nan(stimSize, 1);
% means_d = nan(stimSize, 1);
% % means_diff = nan(stimSize, 1);
% 
% for i = -stimBefore+1:stimAfter
%     means(i + stimBefore) = nanmean(zArea(i+stimBefore,:));
%     means_d(i + stimBefore) = nanmean(zDiam(i+stimBefore,:));
%     %     means_diff(i + stimBefore) = nanmean(zDiffDiam(i+stimBefore,:));
%     %
% end
% 
% sterr =  std(zArea') / (nTrials*2)^(.5);
% sterr_d =  std(zDiam') / (nTrials*2)^(.5);
% % sterr_diff =  std(zDiffDiam') / (nTrials*2)^(.5);
% 
% %std_pupilArea=nanstd(meanAreaMat,0,'all');
% 
% 
% 
% %%  prepare data for trial number simulation power analysis
% 
% %define normal vs state transition trials based on model surprise
% %(threshold?)
% 
% thresh=0.7;
% trialType=surprise>thresh; %normal trials: 0
% trialType=trialType(condNum==-1);
% 
% %get mean eye data for all normal trials in OB block
% eyeNormalOB=zArea(stimBefore:stimBefore+1500,trialType==0);
% meanEyeNormalOB=nanmean(eyeNormalOB,'all');
% stdEyeNormalOB=nanstd(eyeNormalOB,0,'all');
% 
% 
% %get mean for OB trials
% eyeStOB=zArea(stimBefore:stimBefore+1500,trialType==1);
% meanEyeStOB=nanmean(eyeStOB,'all');
% stdEyeStOB=nanstd(eyeStOB,0,'all');
% 
% 
% %calculate d
% real_d=meanEyeStOB-meanEyeNormalOB;
% 
% %overall mean
% meanEyeAll=nanmean(zArea,'all');
% 
% 
% 
% 
% %%
% 
% 
% 
% % Ploting average pupil area at every time point
% 
% % figure;
% % hold on;
% % title('Pupil area', 'fontsize', 14);
% % imagesc(rawArea');
% % colorbar
% % plot([before, before], [1, 2*nTrials], '-.');
% % plot( [meanEndStim + before, meanEndStim+ before], [1, 2*nTrials], '-.');
% % plot([cue1mean+ before, cue1mean+ before],[1, 2*nTrials], '-.','color', 'black');
% % ylabel("Trial", 'fontsize', 14);
% % ylim([1,2*nTrials]);
% % xlabel("Milliseconds after start of trial", 'fontsize', 14);
% % legend({ 'Start of stimulus', 'End of stimulus', 'Start of cue 1'}...
% %         , 'Location', 'northeast');
% 
% 
% %% Both regressions plot (area)
% % defaultPlotParameters;
% % num        = 1;
% % wid        = 17; % total width
% % hts        = [6, 6, 6 ]; % height of each row
% % cols       = {1, 1, 1}; % width of columns
% % [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 2, [], ''); % etc
% % set(axs,'Units','normalized');
% % set(gcf,'color','white')
% % fig_.Position = [10 10 30 40];
% % left_color = [0 0 255] / 256;
% % right_color = [255 0 0] / 256;
% % set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % % draw in each panel, one at a time
% %
% %
% % lw=1.5
% % lw2=3
% % exSub=12
% %
% %
% %
% %
% %
% % for xx = 1:length(axs)
% %     axes(axs(xx)); hold on; cla(gca)
% %     if xx==1
% %
% %         hold on;
% %         meanAreaPlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_intercept ,abs(flip(intercept_params, 2)' - bes_intercept'));
% %         set(meanAreaPlot.mainLine,'HandleVisibility','off');
% %         set(meanAreaPlot.patch,'HandleVisibility','off');
% %         set(meanAreaPlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore, stimAfter], [0, 0], '--', 'HandleVisibility','off', 'color', 'black');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         legend({'Start of stimulus', 'End of stimulus', 'Prompted to report cue 1'}...
% %                 , 'Location', 'south');
% %         legend boxoff
% %         xlabel('Milliseconds after start of stimulus');
% %         ylabel('Average pupil area');
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %     elseif xx==2
% %
% %         hold on;
% %         surprisePlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_surprise ,abs(flip(surprise_params, 2)' - bes_surprise'));
% %         set(surprisePlot.patch,'HandleVisibility','off');
% %         set(surprisePlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore+1,stimAfter], [0,0], '--', 'HandleVisibility','off', 'color', 'black');
% %         xlabel('Milliseconds after start of stimulus');
% %         ylabel('Surprise Coefficient value');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %
% %     elseif xx==3
% %         hold on;
% %         klPlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_kl ,abs(flip(kl_params, 2)' - bes_kl'));
% %         set(klPlot.patch,'HandleVisibility','off');
% %         set(klPlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore+1,stimAfter], [0,0], '--', 'HandleVisibility','off', 'color', 'black');
% %         xlabel('Milliseconds after start of stimulus');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         ylabel('KL Coefficient value', 'fontsize', 14);
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %     end
% %
% %
% %     setPLOT_panelLabel(gca, xx);
% % end
% 
% %% Both regressions plot (diameter)
% 
% % defaultPlotParameters;
% % num        = 1;
% % wid        = 17; % total width
% % hts        = [6, 6, 6 ]; % height of each row
% % cols       = {1, 1, 1}; % width of columns
% % [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 2, [], ''); % etc
% % set(axs,'Units','normalized');
% % set(gcf,'color','white')
% % fig_.Position = [10 10 30 40];
% % left_color = [0 0 255] / 256;
% % right_color = [255 0 0] / 256;
% % set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % % draw in each panel, one at a time
% %
% %
% % lw=1.5
% % lw2=3
% % exSub=12
% %
% %
% %
% %
% %
% % for xx = 1:length(axs)
% %     axes(axs(xx)); hold on; cla(gca)
% %     if xx==1
% %
% %         hold on;
% %         meanAreaPlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_intercept_d ,abs(flip(intercept_params_d, 2)' - bes_intercept_d'));
% %         set(meanAreaPlot.mainLine,'HandleVisibility','off');
% %         set(meanAreaPlot.patch,'HandleVisibility','off');
% %         set(meanAreaPlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore, stimAfter], [0, 0], '--', 'HandleVisibility','off', 'color', 'black');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         legend({'Start of stimulus', 'End of stimulus', 'Prompted to report cue 1'}...
% %                 , 'Location', 'south');
% %         legend boxoff
% %         xlabel('Milliseconds after start of stimulus');
% %         ylabel('Average pupil area');
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %     elseif xx==2
% %
% %         hold on;
% %         surprisePlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_surprise_d ,abs(flip(surprise_params_d, 2)' - bes_surprise_d'));
% %         set(surprisePlot.patch,'HandleVisibility','off');
% %         set(surprisePlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore+1,stimAfter], [0,0], '--', 'HandleVisibility','off', 'color', 'black');
% %         xlabel('Milliseconds after start of stimulus');
% %         ylabel('Surprise Coefficient value');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %
% %     elseif xx==3
% %         hold on;
% %         klPlot = shadedErrorBar([-stimBefore+1:stimAfter],bes_kl_d ,abs(flip(kl_params_d, 2)' - bes_kl_d'));
% %         set(klPlot.patch,'HandleVisibility','off');
% %         set(klPlot.edge,'HandleVisibility','off');
% %         plot([-stimBefore+1,stimAfter], [0,0], '--', 'HandleVisibility','off', 'color', 'black');
% %         xlabel('Milliseconds after start of stimulus');
% %         plot([0, 0], [-1, 1], '--', 'color', 'r');
% %         plot([meanEndStim, meanEndStim], [-1,1], '--', 'color', 'blue');
% %         plot([cue1mean, cue1mean], [-1,1], '--', 'color', 'green');
% %         ylabel('KL Coefficient value', 'fontsize', 14);
% %         set(gca, 'box', 'off')
% %         xlim([-stimBefore, stimAfter]);
% %         ylim([-.5, .5]);
% %
% %     end
% %
% %
% %     setPLOT_panelLabel(gca, xx);
% % end
% %
% % % saveas(gcf,  "~/Dropbox (Brown)/arousalLearningPerception/registered report/figures/.fig", 'fig')
% % % % saveas(gcf,  'figure1new.eps', 'epsc2')
% % % close(gcf)
% %
% % %% Changepoint condition plots
% % % figure;
% % % hold on;
% % % plot(-stimBefore+1:stimAfter, cp_bes_kl, 'HandleVisibility','off');
% % % plot(-stimBefore+1:stimAfter, cp_bes_kl, '.', 'markersize', 12);
% % % plot(-stimBefore+1:stimAfter, cp_bes_surprise, 'HandleVisibility','off');
% % % plot(-stimBefore+1:stimAfter, cp_bes_surprise, '.', 'markersize', 12);
% % % xlabel('Time after start of stimulus', 'fontsize', 14);
% % % ylabel('Coefficient value', 'fontsize', 14);
% % % title('Coefficient value for multiple times after stimulus start', 'fontsize', 14);
% % % legend({'KL','Surprise'},'Location','northeast');
% % % title('Changepoint condition', 'fontsize', 18);
% %
% % % Confidence intervals
% % figure;
% % subplot(2,1,1);
% % title('Changepoint condition', 'fontsize', 18);
% % hold on;
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], cp_surprise_params', 'HandleVisibility','off');
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], cp_surprise_params', '.', 'HandleVisibility','off');
% % plot([-stimBefore+1,stimAfter], [0,0], '-', 'HandleVisibility','off');
% % xlabel('Time after start of stimulus', 'fontsize', 14);
% % ylabel('Surprise Coefficient value', 'fontsize', 14);
% % plot([0, 0], [-1, 1], '-.');
% % plot([meanEndStim, meanEndStim], [-1,1], '-.');
% % plot([cue1mean, cue1mean], [-1,1], '-.');
% % legend({ 'Start of stimulus', 'End of stimulus', 'Start of cue 1'}...
% %         , 'Location', 'northeast', 'Position',[0.2 0.6 0.1 0.2]);
% %
% % subplot(2,1,2);
% % hold on;
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], cp_kl_params', 'HandleVisibility','off');
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], cp_kl_params', '.', 'HandleVisibility','off');
% % plot([-stimBefore+1,stimAfter], [0,0], '-', 'HandleVisibility','off');
% % plot([0, 0], [-1, 1], '-.');
% % plot([meanEndStim, meanEndStim], [-1,1], '-.');
% % plot([cue1mean, cue1mean], [-1,1], '-.');
% % legend({ 'Start of stimulus', 'End of stimulus', 'Start of cue 1'}...
% %         , 'Location', 'northeast');
% % xlabel('Time after start of stimulus', 'fontsize', 14);
% % ylabel('KL Coefficient value', 'fontsize', 14);
% 
% %% Oddball condition plots
% 
% %
% % % Only exact parameter
% % % figure;
% % % hold on;
% % % plot(-stimBefore+1:stimAfter, ob_bes_kl, 'HandleVisibility','off');
% % % plot(-stimBefore+1:stimAfter, ob_bes_kl, '.', 'markersize', 12);
% % % plot(-stimBefore+1:stimAfter, ob_bes_surprise, 'HandleVisibility','off');
% % % plot(-stimBefore+1:stimAfter, ob_bes_surprise, '.', 'markersize', 12);
% % % xlabel('Time after start of stimulus', 'fontsize', 14);
% % % ylabel('Coefficient value', 'fontsize', 14);
% % % title('Coefficient value for multiple times after stimulus start', 'fontsize', 14);
% % % legend({'KL','Surprise'},'Location','northeast');
% % % title('Oddball condition', 'fontsize', 18);
% %
% % % Confidence intervals
% % figure;
% % subplot(2,1,1);
% % hold on;
% % title('Oddball condition', 'fontsize', 18);
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], ob_surprise_params',  'HandleVisibility','off');
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], ob_surprise_params', '.',  'HandleVisibility','off');
% % plot([-stimBefore+1,stimAfter], [0,0], '-', 'HandleVisibility','off');
% % xlabel('Time after start of stimulus', 'fontsize', 14);
% % ylabel('Surprise Coefficient value', 'fontsize', 14);
% % plot([0, 0], [-1, 1], '-.');
% % plot([meanEndStim, meanEndStim], [-1,1], '-.');
% % plot([cue1mean, cue1mean], [-1,1], '-.');
% % legend({ 'Start of stimulus', 'End of stimulus', 'Start of cue 1'}...
% %         , 'Location', 'northeast');
% % subplot(2,1,2);
% % hold on;
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], ob_kl_params',  'HandleVisibility','off');
% % plot([-stimBefore+1:stimAfter; -stimBefore+1:stimAfter], ob_kl_params', '.',  'HandleVisibility','off');
% % plot([-stimBefore+1,stimAfter], [0,0], '-', 'HandleVisibility','off');
% % plot([0, 0], [-1, 1], '-.');
% % plot([meanEndStim, meanEndStim], [-1,1], '-.');
% % plot([cue1mean, cue1mean], [-1,1], '-.');
% % legend({ 'Start of stimulus', 'End of stimulus', 'Start of cue 1'}...
% %         , 'Location', 'northeast');
% % xlabel('Time after start of stimulus', 'fontsize', 14);
% % ylabel('KL Coefficient value', 'fontsize', 14);
% 
% %% Returning results
% 
% result.nTrials = nTrials;
% 
% result.areaMatrix = zArea;
% result.learning_confidence = learning_params;
% result.st_confidence = st_params;
% result.bes_learning = bes_learning;
% result.bes_st = bes_st;
% result.bes_intercept=bes_intercept;
% result.intercept_confidence=intercept_params;
% result.bes_cond=bes_cond;
% result.cond_confidence=cond_params;
% result.means=means;
% result.sterr=sterr;
% 
% result.diameterMatrix = zDiam;
% result.learning_confidence_d = kl_params_d;
% result.st_confidence_d = surprise_params_d;
% result.bes_learning_d = bes_kl_d;
% result.bes_st_d = bes_surprise_d;
% result.bes_intercept_d = bes_intercept_d;
% result.intercept_confidence_d = intercept_params_d;
% result.bes_cond_d=bes_cond_d;
% result.cond_confidence_d=cond_params_d;
% result.means_d = means_d;
% result.sterr_d = sterr_d;
% %
% % result.diffMatrix = zDiffDiam;
% % result.kl_confidence_diff = kl_params_diff;
% % result.surprise_confidence_diff = surprise_params_diff;
% % result.bes_kl_diff = bes_kl_diff;
% % result.bes_surprise_diff = bes_surprise_diff;
% % result.bes_intercept_diff = bes_intercept_diff;
% % result.intercept_confidence_diff = intercept_params_diff;
% % result.bes_cond_diff=bes_cond_diff;
% % result.cond_confidence_diff=cond_params_diff;
% % result.means_diff = means_diff;
% % result.sterr_diff = sterr_diff;
% 
% result.stimBefore = before;
% result.stimAfter = after;
% result.meanEndStim=meanEndStim;
if plotGaze
    % result.gazeLRxStimOn=gazeLRxStimOn;
    % result.gazeLRyStimOn=gazeLRxStimOn;
    
%     result.gazeLRxFix=gazeLRxFix;
%     result.gazeLRyFix=gazeLRxFix;
%     
    % result.centerScreenX=screenX*0.5;
    % result.centerScreenY=screenY*0.5;
    
    result.gazeAttention=gazeAttention;
%     result.towardsEstErr=towardsEstErr;
% result.awayEstErr=awayEstErr;
end
% result.meanEyeNormalOB=meanEyeNormalOB;
% result.stdEyeNormalOB=stdEyeNormalOB;
% 
% result.meanEyeStOB=meanEyeStOB;
% result.stdEyeStOB=stdEyeStOB;
% 
% result.real_d=real_d;
% result.meanEyeAll=meanEyeAll;
% if checkColorPupil
% %result.fValues=fValues;
% %result.fValueSinCos=fValueSinCos;
% result.F2_stored=F2_stored;
% result.fValueDiff=fValueDiff;

% end


end

