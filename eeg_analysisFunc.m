function result=eeg_analysisFunc(subNum,saveText,dirs)
showPlots=true;
%set paths
baseDir = dirs.basePath;
sharedFuncDir =   dirs.smuPath;
eegDataDir = dirs.eegDir;
behaveDataDir = dirs.behaveDir;
modelDataDir = [dirs.behaveDir,'allModelData',saveText];
% if whichComp>2
%     datapath=fullfile(baseDir,'EEG_data'); %eeg data folder
% else
%     datapath=fullfile(baseDir,'eeg_data');
% end


% Add entire set of subdirectories to the path
addpath(genpath(baseDir))
addpath(genpath(sharedFuncDir))


% parameters
nTrials=120;
baseTime = 250;
stimOn=0;
stimOff=stimOn+200; % new paradigm stim on for 200ms
maskOff=stimOff+0; % new paradigm doesn't have mask
cue1on=maskOff+1800; % retention interval for 1800


cd(baseDir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     1) Load EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subno = str2num(subNum);


%load eeg data
eegDat = load(fullfile([eegDataDir,subNum,'_ALP_FILT_STIM.mat']));


% 1) LOAD BEHAVIORAL DATA


% behaveFiles=dir(fullfile(behaveDataDir, ['vwm_multiStim_color_', subNum, '*']));
% behavFNs={behaveFiles.name};


allData=load(fullfile(behaveDataDir,['allsubCombined/', subNum, '_allBlockData.mat']));
allData=allData.alldata;
behaveMat=sprintf('%s_3and4BlockData.mat',subNum);


%find end of practice blocks
practiceEndTrial=find(allData.block==3,1)-1;


condition=allData.condition(1,1);


%add nan values for the last trials of OB/CP blocks to make sure toPredict
%also has the same trial number as other fields in the allData struct
allData.toPredictOld=allData.toPredict;
allData.toPredict=nan(size(allData.chosenTargPredict));
allData.toPredict(1:practiceEndTrial)=allData.toPredictOld(1:practiceEndTrial);
allData.toPredict(practiceEndTrial+2:practiceEndTrial+nTrials,:)=allData.toPredictOld(practiceEndTrial+1:practiceEndTrial+nTrials-1,:);
allData.toPredict(practiceEndTrial+nTrials+2:end,:)=allData.toPredictOld(practiceEndTrial+nTrials:end,:);


%save a condition variable
condition = allData.condition(1);


%remove everything that doesn't have the same number of trials
allData=rmfield(allData,'toPredictOld');
allData=rmfield(allData,'allScore');
allData=rmfield(allData,'totScore');
allData=rmfield(allData,'sumScore');
allData=rmfield(allData,'condition');


%selecting data for just the last two blocks out
selInd=[practiceEndTrial+1:length(allData.toPredict)];
allData=selBehav(allData,selInd);


allDataStruct=allData;


dataPath=[modelDataDir,'/',subNum,'_allBlockData.mat'];


DP = fullfile(dataPath);
vars = shared_variables(DP);
vars.H = .15;


%getting learning/state transition probability from behavioral data


allModelData = load(fullfile(behaveDataDir,['allModelData',saveText,'/', subNum, '_allBlockData.mat'])); 
allModelData = allModelData.allDataStruct;


surpriseCP=(allModelData.surpriseCP);
surpriseOB=(allModelData.surpriseOB);


entropyOB=allModelData.entropyOB;
entropyCP=allModelData.entropyCP;


condCP=length(surpriseCP)/2;
condOB=length(surpriseOB)/2;


surpriseCP2col=[surpriseCP(1:condCP),surpriseCP(condCP+1:end)];
surpriseCP=mean(surpriseCP2col,2);
maxSurpriseCP=max(surpriseCP2col,[],2);
surpriseCPpoint5=mean(-abs(surpriseCP2col-0.5),2);
surpriseOB2col=[surpriseOB(1:condOB),surpriseOB(condOB+1:end)];
surpriseOB=mean(surpriseOB2col,2);
maxSurpriseOB=max(surpriseOB2col,[],2);
surpriseOBpoint5=mean(-abs(surpriseOB2col-0.5),2);


entropyCP2col=[entropyCP(1:condCP),entropyCP(condCP+1:end)];
entropyCP=mean(entropyCP2col,2);
maxEntropyCP=max(entropyCP2col,[],2);
entropyCPpoint5=mean(-abs(entropyCP2col-0.5),2);
entropyOB2col=[entropyOB(1:condOB),entropyOB(condOB+1:end)];
entropyOB=mean(entropyOB2col,2);
maxEntropyOB=max(entropyOB2col,[],2);
entropyOBpoint5=mean(-abs(entropyOB2col-0.5),2);




ambigCP = -abs(surpriseCP2col-0.33);
ambigCP = mean(ambigCP,2);
ambigOB = -abs(surpriseOB2col-0.33);
ambigOB = mean(ambigOB,2);


cpTrial=ones(nTrials,1);
obTrial=ones(nTrials,1).*-1;


if condition==1
    condNum=[cpTrial;obTrial];
    surprise=[(surpriseCP);(surpriseOB)];
    ambig=[ambigCP;ambigOB];
    maxSurprise=[maxSurpriseCP;maxSurpriseOB];
    maxEntropy=[maxEntropyCP;maxEntropyOB];
    surprisepoint5=[surpriseCPpoint5;surpriseOBpoint5];
    surprise2col=[surpriseCP2col;surpriseOB2col];
    entropy=[(entropyCP);(entropyOB)];
else
    condNum=[obTrial;cpTrial];
    surprise=[zscore(surpriseOB);zscore(surpriseCP)];
    ambig=[ambigOB;ambigCP];
    maxEntropy=[maxEntropyOB;maxEntropyCP];
    maxSurprise=[maxSurpriseOB;maxSurpriseCP];
    surprisepoint5=[surpriseOBpoint5;surpriseCPpoint5];
    surprise2col=[surpriseOB2col;surpriseCP2col];
    entropy=[zscore(entropyOB);zscore(entropyCP)];
end


allDataStruct.maxEntropy=maxEntropy;
allDataStruct.entropy=entropy;
allDataStruct.STP=surprise;
allDataStruct.ambig=ambig;
allDataStruct.maxSTP=maxSurprise;
allDataStruct.STPpoint5=surprisepoint5;
allDataStruct.STPMinus1 = [0;allDataStruct.STP(1:end-1,:)];
allDataStruct.STPL=surprise2col(:,1);
allDataStruct.STPR=surprise2col(:,2);
allDataStruct.condNum=condNum;


%add sin/cos of left/right color


colors = allDataStruct.colorArray;


allDataStruct.sinLeftColor = sin(deg2rad(colors(:,1)));
allDataStruct.sinRightColor = sin(deg2rad(colors(:,2)));
allDataStruct.cosLeftColor = cos(deg2rad(colors(:,1)));
allDataStruct.cosRightColor = cos(deg2rad(colors(:,2)));


  
% 2)
% Find trials where EEG data was good:
isGoodEEG=false(length(allDataStruct.chosenTargEst), 1);
subno = str2double(subNum);


if ismember(subno,[2046,2047,2063,2086,2071,2073,2090,2094,2099,2104])
    eegDat.epochNumbers = eegDat.epochNumbers - 1;
    disp(subno)
end


ind_OBCPstart=find(eegDat.epochNumbers>practiceEndTrial,1); %practice 5, random 20
%epoch_OBCP=eegDat.epochNumbers(ind_OBCPstart:end);


epoch_OBCP=eegDat.epochNumbers(ind_OBCPstart:end);


epoch_OBCP=(epoch_OBCP-practiceEndTrial)';




isGoodEEG(epoch_OBCP)=true;


firstTrials = zeros(length(epoch_OBCP),1);


if ismember(1,epoch_OBCP)
    firstTrials(1) = 1;
end


if ismember(nTrials+1,epoch_OBCP)
    firstTrials(find(epoch_OBCP == nTrials+1)) = 1;
end
goodData=selBehav(allDataStruct, epoch_OBCP);




% Run a subject level regression to see whether any electrode at any
% timepont has trial-to-trial responses that relate to
t1=find(eegDat.EEG.times>-4000, 1,'first');
t2=find(eegDat.EEG.times<4000, 1,'last');


%downSampData=eegDat.EEG.data(:, 1:sampNum:end,:);
downSampData=eegDat.EEG.data(:, t1:t2,ind_OBCPstart:end);
%downSampTimes=eegDat.EEG.times(1:sampNum:end);
downSampTimes=eegDat.EEG.times(t1:t2);


x0=ones(length(goodData.colorArray), 1);
x1=zscore(goodData.STP);
x2=zscore(goodData.STP).*goodData.condNum;
x3=zscore(goodData.condNum);
x4=zscore(goodData.entropy);
x1max=zscore(goodData.maxSTP);
x2max=zscore(goodData.maxSTP).*goodData.condNum;
x4max=zscore(goodData.maxEntropy);
xSins=zscore(goodData.sinLeftColor+goodData.sinRightColor);
xCoss=zscore(goodData.cosLeftColor+goodData.cosRightColor);


% control for baseline per channel
baselineData=eegDat.EEG.data(:,eegDat.EEG.times<0&eegDat.EEG.times>-baseTime,ind_OBCPstart:end);
channelBaseline=squeeze(mean(baselineData,2));


%channelBaseline(:,isnan(xPredErr)) = [];


%    noBaselineX=[x0,x1,x2,x3,x4,xSins,xCoss];
   noBaselineX=[x0,x1,x2,x3,x4];
%    noBaselineX=[x0,x1max,x2max,x3,x4max];


%preallocate matrix for each subject
mat(:,:,:)=nan(length({eegDat.EEG.chanlocs.labels}), length(downSampTimes), size(noBaselineX, 2)+1);


for chan=1:length({eegDat.EEG.chanlocs.labels})
    X=[noBaselineX,zscore(channelBaseline(chan,:)')];
    %X = [noBaselineX];
    for t=1:length(downSampTimes)
        Y=squeeze(downSampData(chan,t,:));            
        [coef]=regress(Y(~firstTrials),X(~firstTrials,:));
        %dumbData(t)=coef(1);
        mat(chan,t,:) = coef;
    end
   
end
    
%% save results useful for plot 


result.downSampTimes=downSampTimes;
result.mat=mat;
result.xes=noBaselineX;
result.baseline=channelBaseline;
