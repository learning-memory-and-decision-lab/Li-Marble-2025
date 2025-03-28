 
clear
%Functions Called: eeg_analysisFunc eyeRegressionFunc eyeRegressionPredResp computeLearningRate getEEGClusterSize getPupilClusterSize modelCP modelOB eyeAnalysis shared_variables
%Data Needed: behave data all blocks; behave data 3 and 4; eye data; eeg data; chanlocs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       %%%
%%%  SETUP INSTRUCTIONS:  %%%
%%%                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. you're probably going to need 70 GB worth of space on your computer to run this with every subject
% 1. set up some folder to put data and results in, and set that folder as the basePath
% 2. place this script and the functions "eyeRegressionFunc", "eyeRegressionPredResp", "eeg_analysisFunc", "modelOB", "modelCP", "getEEGClusterSize", "computeLearningRate", "eyeAnalysis", "shared_variables", and "chanlocs.mat" in the basePath folder
% 3. in basePath make subfolders for behaveData, eyeData and eegData, and set them as behaveDir, eyeDir, and eegDir respectively
% 4. download CLEANED eeg data (files titled XXXX_ALP_FILT_STIM.mat) from ALP_Summer2022\vwm_task_Harry\eegData\ and place in eegData subfolder
% 5. download eye data (files titled XXXX.mat) from ALP_Summer2022\vwm_task_Harry\ET_data\ and place in eyeData subfolder
% 6. inside the behaveData subfolder, create folders titled "allSubCombined" and "subCombined"
% 7. download the full behavioral data (files titled XXXX_allBlockData.mat) from ALP_Summer2022\vwm_task_Harry\behaveData\allSubCombined and place in allSubCombined subfolder
% 8. download the prediction phase's behavioral data (files titled XXXX_3and4BlockData.mat) from ALP_Summer2022\vwm_task_Harry\behaveData\subCombined and place in subCombined subfolder
% 9. create a final subfolder for figures that the script will generate figures and set that folder as figDir
% 10. download the version of sharedMatlabUtilities that exists in the ALP_Summer2022 folder and place that as a subfolder in basePath. Set smuPath to the helperFuncs location
% 11. review directories, it should look something like
%
%     basePath
%         *analysis functions go here*
%         behaveData
%             allSubCombined
%                 XXXX_allBlockData.mat
%             subCombined
%                 XXXX_3and4BlockData.mat
%         eyeData
%             XXXX.mat
%         eegData
%             XXXX_ALP_FILT_STIM.mat
%         figDir
%         
%
% 12. review parameters, note that the first time you run the script, runEEGRegression and runModel should both be true
% 13. verify that the EEGSubs, eyeSubs, and behaveSubs lists in the script match up with the data that has been downloaded to your computer
%     (script should stop and inform you if data isn't on your path)

%% Step 1: Set parameters/paths

whichComp = 1;

if whichComp == 1
    basePath =  '~/Documents/GitHub/Li-Marble-2025/analysis/';
    eyeDir =    '~/Documents/GitHub/Li-Marble-2025/data/ET_data/';
    behaveDir = '~/Documents/GitHub/Li-Marble-2025/data/behaveData/';
    smuPath =   '~/Documents/GitHub/Li-Marble-2025/analysis/helperFuncs/';
    subFunc =   '~/Documents/GitHub/Li-Marble-2025/analysis/subFunctions/';
    eegDir =    '~/Documents/GitHub/Li-Marble-2025/data/eegDataSmall/';
    figDir =    '~/Documents/GitHub/Li-Marble-2025/data/figures/generatedFigs/';
else
    
end

addpath(genpath(basePath))
addpath(genpath(subFunc))

dirs.basePath = basePath;
dirs.eyeDir = eyeDir;
dirs.behaveDir = behaveDir;
% dirs.smuPath = smuPath;
dirs.eegDir = eegDir;
dirs.figDir = figDir;

%colors for figures
cpColor = [246 146 30]./255;
obColor = [0 173 238]./255;
cpLColor = [251 200 143]./255;
obLColor = [128 214 247]./255;
cbColors=[0 0 0; 230 159 0; 86 180 233; 200 50 200; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167; 256 256 256]./256;

% Set parameters
nAllTrials = 300;
nPracticeTrials = 60;
sigThresh = 0.025;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       %%%
%%%  set subject IDs      %%%
%%%                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

realData=input('are you using the whole dataset?  Yes(1)/No(0):');


if realData
    EEGSubs = [20280 2030 2033 2034 2035 2036 2037 20380 2039 2040 2041 2042 2043 2044 2045 2046 2047 2048 2050 2053 2055 2056 2058 2060 2062 2063 2065 2066 ...
        2068 2069 2071 2073 2074 2075 2080 2081 2082 2083 2084 2085 2086 2087 2088 2090 2092 2093 2094 2095 2096 2097 2098 2099 2100 2101 2104 2105 2106];
    eyeSubs = [20280 2030 2033 2034 2035 2036 2037 20380 2039 2040 2041 2042 2043 2044 2045 2046 2047 2048 2050 2053 2055 2056 2057 2058 2060 2062 2063 2065 2066 2068 ...
        2069 2070 2071 2073 2074 2075 2079 2080 2081 2082 2083 2084 2086 2087 2088 2090 2092 2093 2094 2095 2096 2097 2098 2099 2100 2101 2103 2104 2105 2106];
    behaveSubs = [20280 2030 2031 2033 2034 2035 2036 2037 20380 2039 2040 2041 2042 2043 2044 2045 2046 2047 2048 2050 2053 2055 2056 2057 2058 2060 2062 2063 2065 ...
        2066 2068 2069 2070 2071 2073 2074 2075 2079 2080 2081 2082 2083 2084 2085 2086 2087 2088 2090 2092 2093 2094 2095 2096 2097 2098 2099 2100 2101 2102 2103 2104 2105 2106];
else
    %Use this data if running from github
    EEGSubs=[1001,1002,1003];
    eyeSubs=[1001,1002,1003];
    behaveSubs=[1001,1002,1003];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       %%%
%%%  Analysis parameters  %%%
%%%                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runModel = 0; %if 0, load previous behaveAll and reuses previous allModelDataCombScript
numModelReps = 40; %behavioral model repetitions, for testing, set this low (1-5) for real analysis, set to 401
doSTPResiduals = 0; %regresses stp out of trial-by-trial signal, learning, and bias
runEEGRegression = 0; %if 0, loads a previous b_mat_eeg
trialMeasure = 3; %3 is the right one (dot product of regressed baseline and STP map) but it takes a while (~40-120 min) to run, 2 doesn't regress baseline but is quicker
                  %need to do 3 if you're doing stp residuals
eegTimestepMode = 0; %if 0, looks at clusters, if 1 looks at bins of timepoints in specified channels
oddballFigNewVer = 1;
sigma_y = 31;

%Stimulus eye params
if realData
    timeBeforeEye = 1000;
    timeAfterEye = 4000;
    baselineTimeEye=1000;
    blinkWindow = 150;
else
    timeBeforeEye = 250;
    timeAfterEye = 1000;
    baselineTimeEye= 250;
    blinkWindow = 35;    
end

blinkThresh = 0.2;
numPermEye=2500; % permutation test number
leftArea = 4;
rightArea = 7;
clustThreshEye = 0.025;
combBlinkThresh = 0.5;

%EEG Params
timeBeforeEEG = 1999;
timeAfterEEG = 2000;
EEGTimes = -timeBeforeEEG:timeAfterEEG;
numPermEEG = 1000;
connectThresh=.40;
clustThreshEEG = 0.01;
recreateConnectionMat = 1;

%Prediction eye params
if realData
    timeBeforePred = 2000;
    timeAfterPred = 4000;
    baselineTimeStart = 2000;
    baselineTimeEnd   = 500;
else
    timeBeforePred = 500;
    timeAfterPred = 1000;
    baselineTimeStart = 500;
    baselineTimeEnd   = 125;
end
noPredTrials = [239,240,479,480];
clustThreshPred = 0.025;

%Behavioral Regression params
conditionTogether=1;
useModel = 0;
splitHalf = true;
nStart = 100;
nr = 1;
narrowWidth = .05;
betaWidth = narrowWidth;
LB = -1.5;
UB = 1.5;

lastTrials = [120,240];

% file names will contain sigma_y, date, and hour you ran script
% e.g. if sigma_y = 33 and you start the script at 5:29 on November 1, files would have the tag 'CombScript_33_11-01-24-17'
% you can also set this manually if you want to include more info about the run

saveText = ['CombScript_',num2str(sigma_y),'_',datestr(now,'mm-dd-yy-hh'),'_test'];
% saveText = ['CombScript_33_10-29-24-13'];

% if you want to load some results from a previous run but keep the new figures with the new time, don't change figTime
figTime = datestr(now,'mm-dd-yy-hh');

if runEEGRegression ~=1
    % if you're skipping the eeg regression, load a b_mat_eeg file of your choice
    load("b_mat_eegCombScript_31_03-17-25-20_test.mat")
end

if runModel ~=1
    % if you're skipping the behavioral regression, load results from a previous run and resave them with the new run's saveText tag so 
    % the script can load the right allModelData files later
    load('behaveAllCombScript_31_03-17-25-20_test.mat')

    for s=1:length(behaveSubs)
        subStr = num2str(behaveSubs(s));
        fileName=sprintf('behaveData/subCombined/%s_3and4BlockData.mat',subStr);
        DP = fullfile(basePath, fileName);
        vars = shared_variables(DP);
       
        allModelData = load(fullfile(behaveDir,['allModelDataCombScript_31_03-17-25-20_test/', subStr, '_allBlockData.mat']));
        allDataStruct = allModelData.allDataStruct;
    
        saveDir=[behaveDir,'allModelData',saveText,'/'];
        if s == 1
            mkdir(saveDir)
        end
        fn=fullfile(saveDir,[subStr,'_allBlockData.mat']);
        save(fn,'allDataStruct')
    end
end


% verify that all subject data (behavioral, eye, and EEG) exists on computer
for s = 1:length(behaveSubs)
    subStr = num2str(behaveSubs(s));
    if exist(['allSubCombined\',subStr,'_allBlockData.mat']) ~= 2
        disp(['Missing data for subject ',subStr,' in allSubCombined folder'])
        keyboard
    end
    if exist(['subCombined\',subStr,'_3and4BlockData.mat']) ~=2
        disp(['Missing data for subject ',subStr,' in subCombined folder'])
        keyboard
    end
end

for s = 1:length(eyeSubs)
    subStr = num2str(eyeSubs(s));
    if exist([subStr,'.mat']) ~= 2
        disp(['Missing data for subject ',subStr,' in eyeData folder'])
        keyboard
    end  
end

for s = 1:length(EEGSubs)
    subStr = num2str(EEGSubs(s));
    if exist([subStr,'_ALP_FILT_STIM.mat']) ~= 2
        disp(['Missing data for subject ',subStr,' in eegData folder'])
        keyboard
    end  
end
%% Behavioral Model stuff
%determine which subjects have eye/eeg data
noEyeSubs = behaveSubs(~ismember(behaveSubs,eyeSubs));
noEEGSubs = behaveSubs(~ismember(behaveSubs,EEGSubs));
missingOneSubs = unique([noEyeSubs,noEEGSubs]);

%if you don't want to run model, you can just load results of a previous run

if runModel == 1
behaveAll=[];   
for s=1:length(behaveSubs)  %loop through all subs
    tic
    subNum=behaveSubs(s);
    subStr=num2str(subNum);
    disp(subNum) 
    fileName=sprintf('subCombined/%s_3and4BlockData.mat',subStr);
    allBlockFileName = sprintf('allSubCombined/%s_allBlockData.mat',subStr);
     
    % load shared variables
    DP = fullfile(behaveDir, fileName);
    vars = shared_variables(DP);
    vars.H = .15;
    dataPath=DP;
     
    % load behavioral data
    data = load(vars.path);
    allBlockData = load([behaveDir,allBlockFileName],'alldata');
    allBlockData = allBlockData.alldata;
    allDataStruct = data.alldata;% (data.alldata.block == 3 | data.alldata.block == 4);
    allDataStruct = straightStruct(allDataStruct);
    alldata = allDataStruct;
    vars.sigma_y = sigma_y;
    disp(vars.sigma_y)
    nTrials = 120;
    % Steps:
    % 1 -- use colorArray in modeling script instead of presented color
    % 2 -- run modeling script lots of times, take "average" of everything
    %      you will use.
    % 3 -- create one structure to store data could be allDataStruct. Add
    %       fields for each variable you will extract from model. Plug in
    %        "averaged" model values to those fields for appropriate condition.
    % 4 -- rerun behavioral analyses with averaged model data. Make sure to
    %      do a check of whether new "surprise" measures make sense. For
    %      example plot update as a function of prediction error... plot
    %      trials with "high surprise" in a different color. Are they the
    %      big prediction error trials? I hope so.
    
    
    % Once we've submitted -- goal: streamlining code so you can hit play,
    % load arbitrary number of datasets, perform all analyses and generate
    % all figures.
    
    %run eye analysis for gazeAttention variable of bias regression
    if ~ismember(subNum,[2031,2085,2102,7777])
        eyePath=sprintf('%s.mat',subStr);
        eyePath=[eyeDir,eyePath];
        resultEye = eyeAnalysis(eyePath,timeBefore, timeAfter, blinkWindow,subStr,basePath,DP,smuPath); 
        if allDataStruct.condition(1)==1
            CPgazeAttention=[resultEye.gazeAttention(1:nTrials,1);resultEye.gazeAttention(1:nTrials,2)];
            OBgazeAttention=[resultEye.gazeAttention(nTrials+1:end,1);resultEye.gazeAttention(nTrials+1:end,2)];
        else
            OBgazeAttention=[resultEye.gazeAttention(1:nTrials,1);resultEye.gazeAttention(1:nTrials,2)];
            CPgazeAttention=[resultEye.gazeAttention(nTrials+1:end,1);resultEye.gazeAttention(nTrials+1:end,2)];
        end
    else
        CPgazeAttention=zeros(2*nTrials,1);
        OBgazeAttention=zeros(2*nTrials,1);
    end

    %preallocate variables for model loop
    allCPEntropy = [];
    allCPSurprise = [];
    allCPobjPE = [];
    allCPsubPE = [];
    allCPestErr = [];
    allCPmodelPred = [];

    allOBEntropy = [];
    allOBSurprise =[];
    allOBobjPE =[];
    allOBsubPE =[];
    allOBestErr = [];
    allOBmodelPred = [];

    for i = 1:numModelReps
        % Get surprise and entropy for changepoint trials:
        dataCP = modelCP(vars);
        allCPEntropy=[allCPEntropy; dataCP.entropyCP'];
        allCPSurprise=[allCPSurprise; dataCP.surpriseCP'];
        allCPobjPE=[allCPobjPE;dataCP.predictionErrorOnX'];
        allCPestErr=[allCPestErr;dataCP.perceptualErrorOnX'];
        allCPsubPE=[allCPsubPE;dataCP.predErrorCP'];
        allCPmodelPred=[allCPmodelPred;dataCP.maxLikePostMu'];

        % Get surprise and entropy for oddball trials:
        dataOB = modelOB(vars);
        allOBEntropy=[allOBEntropy; dataOB.entropyOB'];
        allOBSurprise=[allOBSurprise; dataOB.surpriseOB'];
        allOBobjPE=[allOBobjPE;dataOB.predictionErrorOnB'];
        allOBestErr=[allOBestErr;dataOB.perceptualErrorOnB'];
        allOBsubPE=[allOBsubPE;dataOB.predErrorOB'];
        allOBmodelPred=[allOBmodelPred;dataOB.maxLikePostC'];

        if mod(i,10)==0
            disp(i)
        end
    end

    predictions = allDataStruct.pred;
    outcomes = (allDataStruct.est);
    newBlock = 121;
    
    %run CLR function (done twice, one for left and right stimulus
    [LR1,UP1,subPE1] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polarHalfCorrect');
    [LR2,UP2,subPE2] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polarHalfCorrect');
    adjPE = [subPE1,subPE2];
    UP = [UP1,UP2;nan,nan];
    LR = [LR1,LR2;nan,nan];
    bias = allDataStruct.estErr(:,:)./allDataStruct.predictErr(:,:);

    % order learning and bias based on trial order
    if allDataStruct.condition(1) == 1
        LRCP = [LR(1:nTrials,1);LR(1:nTrials,2)];
        LROB = [LR(nTrials+1:end,1);LR(nTrials+1:end,2)];
        biasCP = [bias(1:nTrials,1);bias(1:nTrials,2)];
        biasOB = [bias(nTrials+1:end,1);bias(nTrials+1:end,2)];
    else
        LRCP = [LR(nTrials+1:end,1);LR(nTrials+1:end,2)];
        LROB = [LR(1:nTrials,1);LR(1:nTrials,2)];
        biasCP = [bias(1:nTrials,1);bias(1:nTrials,2)];
        biasOB = [bias(nTrials+1:end,1);bias(nTrials+1:end,2)];
    end

    %save variables in allDataStruct
    allDataStruct.entropyCP=mean(allCPEntropy, 1);
    allDataStruct.surpriseCP=mean(allCPSurprise,1);
    allDataStruct.entropyOB=mean(allOBEntropy,1);
    allDataStruct.surpriseOB=mean(allOBSurprise,1);
    allDataStruct.subNum=nan(size(allDataStruct.isRand));
    allDataStruct.subNum(:)=subNum;
    allDataStruct.blockCond=nan(size(allDataStruct.isRand));
    allDataStruct.OBgazeAttention = OBgazeAttention;
    allDataStruct.CPgazeAttention = CPgazeAttention;
    allDataStruct.perceptualErrorOnX=circ_mean(allCPestErr);
    allDataStruct.perceptualErrorOnB=circ_mean(allOBestErr);
    allDataStruct.predictionErrorOnX=circ_mean(allCPobjPE);
    allDataStruct.predictionErrorOnB=circ_mean(allOBobjPE);
    allDataStruct.subPredErrorCP=circ_mean(allCPsubPE);
    allDataStruct.subPredErrorOB=circ_mean(allOBsubPE);
    allDataStruct.maxLikePostX=dataCP.maxLikePostX;
    allDataStruct.maxLikePostMu=circ_mean(deg2rad(allCPmodelPred));
    allDataStruct.maxLikePostB=dataOB.maxLikePostB;
    allDataStruct.maxLikePostC=circ_mean(deg2rad(allOBmodelPred));
    allDataStruct.LRCP = LRCP;
    allDataStruct.LROB = LROB;
    allDataStruct.biasCP = biasCP;
    allDataStruct.biasOB = biasOB;  
    allDataStruct = straightStruct(allDataStruct);

    % add blockCond (tells you what trials were in which condition) variable to aDS
    t=length(allDataStruct.isRand);
    if allDataStruct.condition(1)==1
        allDataStruct.blockCond(1:t/2)=1;
        allDataStruct.blockCond(t/2+1:end)=-1;
    else
        allDataStruct.blockCond(1:t/2)=-1;
        allDataStruct.blockCond(t/2+1:end)=1;
    end

    %Save subject's model calculated surprise,entropy,etc.
    saveDir=[behaveDir,'allModelData',saveText,'/'];
    if s == 1
        mkdir(saveDir)
    end
    fn=fullfile(saveDir,[subStr,'_allBlockData.mat']);
    save(fn,'allDataStruct')

    if isempty(behaveAll)
        behaveAll=allDataStruct;
    elseif  exist('behaveAll')&& ~isempty(behaveAll)
        behaveAll=catBehav(allDataStruct,behaveAll);
    end
    toc
end
beep
end
% straighten behave all structure and add not move (variable of distance between estimation and start point)
behaveAll=straightStruct(behaveAll);
behaveAll.notMove=circ_dist(deg2rad(behaveAll.est),deg2rad(behaveAll.startPoint));

%save all subjects' model calculated surprise,entropy,etc.
saveDir=[basePath];
fn=fullfile(saveDir,['behaveAll',saveText,'.mat']);
save(fn,'behaveAll')


% subStr='2035';
% fileName=sprintf('behaveData/subCombined/%s_3and4BlockData.mat',subStr);
% DP = fullfile(basePath, fileName);
% vars = shared_variables(DP);
%% Raw Behavior
% this code saves results for panels B and C of figures 3 and 4
% preallocate variables for loop (specify allerrors/updates variables)
oddballPredErrorAll = [];
oddballPredUpdateAll = [];
changepointPredErrorAll = [];
changepointPredUpdateAll = [];
oddballEstErrAll = [];
oddballobjPredErrorAll = [];
changepointEstErrAll = [];
changepointobjPredErrorAll = [];
allCPSurprise = [];
allOBSurprise = [];
nBlockTrials = 120;
nTrials = 240;

%specify which trials have no bias(first trial) and no learning (first & last
nanTrialsBias = [1,121];
nanTrialsLearning = [1,120,121,240];
for subno = behaveSubs
    disp(subno)
    subnoStr = num2str(subno);
    s = find(behaveSubs==subno);

    %load behavior data
    behaveMat = sprintf('allsubCombined/%s_allBlockData.mat',subnoStr);
    file2 = fullfile(behaveDir,behaveMat);
    load(file2)

    allDataStruct = alldata;
    
    %load model generated parameters
    allModelData = load(fullfile(behaveDir,['allModelData',saveText,'/', subnoStr, '_allBlockData.mat']));
    allModelData = allModelData.allDataStruct;

    surpriseCP = reshape(allModelData.surpriseCP',[nBlockTrials,2]);
    surpriseOB = reshape(allModelData.surpriseOB',[nBlockTrials,2]);
    
    if allDataStruct.condition(1) == 1
        surprise = [surpriseCP;surpriseOB];
    else
        surprise = [surpriseOB;surpriseCP];
    end

    %separate reproduction error for blcoks 3 and 4
    estErr3 = [allDataStruct.estErr(61:180,1);allDataStruct.estErr(61:180,2)];
    estErr4 = [allDataStruct.estErr(181:300,1);allDataStruct.estErr(181:300,2)];

    %separate subjective prediciton error for blcoks 3 and 4
    subPredError3 = [allDataStruct.subPredErr(61:180,1);allDataStruct.subPredErr(61:180,2)];
    subPredError4 = [allDataStruct.subPredErr(181:300,1);allDataStruct.subPredErr(181:300,2)];
    
    %separate objective prediction error for blcoks 3 and 4    
    objPredError3 = [allDataStruct.predictErr(61:180,1);allDataStruct.predictErr(61:180,2)];
    objPredError4 = [allDataStruct.predictErr(181:300,1);allDataStruct.predictErr(181:300,2)];

    %load shortened data (only blocks 3 and 4) 
    behaveMat = sprintf('subCombined/%s_3and4BlockData.mat',subnoStr);
    file2 = fullfile(behaveDir,behaveMat);
    alldataShort = load(file2);

    %define predictions and outcomes(estimations for subjective)
    predictions = alldataShort.alldata.pred;
    outcomes = (alldataShort.alldata.est);
    newBlock = nBlockTrials + 1;
    
    %run CLR function (done twice, one for left and right stimulus
    [LR1,UP1,~] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polarHalfCorrect');
    [LR2,UP2,~] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polarHalfCorrect');
    
    %concatenate LR UP and PE
    LR3 = [LR1(2:nBlockTrials-1);LR2(2:nBlockTrials-1)];
    LR4 = [LR1(nBlockTrials+2:end);LR2(nBlockTrials+2:end)];
    UP3 = [UP1(2:nBlockTrials-1);UP2(2:nBlockTrials-1)];
    UP4 = [UP1(nBlockTrials+2:end);UP2(nBlockTrials+2:end)];

    learningTrials = true(nTrials,1);
    biasTrials = true(nTrials,1);
    learningTrials(nanTrialsLearning) = false;
    biasTrials(nanTrialsBias) = false;

    %specify which block of data is CP/OB
    if allDataStruct.condition(1) == 1
        oddballPredError = subPredError4(learningTrials);
        oddballPredUpdate = UP4;
        oddballobjPredError = objPredError4(biasTrials);
        oddballEstErr = estErr4(biasTrials);
    
        changepointPredError = subPredError3(learningTrials);
        changepointPredUpdate = UP3;
        changepointobjPredError = objPredError3(biasTrials);
        changepointEstErr = estErr3(biasTrials);
    else
        oddballPredError = subPredError3(learningTrials);
        oddballPredUpdate = UP3;
        oddballobjPredError = objPredError3(biasTrials);
        oddballEstErr = estErr3(biasTrials);
    
        changepointPredError = subPredError4(learningTrials);
        changepointPredUpdate = UP4;
        changepointobjPredError = objPredError4(biasTrials);
        changepointEstErr = estErr4(biasTrials);
    end

    %add subject's data to all data variables
    oddballPredErrorAll = [oddballPredErrorAll;oddballPredError];
    oddballPredUpdateAll = [oddballPredUpdateAll;oddballPredUpdate];
    changepointPredErrorAll = [changepointPredErrorAll;changepointPredError];
    changepointPredUpdateAll = [changepointPredUpdateAll;changepointPredUpdate];
    oddballEstErrAll = [oddballEstErrAll;oddballEstErr];
    oddballobjPredErrorAll = [oddballobjPredErrorAll;oddballobjPredError];
    changepointEstErrAll = [changepointEstErrAll;changepointEstErr];
    changepointobjPredErrorAll = [changepointobjPredErrorAll;changepointobjPredError];
    subsEstErr(s) = mean(abs(allDataStruct.estErr),'all');

    %calculate average pred err on non surprise trials
    subsObjPredErr(s) = nanmean(abs(allDataStruct.subPredErr(surprise<0.25)),'all');
    subsSubPredErr(s) = nanmean(abs(allDataStruct.predictErr(surprise<0.25)),'all');    
%     subsObjPredErr(s) = nanmean(abs(allDataStruct.subPredErr),'all');
%     subsSubPredErr(s) = nanmean(abs(allDataStruct.predictErr),'all');
end

%calculate median est/pred err for "good subs" plot threshold, and find list of good subs
medianEstErr=median(subsEstErr);
goodEstSubs = subsEstErr<medianEstErr;

medianObjPredErr=median(subsObjPredErr);
goodPredSubs = subsObjPredErr<medianObjPredErr;

goodEstSel = reshape(repmat(goodEstSubs,nTrials-2,1),[],1);
goodPredSel = reshape(repmat(goodPredSubs,nTrials-4,1),[],1);

% steps for these loops
% 1 - preallocate quantile borders
% 2 - add trials in a quantile to that quantile
% 3 - find mean error/update for that quantile
% 4 - repeat until quantiles are filled
% 5  - repeat whole process for bias and learning for all and good subjects
nBins = 40;
bordersAllCP = quantile(changepointobjPredErrorAll,nBins-1);
bordersAllOB = quantile(oddballobjPredErrorAll,nBins-1);
quantilesAllOB = zeros(size(oddballobjPredErrorAll,1),1);
quantilesAllCP = zeros(size(changepointobjPredErrorAll,1),1);
for q = nBins:-1:1
    if q<nBins
        quantilesAllCP(changepointobjPredErrorAll<=bordersAllCP(q)) = q;
        quantilesAllOB(oddballobjPredErrorAll<=bordersAllOB(q)) = q;
    else
        quantilesAllCP(changepointobjPredErrorAll>bordersAllCP(q-1)) = q;
        quantilesAllOB(oddballobjPredErrorAll>bordersAllOB(q-1)) = q;
    end
end
for q = 1:nBins
    meanQuantilesCPAllErrBias(q) = mean(changepointEstErrAll(quantilesAllCP==q));
    meanQuantilesOBAllErrBias(q) = mean(oddballEstErrAll(quantilesAllOB==q));
    quantileAllCPxesBias(q) = mean(changepointobjPredErrorAll(quantilesAllCP==q));
    quantileAllOBxesBias(q) = mean(oddballobjPredErrorAll(quantilesAllOB==q));
end

bordersAllCP = quantile(changepointPredErrorAll,nBins-1);
bordersAllOB = quantile(oddballPredErrorAll,nBins-1);
quantilesAllOB = zeros(size(oddballPredErrorAll,1),1);
quantilesAllCP = zeros(size(changepointPredErrorAll,1),1);
for q = nBins:-1:1
    if q<nBins
        quantilesAllCP(changepointPredErrorAll<=bordersAllCP(q)) = q;
        quantilesAllOB(oddballPredErrorAll<=bordersAllOB(q)) = q;
    else
        quantilesAllCP(changepointPredErrorAll>bordersAllCP(q-1)) = q;
        quantilesAllOB(oddballPredErrorAll>bordersAllOB(q-1)) = q;
    end
end
for q = 1:nBins
    meanQuantilesCPAllErrLR(q) = mean(changepointPredUpdateAll(quantilesAllCP==q));
    meanQuantilesOBAllErrLR(q) = mean(oddballPredUpdateAll(quantilesAllOB==q));
    quantileAllCPxesLR(q) = mean(changepointPredErrorAll(quantilesAllCP==q));
    quantileAllOBxesLR(q) = mean(oddballPredErrorAll(quantilesAllOB==q));
end

changepointobjPredErrorGood = changepointobjPredErrorAll(goodEstSel);
oddballobjPredErrorGood = oddballobjPredErrorAll(goodEstSel);
changepointEstErrGood = changepointEstErrAll(goodEstSel);
oddballEstErrGood = oddballEstErrAll(goodEstSel);

bordersGoodCP = quantile(changepointobjPredErrorGood,nBins-1);
bordersGoodOB = quantile(oddballobjPredErrorGood,nBins-1);
quantilesGoodOB = zeros(size(oddballobjPredErrorGood,1),1);
quantilesGoodCP = zeros(size(changepointobjPredErrorGood,1),1);
for q = nBins:-1:1
    if q<nBins
        quantilesGoodCP(changepointobjPredErrorGood<=bordersGoodCP(q)) = q;
        quantilesGoodOB(oddballobjPredErrorGood<=bordersGoodOB(q)) = q;
    else
        quantilesGoodCP(changepointobjPredErrorGood>bordersGoodCP(q-1)) = q;
        quantilesGoodOB(oddballobjPredErrorGood>bordersGoodOB(q-1)) = q;
    end
end
for q = 1:nBins
    meanQuantilesCPGoodErrBias(q) = mean(changepointEstErrGood(quantilesGoodCP==q));
    meanQuantilesOBGoodErrBias(q) = mean(oddballEstErrGood(quantilesGoodOB==q));
    quantileGoodCPxesBias(q) = mean(changepointobjPredErrorGood(quantilesGoodCP==q));
    quantileGoodOBxesBias(q) = mean(oddballobjPredErrorGood(quantilesGoodOB==q));
end

changepointPredErrorGood = changepointPredErrorAll(goodPredSel);
oddballPredErrorGood = oddballPredErrorAll(goodPredSel);
changepointPredUpdateGood = changepointPredUpdateAll(goodPredSel);
oddballPredUpdateGood = oddballPredUpdateAll(goodPredSel);

nBins = 40;
bordersGoodCP = quantile(changepointPredErrorGood,nBins-1);
bordersGoodOB = quantile(oddballPredErrorGood,nBins-1);
quantilesGoodOB = zeros(size(oddballPredErrorGood,1),1);
quantilesGoodCP = zeros(size(changepointPredErrorGood,1),1);
for q = nBins:-1:1
    if q<nBins
        quantilesGoodCP(changepointPredErrorGood<=bordersGoodCP(q)) = q;
        quantilesGoodOB(oddballPredErrorGood<=bordersGoodOB(q)) = q;
    else
        quantilesGoodCP(changepointPredErrorGood>bordersGoodCP(q-1)) = q;
        quantilesGoodOB(oddballPredErrorGood>bordersGoodOB(q-1)) = q;
    end
end
for q = 1:nBins
    meanQuantilesCPGoodUpLR(q) = mean(changepointPredUpdateGood(quantilesGoodCP==q));
    meanQuantilesOBGoodUpLR(q) = mean(oddballPredUpdateGood(quantilesGoodOB==q));
    quantileGoodCPxesLR(q) = mean(changepointPredErrorGood(quantilesGoodCP==q));
    quantileGoodOBxesLR(q) = mean(oddballPredErrorGood(quantilesGoodOB==q));
end

%calculate subject's standard deviations of errors and add subject's errors to list of all errors made by all subjects
for s = 1:length(behaveSubs)
    subno = behaveSubs(s);
    subNum = num2str(subno);

    allData=load(fullfile([behaveDir,'allSubCombined/', subNum, '_allBlockData.mat']));
    allData=allData.alldata;

    estErrPractice = allData.estErr(21:60,:);
    estErrStDevs(s) = std(estErrPractice,0,"all");
    estErrPred = allData.estErr(61:300,:);
    predErrPred = allData.predictErr([62:180,182:end],:);

    allEstErrPractice(:,:,s) = estErrPractice;
    allEstErrPred(:,:,s) = estErrPred;
    allPredErrPred(:,:,s) = predErrPred;
    stdSubEstPractice(s) = std(estErrPractice,0,'all');
    stdSubEstPred(s) = std(estErrPred,0,'all');
    stdSubPredPred(s) = std(predErrPred,0,'all');
    meanSubEstPractice(s) = mean(estErrPractice,'all');
    meanSubEstPred(s) = mean(estErrPred,'all');
    meanSubPredPred(s) = mean(predErrPred,'all');
end
numPractice = numel(allEstErrPractice);
numPred = numel(allEstErrPred);
pPracticeErrs = ttest(meanSubEstPractice);
pPredBlockErrs = ttest(meanSubEstPred);

%calculate subject level variance of error
varSubEstPractice = stdSubEstPractice.*stdSubEstPractice;
varSubPredPred = stdSubPredPred.*stdSubPredPred;

%calculate predicted std based on bayesian combination
predictedVar = 1./((1./varSubPredPred)+(1./varSubEstPractice));
predictedStd = sqrt(predictedVar);

%plot figure S1
figure("Position",[250,250,1200,330])
subplot(1,3,1)
scatter(stdSubEstPractice,stdSubEstPred,20,[.5,.5,.5],'filled','MarkerEdgeColor','k');
[h,pEsts] = ttest(stdSubEstPractice-stdSubEstPred);
hold on
plot([0,1.22],[0,1.22],"--k","LineWidth",0.25);
xlabel("Calibration Estimation σ","FontSize",11)
ylabel("Task Estimation σ","FontSize",11)
ylim([0,1.22])
xlim([0,1.22])
xticks(0:0.5:1.5)
yticks(0:0.5:1.5)
set(gca, 'box', 'off')
   set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)

subplot(1,3,2)
scatter(stdSubPredPred,stdSubEstPred,20,[.5,.5,.5],'filled','MarkerEdgeColor','k');
[h,pEstPred,~,stats] = ttest(stdSubPredPred-stdSubEstPred);
hold on
plot([0,1.6],[0,1.6],"--k","LineWidth",0.25);
xlabel("Task Prediction σ","FontSize",11)
ylabel("Task Estimation σ","FontSize",11)
ylim([0,1.6])
xlim([0,1.6])
xticks(0:0.5:1.5)
yticks(0:0.5:1.5)
   set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)

subplot(1,3,3)
scatter(predictedStd,stdSubEstPred,20,[.5,.5,.5],'filled','MarkerEdgeColor','k');
[h,pModel] = ttest(predictedStd-stdSubEstPred);
hold on
plot([0,1],[0,1],"--k","LineWidth",0.25);
ylabel("Task Estimation σ","FontSize",11)
xlabel("Bayesian Predicted σ","FontSize",11)
ylim([0,1])
xlim([0,1])
xticks(0:0.5:1.5)
yticks(0:0.5:1.5)
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)

% save figure S1
if doSTPResiduals == 0
    fig = gcf;
    figName = append("Figure_S1_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_S1_",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
else
    fig = gcf;
    figName = append("Figure_S1_Residual_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_S1_Residual",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
end

%% Run Behavioral regression
% preallocate variables for loop
paramsCircUpdateAll=[];
paramsCircBiasAll  =[];

paramsCircUpdateAllModel=[];
paramsCircBiasAllModel  =[];

normalizedCoefUpdate=[];
normalizedCoefBias=[];

allCPUpdate = [];
allOBUpdate = [];
allCPsPE = [];
allOBsPE = [];

noTrial = [1,120,121,240,241,360,361,480];
LRTrial = true(480,1);
LRTrial(noTrial) = false;

behaveAll.modelPredCP = behaveAll.maxLikePostMu;
behaveAll.modelPredCP(1:241:end) = [];
behaveAll.modelPredOB = behaveAll.maxLikePostC;
behaveAll.modelPredOB(1:241:end) = [];

for s=1:length(behaveSubs)
    subNum=behaveSubs(s);
    subStr=num2str(subNum);
    disp(subNum)
    
    %get eye data       
    
     % Get useful variables for descriptive analyses:
    decomposeSine=0;
    decomposeCosine=0;
    
    nTrials = 120; %length(behaveAll.block(behaveAll.subNum==subNum))*.5;

    %find subject's model variables
    OBpredErrorSubModel = behaveAll.subPredErrorOB(behaveAll.subNum==subNum);
    CPpredErrorSubModel = behaveAll.subPredErrorCP(behaveAll.subNum==subNum);

    perceptualErrorOnXModel = behaveAll.perceptualErrorOnX(behaveAll.subNum==subNum);
    perceptualErrorOnBModel = behaveAll.perceptualErrorOnB(behaveAll.subNum==subNum);

    CPpredErrorModel = behaveAll.predictionErrorOnX(behaveAll.subNum==subNum);
    OBpredErrorModel = behaveAll.predictionErrorOnB(behaveAll.subNum==subNum);

    CPupdateModel = behaveAll.modelPredCP(behaveAll.subNum==subNum);
    CPupdateModel = [nan;circ_dist(CPupdateModel(2:end),CPupdateModel(1:end-1))];
%     CPupdateModel = deg2rad(CPupdateModel);
    OBupdateModel = behaveAll.modelPredOB(behaveAll.subNum==subNum);
    OBupdateModel = [nan;circ_dist(OBupdateModel(2:end),OBupdateModel(1:end-1))];
%     OBupdateModel = deg2rad(OBupdateModel);

    %find subject's error/update variables
    OBpredError = [behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
    CPpredError = [behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];

    OBpredErrorSub = [behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
    CPpredErrorSub = [behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
    
    OBupdate = [behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
    CPupdate = [behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];

    perceptualErrorOnX = [behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
    perceptualErrorOnB = [behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];

    %find subjects surprise, entropy,  lr, bias, condition
    OBsurprise = behaveAll.surpriseOB(behaveAll.subNum==subNum);
    CPsurprise = behaveAll.surpriseCP(behaveAll.subNum==subNum);
    
    OBentropy = behaveAll.entropyOB(behaveAll.subNum==subNum);
    CPentropy = behaveAll.entropyCP(behaveAll.subNum==subNum);

    CPLR = [behaveAll.LRCP(behaveAll.subNum==subNum)];
    CPBias = [behaveAll.biasCP(behaveAll.subNum==subNum)];

    OBLR = [behaveAll.LROB(behaveAll.subNum==subNum)];
    OBBias = [behaveAll.biasOB(behaveAll.subNum==subNum)];

    conNum=behaveAll.blockCond(behaveAll.subNum==subNum);
    condNum=[conNum(1:nTrials);conNum(1:nTrials);conNum(nTrials+1:end);conNum(nTrials+1:end)];
    
    if conNum(1) == 1
        condition=1;
    else
        condition=2;
    end

    %other variables for regressions and uniform model
    predictionErrorOnX = [behaveAll.predictionErrorOnX(behaveAll.subNum==subNum)];
    predictionErrorOnB = [behaveAll.predictionErrorOnB(behaveAll.subNum==subNum)];
    
    OBsurpriseTrialNum=sum(nansum(behaveAll.surpriseTrial(behaveAll.subNum==subNum & behaveAll.blockCond==-1,:)));
    CPsurpriseTrialNum=sum(nansum(behaveAll.surpriseTrial(behaveAll.subNum==subNum & behaveAll.blockCond==1,:)));
    
    OBmoveDist=[behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
    CPmoveDist=[behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
        
    CPgazeAttention=[behaveAll.CPgazeAttention(behaveAll.subNum==subNum,:)];
    OBgazeAttention=[behaveAll.OBgazeAttention(behaveAll.subNum==subNum,:)];
    
    subjectPEobj=[OBpredError;CPpredError];
    modelPEobj=[predictionErrorOnB;predictionErrorOnX];
    
    subModPEobjDiff=modelPEobj-subjectPEobj;        
    subModPEobjDiffModel= modelPEobj-[OBpredErrorModel;CPpredErrorModel];

    clear data
    data.signedError=subModPEobjDiff;
    data.signedError_allTargs=[];
    
    whichParams=[1 1 0 0];
    data.doFit=true;

    % data.signedError  --> signed errors (or just angles if not for VWM task)
    % data.signedError_allTargs --> errors computed as if subject were
    % estimating all of the colors in the array (ie colorArray - subject
    % response).
    % data.xMat        --> all variables that could affect recall or precision
    % data.doFit  --> fit? if not, just evaluate at startPoint

    %params: maximum  likelihood model parameters,
    % 1= proportion gaussian,
    % 2= concentration (ie 1./sigma^2) of von mises,
    % 3= mean of gaussian (should be zero... but who knows!)
    % 4= proportion binding error (ie propr of gaussian that are evenly distributed across targets)
    
    % simplex order: gaussian, binding error, uniform 
    
    [~,~,~, uniformProb]=fit_VWM_mixtureModelTrial(data,whichParams);

    data.signedError=subModPEobjDiffModel;
    data.signedError_allTargs=[];
    [~,~,~, uniformProbModel]=fit_VWM_mixtureModelTrial(data,whichParams);

    if condition==1
        gazeAttention=[CPgazeAttention;OBgazeAttention];
        entropy = [nanzscore(CPentropy);nanzscore(OBentropy)];
    else
        gazeAttention=[OBgazeAttention;CPgazeAttention];
        entropy = [nanzscore(OBentropy);nanzscore(CPentropy)];
    end

    % Regression for the Update
    % xes and ys for both model and subject based regressions
    if condition==1
        xes = [ ones(4*nTrials,1), ... %This is one that actually gets used
            ([CPpredErrorSub;OBpredErrorSub]),...
            ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise]).* condNum),...
            ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise])),...
            ([CPpredErrorSub;OBpredErrorSub] .* entropy), ...
            ([CPpredErrorSub;OBpredErrorSub] .* condNum), ...
            ([CPpredErrorSub;OBpredErrorSub] .* nanzscore(round(uniformProb,9)))];

        xes = [xes(:,1),nanzscore(xes(:,2:end))];
        Y=[ CPupdate; OBupdate];

        xesMod = [ ones(4*nTrials,1), ... %This is one that actually gets used
            ([CPpredErrorSubModel;OBpredErrorSubModel]),...
            ([CPpredErrorSubModel;OBpredErrorSubModel] .* nanzscore([CPsurprise;OBsurprise]).* condNum),...
            ([CPpredErrorSubModel;OBpredErrorSubModel] .* nanzscore([CPsurprise;OBsurprise])),...
            ([CPpredErrorSubModel;OBpredErrorSubModel] .* entropy), ...
            ([CPpredErrorSubModel;OBpredErrorSubModel] .* condNum), ...
            ([CPpredErrorSubModel;OBpredErrorSubModel] .* nanzscore(round(uniformProb,9)))];

        xesMod = [xesMod(:,1),nanzscore(xesMod(:,2:end))];
        YMod=[ CPupdateModel; OBupdateModel];
    else
        xes = [ ones(4*nTrials,1), ... %This is one that actually gets used
            ([OBpredErrorSub;CPpredErrorSub]),...
            ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise;CPsurprise]) .* condNum),...
            ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise;CPsurprise])),...
            ([OBpredErrorSub;CPpredErrorSub] .* entropy), ...
            ([OBpredErrorSub;CPpredErrorSub] .* condNum), ...
            ([OBpredErrorSub;CPpredErrorSub] .* nanzscore(round(uniformProb,9)))];
        
        xes = [xes(:,1),nanzscore(xes(:,2:end))];
        Y=[OBupdate; CPupdate];
                    
        xesMod = [ ones(4*nTrials,1), ... %This is one that actually gets used
            ([OBpredErrorSubModel;CPpredErrorSubModel]),...
            ([OBpredErrorSubModel;CPpredErrorSubModel] .* nanzscore([OBsurprise;CPsurprise]) .* condNum),...
            ([OBpredErrorSubModel;CPpredErrorSubModel] .* nanzscore([OBsurprise;CPsurprise])),...
            ([OBpredErrorSubModel;CPpredErrorSubModel] .* entropy), ...
            ([OBpredErrorSubModel;CPpredErrorSubModel] .* condNum), ...
            ([OBpredErrorSubModel;CPpredErrorSubModel] .* nanzscore(round(uniformProb,9)))];
        
        xesMod = [xesMod(:,1),nanzscore(xesMod(:,2:end))];
        YMod=[OBupdateModel; CPupdateModel];
    end        

    %set up data for subject circular regression
    % data.Y              = ydata
    % data.X              = xdata
    % data.includeUniform = do you want to include a uniform mixture component?
    % data.whichParams    = which parameters should we fit?
    % data.startPoint     = where should we start parameter search
    % data.lb             = lower bound
    % data.ub             = upper bound
    % data.priorMean      = mean of gaussian parameter priors (should have one for each coefficient [NOT PRECISION OR MIXTURE parameters])
    % data.priorWidth     = width of gaussian priors (same as above).
    % data.nStart         = number of start points to use for optimizer
    regData.X=xes;
    regData.Y=Y;
    nanTrials = noTrial;

    regData.startPoint=[5,zeros(1,size(regData.X,2)),0.05];
    regData.whichParams=logical([ones(1,size(regData.X,2)+1),0]);
    regData.includeUniform=1;
    
    regData.priorMean=[0,0,0,0,0,0,0];
    
    regData.priorWidth=[1,1,narrowWidth,narrowWidth,narrowWidth,narrowWidth,narrowWidth];

    regData.lb=[.0001, ones(1, size(regData.X, 2)).*LB];
    regData.ub=[100, ones(1, size(regData.X, 2)).*UB];
    
    regData.nStart=nStart;
    regData.Y(nanTrials)=[];
    regData.X(nanTrials,:)=[];

    regDataModel = regData;
    regDataModel.X = xesMod;
    regDataModel.Y = YMod;
    regDataModel.Y(nanTrials)=[];
    regDataModel.X(nanTrials,:)=[];

    %runs subject and model-behavior circular error model
    [paramsUpdate, negLogLikeUpdate]=fitLinearModWCircErrs(regData);
    
    bestParamsUpdate=paramsUpdate;
    bestNegLogLikeUpdate=negLogLikeUpdate;

    [paramsUpdateModel, negLogLikeUpdateModel]=fitLinearModWCircErrs(regDataModel);
    
    bestParamsUpdateModel=paramsUpdateModel;
    bestNegLogLikeUpdateModel=negLogLikeUpdateModel;    

    % Regression for Perceptual Error
    % specify xes and ys for human and model
    if condition==1
        xes = [ ones(4*nTrials,1), ...
            ([CPpredError;OBpredError]),...
            ([CPpredError;OBpredError] .*  [nanzscore(CPsurprise); nanzscore(OBsurprise)].* condNum),...
            ([CPpredError;OBpredError] .*  [nanzscore(CPsurprise); nanzscore(OBsurprise)]),...
            ([CPpredError;OBpredError] .* entropy),...
            ([CPpredError;OBpredError] .* condNum), ...
            ([CPpredError;OBpredError] .* nanzscore(round(uniformProb,9))),...
            ([CPpredError;OBpredError] .* nanzscore(gazeAttention))];

        xes = [xes(:,1),nanzscore(xes(:,2:end))];
        Y=[perceptualErrorOnX; perceptualErrorOnB];        
        
        xesMod = [ ones(4*nTrials,1), ...
            ([CPpredErrorModel;OBpredErrorModel]),...
            ([CPpredErrorModel;OBpredErrorModel] .*  [nanzscore(CPsurprise); nanzscore(OBsurprise)].* condNum),...
            ([CPpredErrorModel;OBpredErrorModel] .*  [nanzscore(CPsurprise); nanzscore(OBsurprise)]),...
            ([CPpredErrorModel;OBpredErrorModel] .* entropy),...
            ([CPpredErrorModel;OBpredErrorModel] .* condNum), ...
            ([CPpredErrorModel;OBpredErrorModel] .* nanzscore(round(uniformProb,9))),...
            ([CPpredErrorModel;OBpredErrorModel] .* nanzscore(gazeAttention))];

        xesMod = [xesMod(:,1),nanzscore(xesMod(:,2:end))];
        YMod=[perceptualErrorOnXModel; perceptualErrorOnBModel];
        
    else
        xes = [ ones(4*nTrials,1), ...
            ([OBpredError;CPpredError]),...
            ([OBpredError;CPpredError] .*  [nanzscore(OBsurprise); nanzscore(CPsurprise)].* condNum),...
            ([OBpredError;CPpredError] .*  [nanzscore(OBsurprise); nanzscore(CPsurprise)]),...
            ([OBpredError;CPpredError] .* entropy),...
            ([OBpredError;CPpredError] .* condNum), ...
            ([OBpredError;CPpredError] .* nanzscore(round(uniformProb,9))),...
            ([OBpredError;CPpredError] .* nanzscore(gazeAttention))];

        xes = [xes(:,1),nanzscore(xes(:,2:end))];
        Y=[perceptualErrorOnB; perceptualErrorOnX];        
        
        xesMod = [ ones(4*nTrials,1), ...
            ([OBpredErrorModel;CPpredErrorModel]),...
            ([OBpredErrorModel;CPpredErrorModel] .*  [nanzscore(OBsurprise); nanzscore(CPsurprise)].* condNum),...
            ([OBpredErrorModel;CPpredErrorModel] .*  [nanzscore(OBsurprise); nanzscore(CPsurprise)]),...
            ([OBpredErrorModel;CPpredErrorModel] .* entropy),...
            ([OBpredErrorModel;CPpredErrorModel] .* condNum), ...
            ([OBpredErrorModel;CPpredErrorModel] .* nanzscore(round(uniformProb,9))),...
            ([OBpredErrorModel;CPpredErrorModel] .* nanzscore(gazeAttention))];

        xesMod = [xesMod(:,1),nanzscore(xesMod(:,2:end))];
        YMod=[perceptualErrorOnBModel; perceptualErrorOnXModel];
    end
      
    %fill regData for circular models
    regData.X=xes;
    regData.Y=Y;

    regData.startPoint=[5,zeros(1,size(regData.X,2)),0.05];
    regData.whichParams=logical([ones(1,size(regData.X,2)+1),0]);
    regData.includeUniform=1;
    
    regData.priorMean=[0,0,0,0,0,0,0,0];
    regData.priorWidth=[1,1,narrowWidth,narrowWidth,narrowWidth,narrowWidth,narrowWidth,narrowWidth];

    regData.lb=[.0001, ones(1, size(regData.X, 2)).*LB];
    regData.ub=[100, ones(1, size(regData.X, 2)).*UB];
    regData.nStart=nStart;

    regData.Y(nanTrials)=[];
    regData.X(nanTrials,:)=[];

    regDataModel = regData;
    regDataModel.X = xesMod;
    regDataModel.Y = YMod;
    regDataModel.Y(nanTrials)=[];
    regDataModel.X(nanTrials,:)=[];


    %run circular models
    [paramsBias, negLogLikeBias]=fitLinearModWCircErrs(regData);
    
    bestParamsBias=paramsBias;
    bestNegLogLikeBias=negLogLikeBias;       
        
    [paramsBiasModel, negLogLikeBiasModel]=fitLinearModWCircErrs(regDataModel);
    
    bestParamsBiasModel=paramsBiasModel;
    bestNegLogLikeBiasModel=negLogLikeBiasModel;    

    % Saving coefficient values
    paramsCircUpdateAll=cat(1,paramsCircUpdateAll,bestParamsUpdate);
    paramsCircBiasAll=cat(1,paramsCircBiasAll,bestParamsBias);        
    paramsCircUpdateAllModel=cat(1,paramsCircUpdateAllModel,bestParamsUpdateModel);
    paramsCircBiasAllModel=cat(1,paramsCircBiasAllModel,bestParamsBiasModel);
    
    allCPsPE = [allCPsPE;CPpredErrorSubModel];
    allOBsPE = [allOBsPE;OBpredErrorSubModel];
    allCPUpdate = [allCPUpdate;CPupdateModel];
    allOBUpdate = [allOBUpdate;OBupdateModel];
end
for i = 1:size(paramsCircUpdateAll,2)
[~,pUp(i),~,statsUp] = ttest(paramsCircUpdateAll(:,i));
tStatUp(i) =  statsUp.tstat;
end

for i = 1:size(paramsCircBiasAll,2)
[~,pBias(i),~,statsBias] = ttest(paramsCircBiasAll(:,i));
tStatBias(i) =  statsBias.tstat;
end
%% Step 2: Load EEG/eye data & run regression
% Raw EEG data regression
if runEEGRegression == 1
    for s=1:length(EEGSubs)
        subno=EEGSubs(s);
        subNum=num2str(subno);
        disp(subNum)
        resultEEG=eeg_analysisFunc(subNum,saveText,dirs);
        b_mat_eeg(s,:,:,:)=resultEEG.mat; %subject by channel by trial by regressor
    end
    saveDir=[basePath];
    fn=fullfile(saveDir,['b_mat_eeg',saveText,'.mat']);
    save(fn,'b_mat_eeg')
    beep
end
%%
%Eye data Regression
for s=1:length(eyeSubs)
    subno=eyeSubs(s);
    subNum=num2str(subno);
    disp(subNum)
    resultEye=eyeRegressionFunc(subNum,saveText,timeBeforeEye,timeAfterEye,blinkWindow,baselineTimeEye,dirs);
    resultEyePred=eyeRegressionPredResp(subNum,saveText,timeBeforePred,timeAfterPred,blinkWindow,baselineTimeStart,baselineTimeEnd,dirs);
    allBs(:,:,s) = resultEye.B;
    BBaseline(s,:) = resultEye.BBaseline;
    allBsPred(:,:,s) = resultEyePred.B;
    allDiffBsPred(:,:,s) = resultEyePred.diffB;
end
allBsPerm = permute(allBs,[3,2,1]);
allBsPredPerm = permute(allBsPred,[3,2,1]);
allBsDiffPerm = permute(allDiffBsPred,[3,2,1]);

%% Clustering/Permutation Test
clustThreshEEG = 0.01;
clustThreshEye = 0.025;
%load channel locations
load("chanlocs.mat")

% NOTE EEG.dat

% GET CONNECTION MATRIX specifying which electrodes are connected to which:
% if exist('connectionMat.mat')
%     load connectionMat.mat  % if we've already got one made, load it.
% else
    % otherwise, create one from scratch:
    
    % Then get the XYZ coordinates for each channel:
    clear eye
    allLocas=[[chanlocs.X]; [chanlocs.Y]; [chanlocs.Z]]' ;
    % Loop through the channels and get the distance between that channel and
    % all other channels.
    for i = 1:length(allLocas)
        relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
    end
    
    % Set a threshold on distance... and mark channels that fall within that
    % threshold:
    connectionMat=relDist<connectThresh;
    connectionMat=connectionMat-eye(length(connectionMat));
    
    % CHECK OUT CONNECTIONS:
%     figure
%     imagesc(connectionMat);
%     keyboard
    %close all

% end
% STEP 1b: compute t-stats
regDat = b_mat_eeg(:,:,:,:);
EEG_dat=regDat(:,:,:,2);
% get clusters, cluster sizes, cluster masses for positive clusters (ie
% p<.05 in a one tailed positive test):
posClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, clustThreshEEG, 'right');
% and negative effects:  
negClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, clustThreshEEG, 'left');         

%%STEP 1c:
% 1: Randomly flip signs of coefficients for each subject
% 2: redo t-test and cluster formation
% 3: record the cluster mass of the "biggest" cluster
% 4: repeat procedure 1000 or more times
% 5: see where cluster mass values from actual data fall along distribution
%       of permutation cluster mass values


for k=1:numPermEEG
    % create a subject length array containing randomly assigned -1's and
    % 1s:
    permArray=ones(size(EEG_dat, 1), 1);
    permArray(logical(binornd(1, .5, size(EEG_dat, 1), 1)))=-1;
    
    % multiply each subject timeseries by the -1 or 1 assigned randomly on
    % this trial
    sz=size(EEG_dat);
    permMat=repmat(permArray, [1 sz(2:end)]);
    
    % get cluster statistics for permuted dataset:
    permClusterInfo=getEEG_clusterSize(EEG_dat.*permMat, connectionMat, clustThreshEEG, 'right');
    
    % store the maximum of the statistics, for use in null distribution:
    maxSize(k)=(max(permClusterInfo.clustSizeMap(:)));
    maxWt(k)=(max(permClusterInfo.clustWtMap(:)));
    if mod(k,100)==0
        disp(k)
    end
end

% For a two tailed test, find the the minimum cluster statistics necessary
% to beat 97.5% of the null distribution, based on "mass"
massTh = prctile(maxWt,  97.5);
gPosEEG=unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap	>massTh));
gNegEEG=unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap	>massTh));
threshMapEEG=zeros(size(negClusterInfo.tMap));
threshMapEEG(posClusterInfo.clustWtMap	>massTh)=posClusterInfo.tMap(posClusterInfo.clustWtMap	>massTh);
threshMapEEG(negClusterInfo.clustWtMap	>massTh)=negClusterInfo.tMap(negClusterInfo.clustWtMap	>massTh);

clustMasses = [unique(posClusterInfo.clustWtMap(posClusterInfo.clustWtMap>massTh),'stable');unique(negClusterInfo.clustWtMap(negClusterInfo.clustWtMap>massTh),'stable')];
clear pVal
for k = 1:length([gPosEEG;gNegEEG])
    pVal(k) = (sum(maxWt>clustMasses(k))+1)/numPermEEG;
end

%% repeats permutation test for eye data
% 1 - surprise @ clustThreshEye for analysis
% 2 - surprise @ 0.025 for figures
% 3 - entropy @ 0.025 for figures
% 4 - prediction phase STPCPOB @ 0.025 for figures
downSampTimesEye = timeBeforeEye:timeAfterEye;
titles = ["Intercept","MeanSTP","MeanSTPCPOB","MeanEntropy","Condition","SumSinColor","SumCosColor","Baseline 1000 ms"];

for i = 1:size(allBsPerm,3)
    [h,p(i,:)] = ttest(allBsPerm(:,:,i));
end

% same steps as previous test, with possibly more permutations since this goes a lot faster
clusterMaxPerm=zeros(numPermEye,1);
clusterMaxPermSig=zeros(numPermEye,1);
clusterMaxPermEnt=zeros(numPermEye,1);
clusterMaxPermPredSTP=zeros(numPermEye,1);
clusterMaxPermDiffSTP=zeros(numPermEye,1);
for i=1:numPermEye
    allBsVar = allBsPerm(:,:,2);
    allBsEnt = allBsPerm(:,:,4);
    allBsPredSTP= allBsPredPerm(:,:,3);
    allBsDiffSTP= allBsDiffPerm(:,:,3);
    %randomly select 1 and -1
    random_flip=randsample([1, -1], size(allBsPerm,1), true);
    random_flip_mat = random_flip'*ones(1,size(allBsPerm,2));
    data_flippedVar = random_flip_mat.*allBsVar;
    
    %randomly select 1 and -1
    random_flip=randsample([1, -1], size(allBsPerm,1), true);
    random_flip_mat = random_flip'*ones(1,size(allBsPerm,2));
    data_flippedEnt = random_flip_mat.*allBsEnt;

    %randomly select 1 and -1
    random_flip=randsample([1, -1], size(allBsPredPerm,1), true);
    random_flip_mat = random_flip'*ones(1,size(allBsPredPerm,2));
    data_flippedPredSTP = random_flip_mat.*allBsPredSTP;

    %randomly select 1 and -1
    random_flip=randsample([1, -1], size(allBsDiffPerm,1), true);
    random_flip_mat = random_flip'*ones(1,size(allBsDiffPerm,2));
    data_flippedDiffSTP = random_flip_mat.*allBsDiffSTP;

    %calculate 15 t stats for permuted data
    pupilClusterInfo = getPupil_clusterSize(data_flippedVar,0,clustThreshEye,'right');
    clusterMaxPerm(i) = max(pupilClusterInfo.clustWtMap);
    pupilClusterInfo = getPupil_clusterSize(data_flippedVar,0,sigThresh,'right');
    clusterMaxPermSig(i) = max(pupilClusterInfo.clustWtMap);
    pupilClusterInfo = getPupil_clusterSize(data_flippedEnt,0,sigThresh,'right');
    clusterMaxPermEnt(i) = max(pupilClusterInfo.clustWtMap);
    pupilClusterInfo = getPupil_clusterSize(data_flippedPredSTP,0,sigThresh,'right');
    clusterMaxPermPredSTP(i) = max(pupilClusterInfo.clustWtMap);
    pupilClusterInfo = getPupil_clusterSize(data_flippedDiffSTP,0,sigThresh,'right');
    clusterMaxPermDiffSTP(i) = max(pupilClusterInfo.clustWtMap);
    if mod(i,100)==0
        disp(i)
    end
end


% calculate cluster size of real data
pupilClusterInfo = getPupil_clusterSize(allBsVar,0,clustThreshEye*2,'both');
clusterMax = max(pupilClusterInfo.clustWtMap);

pupilClusterInfoSig = getPupil_clusterSize(allBsVar,0,sigThresh*2,'both');
clusterMaxSig = max(pupilClusterInfoSig.clustWtMap);
pValEye = (sum(clusterMaxPermSig>=clusterMaxSig)+1)/numPermEye*2;

pupilClusterInfoEnt = getPupil_clusterSize(allBsEnt,0,sigThresh*2,'both');
clusterMaxEnt = max(pupilClusterInfoEnt.clustWtMap);
pValEnt = (sum(clusterMaxPermEnt>=clusterMaxEnt)+1)/numPermEye*2;

pupilClusterInfoPredSTP = getPupil_clusterSize(allBsPredSTP,0,sigThresh*2,'both');
clusterMaxPredSTP = max(pupilClusterInfoPredSTP.clustWtMap);
pValPredSTP = (sum(clusterMaxPermPredSTP>=clusterMaxPredSTP)+1)/numPermEye*2;

pupilClusterInfoDiffSTP = getPupil_clusterSize(allBsDiffSTP,0,sigThresh*2,'both');
clusterMaxDiffSTP = max(pupilClusterInfoDiffSTP.clustWtMap);
pValDiffSTP = (sum(clusterMaxPermDiffSTP>=clusterMaxDiffSTP)+1)/numPermEye*2;

%find mass threshold and get clusters where data is larger than threshold
%for the 2nd 3rd and 4th, get a map of where cluster size is significantly large
massThEye=prctile(clusterMaxPerm, 97.5);
gSig=unique(pupilClusterInfo.ID_map(pupilClusterInfo.clustWtMap>massThEye));
threshMapEye=zeros(size(pupilClusterInfo.tMap));
threshMapEye(pupilClusterInfo.clustWtMap>massThEye)=pupilClusterInfo.tMap(pupilClusterInfo.clustWtMap>massThEye);

massThSig=prctile(clusterMaxPermSig, 97.5);
gSigSig=unique(pupilClusterInfoSig.ID_map(pupilClusterInfoSig.clustWtMap>massThSig));
threshMapEyeSig=zeros(size(pupilClusterInfoSig.tMap));
threshMapEyeSig(pupilClusterInfoSig.clustWtMap>massThSig)=pupilClusterInfoSig.tMap(pupilClusterInfoSig.clustWtMap>massThSig);
sigTimes = find(threshMapEyeSig~=0);

massThEnt=prctile(clusterMaxPermEnt, 97.5);
gSigEnt=unique(pupilClusterInfoEnt.ID_map(pupilClusterInfoEnt.clustWtMap>massThEnt));
threshMapEyeEnt=zeros(size(pupilClusterInfoEnt.tMap));
threshMapEyeEnt(pupilClusterInfoEnt.clustWtMap>massThEnt)=pupilClusterInfoEnt.tMap(pupilClusterInfoEnt.clustWtMap>massThEnt);
sigTimesEnt = find(threshMapEyeEnt~=0);

massThPredSTP=prctile(clusterMaxPermPredSTP, 97);
gSigPredSTP=unique(pupilClusterInfoPredSTP.ID_map(pupilClusterInfoPredSTP.clustWtMap>massThPredSTP));
threshMapEyePredSTP=zeros(size(pupilClusterInfoPredSTP.tMap));
threshMapEyePredSTP(pupilClusterInfoPredSTP.clustWtMap>massThPredSTP)=pupilClusterInfoPredSTP.tMap(pupilClusterInfoPredSTP.clustWtMap>massThPredSTP);
sigTimesPredSTP = find(threshMapEyePredSTP~=0);

massThDiffSTP=prctile(clusterMaxPermDiffSTP, 97.5);
gSigDiffSTP=unique(pupilClusterInfoDiffSTP.ID_map(pupilClusterInfoDiffSTP.clustWtMap>massThDiffSTP));
threshMapEyeDiffSTP=zeros(size(pupilClusterInfoDiffSTP.tMap));
threshMapEyeDiffSTP(pupilClusterInfoDiffSTP.clustWtMap>massThDiffSTP)=pupilClusterInfoDiffSTP.tMap(pupilClusterInfoDiffSTP.clustWtMap>massThDiffSTP);
sigTimesDiffSTP = find(threshMapEyeDiffSTP~=0);

figure
imagesc(threshMapEEG);

%% Make big summary figure
panel12 = [1:4,13:16,25:28,37:40,49:52,61:64];
panel22a = [6:12,18:24,30:36];
panel22b = [42:48,54:60,66:72];
panel32 = [97:99,109:111,121:123];
panel42 = [100:102,112:114,124:126];
panel52 = [103:105,115:117,127:129];
panel62 = [106:108,118:120,130:132];
panel72 = [133:168];
panel82 = [181:185,193:197];
panel92 = [205:209,217:221];
panel102 = [187:191,199:203];
panel112 = [211:215,223:227];
load('chanlocs.mat')

clusterInfo = getEEG_clusterSize(b_mat_eeg(:,:,:,2), connectionMat, clustThreshEEG*2,'both');

figure("Position",[100,50,700,1000])
set(gcf,'renderer','Painters')
subplot(19,12,panel12);
%Baseline Pupil
    xticklabels(["Entropy","Cond","STP","STP*Cond"])
    xticks(1:4)
    xlim([0,5])
    sem = std(BBaseline)/sqrt(length(eyeSubs));
    hold on
    errorbar(1,mean(BBaseline(:,4)),sem(4),'.',"Color",cbColors(2,:),"MarkerSize",20,"Marker",'.',"LineWidth",2)
    errorbar(2,mean(BBaseline(:,5)),sem(5),'.',"Color",cbColors(3,:),"MarkerSize",20,"Marker",'.',"LineWidth",2)
    errorbar(3,mean(BBaseline(:,2)),sem(2),'.',"Color",cbColors(5,:),"MarkerSize",20,"Marker",'.',"LineWidth",2)
    errorbar(4,mean(BBaseline(:,3)),sem(3),'.',"Color",cbColors(4,:),"MarkerSize",20,"Marker",'.',"LineWidth",2)
    yline(0)
    yticks([0])
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')
    ylabel("Regression Coefficient")

subplot(19,12,panel22b);
%Running pupil regressors
    %plot each regressor with error bars in right color
    if realData
        sem = std(allBsPerm(:,:,2))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./1000,mean(allBsPerm(:,:,2)),[mean(allBsPerm(:,:,2))-sem;mean(allBsPerm(:,:,2))+sem],{'-','color',cbColors(5,:),'markerfacecolor',cbColors(5,:)},1)
        hold on
        sem = std(allBsPerm(:,:,3))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./1000,mean(allBsPerm(:,:,3)),[mean(allBsPerm(:,:,3))-sem;mean(allBsPerm(:,:,3))+sem],{'-','color',cbColors(4,:),'markerfacecolor',cbColors(4,:)},1)
        sem = std(allBsPerm(:,:,4))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./1000,mean(allBsPerm(:,:,4)),[mean(allBsPerm(:,:,4))-sem;mean(allBsPerm(:,:,4))+sem],{'-','color',cbColors(2,:),'markerfacecolor',cbColors(2,:)},1)
        sem = std(allBsPerm(:,:,5))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./1000,mean(allBsPerm(:,:,5)),[mean(allBsPerm(:,:,5))-sem;mean(allBsPerm(:,:,5))+sem],{'-','color',cbColors(3,:),'markerfacecolor',cbColors(3,:)},1)
    else
        sem = std(allBsPerm(:,:,2))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./250,mean(allBsPerm(:,:,2)),[mean(allBsPerm(:,:,2))-sem;mean(allBsPerm(:,:,2))+sem],{'-','color',cbColors(5,:),'markerfacecolor',cbColors(5,:)},1)
        hold on
        sem = std(allBsPerm(:,:,3))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./250,mean(allBsPerm(:,:,3)),[mean(allBsPerm(:,:,3))-sem;mean(allBsPerm(:,:,3))+sem],{'-','color',cbColors(4,:),'markerfacecolor',cbColors(4,:)},1)
        sem = std(allBsPerm(:,:,4))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./250,mean(allBsPerm(:,:,4)),[mean(allBsPerm(:,:,4))-sem;mean(allBsPerm(:,:,4))+sem],{'-','color',cbColors(2,:),'markerfacecolor',cbColors(2,:)},1)
        sem = std(allBsPerm(:,:,5))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforeEye:timeAfterEye)./250,mean(allBsPerm(:,:,5)),[mean(allBsPerm(:,:,5))-sem;mean(allBsPerm(:,:,5))+sem],{'-','color',cbColors(3,:),'markerfacecolor',cbColors(3,:)},1)
    end
    % the part where you add significant clusters to figure
    if ~isempty(sigTimes)
        scatter((sigTimes-timeBeforeEye-1)/1000,mean(allBsPerm(:,sigTimes,2)),5,[0.1621    0.3301    0.1992],"filled")
    end
    if ~isempty(sigTimesEnt)
        scatter((sigTimesEnt-timeBeforeEye-1)/1000,mean(allBsPerm(:,sigTimesEnt,4)),5,[0.6289    0.4348         0],"filled") 
    end

    xlabel('Time Relative to Stimulus Onset (s)')
    ylabel('Coefficients')
    xticks([-1:4])
    xlim([-1,4])
    yticks([0])
    yline(0)
    xline(0) 
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

subplot(19,12,panel22a);
%Running pupil intercept
    sem = std(allBsPerm(:,:,1))./sqrt(length(eyeSubs)-1);
    shadedErrorBar((-timeBeforeEye:timeAfterEye)./1000,mean(allBsPerm(:,:,1)),[mean(allBsPerm(:,:,1))-sem;mean(allBsPerm(:,:,1))+sem],{'-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5]},1)
    ylabel('Intercept')
    xticks([])
    xlim([-1,4])
    yticks([0])
    yline(0)
    xline(0) 
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

subplot(19,12,panel32);
%topoplot 1
    currplot = topoplot(clusterInfo.tMap(:,EEGTimes==319), chanlocs,'Colormap',jet);
    caxis([-6,6])
    
    title('320 ms')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)

subplot(19,12,panel42);
%topoplot 2
    currplot = topoplot(clusterInfo.tMap(:,EEGTimes==449), chanlocs,'Colormap',jet);
    caxis([-6,6])
    
    title('450 ms')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)

subplot(19,12,panel52);
%topoplot 3
    currplot = topoplot(clusterInfo.tMap(:,EEGTimes==579), chanlocs,'Colormap',jet);
    caxis([-6,6])
    title('580 ms')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)

subplot(19,12,panel62);
%topoplot 4
    currplot = topoplot(clusterInfo.tMap(:,EEGTimes==999), chanlocs,'Colormap',jet);
    caxis([-6,6])
    title('1000 ms')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)

traceTimes = -500:2000;
subplot(19,12,panel72);
%eeg heatmap
    imagesc(traceTimes./1000, [],clusterInfo.tMap(:,traceTimes+timeBeforeEEG))
    cb=colorbar;
    ylabel('Channel')
    xlabel('Time (s)')
    yticklabels([])
    caxis([-6,6])
    cb.Ticks = [-6,0,6];
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')


    mat1=b_mat_eeg(:,[63],:,:); %FCz
    mat2=reshape(mean(mean(mat1,2),1),4000,[]); 
    mat3 = mat2';
    mat4 = reshape(mean(mat1,2),size(b_mat_eeg,[1,3,4]));

subplot(19,12,panel82);
%intercept plot (FCz)
    sem = std(mat4(:,:,1))./sqrt(size(b_mat_eeg,1)-1);
    shadedErrorBar(traceTimes./1000,movmean(mat3(1,traceTimes+timeBeforeEEG),40),[movmean(mat3(1,traceTimes+timeBeforeEEG),40)-sem(traceTimes+timeBeforeEEG);movmean(mat3(1,traceTimes+timeBeforeEEG),40)+sem(traceTimes+timeBeforeEEG)],{'b-','markerfacecolor','b'},1)
    yline(0,'--')
    yticks(0)
    xlim([-0.5,1.5])
    xticks([])
    ylabel("Intercept")
    title("FCz")
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

subplot(19,12,panel92);
%stp plot (FCz)
    sem = std(mat4(:,:,2))./sqrt(size(b_mat_eeg,1)-1);
    shadedErrorBar(traceTimes./1000,movmean(mat3(2,traceTimes+timeBeforeEEG),40),[movmean(mat3(2,traceTimes+timeBeforeEEG),40)-sem(traceTimes+timeBeforeEEG);movmean(mat3(2,traceTimes+timeBeforeEEG),40)+sem(traceTimes+timeBeforeEEG)],{'b-','markerfacecolor','b'},1)
    yline(0,'--')
    yticks(0)
    xlim([-0.5,1.5])
    xticks(-0.5:0.5:1.5)
    xtickangle(0)
    xlabel("Time (s)")
    ylabel("STP")
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

    mat1=b_mat_eeg(:,[13],:,:); %Pz
    mat2=reshape(mean(mean(mat1,2),1),4000,[]); 
    mat3 = mat2';
    mat4 = reshape(mean(mat1,2),size(b_mat_eeg,[1,3,4]));

subplot(19,12,panel102);
%intercept plot (Pz)
    sem = std(mat4(:,:,1))./sqrt(size(b_mat_eeg,1)-1);
    shadedErrorBar(traceTimes./1000,movmean(mat3(1,traceTimes+timeBeforeEEG),40),[movmean(mat3(1,traceTimes+timeBeforeEEG),40)-sem(traceTimes+timeBeforeEEG);movmean(mat3(1,traceTimes+timeBeforeEEG),40)+sem(traceTimes+timeBeforeEEG)],{'b-','markerfacecolor','b'},1)
    yline(0,'--')
    yticks(0)
    xlim([-0.5,1.5])
    xticks([])
    ylabel("Intercept")
    title("Pz")
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

subplot(19,12,panel112);
%stp plot (Pz)
    sem = std(mat4(:,:,2))./sqrt(size(b_mat_eeg,1)-1);
    shadedErrorBar(traceTimes./1000,movmean(mat3(2,traceTimes+timeBeforeEEG),40),[movmean(mat3(2,traceTimes+timeBeforeEEG),40)-sem(traceTimes+timeBeforeEEG);movmean(mat3(2,traceTimes+timeBeforeEEG),40)+sem(traceTimes+timeBeforeEEG)],{'b-','markerfacecolor','b'},1)
    yline(0,'--')
    yticks(0)
    xlim([-0.5,1.5])
    xticks(-0.5:0.5:1.5)
    xtickangle(0)
    xlabel("Time (s)")
    ylabel("STP")
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    set(gca, 'box', 'off')

%save figure 2 (summary)
if doSTPResiduals == 0
    fig = gcf;
    figName = append("Figure_2_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_2_",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
else
    fig = gcf;
    figName = append("Figure_2_Residual_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_2_Residual",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
end

%% Make Figure 5 (oddball pupil)
if ~isempty(sigTimesPredSTP)
    if oddballFigNewVer == 1
        figure('Position',[200 250 1300 600])
        subplot(1,2,1)
    else
        figure('Position',[400 250 700 600])
    end
    set(gcf,'renderer','Painters')
    
        sem = std(allBsPredPerm(:,:,3))./sqrt(length(eyeSubs)-1);
        shadedErrorBar((-timeBeforePred:timeAfterPred)./1000,(mean(allBsPredPerm(:,:,3))),[(mean(allBsPredPerm(:,:,3)))-sem;(mean(allBsPredPerm(:,:,3)))+sem],{'-','color',cbColors(4,:),'markerfacecolor',cbColors(4,:)},1)
        hold on
        % the part where you add significant clusters to figure
        scatter((sigTimesPredSTP-timeBeforePred-1)./1000,mean(allBsPredPerm(:,sigTimesPredSTP,3)),5,[0.5004    0.2352    0.5469],"filled") 
        xlabel("Time Relative to Prediction (s)")
        xticks([-2,0,2,4])
        yticks([0])
        yline(0)
        xline(0)
        ylim([-0.021,0.021])
        ylabel("STP*CP/OB Coefficient")
        text(1.7,-0.0155,"Larger Oddball Pupil","FontName","Arial","FontWeight","bold","FontSize",16)
        text(1.1,0.0155,"Larger Changepoint Pupil","FontName","Arial","FontWeight","bold","FontSize",16)
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",16)
        set(gca, 'box', 'off')
    
    if oddballFigNewVer == 1
        subplot(1,2,2)
    %     meanResponse = mean(allBsPredPerm(:,timeBeforePred-1500:timeBeforePred-800,:),2);
        meanResponse = mean(allBsPredPerm(:,sigTimesPredSTP,:),2);
        scatter(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4),meanResponse(:,3),30,[227 152 227]./255,'filled','MarkerEdgeColor',[200 50 200]./255)
        hold on
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
        ylabel("Eye STP*CP/OB Coefficient");
        xlabel("Behavioral PE*STP*CP/OB Coef");
        yticks(0)
        xticks(0)
        maxVal = max(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4));
        minVal = min(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4));
        xlim([minVal-0.05*(maxVal-minVal),maxVal+0.05*(maxVal-minVal)])
        maxVal = max(meanResponse(:,3))*1.05;
        minVal = min(meanResponse(:,3))*1.05;
        ylim([minVal-0.05*(maxVal-minVal),maxVal+0.05*(maxVal-minVal)])
        yline(0,'--','LineWidth',1)
        xline(0,'--','LineWidth',1)
        mOdd = polyfit(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4),meanResponse(:,3),1);
        plot((linspace(min(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4)),max(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4)),2)),mOdd(1)*(linspace(min(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4)),max(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4)),2))+mOdd(2),'Color',[200 50 200]./255);
    end
    [eyeCorr,eyeCorrP] = corr(paramsCircUpdateAll(ismember(behaveSubs,eyeSubs),4),meanResponse(:,3),"Type","spearman");
    disp(eyeCorr)
    disp(eyeCorrP)
    %save Figure 5 (oddball 2nd transition in pupil figure)
    if doSTPResiduals == 0
        fig = gcf;
        figName = append("Figure_5_",figTime,'.eps');
        figLoc = append(figDir,figName);
        exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
        figName = append("Figure_5_",figTime,'.png');
        figLoc = append(figDir,figName);
        saveas(fig,figLoc)
    else
        fig = gcf;
        figName = append("Figure_5_Residual_",figTime,'.eps');
        figLoc = append(figDir,figName);
        exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
        figName = append("Figure_5_Residual",figTime,'.png');
        figLoc = append(figDir,figName);
        saveas(fig,figLoc)
    end
    beep
end
%% Make Figure S3 (oddball pupil Derivative)
figure('Position',[400 250 700 600])
set(gcf,'renderer','Painters')

    sem = std(allBsDiffPerm(:,:,3))./sqrt(length(eyeSubs)-1);
    shadedErrorBar((-timeBeforePred:timeAfterPred)./1000,(mean(allBsDiffPerm(:,:,3))),[(mean(allBsDiffPerm(:,:,3)))-sem;(mean(allBsDiffPerm(:,:,3)))+sem],{'-','color',cbColors(4,:),'markerfacecolor',cbColors(4,:)},1)
    hold on
    % the part where you add significant clusters to figure
    scatter((sigTimesDiffSTP-timeBeforePred-1)./1000,mean(allBsDiffPerm(:,sigTimesDiffSTP,3)),5,[0.5004    0.2352    0.5469],"filled") 
    xlabel("Time Relative to Prediction (s)")
    xticks([-2,0,2,4])
    yticks([0])
    yline(0)
    xline(0)
    ylim([-0.002,0.002])
    ylabel("STP*CP/OB Coefficient")
    text(1.1,-0.00155,"Larger Oddball Derivative","FontName","Arial","FontWeight","bold","FontSize",16)
    text(.65,0.00155,"Larger Changepoint Derivative","FontName","Arial","FontWeight","bold","FontSize",16)

    set(gca,"FontName","Arial","FontWeight","bold","FontSize",16)
    set(gca, 'box', 'off')

%save Figure 5 (oddball 2nd transition in pupil figure)
if doSTPResiduals == 0
    fig = gcf;
    figName = append("Figure_S3_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_S3_",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
else
    fig = gcf;
    figName = append("Figure_S3_Residual_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_S3_Residual",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
end
%% Trial-by-trial Effect
clear relROIEEG relROIEye
% setup, specify where the clusters we're looking at are
% make a matrix for each that is the cluster mass at the cluster and 0 everywhere else
if ~isempty(gSig)
    k=1;
    while k <=length(gSig)
        relROIEye.maps(:,:,k)=pupilClusterInfo.ID_map==gSig(k);
        k=k+1;
    end 
     
    for k = 1:length(gSig)
        clustTs=relROIEye.maps(:,:,k).*abs(pupilClusterInfo.tMap);
        [relROIEye.peakChannel(k),J] = find(clustTs==max(clustTs(:)));
    end
    relROIEye.fullTMap=pupilClusterInfo.tMap;
    
    
    relROIEEG.sign=[ones(size(gPosEEG)); -ones(size(gNegEEG))];
    k =1;
    while k <=length(gPosEEG)
        relROIEEG.maps(:,:,k)=posClusterInfo.ID_map==gPosEEG(k);
        k=k+1;
    end
    
    while k <= length(gPosEEG) +length(gNegEEG)
        relROIEEG.maps(:,:,k)=negClusterInfo.ID_map==gNegEEG(k-length(gPosEEG));
        k=k+1;
    end
     
    for k = 1:length(relROIEEG.sign)
        clustTs=relROIEEG.maps(:,:,k).*abs(posClusterInfo.tMap);
        [relROIEEG.peakChannel(k),J] = find(clustTs==max(clustTs(:)));
    end
    %
    relROIEEG.fullTMap=posClusterInfo.tMap;
    
    if eegTimestepMode == 1
        % Instead of using the clusters, you can use regions of clusters or regions of timepoints
        
        relROIEEG.maps = zeros(size(relROIEEG.fullTMap,1),size(relROIEEG.fullTMap,2));
    %     relROIEEG.maps([2,3,6,7,27,28,29,31,33,34,35,36,57,59,60,61,62,63],:) = 1; %Frontals
    %     relROIEEG.maps([12 13 14 19 23 42 43 44 46 47 48 50 51 52 53],:) = -1; %Parietals
    %     relROIEEG.maps([3,4,5,6,8,9,10,11,14,15,32,33,36,37,38,40,41,42,44,45,46],:) = 1; %Lefts
    %     relROIEEG.maps([19,20,21,22,24,25,26,27,29,30,48,49,50,53,54,55,57,58,59,60,61],:) = -1; %Rights
    %     relROIEEG.maps = ones(size(relROIEEG.fullTMap,1),size(relROIEEG.fullTMap,2)); %All
        relROIEEG.maps([2,3,6,7,27,28,29,31,33,34,35,36,57,59,60,61,62,63,12 13 14 19 23 42 43 44 46 47 48 50 51 52 53],:) = 1; %Frontal & Parietal
        relROIEEG.maps = sum(relROIEEG.maps,3);
        %relROI.maps = relROI.fullTMap;
        
        increment = 100;
        timeOI = [1500,4000];
        borderTimes = timeOI(1):increment:timeOI(end);
        trueMap = relROIEEG.maps(:,:,1);
        
        relROIEEG.maps =  false(size(relROIEEG.maps,1),size(relROIEEG.maps,2),size(relROIEEG.maps,3));
        for i = 1:length(borderTimes)-1
            relROIEEG.maps(:,borderTimes(i):borderTimes(i+1),i) = trueMap(:,borderTimes(i):borderTimes(i+1));
        end
    
        %plots all clusters you will be using, each cluster is a different color
        figMap = zeros(size(relROIEEG.maps,1),size(relROIEEG.maps,2));
        for m = 1:size(relROIEEG.maps,3)
            figMap(relROIEEG.maps(:,:,m)) = m;
        end
        
        figure
        imagesc(figMap)
    %     relROIEEG.fullTMap = ones(size(relROIEEG.maps,1),size(relROIEEG.maps,2));
    end
    
    
    
    % Instead of using the whole period of pupil dilation, you can use bins of timepoints
    
    % relROIEye.maps = ones(1,size(relROIEye.fullTMap,2));
    % relROIEye.maps = sum(relROIEye.maps,3);
    % 
    % increment = 200;
    % timeOI = [2001,10001];
    % borderTimes = timeOI(1):increment:timeOI(end);
    % trueMap = relROIEye.maps(:,:,1);
    % 
    % relROIEye.maps = false(size(relROIEye.maps,1),size(relROIEye.maps,2),size(relROIEye.maps,3));
    % 
    % for i = 1:length(borderTimes)-1
    %     relROIEye.maps(:,borderTimes(i):borderTimes(i+1),i) = trueMap(:,borderTimes(i):borderTimes(i+1));
    % end
    % 
    nTrials = 240;
    eegEyeNumSubs = sum(ismember(EEGSubs,eyeSubs));
    
    %preallocate for loop
    allContext = nan(nTrials,1,length(behaveSubs));
    allRegLRs = nan(nTrials,1,length(behaveSubs));
    allRegBias = nan(nTrials,1,length(behaveSubs));
    allContextAll = nan(nTrials,1,length(behaveSubs));
    allRegLRsAll = nan(nTrials,1,length(behaveSubs));
    allRegBiasAll = nan(nTrials,1,length(behaveSubs));
    allIndivRegLRs = nan(nTrials,1,length(behaveSubs));
    allIndivRegBias = nan(nTrials,1,length(behaveSubs));
    allIndivRegLRsAll = nan(nTrials,1,length(behaveSubs));
    allIndivRegBiasAll = nan(nTrials,1,length(behaveSubs));
    allDoEye = zeros(length(behaveSubs),1);
    allDoEEG = zeros(length(behaveSubs),1);
    allDoEyeEEG = zeros(length(behaveSubs),1);
    meanTrialEffectEEG = nan(nTrials,size(relROIEEG.maps,3),length(EEGSubs));
    meanEffectEpochNumbers = nan(nTrials,size(relROIEEG.maps,3),length(EEGSubs));
    meanTrialEffectPupil = nan(nTrials,size(relROIEye.maps,3),length(eyeSubs));
    meanTrialEffectPGoodEEG = nan(nTrials,size(relROIEye.maps,3),eegEyeNumSubs);
    allBadBlinksEEG = nan(nTrials,1,eegEyeNumSubs);
    
    for s = 1:length(behaveSubs)
        subNum = num2str(behaveSubs(s));
        subno = behaveSubs(s);
    
        %specify whether subject is eeg/pupil/both/neither
        %a lot of indexing will use the sums of allDoEEG/Eye, as that is essentially an index for what subject we're on
        %except it skips the subjects without that type of data
        doEEG = ismember(subno,EEGSubs);
        allDoEEG(s) = doEEG;
        doEye = ismember(subno,eyeSubs);
        allDoEye(s) = doEye;
        if doEye && doEEG
            doEyeEEG = 1;
        else
            doEyeEEG = 0;
        end
        allDoEyeEEG(s) = doEyeEEG;
    
        if doEEG == 1
            %load eeg data
            eegDat = load(fullfile([eegDir,subNum,'_ALP_FILT_STIM.mat']));
            if ismember(subno,[2046,2047,2063,2086,2071,2073,2090,2094,2099,2104])
            eegDat.epochNumbers = eegDat.epochNumbers - 1;
            disp(subno)
            end
    
            %define "good" trials for EEG (after practice and non-rejected)
            isGoodEEG = false(size(eegDat.EEG.data,3),1);
            validEEG = false(size(eegDat.EEG.data,3),1);
            ind_OBCPstart=find(eegDat.epochNumbers>nPracticeTrials,1); %practice 20, random 40
            validEEG(ind_OBCPstart:end) = true;
            epochNumbers = eegDat.epochNumbers;
            epoch_OBCP=epochNumbers(ind_OBCPstart:end);
        end
        if doEye == 1
            load(fullfile([eyeDir, subNum, '.mat']));
            
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
    
            % identify blinks change from 0 to NaN
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
    
            %Find start of prediction phases for each trial
            if ~ismember(subno,[20280 2030 2079 2103])
                predStartAll = [0; find(eyeData(:,8) == 10)]; %predCueOn1
            else
                pred1StartAll = [0; find(eyeData(:,8) == 10)]; %predCueOn1
                pred2StartAll = [0; find(eyeData(:,8) == 12)]; %predCueOn2
                predStartAll = sort([pred1StartAll;pred2StartAll]); %put them all in order
            end
        
            predStart = [];
            for i = 2:(length(predStartAll))
                diffStart = predStartAll(i,1) - predStartAll(i-1,1);
                if diffStart > 5
                    predStart = [predStart predStartAll(i)];
                end
            end
            predStart = predStart';
        
            %establish time window in terms of eyeData timesteps
            timeWindow = length([-timeBeforeEye:timeAfterEye]);
            allTrialAreas = zeros([length(trialStart),timeWindow]);
            allTrialBlinks = zeros([length(trialStart),timeWindow]);
    
            %fill trial areas based on trial start times
            for i = 1:height(allTrialAreas)
                if trialStart(i) + timeAfterEye <= length(eyeData(:,1))
                    allTrialAreas(i,:) = eyeData(trialStart(i)-timeBeforeEye:trialStart(i)+timeAfterEye,9);
                    allTrialBlinks(i,:) = eyeData(trialStart(i)-timeBeforeEye:trialStart(i)+timeAfterEye,10);
                else
                    allTrialAreas(i,:) = nan;
                    allTrialBlinks(i,:) = 1;
                end
            end
    
            %subjects 20280 and 20380 have very shortened first 2 blocks for eye data so that's why they show up a lot
            %finding bad blink trials
            if ismember(subno,[20280,20380])
                allTrialAreas(1:4,:) = [];
                allTrialBlinks(1:4,:) = [];
                badBlinks = mean(allTrialBlinks,2)>blinkThresh;
                badBlinksComb = mean(allTrialBlinks,2)>combBlinkThresh;
            else
                badBlinks = mean(allTrialBlinks,2)>blinkThresh;
                badBlinks(1:60)=[];
                badBlinksComb = mean(allTrialBlinks,2)>combBlinkThresh;
                badBlinksComb(1:60)=[];
            end
        
            %adding bad blinks to list of bad blink trials (bad blinks comb is for when EEG and eye data are being combined,
            %i use a less stringent threshold there since trials are being removed for bad eeg and eye data, and i want to 
            %keep as many trials as possible
            allBadBlinks(:,:,sum(allDoEye)) = badBlinks;
            allBadBlinksComb(:,:,sum(allDoEye)) = badBlinksComb;
            %eliminate first 2 blocks and normalize trial areas
            if ~ismember(subno,[20280,20380])
                allTrialAreas(1:60,:) = [];
            end
        
            allTrialAreasNorm = reshape(nanzscore(allTrialAreas(:)),size(allTrialAreas));
        
            %calculate mean baseline for each trial
            allBaselineMeansNorm(sum(allDoEye),:) = mean(allTrialAreasNorm(:,(timeBeforeEye-baselineTimeEye+1:timeBeforeEye+1)),2);
            allBaselineMeans(sum(allDoEye),:) = mean(allTrialAreas(:,(timeBeforeEye-baselineTimeEye+1:timeBeforeEye+1)),2);
        
            subBaselineMeansNorm = allBaselineMeansNorm(sum(allDoEye),:);
        end
    
        if doEye || doEEG
            %This is where we calculate all the trial by trial parameters (pupil, eeg, LR, bias, etc)
            %load behavioral data
            allData=load(fullfile([behaveDir,'subCombined/', subNum, '_3and4BlockData.mat']));
            allData=allData.alldata;
        
            %load surprise previously calculated using model
            allModelData = load(fullfile(behaveDir,['allModelData',saveText,'/', subNum, '_allBlockData.mat']));
            allModelData = allModelData.allDataStruct;
        
            meanSurpriseCP = mean(reshape(allModelData.surpriseCP',[nTrials/2,2]),2);
            maxSurpriseCP = max(reshape(allModelData.surpriseCP',[nTrials/2,2]),[],2);
        
            meanSurpriseOB = mean(reshape(allModelData.surpriseOB',[nTrials/2,2]),2);
            maxSurpriseOB = max(reshape(allModelData.surpriseOB',[nTrials/2,2]),[],2);
        
            entropyOB=allModelData.entropyOB;
            entropyCP=allModelData.entropyCP;
        
            entropyOB=[entropyOB(1:nTrials/2),entropyOB(nTrials/2+1:end)];
            entropyOB=mean(entropyOB,2);
            entropyCP=[entropyCP(1:nTrials/2),entropyCP(nTrials/2+1:end)];
            entropyCP=mean(entropyCP,2);
        
            %add surprise and entropy to regular data structure
            if allData.condition(1) == 1
            allData.meanSurprise = [zscore(meanSurpriseCP);zscore(meanSurpriseOB)];
            allData.maxSurprise = [(maxSurpriseCP);(maxSurpriseOB)];
            allData.entropy = [zscore(entropyCP);zscore(entropyOB)];
            else
            allData.meanSurprise = [zscore(meanSurpriseOB);zscore(meanSurpriseCP)];
            allData.maxSurprise = [(maxSurpriseOB);(maxSurpriseCP)];
            allData.entropy = [zscore(entropyOB);zscore(entropyCP)];
            end
        
            %replace condition with "context" which is just nTrials long
            context = zeros(nTrials,1);
        
            if allData.condition(1)==1
                context(1:nTrials/2) = 1;
                context((nTrials/2+1):nTrials) = -1;
            elseif allData.condition(1)==2
                context(1:nTrials/2) = -1;
                context((nTrials/2+1):nTrials) = 1;
            end
            allData.context = context;
        
            %remove non nTrial long fields
            allData=rmfield(allData,'toPredict');
            allData=rmfield(allData,'allScore');
            allData=rmfield(allData,'totScore');
            allData=rmfield(allData,'sumScore');
            allData=rmfield(allData,'condition');
            allData=rmfield(allData,'predRT');
            %define parameters to compute learning rate
            predictions = allData.pred;
            outcomes = (allData.est);
            newBlock = 121;
            
            %run CLR function (done twice, one for left and right stimulus
            [LR1,UP1,subPE1] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polarHalfCorrect');
            [LR2,UP2,subPE2] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polarHalfCorrect');
            subPE = [subPE1,subPE2];
            UP = [UP1,UP2;nan,nan];
            LR = [LR1,LR2;nan,nan];
        
            %define parameters to compute learning rate (for bias/objective prediction error)
            predictions = allData.pred;
            outcomes = (allData.colorArray);
            newBlock = 121;
            
            %run CLR function (done twice, one for left and right stimulus)
            %this gives corrected prediction update so you don't have that 
            %problem from going around the circle, these are used for 
            %calculating avg trial learning/bias
            [~,~,objPE1] = computeLearningRate(outcomes(:,1),predictions(:,1),newBlock,'polarHalfCorrect');
            [~,~,objPE2] = computeLearningRate(outcomes(:,2),predictions(:,2),newBlock,'polarHalfCorrect');
            objPE = [objPE1,objPE2];
           
            %define parameters to compute learning rate
            predictions = allData.pred;
            outcomes = (allData.est);
            newBlock = 121;
            
            newUpdate = UP;
        
            allData.newLR = newUpdate./subPE;
            allData.newLR(allData.newLR > 1) = 1;
            allData.newLR(allData.newLR < 0) = 0;
        
            %Truncate learning rate to 0-1 also add it to allData
            allData.LR = LR;
            allData.LR(allData.LR > 1) = 1;
            allData.LR(allData.LR < 0) = 0;
            
            %new regressed way of calculating learning
            updateX= [zeros(nTrials,1),subPE(:,:)];
            updateY = [zeros(nTrials,1),allData.predUpdate(:,:)];
            updateXes = zeros(nTrials,nTrials+1);
            updateXes(:,1) = 1;
            
            for i = 1:nTrials
                indivRegLR(i,:) = regress(updateY(i,:)',updateX(i,:)');
                updateXes(2*i-1,i+1) = subPE(i,1);
                updateXes(2*i,i+1) = subPE(i,2);
            end
            
            updateYcol = reshape(allData.predUpdate',[],1);
            regLRs = regress(updateYcol,updateXes);
            allData.regLRs = regLRs(2:end);
            allData.indivRegLRs = indivRegLR;
        
            %Calculate regressed bias and add to allData
            allData.bias = allData.estErr(:,:)./allData.predictErr(:,:);
            biasX = [zeros(nTrials,1),allData.predictErr(:,:)];
            biasY = [zeros(nTrials,1),allData.estErr(:,:)];
            biasXes = zeros(nTrials,nTrials+1);
            biasXes(:,1) = 1;
            for i = 1:nTrials
                biasXes(2*i-1,i+1) = allData.predictErr(i,1);
                biasXes(2*i,i+1) = allData.predictErr(i,2);
                indivRegBias(i,:) = regress(biasY(i,:)',biasX(i,:)');
            end
            biasYcol = reshape(allData.estErr',[],1);
            regBias = regress(biasYcol,biasXes);
            allData.regBias = regBias(2:end);
            allData.indivRegBias = indivRegBias;
            
            %if you're regressing STP, regress STP out of learning and bias now
            if doSTPResiduals == 1
                [~,~,allData.regBias] = regress(allData.regBias,allData.meanSurprise);
                [~,~,allData.indivRegBias] = regress(allData.indivRegBias,allData.meanSurprise);
                [~,~,allData.regLRs] = regress(allData.regLRs,allData.meanSurprise);
                [~,~,allData.indivRegLRs] = regress(allData.indivRegLRs,allData.meanSurprise);
            end
    
            %add error vars to allData so good EEG trials can be selected
            allData.subPE = subPE;
            allData.objPE = objPE;
            
            if doEEG == 1
                %select out trials with good EEG data for behavioral stuff
                epoch_OBCP=epochNumbers(ind_OBCPstart:end);
                goodData=selBehav(allData, epoch_OBCP-nPracticeTrials);
            else
                goodData=allData;
            end
            %add LR, Bias and Context to "all" variables for EEG
            if doEEG
                allContext(1:size(goodData.context,1),:,s) = goodData.context;
                allRegLRs(1:size(goodData.regLRs,1),:,s) = goodData.regLRs;
                allRegBias(1:size(goodData.regBias,1),:,s) = goodData.regBias;
                allIndivRegLRs(1:size(goodData.indivRegLRs,1),:,s) = goodData.indivRegLRs;
                allIndivRegBias(1:size(goodData.indivRegBias,1),:,s) = goodData.indivRegBias;
            end
            %and eye data
            if doEye
                allContextAll(1:size(allData.context,1),:,s) = allData.context;
                allRegLRsAll(1:size(allData.regLRs,1),:,s) = allData.regLRs;
                allRegBiasAll(1:size(allData.regBias,1),:,s) = allData.regBias;
                allIndivRegLRsAll(1:size(allData.indivRegLRs,1),:,s) = allData.indivRegLRs;
                allIndivRegBiasAll(1:size(allData.indivRegBias,1),:,s) = allData.indivRegBias;
            end
    
            %regress baseline (and STP if set to) out of eye and eeg data
            if doEEG == 1
                relDataEEG = eegDat.EEG.data(:,:,validEEG);
                baseTimes_roi= [1800:1999]; % figure out what times are in "baseline" period
                if trialMeasure==3 % pull out variance that can be explained by baseline before taking dot product: 
                    resDataEEG=nan(size(relDataEEG)); % preallocate space for residual data.
                    base_roi=(squeeze(nanmean(relDataEEG(:, baseTimes_roi,:), 2)));
                    
                    for CH=1:size(base_roi, 1)
                        ch_base=base_roi(CH,:)';
                        ch_dat =squeeze(relDataEEG(CH,:,:)); % zscore original...
                        
                        for TT=1:size(ch_dat, 1) % loop through time within trial
                            if doSTPResiduals == 1
                                % include STP in xes
                                xMat=[ones(size(ch_base)),goodData.meanSurprise,ch_base];
                            else
                                % exclude STP from xes
                                xMat=[ones(size(ch_base)),ch_base];
                            end
                            [~,~,resDataEEG(CH,TT,:)] = regress(zscore(ch_dat(TT,:)'),xMat);                 
                        end
                    end
                end
            end
            if doEye == 1
                relDataEye = allTrialAreasNorm;
                if trialMeasure==3 % pull out variance that can be explained by baseline before taking dot product: 
                    resDataEye=nan(size(relDataEye)); % preallocate space for residual data.
                    for TT=1:size(relDataEye, 2) % loop through time within trial
                        if doSTPResiduals == 1
                            % include STP in xes
                            xMat=[ones(size(subBaselineMeansNorm,2),1),allData.meanSurprise,subBaselineMeansNorm']; 
                        else
                            % exclude STP from xes
                            xMat=[ones(size(subBaselineMeansNorm,2),1),subBaselineMeansNorm']; 
                        end
                        y = allTrialAreasNorm(:,TT);
                        [~,~,resDataEye(:,TT)] = regress(y,xMat); 
                    end
                end
                eye.resData = resDataEye;
            end
            
            %if subject has both good EEG and eye data, then for combined analyses you must remove bad eeg trials from 
            %eye data and then calculate mean trial effect again separately
            if doEEG && doEye
                goodData=selBehav(allData, epoch_OBCP-nPracticeTrials);
                goodEyeData=selBehav(eye, epoch_OBCP-nPracticeTrials);
                goodEyeData = goodEyeData.resData;
                epochNumbersSel = epochNumbers-nPracticeTrials;
                epochNumbersSel(epochNumbersSel<=0)=[];
                allBadBlinksEEG(1:length(epochNumbersSel),:,sum(allDoEyeEEG))=allBadBlinksComb(epochNumbersSel,:,sum(allDoEye));
            end
        
            %run loop to calculate trial by trial dot product for eeg...
            if doEEG
                for rr=1:size(relROIEEG.maps,3)
                    for t=1:size(relDataEEG,3)
                        tDataEEG=relDataEEG(:,:,t); % loop through trials.
                        if trialMeasure==1
                            meanTrialEffectEEG(t,rr,sum(allDoEEG))=nanmean(tDataEEG(relROIEEG.maps(:,:,rr)));
                        elseif trialMeasure==2   
                            % use Anne Collins dot product method:
                            meanTrialEffectEEG(t,rr,sum(allDoEEG))=tData(relROIEEG.maps(:,:,rr))'*relROIEEG.fullTMap(relROIEEG.maps(:,:,rr));
                        elseif trialMeasure==3  
                        % use dot product method on baseline-regressed residuals:
                            tDataEEG=resDataEEG(:,:,t); % loop through trials. 
                            meanTrialEffectEEG(t,rr,sum(allDoEEG))=dot(tDataEEG(relROIEEG.maps(:,:,rr)),relROIEEG.fullTMap(relROIEEG.maps(:,:,rr)));
                        end
                    end
                end
            end
            % ... and eye ...
            if doEye
                for rr=1:size(relROIEye.maps,3)
                    for t=1:size(relDataEye,1)
                        tDataEye=relDataEye(t,:); % loop through trials.
                        if trialMeasure==1
                            meanTrialEffectPupil(t,rr,sum(allDoEye))=nanmean(tDataEye(relROIEye.maps(:,:,rr)));
                        elseif trialMeasure==2   
                            % use Anne Collins dot product method:
                            meanTrialEffectPupil(t,rr,sum(allDoEye))=tDataEye(relROIEye.maps(:,:,rr))'*relROIEye.fullTMap(relROIEye.maps(:,:,rr));
                        elseif trialMeasure==3  
                        % use dot product method on baseline-regressed residuals:
                            tDataEye=resDataEye(t,:); % loop through trials. 
                            meanTrialEffectPupil(t,rr,sum(allDoEye))=dot(tDataEye(relROIEye.maps(:,:,rr)),relROIEye.fullTMap(relROIEye.maps(:,:,rr)));                        
                        end
                    end
                end
            end
            if doEEG
                %make a version of mean trial effect with trials in true order rather than EEG epoch order
                meanEffectEpochNumbers(epoch_OBCP-nPracticeTrials,:,sum(allDoEEG)) = meanTrialEffectEEG(1:length(epoch_OBCP),:,sum(allDoEEG));     
            end
            % ... and eye data with EEG removed
            if doEEG && doEye
                for rr = 1:size(relROIEye.maps,3)
                    for t=1:size(goodEyeData,1)
                        tGoodEEGData = goodEyeData(t,:);            
                        meanTrialEffectPGoodEEG(t,rr,sum(allDoEyeEEG))=dot(tGoodEEGData(relROIEye.maps(:,:,rr)),relROIEye.fullTMap(relROIEye.maps(:,:,rr)));
                    end
                end
            end
        end
        disp(subNum) 
    end
    
    beep
    %% Step 5: Quantify Relationship Between Learning/Bias and EEG/Pupil
    
    % what follows is everything that needs to be done so mean trial effects can be combined in different ways
    % eeg and eye data are different lengths, so the right trials have to be removed from both for them to be combined effectively
    % this goes for trial average leaning and bias measurements as well
    
    %sum/average meanTrialEffects
    sumTrialEffectEEG = sum(meanTrialEffectEEG(:,:,:),2);                                        %%%
    meanTrialEffectEEGforEye = meanTrialEffectEEG(:,:,ismember(EEGSubs,eyeSubs));                %%%
    sumTrialEffectEEGforEye = sumTrialEffectEEG(:,:,ismember(EEGSubs,eyeSubs));
    eegEyeNumSubs = sum(ismember(EEGSubs,eyeSubs));
    
    zMeanTrialEffectEEGforEye = [];
    zMeanTrialEffectPGoodEEG = [];
    zSumTrialEffectEEGforEye = [];
    zMeanTrialEffectEEG = [];
    slopes = [];
    
    %zscore mean trial effects
    for s = 1:eegEyeNumSubs
        zSumTrialEffectEEGforEye(:,:,s) = nanzscore(sumTrialEffectEEGforEye(:,:,s));
        zMeanTrialEffectEEGforEye(:,:,s) = nanzscore(meanTrialEffectEEGforEye(:,:,s));
        zMeanTrialEffectPGoodEEG(:,:,s) = nanzscore(meanTrialEffectPGoodEEG(:,:,s));
    end
    
    for s = 1:length(EEGSubs)
        zMeanTrialEffectEEG(:,:,s) = nanzscore(meanTrialEffectEEG(:,:,s));                       %%%
    end
    sumTrialEffectEEG = sum(zMeanTrialEffectEEG(:,:,:),2);
    
    %remove blinks from eye mean trial effects
    sumTrialEffectEye = sum(zscore(meanTrialEffectPupil(:,:,:)),2);
    sumTrialEffectEye(allBadBlinks) = nan;
    
    %sum mean trial effect for eye
    sumTrialEffectPGoodEEG = sum(zMeanTrialEffectPGoodEEG(:,:,:),2);
    
    %preallocate variables for mini loop to combine eeg and eye data and remove the right trials from behavioral variables
    allContextEEGEye = nan(nTrials,1,eegEyeNumSubs);
    allRegBiasEEGEye = nan(nTrials,1,eegEyeNumSubs);
    allRegLRsEEGEye = nan(nTrials,1,eegEyeNumSubs);
    allSurpriseEEGEye = nan(nTrials,1,eegEyeNumSubs);
    zMeanTrialEffectEEGEye = nan(nTrials,2,eegEyeNumSubs);
    zMeanTrialEffectEEGEyeSep = nan(nTrials,size(meanTrialEffectEEGforEye,2)+size(meanTrialEffectPGoodEEG,2),eegEyeNumSubs);
    
    %concatenate eeg and eye effects, one with clusters separated, one with them together
    meanTrialEffectEEGEye = [zSumTrialEffectEEGforEye,sumTrialEffectPGoodEEG];
    meanTrialEffectEEGEyeSep = [zMeanTrialEffectEEGforEye,zMeanTrialEffectPGoodEEG];
    
    clear allDoEye allDoEEG allDoEyeEEG
    
    for s = 1:length(behaveSubs)
        subno = behaveSubs(s);
    
        %find whether subject has eye/eeg data
        doEEG = ismember(subno,EEGSubs);
        allDoEEG(s) = doEEG;
        doEye = ismember(subno,eyeSubs);
        allDoEye(s) = doEye;
    
        %find if subject has both eeg and eye data and increment counter if so
        if doEye && doEEG
            doEyeEEG = 1;
        else
            doEyeEEG = 0;
        end
        allDoEyeEEG(s) = doEyeEEG;
        if doEEG == 1
        %remove bad blinks with more permissive threshold to compare with combined EEG + Eye data
        zMeanTrialEffectEEGEye(1:sum(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,sum(allDoEyeEEG)) = nanzscore(meanTrialEffectEEGEye(find(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,sum(allDoEyeEEG)));
        zMeanTrialEffectEEGEyeSep(1:sum(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,sum(allDoEyeEEG)) = nanzscore(meanTrialEffectEEGEyeSep(find(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,sum(allDoEyeEEG)));
        allContextEEGEye(1:sum(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s) = allContext(find(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s);
        allRegBiasEEGEye(1:sum(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s) = allIndivRegBias(find(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s);
        allRegLRsEEGEye(1:sum(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s) = allIndivRegLRs(find(allBadBlinksEEG(:,:,sum(allDoEyeEEG))==0),:,s);
        end
    end
    
    sumTrialEffectEEGEye = sum(zMeanTrialEffectEEGEye,2);
    
    numQuant = 5;
    
    %preallocate variables for loop, mostly just LR and Bias for different data lengths
    meanQuantValsEEG = nan(length(EEGSubs),numQuant);
    meanQuantValsEye = nan(length(eyeSubs),numQuant);
    meanAllLRs = mean(allIndivRegLRs,2);
        meanAllLRs(meanAllLRs>1) = 1;
        meanAllLRs(meanAllLRs<0) = 0;
    meanAllBias = mean(allIndivRegBias,2);
        meanAllBias(meanAllBias>1) = 1;
        meanAllBias(meanAllBias<0) = 0;
    
    meanAllLRsEye = mean(allIndivRegLRsAll,2);
        meanAllLRsEye(meanAllLRsEye>1) = 1;
        meanAllLRsEye(meanAllLRsEye<0) = 0;
    meanAllBiasEye = mean(allIndivRegBiasAll,2);
        meanAllBiasEye(meanAllBiasEye>1) = 1;
        meanAllBiasEye(meanAllBiasEye<0) = 0;
    
    meanAllLRsEEGEye = mean(allRegLRsEEGEye,2);
        meanAllLRsEEGEye(meanAllLRsEEGEye>1) = 1;
        meanAllLRsEEGEye(meanAllLRsEEGEye<0) = 0;
    meanAllBiasEEGEye = mean(allRegBiasEEGEye,2);
        meanAllBiasEEGEye(meanAllBiasEEGEye>1) = 1;
        meanAllBiasEEGEye(meanAllBiasEEGEye<0) = 0;
        
    quantSubsCPEEG = nan(length(EEGSubs),nTrials/2);
    quantSubsOBEEG = nan(length(EEGSubs),nTrials/2);
    quantSubsCPEye = nan(length(eyeSubs),nTrials/2);
    quantSubsOBEye = nan(length(eyeSubs),nTrials/2);
    
    LRCPslopeEEG = nan(length(EEGSubs),2);
    biasCPslopeEEG = nan(length(EEGSubs),2);
    LROBslopeEEG = nan(length(EEGSubs),2);
    biasOBslopeEEG = nan(length(EEGSubs),2);
    biasAllslopeEEG = nan(length(EEGSubs),2);
    
    LRCPslopeEye = nan(length(eyeSubs),2);
    biasCPslopeEye = nan(length(eyeSubs),2);
    LROBslopeEye = nan(length(eyeSubs),2);
    biasOBslopeEye = nan(length(eyeSubs),2);
    biasAllslopeEye = nan(length(eyeSubs),2);
    
    LRCPslopeEEGSep = nan(length(EEGSubs),2);
    biasCPslopeEEGSep = nan(length(EEGSubs),2);
    LROBslopeEEGSep = nan(length(EEGSubs),2);
    biasOBslopeEEGSep = nan(length(EEGSubs),2);
    biasAllslopeEEGSep = nan(length(EEGSubs),2);
    
    LRCPslopeEEGEye = nan(eegEyeNumSubs,2);
    biasCPslopeEEGEye = nan(eegEyeNumSubs,2);
    LROBslopeEEGEye = nan(eegEyeNumSubs,2);
    biasOBslopeEEGEye = nan(eegEyeNumSubs,2);
    biasAllslopeEEGEye = nan(eegEyeNumSubs,2);
    
    biasAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    LRAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    biasCPAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    LRCPAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    biasOBAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    LROBAllSlopes = nan(size(zMeanTrialEffectEEG,2),1,length(EEGSubs));
    
    varSlopes = nan(size(zMeanTrialEffectEEG,2),6,length(EEGSubs));
    clear allDoEye allDoEEG allDoEyeEEG
    
    for s = [1:length(behaveSubs)]
        disp(s);
        subNum = num2str(behaveSubs(s));
        subno = behaveSubs(s);
    
        %find whether subject has good eeg/eye/both data, and increment counters (allDo variables) if true
        doEEG = ismember(subno,EEGSubs);
        allDoEEG(s) = doEEG;
        doEye = ismember(subno,eyeSubs);
        allDoEye(s) = doEye;
        if doEye && doEEG
            doEyeEEG = 1;
        else
            doEyeEEG = 0;
        end
        allDoEyeEEG(s) = doEyeEEG;
    
        if doEEG
            %quantiling and stuff for EEG data only
            numTrials = length(find(isfinite(sumTrialEffectEEG(:,1,sum(allDoEEG)))));
            allData=load(fullfile([behaveDir,'subCombined/', subNum, '_3and4BlockData.mat']));
            allData = allData.alldata;
            condition = allData.condition(1);
        
            %calculate quantiles for sum of mean trial effects
            borders = quantile(sumTrialEffectEEG(1:numTrials,:,sum(allDoEEG)),numQuant-1);
            bordersCP = quantile(sumTrialEffectEEG(allContext(:,:,s)==1,:,sum(allDoEEG)),numQuant-1);
            bordersOB = quantile(sumTrialEffectEEG(allContext(:,:,s)==-1,:,sum(allDoEEG)),numQuant-1);
            quantiles = zeros(numTrials,1);
        
            quantilesCP = [];
            quantilesOB = [];
        
            %assign each trial to a quantile (1-numQuant) based on sum (lowest = 1)
            for q = numQuant:-1:1
                if q<numQuant
                quantilesCP(sumTrialEffectEEG(allContext(:,:,s)==1,:,sum(allDoEEG))<=bordersCP(q)) = q;
                quantilesOB(sumTrialEffectEEG(allContext(:,:,s)==-1,:,sum(allDoEEG))<=bordersOB(q)) = q;
                quantiles(sumTrialEffectEEG(1:numTrials,:,sum(allDoEEG))<=borders(q)) = q;
                else
                quantilesCP(sumTrialEffectEEG(allContext(:,:,s)==1,:,sum(allDoEEG))>bordersCP(q-1)) = q;
                quantilesOB(sumTrialEffectEEG(allContext(:,:,s)==-1,:,sum(allDoEEG))>bordersOB(q-1)) = q;
                quantiles(sumTrialEffectEEG(1:numTrials,:,sum(allDoEEG))>borders(q-1)) = q;
                end
            end
        
            % add number of zeros of trials in first block (CP/OB) to beginning of second block (OB/CP) so everything lines up right
            if allContext(1,:,s) == 1
                quantilesOB = [zeros(length(quantilesCP),1);quantilesOB'];
                quantilesCP = quantilesCP';
                disp(1);
            else
                quantilesCP = [zeros(length(quantilesOB),1);quantilesCP'];
                quantilesOB = quantilesOB';
                disp(2);
            end
           
            for q = numQuant:-1:1
                % find which CP and OB trials are within quantile
                quantCP = quantilesCP==q;
                quantOB = quantilesOB==q;
    
                % find LR & bias for trials within the quantile
                quantCPLRs = meanAllLRs(quantCP,:,s);
                quantCPBias= meanAllBias(quantCP,:,s);
                quantCPVals= sumTrialEffectEEG(quantCP,:,sum(allDoEEG));
                quantOBLRs = meanAllLRs(quantOB,:,s);
                quantOBBias= meanAllBias(quantOB,:,s);
                quantOBVals= sumTrialEffectEEG(quantOB,:,sum(allDoEEG));
                
                %take mean LR and Bias for each quantile
                meanQuantCPLRsEEG(sum(allDoEEG),q) = nanmean(quantCPLRs);
                meanQuantCPBiasEEG(sum(allDoEEG),q) = nanmean(quantCPBias);
                meanQuantCPValsEEG(sum(allDoEEG),q) = mean(quantCPVals);
                meanQuantOBLRsEEG(sum(allDoEEG),q) = nanmean(quantOBLRs);
                meanQuantOBBiasEEG(sum(allDoEEG),q) = nanmean(quantOBBias);
                meanQuantOBValsEEG(sum(allDoEEG),q) = mean(quantOBVals); 
            end
            
            quantilesCP(isnan(quantilesCP)) = [];
            quantilesOB(isnan(quantilesOB)) = [];
            quantilesCP(quantilesCP==0) = [];
            quantilesOB(quantilesOB==0) = [];
        
            %save info of what trial is in which quantile for each subject
            quantSubsCPEEG(sum(allDoEEG),1:length(quantilesCP)) = quantilesCP;
            quantSubsOBEEG(sum(allDoEEG),1:length(quantilesOB)) = quantilesOB;
    
            %testing hypothesis 5
                
            %Rank order column 1 physio effect, col 2 LR, col3 bias
            rankOrderCP = sumTrialEffectEEG(allContext(:,:,s)==1,:,sum(allDoEEG));
            rankOrderCP(:,2) = meanAllLRs(allContext(:,:,s)==1,:,s);
            rankOrderCP(:,3) = meanAllBias(allContext(:,:,s)==1,:,s);
        
            %sort rank order based on physio effect
            rankOrderCPSorted = sortrows(rankOrderCP,1,"ascend");
        
            %repeat for oddball
            rankOrderOB = sumTrialEffectEEG(allContext(:,:,s)==-1,:,sum(allDoEEG));
            rankOrderOB(:,2) = meanAllLRs(allContext(:,:,s)==-1,:,s);
            rankOrderOB(:,3) = meanAllBias(allContext(:,:,s)==-1,:,s);
            rankOrderOBSorted = sortrows(rankOrderOB,1,"ascend");
        
            xIntCP = ones(length(find(allContext(:,:,s)==1)),1);
            xIntOB = ones(length(find(allContext(:,:,s)==-1)),1);
        
            %regress bias against each trial's physiologial effect (how much eeg cluster) rank in block 
            LRCPslopeEEG(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,(1:length(xIntCP))']);
            biasCPslopeEEG(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,(1:length(xIntCP))']);
            LROBslopeEEG(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,(1:length(xIntOB))']);
            biasOBslopeEEG(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,(1:length(xIntOB))']);
            biasAllslopeEEG(sum(allDoEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[(1:length(xIntCP))';(1:length(xIntOB))']]);
         
            %OR regress against raw signal strength rather than rank order (gives better results probably because of distribution
            LRCPslopeEEG(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,rankOrderCPSorted(:,1)]);
            err = rankOrderCPSorted;
            biasCPslopeEEG(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,rankOrderCPSorted(:,1)]);
            LROBslopeEEG(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,rankOrderOBSorted(:,1)]);
            biasOBslopeEEG(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,rankOrderOBSorted(:,1)]);
            biasAllslopeEEG(sum(allDoEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[rankOrderCPSorted(:,1);rankOrderOBSorted(:,1)]]);
            for rr = 1:size(zMeanTrialEffectEEG,2)
                %calculate regression results for each eeg cluster and pupil data separately
                %Rank order column 1 physio effect, col 2 LR, col3 bias
                rankOrderCP = zMeanTrialEffectEEG(allContext(:,:,s)==1,rr,sum(allDoEEG));
                rankOrderCP(:,2) = meanAllLRs(allContext(:,:,s)==1,:,s);
                rankOrderCP(:,3) = meanAllBias(allContext(:,:,s)==1,:,s);
            
                %sort rank order based on physio effect 
                rankOrderCPSorted = sortrows(rankOrderCP,1,"ascend");
            
                %repeat for oddball
                rankOrderOB = zMeanTrialEffectEEG(allContext(:,:,s)==-1,rr,sum(allDoEEG));
                rankOrderOB(:,2) = meanAllLRs(allContext(:,:,s)==-1,:,s);
                rankOrderOB(:,3) = meanAllBias(allContext(:,:,s)==-1,:,s);
                rankOrderOBSorted = sortrows(rankOrderOB,1,"ascend");
            
                xIntCP = ones(length(find(allContext(:,:,s)==1)),1);
                xIntOB = ones(length(find(allContext(:,:,s)==-1)),1);
            
                %regress bias against each trial's physiologial effect (how much eeg cluster) rank in block 
                LRCPslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,(1:length(xIntCP))']);
                biasCPslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,(1:length(xIntCP))']);
                LROBslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,(1:length(xIntOB))']);
                biasOBslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,(1:length(xIntOB))']);
                biasAllslopeEEGSep(sum(allDoEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[(1:length(xIntCP))';(1:length(xIntOB))']]);
             
                %OR regress against raw signal strength rather than rank order (gives better results probably because of distribution
                LRCPslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,rankOrderCPSorted(:,1)]);
                biasCPslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,rankOrderCPSorted(:,1)]);
                LROBslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,rankOrderOBSorted(:,1)]);
                biasOBslopeEEGSep(sum(allDoEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,rankOrderOBSorted(:,1)]);
                biasAllslopeEEGSep(sum(allDoEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[rankOrderCPSorted(:,1);rankOrderOBSorted(:,1)]]);
                
                varSlopes(rr,:,sum(allDoEEG)) = [biasCPslopeEEGSep(sum(allDoEEG),2),biasOBslopeEEGSep(sum(allDoEEG),2),biasCPslopeEEGSep(sum(allDoEEG),2)+biasOBslopeEEGSep(sum(allDoEEG),2),LRCPslopeEEGSep(sum(allDoEEG),2),LROBslopeEEGSep(sum(allDoEEG),2),LRCPslopeEEGSep(sum(allDoEEG),2)-LROBslopeEEGSep(sum(allDoEEG),2)];
                biasAllSlopes(rr,:,sum(allDoEEG)) = biasCPslopeEEGSep(sum(allDoEEG),2)+biasOBslopeEEGSep(sum(allDoEEG),2);
                LRAllSlopes(rr,:,sum(allDoEEG)) = LRCPslopeEEGSep(sum(allDoEEG),2)-LROBslopeEEGSep(sum(allDoEEG),2);
                biasCPAllSlopes(rr,:,sum(allDoEEG)) = biasCPslopeEEGSep(sum(allDoEEG),2);
                LRCPAllSlopes(rr,:,sum(allDoEEG)) = LRCPslopeEEGSep(sum(allDoEEG),2);
                biasOBAllSlopes(rr,:,sum(allDoEEG)) = biasOBslopeEEGSep(sum(allDoEEG),2);
                LROBAllSlopes(rr,:,sum(allDoEEG)) = LROBslopeEEGSep(sum(allDoEEG),2);
            end
        end
        if doEye
            %quantiling and stuff for eye data only
            numTrials = nTrials;
            allData=load(fullfile([behaveDir,'subCombined/', subNum, '_3and4BlockData.mat']));
            allData = allData.alldata;
            condition = allData.condition(1);
        
            %calculate quantiles for sum of mean trial effects
            borders = quantile(sumTrialEffectEye(1:numTrials,:,sum(allDoEye)),numQuant-1);
            bordersCP = quantile(sumTrialEffectEye(allContextAll(:,:,s)==1,:,sum(allDoEye)),numQuant-1);
            bordersOB = quantile(sumTrialEffectEye(allContextAll(:,:,s)==-1,:,sum(allDoEye)),numQuant-1);
            quantiles = zeros(numTrials,1);
    
            quantilesCP = zeros(numTrials/2,1);
            quantilesOB = zeros(numTrials/2,1);
        
            %assign each trial to a quantile (1-numQuant) based on sum (lowest = 1)
            for q = numQuant:-1:1
                if q<numQuant
                quantilesCP(sumTrialEffectEye(allContextAll(:,:,s)==1,:,sum(allDoEye))<=bordersCP(q)) = q;
                quantilesOB(sumTrialEffectEye(allContextAll(:,:,s)==-1,:,sum(allDoEye))<=bordersOB(q)) = q;
                quantiles(sumTrialEffectEye(1:numTrials,:,sum(allDoEye))<=borders(q)) = q;
                else
                quantilesCP(sumTrialEffectEye(allContextAll(:,:,s)==1,:,sum(allDoEye))>bordersCP(q-1)) = q;
                quantilesOB(sumTrialEffectEye(allContextAll(:,:,s)==-1,:,sum(allDoEye))>bordersOB(q-1)) = q;
                quantiles(sumTrialEffectEye(1:numTrials,:,sum(allDoEye))>borders(q-1)) = q;
                end
            end
        
            % add number of zeros of trials in first block (CP/OB) to beginning of second block (OB/CP) so everything lines up right
    
            if allContextAll(1,:,s) == 1
                quantilesOB = [zeros(length(quantilesCP),1);quantilesOB];
                quantilesCP = quantilesCP;
                disp(1);
            else
                quantilesCP = [zeros(length(quantilesOB),1);quantilesCP];
                quantilesOB = quantilesOB;
                disp(2);
            end
    
            for q = numQuant:-1:1
                % find which CP and OB trials are within quantile
                quantCP = quantilesCP==q;
                quantOB = quantilesOB==q;
                % find LR & bias for trials within the quantile
                quantCPLRs = meanAllLRsEye(quantCP,:,s);
                quantCPBias= meanAllBiasEye(quantCP,:,s);
                quantCPVals= sumTrialEffectEye(quantCP,:,sum(allDoEye));
                quantOBLRs = meanAllLRsEye(quantOB,:,s);
                quantOBBias= meanAllBiasEye(quantOB,:,s);
                quantOBVals= sumTrialEffectEye(quantOB,:,sum(allDoEye));
                
                %take mean LR and Bias for each quantile
                meanQuantCPLRsEye(sum(allDoEye),q) = nanmean(quantCPLRs);
                meanQuantCPBiasEye(sum(allDoEye),q) = nanmean(quantCPBias);
                meanQuantCPValsEye(sum(allDoEye),q) = mean(quantCPVals);
                meanQuantOBLRsEye(sum(allDoEye),q) = nanmean(quantOBLRs);
                meanQuantOBBiasEye(sum(allDoEye),q) = nanmean(quantOBBias);
                meanQuantOBValsEye(sum(allDoEye),q) = mean(quantOBVals); 
            end
            
            quantilesCP(isnan(quantilesCP)) = [];
            quantilesOB(isnan(quantilesOB)) = [];
            quantilesCP(quantilesCP==0) = [];
            quantilesOB(quantilesOB==0) = [];
        
            %save info of what trial is in which quantile for each subject
        %     if condition == 1
                quantSubsCPEye(sum(allDoEye),1:length(quantilesCP)) = quantilesCP;
                quantSubsOBEye(sum(allDoEye),1:length(quantilesOB)) = quantilesOB;
    
            %testing hypothesis 5
            %Rank order column 1 physio effect, col 2 LR, col3 bias
            rankOrderCP = sumTrialEffectEye(allContextAll(:,:,s)==1,:,sum(allDoEye));
            rankOrderCP(:,2) = meanAllLRsEye(allContextAll(:,:,s)==1,:,s);
            rankOrderCP(:,3) = meanAllBiasEye(allContextAll(:,:,s)==1,:,s);
        
            %sort rank order based on physio effect
            rankOrderCPSorted = sortrows(rankOrderCP,1,"ascend");
        
            %repeat for oddball
            rankOrderOB = sumTrialEffectEye(allContextAll(:,:,s)==-1,:,sum(allDoEye));
            rankOrderOB(:,2) = meanAllLRsEye(allContextAll(:,:,s)==-1,:,s);
            rankOrderOB(:,3) = meanAllBiasEye(allContextAll(:,:,s)==-1,:,s);
            rankOrderOBSorted = sortrows(rankOrderOB,1,"ascend");
        
            xIntCP = ones(length(find(allContextAll(:,:,s)==1)),1);
            xIntOB = ones(length(find(allContextAll(:,:,s)==-1)),1);
        
            %regress bias against each trial's physiologial effect (how much eeg cluster) rank in block 
            LRCPslopeEye(sum(allDoEye),:) = regress(rankOrderCPSorted(:,2),[xIntCP,(1:length(xIntCP))']);
            biasCPslopeEye(sum(allDoEye),:) = regress(rankOrderCPSorted(:,3),[xIntCP,(1:length(xIntCP))']);
            LROBslopeEye(sum(allDoEye),:) = regress(rankOrderOBSorted(:,2),[xIntOB,(1:length(xIntOB))']);
            biasOBslopeEye(sum(allDoEye),:) = regress(rankOrderOBSorted(:,3),[xIntOB,(1:length(xIntOB))']);
            biasAllslopeEye(sum(allDoEye),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[(1:length(xIntCP))';(1:length(xIntOB))']]);
         
            %OR regress against raw signal strength rather than rank order (gives better results probably because of distribution
            LRCPslopeEye(sum(allDoEye),:) = regress(rankOrderCPSorted(:,2),[xIntCP,rankOrderCPSorted(:,1)]);
            biasCPslopeEye(sum(allDoEye),:) = regress(rankOrderCPSorted(:,3),[xIntCP,rankOrderCPSorted(:,1)]);
            LROBslopeEye(sum(allDoEye),:) = regress(rankOrderOBSorted(:,2),[xIntOB,rankOrderOBSorted(:,1)]);
            biasOBslopeEye(sum(allDoEye),:) = regress(rankOrderOBSorted(:,3),[xIntOB,rankOrderOBSorted(:,1)]);
            biasAllslopeEye(sum(allDoEye),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[rankOrderCPSorted(:,1);rankOrderOBSorted(:,1)]]);
        end
        if doEyeEEG
            %quantiling and stuff for eeg + eye data where eeg clusters have been summed and equally weighted with pupil 
            numTrials = length(find(isfinite(sumTrialEffectEEGEye(:,1,sum(allDoEyeEEG)))));
            allData=load(fullfile([behaveDir,'subCombined/', subNum, '_3and4BlockData.mat']));
            allData = allData.alldata;
            condition = allData.condition(1);
        
            %calculate quantiles for sum of mean trial effects
            borders = quantile(sumTrialEffectEEGEye(1:numTrials,:,sum(allDoEyeEEG)),numQuant-1);
            bordersCP = quantile(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==1,:,sum(allDoEyeEEG)),numQuant-1);
            bordersOB = quantile(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==-1,:,sum(allDoEyeEEG)),numQuant-1);
            quantiles = zeros(numTrials,1);
        
            quantilesCP = [];
            quantilesOB = [];
        
            %assign each trial to a quantile (1-numQuant) based on sum (lowest = 1)
            for q = numQuant:-1:1
                if q<numQuant
                quantilesCP(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==1,:,sum(allDoEyeEEG))<=bordersCP(q)) = q;
                quantilesOB(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==-1,:,sum(allDoEyeEEG))<=bordersOB(q)) = q;
                quantiles(sumTrialEffectEEGEye(1:numTrials,:,sum(allDoEyeEEG))<=borders(q)) = q;
                else
                quantilesCP(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==1,:,sum(allDoEyeEEG))>bordersCP(q-1)) = q;
                quantilesOB(sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==-1,:,sum(allDoEyeEEG))>bordersOB(q-1)) = q;
                quantiles(sumTrialEffectEEGEye(1:numTrials,:,sum(allDoEyeEEG))>borders(q-1)) = q;
                end
            end
        
            % add number of zeros of trials in first block (CP/OB) to beginning of second block (OB/CP) so everything lines up right
            if allContextEEGEye(1,:,s) == 1
                quantilesOB = [zeros(length(quantilesCP),1);quantilesOB'];
                quantilesCP = quantilesCP';
                disp(1);
            else
                quantilesCP = [zeros(length(quantilesOB),1);quantilesCP'];
                quantilesOB = quantilesOB';
                disp(2);
            end
           
            for q = numQuant:-1:1
                % find which CP and OB trials are within quantile
                quantCP = quantilesCP==q;
                quantOB = quantilesOB==q;
    
                % find LR & bias for trials within the quantile
                quantCPLRs = meanAllLRsEEGEye(quantCP,:,s);
                quantCPBias= meanAllBiasEEGEye(quantCP,:,s);
                quantCPVals= sumTrialEffectEEGEye(quantCP,:,sum(allDoEyeEEG));
                quantOBLRs = meanAllLRsEEGEye(quantOB,:,s);
                quantOBBias= meanAllBiasEEGEye(quantOB,:,s);
                quantOBVals= sumTrialEffectEEGEye(quantOB,:,sum(allDoEyeEEG));
                
                %take mean LR and Bias for each quantile
                meanQuantCPLRsEEGEye(sum(allDoEyeEEG),q) = nanmean(quantCPLRs);
                meanQuantCPBiasEEGEye(sum(allDoEyeEEG),q) = nanmean(quantCPBias);
                meanQuantCPValsEEGEye(sum(allDoEyeEEG),q) = mean(quantCPVals);
                meanQuantOBLRsEEGEye(sum(allDoEyeEEG),q) = nanmean(quantOBLRs);
                meanQuantOBBiasEEGEye(sum(allDoEyeEEG),q) = nanmean(quantOBBias);
                meanQuantOBValsEEGEye(sum(allDoEyeEEG),q) = mean(quantOBVals); 
            end
            
            quantilesCP(isnan(quantilesCP)) = [];
            quantilesOB(isnan(quantilesOB)) = [];
            quantilesCP(quantilesCP==0) = [];
            quantilesOB(quantilesOB==0) = [];
        
            %save info of what trial is in which quantile for each subject
            quantSubsCPEEGEye(sum(allDoEyeEEG),1:length(quantilesCP)) = quantilesCP;
            quantSubsOBEEGEye(sum(allDoEyeEEG),1:length(quantilesOB)) = quantilesOB;
    
            %testing hypothesis 5
                
            %Rank order column 1 physio effect, col 2 LR, col3 bias
            rankOrderCP = sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==1,:,sum(allDoEyeEEG));
            rankOrderCP(:,2) = meanAllLRsEEGEye(allContextEEGEye(:,:,s)==1,:,s);
            rankOrderCP(:,3) = meanAllBiasEEGEye(allContextEEGEye(:,:,s)==1,:,s);
        
            %sort rank order based on physio effect
            rankOrderCPSorted = sortrows(rankOrderCP,1,"ascend");
        
            %repeat for oddball
            rankOrderOB = sumTrialEffectEEGEye(allContextEEGEye(:,:,s)==-1,:,sum(allDoEyeEEG));
            rankOrderOB(:,2) = meanAllLRsEEGEye(allContextEEGEye(:,:,s)==-1,:,s);
            rankOrderOB(:,3) = meanAllBiasEEGEye(allContextEEGEye(:,:,s)==-1,:,s);
            rankOrderOBSorted = sortrows(rankOrderOB,1,"ascend");
        
            xIntCP = ones(length(find(allContextEEGEye(:,:,s)==1)),1);
            xIntOB = ones(length(find(allContextEEGEye(:,:,s)==-1)),1);
        
            %regress bias against each trial's physiologial effect (how much eeg cluster) rank in block 
            LRCPslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,(1:length(xIntCP))']);
            biasCPslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,(1:length(xIntCP))']);
            LROBslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,(1:length(xIntOB))']);
            biasOBslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,(1:length(xIntOB))']);
            biasAllslopeEEGEye(sum(allDoEyeEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[(1:length(xIntCP))';(1:length(xIntOB))']]);
         
            %OR regress against raw signal strength rather than rank order (gives better results probably because of distribution
            LRCPslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderCPSorted(:,2),[xIntCP,rankOrderCPSorted(:,1)]);
            biasCPslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderCPSorted(:,3),[xIntCP,rankOrderCPSorted(:,1)]);
            LROBslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderOBSorted(:,2),[xIntOB,rankOrderOBSorted(:,1)]);
            biasOBslopeEEGEye(sum(allDoEyeEEG),:) = regress(rankOrderOBSorted(:,3),[xIntOB,rankOrderOBSorted(:,1)]);
            biasAllslopeEEGEye(sum(allDoEyeEEG),:) = regress([rankOrderCPSorted(:,3);rankOrderOBSorted(:,3)],[[xIntCP;xIntOB],[rankOrderCPSorted(:,1);rankOrderOBSorted(:,1)]]);
        end
    end
    
    % calculate p values of bias/lr cp/ob slopes for eeg ...
    [~,pLRCPEEG] = ttest(LRCPslopeEEG(:,2));
    [~,pLROBEEG] = ttest(LROBslopeEEG(:,2));
    [~,pBiasCPEEG] = ttest(biasCPslopeEEG(:,2));
    [~,pBiasOBEEG] = ttest(biasOBslopeEEG(:,2));
    [~,pBiasAllEEG] = ttest(biasAllslopeEEG(:,2));
    [~,pBiasSumEEG] = ttest(biasCPslopeEEG(:,2)+biasOBslopeEEG(:,2));
    [~,pLRDiffEEG] = ttest(LRCPslopeEEG(:,2)-LROBslopeEEG(:,2));
    
    [~,pBiasSumEEGSep] = ttest(squeeze(biasAllSlopes)');
    [~,pLRDiffEEGSep] = ttest(squeeze(LRAllSlopes)');
    [~,pBiasCPEEGSep] = ttest(squeeze(biasCPAllSlopes)');
    [~,pLRCPEEGSep] = ttest(squeeze(LRCPAllSlopes)');
    [~,pBiasOBEEGSep] = ttest(squeeze(biasOBAllSlopes)');
    [~,pLROBEEGSep] = ttest(squeeze(LROBAllSlopes)');
    meanVarSlopes = mean(varSlopes,3);
    varP = [pBiasCPEEGSep',pBiasOBEEGSep',pBiasSumEEGSep',pLRCPEEGSep',pLROBEEGSep',pLRDiffEEGSep'];
    
    % ... eye data ...
    [~,pLRCPEye] = ttest(LRCPslopeEye(:,2));
    [~,pLROBEye] = ttest(LROBslopeEye(:,2));
    [~,pBiasCPEye] = ttest(biasCPslopeEye(:,2));
    [~,pBiasOBEye] = ttest(biasOBslopeEye(:,2));
    [~,pBiasAllEye] = ttest(biasAllslopeEye(:,2));
    [~,pBiasSumEye] = ttest(biasCPslopeEye(:,2)+biasOBslopeEye(:,2));
    [~,pLRDiffEye] = ttest(LRCPslopeEye(:,2)-LROBslopeEye(:,2));
    
    % ... and both forms of data combined
    [~,pLRCPEEGEye] = ttest(LRCPslopeEEGEye(:,2));
    [~,pLROBEEGEye] = ttest(LROBslopeEEGEye(:,2));
    [~,pBiasCPEEGEye] = ttest(biasCPslopeEEGEye(:,2));
    [~,pBiasOBEEGEye] = ttest(biasOBslopeEEGEye(:,2));
    [~,pBiasAllEEGEye] = ttest(biasAllslopeEEGEye(:,2));
    [~,pBiasSumEEGEye] = ttest(biasCPslopeEEGEye(:,2)+biasOBslopeEEGEye(:,2));
    [~,pLRDiffEEGEye] = ttest(LRCPslopeEEGEye(:,2)-LROBslopeEEGEye(:,2));
    
    %set up regression for testing hypothesis 5 
    subContext = [ones(eegEyeNumSubs,1);-1*ones(eegEyeNumSubs,1)];
    allLRSlopes = [LRCPslopeEEGEye(:,2);LROBslopeEEGEye(:,2)];
    allBiasSlopes = [biasAllslopeEEGEye(:,2);biasAllslopeEEGEye(:,2)];
    
    X = [ones(length(subContext),1),allBiasSlopes,subContext,allBiasSlopes.*subContext];
    useDummy = 0;
    if useDummy == 1
        dummyMat = [eye(eegEyeNumSubs);eye(eegEyeNumSubs)];
        zX = [X(:,1),zscore(X(:,2:end)),dummyMat];
    else
        zX = [X(:,1),zscore(X(:,2:end))];
    end
    %regressing Bias, condition, bias*condition against learning
    [slopeBsStats] = regstats(allLRSlopes,zX(:,2:end));
    
    disp(slopeBsStats.tstat.t)
    disp(slopeBsStats.tstat.pval)
    
    posClustLabels = ["1st +","2nd +","3rd +","4th +","5th +","6th +","7th +","8th +"];
    negClustLabels = ["1st -","2nd -","3rd -","4th -","5th -","6th -","7th -","8th -"];
    
    clusterLabels = [posClustLabels(1:length(gPosEEG)),negClustLabels(1:length(gNegEEG)),"Pupil"]; 
end

%%  BONUS FIGURES
% find correlation between eeg clusters and pupil dilation
if ~isempty(gSig)
    clear eegCorrs eegEyeCorrs
    for s = 1:length(EEGSubs)
    eegCorrs(:,:,s) = corr(meanTrialEffectEEG(:,:,s),'rows','complete');
    end
    meanEEGCorrs = mean(eegCorrs,3);
    for s = 1:eegEyeNumSubs
    eegEyeCorrs(:,:,s) = corr(meanTrialEffectEEGEyeSep(:,:,s),'rows','complete');
    end
    meanEEGEyeCorrs = mean(eegEyeCorrs,3);
    allMeanCorrs = [[meanEEGCorrs;meanEEGEyeCorrs(end,1:end-1)],meanEEGEyeCorrs(:,end)];
    allCorrsSig = [[ttest(eegCorrs,0,'dim',3);ttest(eegEyeCorrs(end,1:end-1,:),0,'dim',3)],ttest(eegEyeCorrs(:,end,:),0,'dim',3)];
    allCorrsSig01 = [[ttest(eegCorrs,0,'dim',3,'alpha',0.01);ttest(eegEyeCorrs(end,1:end-1,:),0,'dim',3,'alpha',0.01)],ttest(eegEyeCorrs(:,end,:),0,'dim',3,'alpha',0.01)];
    allCorrsSig001 = [[ttest(eegCorrs,0,'dim',3,'alpha',0.001);ttest(eegEyeCorrs(end,1:end-1,:),0,'dim',3,'alpha',0.001)],ttest(eegEyeCorrs(:,end,:),0,'dim',3,'alpha',0.001)];
    allCorrsSig0001 = [[ttest(eegCorrs,0,'dim',3,'alpha',0.0001);ttest(eegEyeCorrs(end,1:end-1,:),0,'dim',3,'alpha',0.0001)],ttest(eegEyeCorrs(:,end,:),0,'dim',3,'alpha',0.0001)];
    
    %image plot for correlation between eeg clusters and pupil dilation
    % figure
    % imagesc(allMeanCorrs)
    % yticks(1:length(allMeanCorrs))
    % yticklabels([clusterLabels])
    % xticks(1:length(allMeanCorrs))
    % xticklabels([clusterLabels])
    % 
    % colorbar
    % colormap parula
    % caxis([-1,1]);
    % [R, C] = ndgrid(1:size(allMeanCorrs,1), 1:size(allMeanCorrs,1));
    % R = R(:); C = C(:) - 1/10;
    % 
    % clear mask
    % vals = round(allMeanCorrs,2);
    % mask = vals <= 0;
    % text(C(mask), R(mask), string(vals(mask)), 'color', 'w',"FontName","Arial","FontWeight","bold","FontSize",14)
    % text(C(~mask), R(~mask), string(vals(~mask)), 'color', 'k',"FontName","Arial","FontWeight","bold","FontSize",14)
    %         set(gca, 'box', 'off')
    %         set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
    % 
    % if doSTPResiduals == 0
    %     fig = gcf;
    %     figName = append("Figure_S4_",figTime,'.eps');
    %     figLoc = append(figDir,figName);
    %     exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    %     figName = append("Figure_S4_",figTime,'.png');
    %     figLoc = append(figDir,figName);
    %     saveas(fig,figLoc)
    % else
    %     fig = gcf;
    %     figName = append("Figure_S4_Residual_",figTime,'.eps');
    %     figLoc = append(figDir,figName);
    %     exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    %     figName = append("Figure_S4_Residual",figTime,'.png');
    %     figLoc = append(figDir,figName);
    %     saveas(fig,figLoc)
    % end
    % figure
    % imagesc(allMeanCorrs)
    % yticks(1:length(allMeanCorrs))
    % xticks(1:length(allMeanCorrs))
    % yticklabels([clusterLabels])
    % xticklabels([clusterLabels])
    % 
    % colorbar
    % colormap parula
    % caxis([-1,1]);
    % [i,j]=find(allCorrsSig==1);
    % ms=10;
    % hold on
    % for k=1:length(i)
    %     plot(i(k),j(k), 'o', 'markerFaceColor', 'k', 'markerSize', ms, 'markerEdgeColor', 'none')
    % end
    %         set(gca, 'box', 'off')
    %         set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
    
    if eegTimestepMode == 1
        figure
        imagesc(meanVarSlopes)
        ax = gca;
        ax.FontSize = 10;
        
        [i,j]=find(varP<.05);
        ms=3;
        hold on
        for k=1:length(i)
            plot(j(k),i(k), 'o', 'markerFaceColor', 'k', 'markerSize', ms, 'markerEdgeColor', 'none')
        end
        ylabel("Pupil Data")
        yticks(1:size(varP,1))
        xlabel("Parameter")
        tickLabels = [];
        for i = 1:length(borderTimes)-1
            tickLabels = string([tickLabels;[num2str(borderTimes(i)-timeBeforeEEG),'ms to ',num2str(borderTimes(i+1)-timeBeforeEEG),'ms']]);
        end
        yticklabels(tickLabels')
        xticklabels(["BiasCP","BiasOB","BiasSum","LRCP","LROB","LRDiff"])
        colorbar
        
        %
        figure 
        imagesc(log10(varP))
        ax = gca;
        ax.FontSize = 10;
        [R, C] = ndgrid(1:size(varP,1), 1:size(varP,2));
        R = R(:);
        C = C(:) - 1/4;
        
        clear mask
        vals = round(varP,4);
        mask = vals <= 0;
        text(C(mask), R(mask), string(vals(mask)), 'color', 'w')
        text(C(~mask), R(~mask), string(vals(~mask)), 'color', 'k')
        ylabel("Pupil Data")
        yticks(1:1:size(varP,1))
        xlabel("Parameter")
        yticklabels(tickLabels')
        xticklabels(["BiasCP","BiasOB","BiasSum","LRCP","LROB","LRDiff"])
        colorbar
        
        
        %find correlation between eeg clusters and pupil dilations
        stackedEEGCorrs = corr(reshape(permute(meanTrialEffectEEG,[1,3,2]),[],size(meanTrialEffectEEG,2),1),'rows','complete');
        stackedEEGEyeCorrs = corr(reshape(permute(meanTrialEffectEEGEyeSep,[1,3,2]),[],size(meanTrialEffectEEGEyeSep,2),1),'rows','complete');
        allStackedCorrs = [[stackedEEGCorrs;stackedEEGEyeCorrs(end,1:end-1)],stackedEEGEyeCorrs(:,end)];
        
        %image plot for correlation between eeg clusters and pupil dilation
        figure
        imagesc(allStackedCorrs)
        yticks(1:length(allStackedCorrs))
        yticklabels([clusterLabels])
        xticks(1:length(allStackedCorrs))
        xticklabels([clusterLabels])
        
        colorbar
        colormap parula
        caxis([-1,1]);
        [R, C] = ndgrid(1:size(allStackedCorrs,1), 1:size(allStackedCorrs,1));
        R = R(:); C = C(:) - 1/10;
        
        clear mask
        vals = round(allStackedCorrs,2);
        mask = vals <= 0;
        text(C(mask), R(mask), string(vals(mask)), 'color', 'w',"FontName","Arial","FontWeight","bold","FontSize",14)
        text(C(~mask), R(~mask), string(vals(~mask)), 'color', 'k',"FontName","Arial","FontWeight","bold","FontSize",14)
                set(gca, 'box', 'off')
                set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
    end
    
    %%
    % make figure 6 if not doing stp residuals, otherwise make figure S2
    if doSTPResiduals == 1
        %figure s2
        %scattered subject bias CP vs OB
        figure("Position",[100,100,900,900])
        subplot(2,8,1:4)
            scatter(biasOBslopeEEGEye(:,2),biasCPslopeEEGEye(:,2),25,[.5,.5,.5],'filled','MarkerEdgeColor','k')
            hold on
            xlabel("OB Bias vs Arousal Slope")
            ylabel("CP Bias vs Arousal Slope")
            set(gca, 'box', 'off')
            set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
            yline(0)
            xline(0)
            xticks(0)
            yticks(0)
            maxBias = max([biasOBslopeEEGEye(:,2);biasCPslopeEEGEye(:,2)])*1.05;
            ylim([-maxBias,maxBias])
            xlim([-maxBias,maxBias])
            plot([-1,1],[1,-1],'--','Color',[0.5,0.5,0.5])
        
        %scattered subject learning cp vs ob
        subplot(2,8,5:8)
            scatter(LROBslopeEEGEye(:,2),LRCPslopeEEGEye(:,2),25,[.5,.5,.5],'filled','MarkerEdgeColor','k')
            hold on
            xlabel("OB LR vs Arousal Slope")
            ylabel("CP LR vs Arousal Slope")
            set(gca, 'box', 'off')
            set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
            yline(0)
            xline(0)
            xticks(0)
            yticks(0)
            maxLR = max([LROBslopeEEGEye(:,2);LRCPslopeEEGEye(:,2)])*1.05;
            ylim([-maxLR,maxLR])
            xlim([-maxLR,maxLR])
            plot([-1,1],[-1,1],'--','Color',[0.5,0.5,0.5])
        
        %hypothesis 5 main plot with residuals
        subplot(2,8,[9:12])
            scatter(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2),25,[251 200 143]./255,'filled','MarkerEdgeColor',[246 146 30]./255)
            hold on
            scatter(biasAllslopeEEGEye(:,2),LROBslopeEEGEye(:,2),25,[128 214 247]./255,'filled','MarkerEdgeColor',[0 173 238]./255)
            mCP = polyfit(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2),1);
            plot((linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2)),mCP(1)*(linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2))+mCP(2),'Color',[246 146 30]./255);
            mOB = polyfit(biasAllslopeEEGEye(:,2),LROBslopeEEGEye(:,2),1);
            mCPOB = polyfit(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2)-LROBslopeEEGEye(:,2),1);
            plot((linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2)),mOB(1)*(linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2))+mOB(2),'Color',[0 173 238]./255);
            hold on
            axLim = 1.05*max(abs([biasAllslopeEEGEye(:,2);LRCPslopeEEGEye(:,2);LROBslopeEEGEye(:,2)]));
            ylim([-axLim,axLim])
            xlim([-axLim,axLim])
            xline(0)
            yline(0)
            yticks(0)
            xticks(0)
            xlabel("Bias vs Arousal Slope")
            ylabel("Learning vs Arousal Slope")
            legend('Changepoint','Oddball')
            set(gca, 'box', 'off')
            set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
        
    
        
        %hypotheses 5 regression result plot
        subplot(2,8,14:16)    
        [slopeBs,slopeBsInt] = regress(allLRSlopes,zX);
        %testing hypothesis 5 by regressing learning slopes against bias slopes, condition, etc
            bar(slopeBs(1:4))
            hold on
            errorbar(slopeBs(1:4),slopeBsInt(1:4,2)-slopeBs(1:4),'.','Color','k');
            xticklabels(["Intercept","Bias Slope","Cond","Bias Slope*Cond"])
            ylabel("Regression Coefficient")
            yticks(0)
            set(gca, 'box', 'off')
            set(gca,"FontName","Arial","FontWeight","bold","FontSize",12)
    
        %save figure S2
        if doSTPResiduals == 0
            fig = gcf;
            figName = append("Figure_S2_",figTime,'.eps');
            figLoc = append(figDir,figName);
            exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
            figName = append("Figure_S2_",figTime,'.png');
            figLoc = append(figDir,figName);
            saveas(fig,figLoc)
        else
            fig = gcf;
            figName = append("Figure_S2_Residual_",figTime,'.eps');
            figLoc = append(figDir,figName);
            exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
            figName = append("Figure_S2_Residual",figTime,'.png');
            figLoc = append(figDir,figName);
            saveas(fig,figLoc)
        end
    else
        % make figure 6 with a panel open for illustrator creation
        figure('Position',[0 250 1800 600])
        %main hypothesis 5 panel
        subplot(1,5,3:4)
        scatter(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2),25,[251 200 143]./255,'filled','MarkerEdgeColor',[246 146 30]./255)
        hold on
        scatter(biasAllslopeEEGEye(:,2),LROBslopeEEGEye(:,2),25,[128 214 247]./255,'filled','MarkerEdgeColor',[0 173 238]./255)
        mCP = polyfit(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2),1);
        plot((linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2)),mCP(1)*(linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2))+mCP(2),'Color',[246 146 30]./255);
        mOB = polyfit(biasAllslopeEEGEye(:,2),LROBslopeEEGEye(:,2),1);
        mCPOB = polyfit(biasAllslopeEEGEye(:,2),LRCPslopeEEGEye(:,2)-LROBslopeEEGEye(:,2),1);
        plot((linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2)),mOB(1)*(linspace(min(biasAllslopeEEGEye(:,2)),max(biasAllslopeEEGEye(:,2)),2))+mOB(2),'Color',[0 173 238]./255);
        hold on
        axLim = 1.05*max(abs([biasAllslopeEEGEye(:,2);LRCPslopeEEGEye(:,2);LROBslopeEEGEye(:,2)]));
        ylim([-axLim,axLim])
        xlim([-axLim,axLim])
        xline(0)
        yline(0)
        yticks(0)
        xticks(0)
        xlabel("Bias vs Arousal Slope")
        ylabel("Learning vs Arousal Slope")
        legend('Changepoint','Oddball')
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
    
        %set up regression for testing hypothesis 5 
        subContext = [ones(eegEyeNumSubs,1);-1*ones(eegEyeNumSubs,1)];
        allLRSlopes = [LRCPslopeEEGEye(:,2);LROBslopeEEGEye(:,2)];
        allBiasSlopes = [biasAllslopeEEGEye(:,2);biasAllslopeEEGEye(:,2)];
        
        X = [ones(length(subContext),1),allBiasSlopes,subContext,allBiasSlopes.*subContext];
        useDummy = 0;
        if useDummy == 1
            dummyMat = [eye(eegEyeNumSubs);eye(eegEyeNumSubs)];
            zX = [X(:,1),zscore(X(:,2:end)),dummyMat];
        else
            zX = [X(:,1),zscore(X(:,2:end))];
        end
        %regressing Bias, condition, bias*condition against learning
        [slopeBs,slopeBsInt] = regress(allLRSlopes,zX);
        
        %hypothesis 5 regression plot
        subplot(1,5,5)
            bar(slopeBs(1:4))
            hold on
            errorbar(slopeBs(1:4),slopeBsInt(1:4,2)-slopeBs(1:4),'.','Color','k');
            xticklabels(["Intercept","Bias Slope","Condition","Bias Slope*Condition"])
            ylabel("Regression Coefficient")
            yticks(0)
            set(gca, 'box', 'off')
            set(gca,"FontName","Arial","FontWeight","bold","FontSize",14)
    
        %save figure 6 (hypothesis 5)
        if doSTPResiduals == 0
            fig = gcf;
            figName = append("Figure_6_",figTime,'.eps');
            figLoc = append(figDir,figName);
            exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
            figName = append("Figure_6_",figTime,'.png');
            figLoc = append(figDir,figName);
            saveas(fig,figLoc)
        else
            fig = gcf;
            figName = append("Figure_6_Residual_",figTime,'.eps');
            figLoc = append(figDir,figName);
            exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
            figName = append("Figure_6_Residual",figTime,'.png');
            figLoc = append(figDir,figName);
            saveas(fig,figLoc)
        end
    end
end
%% MAKE bias and learning FIGURES FOR ILLUSTRATOR
%set panels for learning/bias figures
panel1 = [1:4,15:18,29:32];
panel2 = [6:9,20:23,34:37];
panel3 = [11:14,25:28,39:42];
panel4 = [57:59,71:73,85:87];
panel5 = [61:63,75:77,89:91];
panel6 = [65:67,79:81,93:95];
panel7 = [69,70,83,84,97,98];

% Bias Figure
figure("Position",[100,100,1600,850])
set(gcf,'renderer','Painters')
subplot(7,14,panel1)
%raw estimation error quantiled
    scatter(quantileAllCPxesBias,meanQuantilesCPAllErrBias,30,'o',"MarkerEdgeColor",cpColor,"MarkerFaceColor",cpLColor)
    hold on
    scatter(quantileAllOBxesBias,meanQuantilesOBAllErrBias,30,'o',"MarkerEdgeColor",obColor,"MarkerFaceColor",obLColor)
    plot([-3,3],[-3,3],'--','Color',[0.5,0.5,0.5])
    hold off
    ylabel("Perceptual Error")
    xlabel("Objective Prediction Error")
    ylim([-.5,.5])
    yline(0)
    xline(0)
    xlim([-3,3])
    xticks([-3,3])
    legend("Changepoint","Oddball","Location","northwest")
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

subplot(7,14,panel2)
%raw estimation error quantiled for good subjects
    scatter(quantileAllCPxesBias,meanQuantilesCPAllErrBias,30,'o',"MarkerEdgeColor",cpColor,"MarkerFaceColor",cpLColor)
    hold on
    plot([-3,3],[-3,3],'--','Color',[0.5,0.5,0.5])
    scatter(quantileAllOBxesBias,meanQuantilesOBAllErrBias,30,'o',"MarkerEdgeColor",obColor,"MarkerFaceColor",obLColor)
    hold off
    ylabel("Perceptual Error")
    xlabel("Objective Prediction Error")
    ylim([-.5,.5])
    yline(0)
    xline(0)
    xlim([-3,3])
    xticks([-3,3])
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

subplot(7,14,panel3)
%bias behavioral regression
    semBias=std(paramsCircBiasAll,1)./sqrt(size(paramsCircBiasAll,2));
    x=[0.01:.01:0.01*length(behaveSubs)];
    stdBias=std(paramsCircBiasAll);
    normalizedCoef=[];
    l = 0;
    cAll =[];
    for c=3:6    
        if c==9
            normalizedCoef=paramsCircBiasAll(:,end)./stdBias(end);
            normalizedCoefBias(:,c-1)=normalizedCoef;
            meanToPlot=mean(normalizedCoef);
            semToPlot=std(normalizedCoef)./sqrt(length(normalizedCoef));   
            scatter(x+c,normalizedCoef,30,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:) )
            hold on
            errorbar(mean(x+c),meanToPlot,1.96*semToPlot,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', .75, 'markerSize', 8 )
            scatter(mean(x+c),mean(paramsCircUpdateAllModel(:,c+1))./stdUpdate(c),100,'xk','LineWidth',3);
            xtickVal(c)=mean(x+c);
            l = l+1;
            cAll(l) = c;
        else
            normalizedCoef=paramsCircBiasAll(:,c)./stdBias(c);
            normalizedCoefBias(:,c-1)=normalizedCoef;
            meanToPlot=mean(normalizedCoef);
            semToPlot=std(normalizedCoef)./sqrt(length(normalizedCoef));
            scatter(x+c,normalizedCoef,30,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:))
            hold on
            errorbar(mean(x+c),meanToPlot,1.96*semToPlot,'^k','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 8)
            scatter(mean(x+c),mean(paramsCircBiasAllModel(:,c))./stdBias(c),100,'xk','LineWidth',2);
            xtickVal(c)=mean(x+c);
            l = l+1;
            cAll(l) = c;
        end
        %plot(x+c,normalizedCoef,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 10 )
    end
    %xlim([0.1 0.6])
    %ylim([-0.1 0.1])
    set(gca, 'box', 'off')
    ylabel('Normalized Coefficient')
    xticks(xtickVal(cAll))
    xticklabels({'PE','PE*STP*CP/OB','PE*STP','PE*entropy','PE x CP/OB','PE x uniform','Gaze attention'})
    xtickangle(45)
    yline(0);
    xlim([min(cAll)-0.5 c+1])
    yticks(0)
    %ylim([-4 4])
    hold off
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

if ~isempty(gSig)
    subplot(7,14,panel4)
    %eeg bias quantiles
        semCP = std(meanQuantCPBiasEEG(:,:))/sqrt(length(EEGSubs));
        semOB = std(meanQuantOBBiasEEG(:,:))/sqrt(length(EEGSubs));
        errorbar(mean(meanQuantCPBiasEEG(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBBiasEEG(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.2,0.35]);
        yticks([0.2,0.35]) 
        xticks(1:numQuant)
        xlabel("EEG Quantile")
        ylabel("Mean Bias")
        fitCPBias = polyfit(1:numQuant,mean(meanQuantCPBiasEEG), 1);
        fitOBBias = polyfit(1:numQuant,mean(meanQuantOBBiasEEG), 1);
        x = 1:numQuant;
        yCPBias = polyval(fitCPBias , x);
        yOBBias = polyval(fitOBBias , x);
        plot(x,yCPBias,'Color',cpColor)
        plot(x,yOBBias,'Color',obColor)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel5)
    %pupil bias quantiles
        semCP = std(meanQuantCPBiasEye(:,:))/sqrt(length(eyeSubs));
        semOB = std(meanQuantOBBiasEye(:,:))/sqrt(length(eyeSubs));
        errorbar(mean(meanQuantCPBiasEye(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBBiasEye(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.2,0.35]);
        yticks([0.2,0.35]) 
        xticks(1:numQuant)
        xlabel("Pupil Quantile")
        ylabel("Mean Bias")
        fitCPBias = polyfit(1:numQuant,mean(meanQuantCPBiasEye), 1);
        fitOBBias = polyfit(1:numQuant,mean(meanQuantOBBiasEye), 1);
        x = 1:numQuant;
        yCPBias = polyval(fitCPBias , x);
        yOBBias = polyval(fitOBBias , x);
        plot(x,yCPBias,'Color',cpColor)
        plot(x,yOBBias,'Color',obColor)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel6)
    %eeg+pupil bias quantiles
        semCP = std(meanQuantCPBiasEEGEye(:,:))/sqrt(eegEyeNumSubs);
        semOB = std(meanQuantOBBiasEEGEye(:,:))/sqrt(eegEyeNumSubs);
        errorbar(mean(meanQuantCPBiasEEGEye(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBBiasEEGEye(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.2,0.35]);
        yticks([0.2,0.35]) 
        xticks(1:numQuant)
        xlabel("EEG + Pupil Quantile")
        ylabel("Mean Bias")
        fitCPBias = polyfit(1:numQuant,mean(meanQuantCPBiasEEGEye), 1);
        fitOBBias = polyfit(1:numQuant,mean(meanQuantOBBiasEEGEye), 1);
        x = 1:numQuant;
        yCPBias = polyval(fitCPBias , x);
        yOBBias = polyval(fitOBBias , x);
        plot(x,yCPBias,'Color',cpColor)
        plot(x,yOBBias,'Color',obColor)
        legend("","","Changepoint","Oddball")
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel7)
    %separate clusters bias slopes
     bar([mean(varSlopes(:,3,:),3);mean(biasCPslopeEye(:,2)'+biasOBslopeEye(:,2)')])
        hold on
        EEGSlopes = squeeze(varSlopes(:,3,:))';
        eyeSlopes = biasCPslopeEye(:,2)'+biasOBslopeEye(:,2)';
        sem = [std(EEGSlopes)./sqrt(length(EEGSubs)),std(eyeSlopes)./sqrt(length(eyeSubs))];
        errorbar(1:size(varSlopes,1)+1,[mean(EEGSlopes),mean(eyeSlopes')],sem,"Color",'k','LineStyle','none');
        xticklabels(clusterLabels)
        %FIX XTICK LABELS
        ylabel("Average Bias Slope")
        yticks(0)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
end

    %save bias figure (figure 3)
if doSTPResiduals == 0
    fig = gcf;
    figName = append("Figure_3_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_3_",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
else
    fig = gcf;
    figName = append("Figure_3_Residual_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_3_Residual",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
end

%%
% Learning Rate Figure
figure("Position",[100,100,1600,850])
set(gcf,'renderer','Painters')
subplot(7,14,panel1)
%raw prediction update quantiled
    scatter(quantileAllCPxesLR,meanQuantilesCPAllErrLR,30,'o',"MarkerEdgeColor",cpColor,"MarkerFaceColor",cpLColor)
    hold on
    scatter(quantileAllOBxesLR,meanQuantilesOBAllErrLR,30,'o',"MarkerEdgeColor",obColor,"MarkerFaceColor",obLColor)
    plot([-3,3],[-3,3],'--','Color',[0.5,0.5,0.5])
    hold off
    ylabel("Prediction Update")
    xlabel("Subjective Prediction Error")
    ylim([-3,3])
    yticks([-3,3])
    yline(0)
    xline(0)
    xlim([-3,3])
    xticks([-3,3])
    legend("Changepoint","Oddball","Location","northwest")
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

subplot(7,14,panel2)
%raw prediction update quantiled for good subjects
    scatter(quantileGoodCPxesLR,meanQuantilesCPGoodUpLR,30,'o',"MarkerEdgeColor",cpColor,"MarkerFaceColor",cpLColor)
    hold on
    scatter(quantileGoodOBxesLR,meanQuantilesOBGoodUpLR,30,'o',"MarkerEdgeColor",obColor,"MarkerFaceColor",obLColor)
    plot([-3,3],[-3,3],'--','Color',[0.5,0.5,0.5])
    hold off
    ylabel("Prediction Update")
    xlabel("Subjective Prediction Error")
    ylim([-3,3])
    yticks([-3,3])
    yline(0)
    xline(0)
    xlim([-3,3])
    xticks([-3,3])
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

subplot(7,14,panel3)
%learning behavioral regression
    semUpdate=std(paramsCircUpdateAll,1)./sqrt(size(paramsCircUpdateAll,2));
    x=[0.01:.01:0.01*length(behaveSubs)];
    %plotting 5 coefficients
    stdUpdate=std(paramsCircUpdateAll);
    normalizedCoef=[];
    l = 0;
    cAll =[];
    for c=3:6
        normalizedCoef=paramsCircUpdateAll(:,c)./stdUpdate(c);
        normalizedCoefUpdate(:,c-1)=normalizedCoef;
        meanToPlot=mean(normalizedCoef);
        semToPlot=std(normalizedCoef)./sqrt(length(normalizedCoef));
        %keyboard
        scatter(x+c,normalizedCoef,30,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:))
        hold on
        errorbar(mean(x+c),meanToPlot,1.96*semToPlot,'^k','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 8 )
        scatter(mean(x+c),mean(paramsCircUpdateAllModel(:,c))./stdUpdate(c),100,'xk','LineWidth',2);
        xtickVal(c)=mean(x+c);
        l = l+1;
        cAll(l) = c;
    end
    %xlim([0.1 0.6])
    %ylim([-1 1])
    % title('Learning Regression Circ','color',cbColors(4,:))
    ylabel('Normalized Coefficient')
    xticks(xtickVal(cAll))
    xticklabels({'PE','PE*STP*Cond','PE*STP','PE*entropy','PE x condition','PE x uniform'})
    xtickangle(45)
    yline(0);
    %ylim([-4,8])
    xlim([min(cAll)-0.5 c+1])
    yticks(0)
    %ylim([-4.5 4.5])
    hold off
    set(gca, 'box', 'off')
    set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)

if ~isempty(gSig)
    subplot(7,14,panel4)
    %eeg learning quantiles
        semCP = std(meanQuantCPLRsEEG(:,:))/sqrt(length(EEGSubs));
        semOB = std(meanQuantOBLRsEEG(:,:))/sqrt(length(EEGSubs));
        errorbar(mean(meanQuantCPLRsEEG(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBLRsEEG(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.45,0.67]);
        yticks([0.45,0.67]) 
        xticks(1:numQuant)
        xlabel("EEG Quantile")
        ylabel("Mean LR")
        fitCPLR = polyfit(1:numQuant,mean(meanQuantCPLRsEEG), 1);
        fitOBLR = polyfit(1:numQuant,mean(meanQuantOBLRsEEG), 1);
        x = 1:numQuant;
        yCPLR = polyval(fitCPLR , x);
        yOBLR = polyval(fitOBLR , x);
        plot(x,yCPLR,'Color',cpColor)
        plot(x,yOBLR,'Color',obColor)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel5)
    %pupil learning quantiles
        semCP = std(meanQuantCPLRsEye(:,:))/sqrt(length(eyeSubs));
        semOB = std(meanQuantOBLRsEye(:,:))/sqrt(length(eyeSubs));
        errorbar(mean(meanQuantCPLRsEye(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBLRsEye(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.45,0.67]);
        yticks([0.45,0.67]) 
        xticks(1:numQuant)
        xlabel("Pupil Quantile")
        ylabel("Mean LR")
        fitCPLR = polyfit(1:numQuant,mean(meanQuantCPLRsEye), 1);
        fitOBLR = polyfit(1:numQuant,mean(meanQuantOBLRsEye), 1);
        x = 1:numQuant;
        yCPLR = polyval(fitCPLR , x);
        yOBLR = polyval(fitOBLR , x);
        plot(x,yCPLR,'Color',cpColor)
        plot(x,yOBLR,'Color',obColor)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel6)
    %eeg+pupil learning quantiles
        semCP = std(meanQuantCPLRsEEGEye(:,:))/sqrt(eegEyeNumSubs);
        semOB = std(meanQuantOBLRsEEGEye(:,:))/sqrt(eegEyeNumSubs);
        errorbar(mean(meanQuantCPLRsEEGEye(:,:)),semCP,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",cpColor,"CapSize",0,"Color",cpColor,"MarkerEdgeColor",cpColor)
        hold on
        errorbar(mean(meanQuantOBLRsEEGEye(:,:)),semOB,'LineStyle','none','LineWidth',1,'Marker','o',"MarkerFaceColor",obColor,"CapSize",0,"Color",obColor,"MarkerEdgeColor",obColor)
        xlim([0,numQuant+1])
        ylim([0.45,0.67]);
        yticks([0.45,0.67]) 
        xticks(1:numQuant)
        xlabel("EEG + Pupil Quantile")
        ylabel("Mean LR")
        fitCPLR = polyfit(1:numQuant,mean(meanQuantCPLRsEEGEye), 1);
        fitOBLR = polyfit(1:numQuant,mean(meanQuantOBLRsEEGEye), 1);
        x = 1:numQuant;
        yCPLR = polyval(fitCPLR , x);
        yOBLR = polyval(fitOBLR , x);
        plot(x,yCPLR,'Color',cpColor)
        plot(x,yOBLR,'Color',obColor)
        legend("","","Changepoint","Oddball","Location","east")
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)
    
    subplot(7,14,panel7)
    %separate clusters learning slopes
    bar([mean(varSlopes(:,6,:),3);mean(LRCPslopeEye(:,2)'-LROBslopeEye(:,2)')])
        hold on
        EEGSlopes = squeeze(varSlopes(:,6,:))';
        eyeSlopes = LRCPslopeEye(:,2)'-LROBslopeEye(:,2)';
        sem = [std(EEGSlopes)./sqrt(length(EEGSubs)),std(eyeSlopes)./sqrt(length(eyeSubs))];
        errorbar(1:size(varSlopes,1)+1,[mean(EEGSlopes),mean(eyeSlopes')],sem,"Color",'k','LineStyle','none');
        xticklabels(clusterLabels)
        %FIX XTICK LABELS
        ylabel("Average LR Slope")
        yticks(0)
        set(gca, 'box', 'off')
        set(gca,"FontName","Arial","FontWeight","bold","FontSize",18)        
end

%save lr figure (figure 4)
if doSTPResiduals == 0
    fig = gcf;
    figName = append("Figure_4_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_4_",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
else
    fig = gcf;
    figName = append("Figure_4_Residual_",figTime,'.eps');
    figLoc = append(figDir,figName);
    exportgraphics(fig,figLoc,'BackgroundColor','none','ContentType','vector')
    figName = append("Figure_4_Residual",figTime,'.png');
    figLoc = append(figDir,figName);
    saveas(fig,figLoc)
end
