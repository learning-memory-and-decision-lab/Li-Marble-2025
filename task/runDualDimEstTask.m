% run VSTM task:
% shell script to run VSTM task

clear

%set up paths for task folder and colorTools and data paths
whichComp=1;

if whichComp==1
    basePath='~/Documents/GitHub/Li-Marble_2025/task/';
    colorPath='~/Documents/GitHub/Li-Marble_2025/task/colorTools/';
    eyePath='EyeDataLocal/';
    behavePath='/behaviorDataLocal/';
    
else

end

cd(basePath)
addpath(genpath('./'))
addpath(genpath(colorPath))
addpath(genpath(eyePath))
addpath(genpath(behavePath))



calibrationFile=[colorPath,'labMonitorData_2022-03-16_426C.mat'];


load(calibrationFile)

capFit=input('Are you doing cap fitting? Yes(1)/No(0):');
if capFit
    doEEG=0;
else
    doEEG=input('Are you collecting EEG data? Yes(1)/No(0):');
end



%prefs.taskType = 0;          %: Pick one of several task types, previously 3
%%
if capFit
    prefs.numBlocks=1;
else
    prefs.numBlocks=4;
end

%prefs.nTrials = 4;
prefs.nItems = 2;
prefs.stimulusDuration = .2;%.5
prefs.retentionInterval = 1.8;
prefs.squareSize = 40;       % in pixels
prefs.radius = 160;
prefs.fixationSize = 3;
prefs.colorWheelRadius = 350;
prefs.propFixedSpace = 0;
prefs.fixedSpacings  = [360./prefs.nItems: 360./prefs.nItems:359];
prefs.fixedOSpacings = [360./prefs.nItems: 360./prefs.nItems:359];
prefs.useCieLab=true;
prefs.monitorData=monitorData;
prefs.rectSize1= 120;%80; %120;         %: length of stimulus rectangle
prefs.rectSize2= 120; %80; %120;         %: width of stimulus rectangle
prefs.debug = 0;         %: Are you debugging the task? true = small screen
prefs.preEstPause=1;
%prefs.doFeedback=0;          %  Display + or - feedback after each trial?
prefs.fdbkProp=0;            %  Prop b feedback.
prefs.fdbk_oThresh=pi./6;    %  orientation error threshold
prefs.fdbk_cThresh=pi./6;    %  color threshold
prefs.fdbk_duration= 1;
prefs.ITI= .5;
prefs.pdw_prop=0;          %  Proportion of trials that will have a post decision wager
prefs.pdw_vals=[0 2];        %  Possible bet values [array of two integers]
prefs.openWindow=true;
prefs.closeWindow=false;
prefs.initRTthresh=5;
prefs.responseRTthresh=30;
prefs.skipSyncTests=true;
prefs.jittered = false;
%prefs.doPredict =true;t+
prefs.showMask=0;
prefs.testGray = round(LabTo_calRGB([prefs.monitorData.lVal, 0,0], prefs.monitorData));
prefs.pointsPerDollar=30;
prefs.doPredInstructions=0;  % whether task will run with instructions during CP/OB blocks
if doEEG
    prefs.doEEG=1;
else
    prefs.doEEG=0;
end
prefs.doET=1;

if whichComp==3
    prefs.fontSize=40;
    prefs.imageSizeX=1500;
    prefs.imageSizeY=1000;
else
    prefs.fontSize=20;
    prefs.imageSizeX=1200;
    prefs.imageSizeY=800;
end


SN=input('Enter subject number:', 's');

ETname=input('Enter name for Eyetracking data file:','s');

edfFile=[eyePath ETname '.edf'];
edfFileFit=[eyePath ETname '0.edf'];

%1 is CP first, 2 is OB first; counter-balance across subjects
Cond=input('Enter condition for subject:','s');

CondN=str2double(Cond);

TN='vwm_multiStim_color_';



%%%%%%%%%%%%%%%%% Initialize EEG Calls %%%%%%%%%%%%%%%%%%%%%%%%%%
if prefs.doEEG || prefs.doET
    config_parallel
    ioObject = io64;
    LPT1address = hex2dec('3FB8');
    %'E050'); %standard location of LPT1 port
    status = io64(ioObject);

    % Create structure of trigger numbers:
    triggerNum.start           =15;
    triggerNum.end             =14;
    
    triggerNum.practice           =16;
    %triggerNum.block           =17;
    
    triggerNum.block           =1;
    triggerNum.instructionsOn  =2;
    triggerNum.instructionsOff =3;
    triggerNum.stimOn          =4;
    triggerNum.stimOff         =5;
    triggerNum.cue1on          =6;
    triggerNum.responseMade1   =7;  %16;  %240;
    triggerNum.cue2on          =8;
    triggerNum.responseMade2   =9;  %32;  %224;
    triggerNum.predCue1on      =10;
    triggerNum.predResp1       =11;    %48;  %208;
    triggerNum.predCue2on      =12;
    triggerNum.predResp2       =13;   %64;  %192;
    
    % put structure in prefs:
    
    prefs.triggerNum=triggerNum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Test trigger
% % 
% while true
%     io64(ioObject,LPT1address,prefs.triggerNum.block);
%     pause(0.05)
%     disp('got here')
%     io64(ioObject,LPT1address,0);
% end

%io64(ioObject,LPT1address,0)

%%
% if prefs.doEEG
%          io64(ioObject,LPT1address,prefs.triggerNum.start);
%          WaitSecs(0.01)
%          io64(ioObject,LPT1address,0);
% end

%EyelinkSetup(1)

rng('default')
rng shuffle %try this later
% a = clock;
% randomSeed = rng(a(6));

clear sumScoreAll blockScore
for bb=1:prefs.numBlocks
     % send trigger indicate the start of every block
%      if prefs.doEEG
%          io64(ioObject,LPT1address,prefs.triggerNum.block);
%          WaitSecs(0.01)
%          io64(ioObject,LPT1address,0);
%      end
%      
%      if prefs.doET
%          myMessage = 'block start';
%          Eyelink('SendMessage', myMessage);
%      end
     
    if bb>1
        prefs.openWindow=false;
        
    end
    
    if bb==prefs.numBlocks
            prefs.closeWindow=true;
    end
    
    fn=[behavePath TN SN num2str(bb) '.mat'];
    
%     if CondN==2 && bb==3
%         bb=4;
%     elseif CondN==2 && bb==4
%         bb=3;
%     end

    if bb==1 || bb==2
        sumScore=0;
    else
        sumScore=sumScoreAll(bb-1);
    end
    
    if capFit
        SNfit=[SN,'0'];
        prefs.capFit=1;
        [data, prefs]=colorOrient_estTask_colorPupil(prefs, fn, eyePath,SN,bb,CondN);
        %pupilCapFitting(prefs,SNfit)
    else
        if doEEG
            [data, prefs]=colorOrient_estTask_test(prefs, fn, eyePath,SN,bb,CondN,sumScore);
            blockScore(bb)=data.totScore;
            sumScoreAll(bb)=data.sumScore;
        else
            [data, prefs]=colorOrient_estTask_noEEG(prefs, fn, eyePath,SN,bb,CondN,sumScore);
            blockScore(bb)=data.totScore;
            sumScoreAll(bb)=data.sumScore;
            
        end
    end 
    
end

%disp(num2str(randomSeed.Seed));

if ~capFit
    disp(sprintf('You Earned %s Points for %s bonus dollars!!!', num2str(sumScoreAll(end)), num2str(ceil(sumScoreAll(end)./prefs.pointsPerDollar))))
end


