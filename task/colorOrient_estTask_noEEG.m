function [data, prefs]=colorOrient_estTask_Harry_noEEG(prefs, saveFileName,edfFile,subStr, block,condition,sumScore)

cd(edfFile)

rng(str2num(datestr(now,"yyddhhmmss")),'twister')
if prefs.doET
    ioObject = io64;
    LPT1address = hex2dec('3FB8');
    status = io64(ioObject);
end



if ~isfield(prefs, 'mouseScale')
    prefs.mouseScale=1;
end

if ~isfield(prefs, 'negPointScale');
    prefs.negPointScale=-1;
end

% Block specific paramters
% pilot: 1 hour w/ block1 5 trials; block2 20 trials; block3&4 70 trials each
%actual task: 20,40,120,120
if block==1
    prefs.taskType = 0;
    prefs.doPredict=false;
    prefs.doFeedback=1;
    prefs.nTrials = 2;
elseif block==2
    prefs.taskType = 0;
    prefs.doPredict=false;
    prefs.doFeedback=0;
    prefs.nTrials = 2;
elseif block==3 && condition==1
    prefs.taskType = 10;
    prefs.doPredict=true;
    prefs.doFeedback=0;
    prefs.nTrials = 2;
elseif block==3 && condition==2
    prefs.taskType = 11;
    prefs.doPredict=true;
    prefs.doFeedback=0;
    prefs.nTrials = 2;
elseif block==4 && condition==1
    prefs.taskType = 11;
    prefs.doPredict=true;
    prefs.doFeedback=0;
    prefs.nTrials = 2;
elseif block==4 && condition==2
    prefs.taskType = 10;
    prefs.doPredict=true;
    prefs.doFeedback=0;
    prefs.nTrials = 2;
end




%%

try

    if prefs.openWindow
        prepareEnvironment(prefs);
        % light gray starts here:
        window = openWindow(prefs);
        if prefs.doET
            [el] = EyelinkSetup(1,window.onScreen) %initialize eyelink setup (gray screen)
            Eyelink('Openfile', subStr); % create a file for eye data
            Eyelink('StartRecording'); % start recording from eyelink
            if prefs.doET
                io64(ioObject,LPT1address,prefs.triggerNum.start); %send a trigger through parallel port that signals start of task
                WaitSecs(0.01) %wait 0.01 seconds
                io64(ioObject,LPT1address,0); %stop trigger
            end
            % record a few samples before we actually start displaying
            WaitSecs(0.1);
        end
        %keyboard
    else
        window = prefs.window;
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.block);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
    end
    
    
    if prefs.doFeedback
        negFdbk=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Feedback Images/Gray - Neg Feedback.jpg');
        posFdbk=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Feedback Images/Gray - Pos Feedback.jpg');
        
        % first texture is negative feedback, second is positive.
        fdbk_textures{1} = Screen('MakeTexture', window.onScreen, negFdbk);
        fdbk_textures{2} = Screen('MakeTexture', window.onScreen, posFdbk);
    end
    
    if prefs.doPredict
        
        figCP=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/CP.jpg');
        figOB=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/OB.jpg');
        
        % first is change point stats figure
        instr_textures{1} = Screen('MakeTexture', window.onScreen, figCP);
        instr_textures{2} = Screen('MakeTexture', window.onScreen, figOB);
    end
    
    %if first block, then show elaborated version of instruction to whole task
    if block==1
        Screen('TextSize',window.onScreen,prefs.fontSize);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        Screen('FillRect', window.onScreen, window.bcolor);
        
        % send trigger indicate participant currently reading instruction
        %         io64(ioObject,LPT1address,prefs.triggerNum.instructions);
        %         WaitSecs(0.01)
        %         io64(ioObject,LPT1address,0);
        
        %elaborate whole task instruction
        %         instructionText1 = 'Welcome to our study! \n\n\n In this task, you will be asked to reproduce colors \n that briefly flash on the computer screen. \n\n\n Please keep your gaze at the center of the screen throughout the task. \n\n\n Left click to continue.';
        %         DrawFormattedText(window.onScreen,instructionText1,'center','center')
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.instructionsOn);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
        fig1=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide1.jpg');
        Screen('FillRect', window.onScreen, window.bcolor);
        Slide1=Screen('MakeTexture', window.onScreen, fig1);
        Screen('DrawTexture', window.onScreen, Slide1, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        fig2=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide2.jpg');
        Slide2=Screen('MakeTexture', window.onScreen, fig2);
        Screen('DrawTexture', window.onScreen, Slide2, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
%         fig3=imread('./Instruction Images/Slide3.jpg');
%         Slide3=Screen('MakeTexture', window.onScreen, fig3);
%         Screen('DrawTexture', window.onScreen, Slide3, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
%         Screen(window.onScreen, 'Flip');
%         pause(1)
%         [x,y,buttons]=GetMouse;
%         while ~any(buttons)
%             [x,y,buttons]=GetMouse
%             WaitSecs(.01)
%         end
        
        fig4=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide4.jpg');
        Slide4=Screen('MakeTexture', window.onScreen, fig4);
        Screen('DrawTexture', window.onScreen, Slide4, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        fig6=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide6.jpg');
        Slide6=Screen('MakeTexture', window.onScreen, fig6);
        Screen('DrawTexture', window.onScreen, Slide6, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        fig7=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide7.jpg');
        Slide7=Screen('MakeTexture', window.onScreen, fig7);
        Screen('DrawTexture', window.onScreen, Slide7, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        fig8=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide8.jpg');
        Slide8=Screen('MakeTexture', window.onScreen, fig8);
        Screen('DrawTexture', window.onScreen, Slide8, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        goalFig=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/goal1.jpg');
        goal1=Screen('MakeTexture', window.onScreen, goalFig);
        Screen('DrawTexture', window.onScreen, goal1, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        paramFig=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/taskParams.jpg');
        param=Screen('MakeTexture', window.onScreen, paramFig);
        Screen('DrawTexture', window.onScreen, param, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,~,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        instructionText = 'We will now get into the practice block.\n You will receive feedback on your estimation of the colors.\n\n\n Try not to blink a lot after stimulus onset.\n You can use the period before stimulus onset to blink.';
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        DrawFormattedText(window.onScreen,instructionText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        instructionText = 'You will only have 10s to make a response.\n\n\n Ready? Left click to continue.';
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        DrawFormattedText(window.onScreen,instructionText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.01)
        end
        
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.instructionsOff);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
        
        
    elseif block==2
        
        %         io64(ioObject,LPT1address,prefs.triggerNum.instructions);
        %         WaitSecs(0.01)
        %         io64(ioObject,LPT1address,0);
        instructionText = 'We will get into the actual task now.\n\n\n In this block, you will be asked to estimate the color you just saw at the prompted stimulus location.\n Remember,you will only have 10s for each response. \n\n\n You will no longer receive feedback on your perception from now on. \n\n\n Left click to continue.';
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);DrawFormattedText(window.onScreen,instructionText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            WaitSecs(.1)
        end
        
    else
        if block==3 
            
            Screen('TextSize',window.onScreen,40);
            instructionText = 'You have successfully completed the first half of the task! \n\n For further instructions for the next blocks, \n please wait for your experimenter. ';
            Screen('TextFont',window.onScreen,'Times');
            Screen('TextStyle',window.onScreen,0);
            Screen('TextColor',window.onScreen, 255);
            DrawFormattedText(window.onScreen,instructionText,'center','center')
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
            
            if prefs.doET
                io64(ioObject,LPT1address,prefs.triggerNum.instructionsOn);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
            end
            fig9=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide9.jpg');
            Slide9=Screen('MakeTexture', window.onScreen, fig9);
            Screen('DrawTexture', window.onScreen, Slide9, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.01)
            end
            
            if prefs.doPredInstructions == 1
                fig10=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide10.jpg');
                Slide10=Screen('MakeTexture', window.onScreen, fig10);
                Screen('DrawTexture', window.onScreen, Slide10, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
                Screen(window.onScreen, 'Flip');
                pause(1)
                [x,y,buttons]=GetMouse;
                while ~any(buttons)
                    [x,y,buttons]=GetMouse
                    WaitSecs(.01)
                end
                
                fig11=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/Slide11.jpg');
                Slide11=Screen('MakeTexture', window.onScreen, fig11);
                Screen('DrawTexture', window.onScreen, Slide11, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
                Screen(window.onScreen, 'Flip');
                pause(1)
                [x,y,buttons]=GetMouse;
                while ~any(buttons)
                    [x,y,buttons]=GetMouse
                    WaitSecs(.01)
                end

            else
                Screen('TextSize',window.onScreen,40);
                instructionText ='While you cannot predict the colors you will see exactly, \n you might notice that some colors are shown more frequently during certain periods of the task. \n\n The way that colors change over time will differ in the next two blocks. \n\n It is your job to figure out how the colors are changing to \n help you make more accurate predictions and score more points.';
                DrawFormattedText(window.onScreen,instructionText,'center','center')
                Screen(window.onScreen, 'Flip');
                pause(1)
                [x,y,buttons]=GetMouse;
                while ~any(buttons)
                    [x,y,buttons]=GetMouse
                    WaitSecs(.1)
                end

                Screen('TextSize',window.onScreen,40);
                instructionText ='Here are a few pointers before we start \n\n 1. Nothing you do will affect the sequence, it is all pregenerated. \n\n 2. The left and right colors are still independent of one another, but they follow similar structures. \n\n 3. Pay the most attention to what happens after a drastic color change';
                DrawFormattedText(window.onScreen,instructionText,'center','center')
                Screen(window.onScreen, 'Flip');
                pause(1)
                [x,y,buttons]=GetMouse;
                while ~any(buttons)
                    [x,y,buttons]=GetMouse
                    WaitSecs(.1)
                end
            end

            goalFig2=imread('/Users/mattLab/Dropbox (Brown)/ALP_taskFolder/Instruction Images/goal2.jpg');
            goal2=Screen('MakeTexture', window.onScreen, goalFig2);
            Screen('DrawTexture', window.onScreen, goal2, [],CenterRectOnPoint([0 0 prefs.imageSizeX prefs.imageSizeY],window.centerX, window.centerY));
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.01)
            end
            Screen('TextSize',window.onScreen,40);
            instructionText ='Remember you will only have 10s per response. \n\n Ready? Left click to enter the next two blocks.';
            DrawFormattedText(window.onScreen,instructionText,'center','center')
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
            
            if prefs.doET
                io64(ioObject,LPT1address,prefs.triggerNum.instructionsOff);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
            end
        
        end
        
        
        if prefs.taskType == 10 && prefs.doPredInstructions == 1
            
            %             io64(ioObject,LPT1address,prefs.triggerNum.instructions);
            %             WaitSecs(0.01)
            %             io64(ioObject,LPT1address,0);
            %             instructionText1 = 'In the next block, you will be asked to: \n\n 1). Reproduce the color you just saw at the prompted stimulus location \n 2). Predict the color you will see in the next trial at the prompted stimulus location. \n\n\n *REMEMBER* \n\n Your goal is to accurately reproduce the colors you saw on each trial and predict what you will see in the next trial based on your belief.';
            %             DrawFormattedText(window.onScreen,instructionText1,'center',window.screenY*0.1);
            Screen('DrawTexture', window.onScreen, instr_textures{1}, [],CenterRectOnPoint([0 0 1000 250],window.centerX, window.screenY*0.3));
            %Screen('DrawTexture', window.onScreen, instr_textures{1}, [], CenterRectOnPoint([0 0 200 200],window.centerX, window.centerY));
            instructionText2 = '\n\n *HINT*: \n\n In this block, the colors of each stimulus will change in a way similar to the situation above across trials.\n\n Each circle in the above figure represents the color of the same stimulus at each trial. \n\n Pay attention to how the color is changing. \n\n\n Left click to continue.';
            %The uncertainty in this block results from a discontinuous changes , \n they reflect a complete change in the way the colors are generated. \n\n\n *HINT* \n\n The most recent information is the most predictive of the future. \n Click to continue.';
            Screen('TextSize',window.onScreen,40);
            Screen('TextFont',window.onScreen,'Times');
            Screen('TextStyle',window.onScreen,0);
            Screen('TextColor',window.onScreen, 255);
            DrawFormattedText(window.onScreen,instructionText2,'center', window.screenY*0.5);
            %DrawFormattedText(window.onScreen,instructionText,'center','center')
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
        elseif prefs.taskType == 11 && prefs.doPredInstructions == 1
            
            %             io64(ioObject,LPT1address,prefs.triggerNum.instructions);
            %             WaitSecs(0.01)
            %             io64(ioObject,LPT1address,0);
            %             instructionText1 = 'In the next block, you will be asked to: \n\n 1). Reproduce the color you just saw at the prompted stimulus location \n 2). Predict the color you will see in the next trial at the prompted stimulus location. \n\n\n *REMEMBER* \n\n Your goal is to accurately reproduce the colors you saw on each trial and predict what you will see in the next trial based on your belief.';
            %             DrawFormattedText(window.onScreen,instructionText1,'center',window.screenY*0.1);
            Screen('DrawTexture', window.onScreen, instr_textures{2}, [], CenterRectOnPoint([0 0 1000 250],window.centerX, window.screenY*0.3));
            instructionText2 = '\n\n *HINT*: \n\n In this block, the colors of each stimulus will change in a way similar to the situation above across trials.\n\n Each circle in the above figure represents the color of the same stimulus at each trial. \n\n Pay attention to how the color is changing. \n\n\n Left click to continue.';
            Screen('TextSize',window.onScreen,40);
            Screen('TextFont',window.onScreen,'Times');
            Screen('TextStyle',window.onScreen,0);
            Screen('TextColor',window.onScreen, 255);
            DrawFormattedText(window.onScreen,instructionText2,'center', window.screenY*0.5);
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
        end
        
        
    end
    
    %
    %         %####################### GIVE SOME INSTRUCTIONS ###########################
    %
    %         instructionText = 'Remember to fixate on the central point at all times! \n\n\n';
    %         DrawFormattedText(window.onScreen,instructionText,'center','center')
    %         Screen(window.onScreen, 'Flip');
    %         [x,y,buttons]=GetMouse;
    %         while ~any(buttons)
    %             [x,y,buttons]=GetMouse
    %             WaitSecs(.1)
    %         end
    
    % put up instructions and wait for keypress
    %drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
    if block==1
        instruct(window);
    else
        %         Screen('TextColor',window.onScreen, 255);
        %     Screen('TextSize',window.onScreen,40);
        %     instructionText = 'Please take a break and close your eyes for a minute.';
        %
        %     DrawFormattedText(window.onScreen,instructionText,'center','center')
        %     %returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize, prefs)
        %
        %     Screen('Flip', window.onScreen);
        %     WaitSecs(60);
        instructBetween(window);
    end
    %returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize, prefs)
    
    %     WaitSecs(1);
    %
    %     if block==1
    %         Screen('TextSize',window.onScreen,50);
    %         Screen('TextColor',window.onScreen, 255);
    %         blockName='0: Practice Block';
    %         DrawFormattedText(window.onScreen,blockName,'center','center')
    %     elseif block==2
    %         Screen('TextSize',window.onScreen,50);
    %         Screen('TextColor',window.onScreen, 255);
    %         blockName='1: Estimation Block, no prediction';
    %         DrawFormattedText(window.onScreen,blockName,'center','center')
    %     elseif block==3
    %         Screen('TextSize',window.onScreen,50);
    %         Screen('TextColor',window.onScreen, 255);
    %         blockName='2: Estimation and Prediction Block A';
    %         DrawFormattedText(window.onScreen,blockName,'center','center')
    %     elseif block==4
    %         Screen('TextSize',window.onScreen,50);
    %         Screen('TextColor',window.onScreen, 255);
    %         blockName='3: Estimation and Prediction Block B';
    %         DrawFormattedText(window.onScreen,blockName,'center','center')
    %     end
    
    Screen('Flip', window.onScreen);
    WaitSecs(1.5);
    
    %drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
    Screen('TextColor',window.onScreen, 255);
    Screen('TextSize',window.onScreen,40);
    instructionText = 'Remember to fixate on the central point at all times!';
    
    DrawFormattedText(window.onScreen,instructionText,'center','center')
    %returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize, prefs)
    
    Screen('Flip', window.onScreen);
    WaitSecs(2);
    
    
    % put up instructions and wait for keypress
    %     drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
    %     instruct(window);
    %     returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize, prefs)
    %
    %     WaitSecs(1);
    
    
    
    % decide colors and orientation.
    if  prefs.propFixedSpace==0 & prefs.taskType ==0% if we're picking at random, do it.
        
        colorsInDegrees = random('Uniform', 0, 360, prefs.nTrials, prefs.nItems);
        colors = colorsInDegrees ./ 360;
        %         anglesInDegrees= random('Uniform', 0, 180, prefs.nTrials, prefs.nItems);
        %         angles=anglesInDegrees./180;
        anglesInDegrees = repmat([90], prefs.nTrials, prefs.nItems);
        angles = anglesInDegrees./360;
        
        isRand=true(length(colorsInDegrees), 1);
        
        
        %----------------------------------------
        %This is the generative structure for the change-point condition
        %----------------------------------------
    elseif prefs.taskType == 10
        %         colorsInDegrees = repmat([180], prefs.nTrials, prefs.nItems);
        %         colors = colorsInDegrees./ 360;
        %         anglesInDegrees = repmat([90], prefs.nTrials, prefs.nItems);
        %         angles = anglesInDegrees./360;
        
        nTrials   = prefs.nTrials;   % number of trials
        haz       = .15;   % Hazard rate
        sigma_x   = 10;   % variance on distribution of stimuli/outcomes around mean
        X= ones(nTrials,prefs.nItems).*-1;
        
        
        
        for t = 1:nTrials
            
            % SAMPLE CHANGEPOINT:
            % S_t    ~B(H)
            for j = 1:prefs.nItems
                S(t,j)=rand<haz;
                
                % SAMPLE MU:
                % if S_t = 0:
                % mu_t   =mu_{t-1}
                % elseif S_t = 1
                % mu_t   ~ U(0-180)
                
                if S(t,j)==1 || t ==1
                    Mu(t,j)=rand*360;
                    change(t,j)=1;
                else
                    Mu(t,j)=Mu(t-1,j);
                    change(t,j)=0;
                end
                
                
                
                % SAMPLE X:
                % X_t    = N(mu_t, sigma_s)
                X(t,j) = normrnd(Mu(t,j), sigma_x);
                if X(t,j)<0|| X(t,j)>360
                    X(t,j)=mod(X(t,j), 360);
                end
            end
        end
        generativeMean=Mu;
        surprise=change;
        
        %data.generativeMean = Mu;
        colorsInDegrees = X;
        colors = colorsInDegrees./ 360;
        anglesInDegrees = repmat([90], prefs.nTrials, prefs.nItems);
        angles = anglesInDegrees./360;
        isRand=false(length(colorsInDegrees), 1);
        
        
        %This is the odd-ball condition
    elseif prefs.taskType == 11
        nTrials = prefs.nTrials;
        H = .15;
        sigma_c = 10;
        sigma_b = 6.5;
        
        B= ones(nTrials, prefs.nItems).*-1;
        C= ones(nTrials, prefs.nItems).*-1;
        S= ones(nTrials, prefs.nItems).*-1;
        for t = 1:nTrials
            for j = 1:prefs.nItems
                %     If it is the first trial we need to create a random point for C to
                %     start
                if t == 1
                    C(t,j) = rand*360;
                    change(t,j)=1;
                else
                    C(t,j) = normrnd(C(t-1,j), sigma_c);
                    if C(t,j)<0|C(t,j)>360
                        C(t,j)=mod(C(t,j), 360);
                    end
                end
                
                %     Generate the S value.
                %        ==>  If S  = 1 with p=H, then the value of B is not
                %             going to be drawn from C but from a uniform distribution
                %        ==> .If S = 0, with probability 1-H then the vaue of B will be
                %        drawn from a normal distribution with mean C and variance of
                %        sigma_b
                if rand < H
                    S(t,j) = 1;
                    B(t,j) = rand*360
                    change(t,j)=1;
                else
                    S(t,j) = 0;
                    B(t,j) = normrnd(C(t,j), sigma_b);
                    if B(t,j)<0|B(t,j)>360
                        B(t,j)=mod(B(t,j), 360);
                    end
                    change(t,j)=0;
                end
            end
            
        end
        generativeMean=C;
        surprise=change;
        
        %data.generativeMean = C;
        colorsInDegrees = B;
        colors = colorsInDegrees./ 360;
        if prefs.jittered
            anglesInDegrees= random('Uniform', 0, 180, prefs.nTrials, prefs.nItems);
        else
            anglesInDegrees = repmat([90], prefs.nTrials, prefs.nItems);
        end
        angles = anglesInDegrees./360;
        isRand=false(length(colorsInDegrees), 1);
        
        
    else % or else, choose which trials will be fixed, choose one color for each fixed trial and space the others appropriately.
        
        % OK. lets only randomize the estimation dimension, lets always
        % keep fixed spacing on the probe dimension
        
        % how many trials of each type?
        numFixed=ceil(prefs.nTrials.*prefs.propFixedSpace);
        numRand =floor(prefs.nTrials.*(1-prefs.propFixedSpace));
        
        
        % create mid colors for fixed trials
        % for each fixed trial, choose a random between 0 and 360
        midColor=random('Uniform', 0, 360, numFixed, 1);
        midOrient=random('Uniform', 0, 360, numFixed, 1);
        oColor=nan(length(midColor), size(prefs.fixedSpacings, 2));
        oOrient=nan(length(midOrient), size(prefs.fixedOSpacings, 2));
        % run through each trial, determine a set of fixed dists and
        % create target locations:
        for tt=1:length(midColor)
            % pick from the rows of fixed spacings
            fixPick=randperm(size(prefs.fixedSpacings, 1), 1);
            fixPickO=randperm(size(prefs.fixedOSpacings, 1), 1);
            % add in all columns from that row:
            oColor(tt, :)=midColor(tt)+prefs.fixedSpacings(fixPick, :);
            oOrient(tt,:)=midOrient(tt)+prefs.fixedOSpacings(fixPickO, :)
            
        end
        
        fSpaceColors=rad2deg(circ_dist(deg2rad([midColor oColor]), 0)+pi);
        fSpaceOrient= (rad2deg(circ_dist(deg2rad([midOrient oOrient]), 0)+pi))./2;
        randColors  =random('Uniform', 0, 360, numRand, prefs.nItems);
        
        
        % if we are estimating color, do NOT randomize orientation, fixed
        % space it.
        if prefs.taskType==4
            %
            midOrient=random('Uniform', 0, 360, numRand, 1);
            oOrient=nan(length(midOrient), size(prefs.fixedOSpacings, 2));
            % run through each trial, determine a set of fixed dists and
            % create target locations:
            
            %keyboard
            for tt=1:length(midColor)
                fixPickO=randperm(size(prefs.fixedOSpacings, 1), 1);
                oOrient(tt,:)=midOrient(tt)+prefs.fixedOSpacings(fixPickO, :)
            end
            randOrient= (rad2deg(circ_dist(deg2rad([midOrient oOrient]), 0)+pi))./2;
        else
            randOrient  =random('Uniform', 0, 360, numRand, prefs.nItems)./2;
        end
        
        
        
        
        %interleave types randomly.:
        colorsInDegrees=[ randColors; fSpaceColors];
        anglesInDegrees=[ randOrient; fSpaceOrient];
        
        isRand=false(length(colorsInDegrees), 1);
        isRand(1:length(randColors))=true;
        
        
        ind=randperm(prefs.nTrials)
        
        % store an array that tells us whether each trial was fixed or
        % random spacing.
        isRand=isRand(ind);
        colorsInDegrees=colorsInDegrees(ind,:);
        
        for k = 1:length(colorsInDegrees)
            colorsInDegrees(k,:) = colorsInDegrees(k,(randperm(prefs.nItems)));
            anglesInDegrees(k,:) = anglesInDegrees(k,(randperm(prefs.nItems)));
        end
        
        colors = colorsInDegrees ./ 360;
    end
    
    
    allStartPoint = nan(prefs.nTrials,prefs.nItems);
    allPredStartPoint = nan((prefs.nTrials-1),prefs.nItems);
    reportedValue = nan(prefs.nTrials,prefs.nItems);
    predictedValue = nan(prefs.nTrials,prefs.nItems);
    allMousePos = cell(prefs.nTrials,prefs.nItems);
    allMousePosPred = cell(prefs.nTrials,prefs.nItems);
    
    %  Each of the actual trials with both the probe and the estimate
    for trialIndex = 1:prefs.nTrials
        if block==1 && trialIndex==2
            instructionText='Continue practice with the next few trials. \n\n Pay attention to the feedback and \n learn how accurate your color choice has to be \n in order to be counted as a correct answer.\n\n Left click to continue.';
            DrawFormattedText(window.onScreen,instructionText,'center','center');
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
        end
        if sum([block==[3,4], trialIndex==[11,51,91]]) == 2
            instructionText='Take a short break and close your eyes for a few seconds. \n\n Use this time to contemplate what you have seen on the last few trials \n and come up with a strategy for predicting future trials in the block. \n\n The sequence of colors will continue from where you left off';
            DrawFormattedText(window.onScreen,instructionText,'center','center');
            Screen(window.onScreen, 'Flip');
            pause(1)
            [x,y,buttons]=GetMouse;
            while ~any(buttons)
                [x,y,buttons]=GetMouse
                WaitSecs(.1)
            end
        end
        drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
        Screen('Flip', window.onScreen);
        WaitSecs(prefs.stimulusDuration);
        %draw larger fixation for attention before stimulus show
        drawFixation(window, window.centerX, window.centerY, 10, 255);
        Screen('Flip', window.onScreen);
        WaitSecs(prefs.stimulusDuration);
        
        % Get new target coordinates (including jitter).
        [coor, locJit] = circularArrayRects([0, 0, prefs.rectSize1, prefs.rectSize2], prefs.nItems, prefs.radius, window.centerX, window.centerY);
        space = (window.screenX/(prefs.nItems + 1));
        if prefs.taskType == 10 | prefs.taskType == 11 |prefs.taskType == 0
            Xpositions = [space:space:(window.screenX-space)];
            Xpositions = Xpositions;
            display(size(Xpositions, 2));
            Ypositions = repmat([(window.screenY/2)], 1, prefs.nItems);
            display(size(Ypositions, 2));
            coor = [Xpositions; Ypositions];
        end
        allLocJit(trialIndex)=locJit;
        
        
        % CHOOSE COLORS TO DISPLAY
        if prefs.useCieLab
            if isfield(prefs.monitorData, 'constantB')
                
                colorsToDisplay=round(hueTo_calRGB(colors(trialIndex, :), prefs.monitorData))';
            else
                colorsToDisplay=round(255.*hueTo_calRGB(colors(trialIndex, :), prefs.monitorData));
            end
        else
            disp('HSV colors NOT enabled!!!!')
        end
        
        
        % draw fixation
        drawFixation(window, window.centerX, window.centerY, prefs.fixationSize, 255);
        
        %             rotate stimuli to appropriate orientation:
        %
        %           open GL CODE TO IMPLEMENT RECTANGLE ROTATION
        
        baseRect=[0 0 prefs.rectSize1, prefs.rectSize2];
        %  SHOWING THE Stimulus TO BE ESTIMATED-----------------------------
        for i = 1:size(coor, 2)
            
            % Get the current squares position and rotation angle
            posX = coor(1, i);
            posY = coor(2, i);
            display(posX);
            display(posY);
            targAngle = anglesInDegrees(trialIndex, i);
            targColor=colorsToDisplay(:,i)
            %  keyboard
            % Translate, rotate, re-tranlate and then draw our square
            Screen('glPushMatrix', window.onScreen)
            Screen('glTranslate', window.onScreen, posX, posY)
            Screen('glRotate', window.onScreen, targAngle, 0, 0);
            Screen('glTranslate', window.onScreen, -posX, -posY)
            %             if prefs.doEEG
            %                 io64(ioObject,LPT1address,prefs.triggerNum.stimOn);
            %                 WaitSecs(0.01)
            %                 io64(ioObject,LPT1address,0);
            %             end
            Screen('FillRect', window.onScreen, targColor, CenterRectOnPoint(baseRect, posX, posY));
            Screen('glPopMatrix', window.onScreen)
        end
        
        % Flip to the screen
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.stimOn);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
        Screen('Flip', window.onScreen);
        WaitSecs(prefs.stimulusDuration);
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.stimOff);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
        
        %       END OF Stimulus---------------------
        %       BEGINING OF PROBE--------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%       Optionally show a gray mask to make sure color after images don't remain      %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isfield(prefs, 'showMask') && prefs.showMask==true
            
            for i = 1:size(coor, 2)
                
                % Get the current squares position ans rotation angle
                posX = coor(1, i);
                posY = coor(2, i);
                targAngle = anglesInDegrees(trialIndex, i);
                targColor=colorsToDisplay(:,i)
                %  keyboard
                % Translate, rotate, re-tranlate and then draw our square
                Screen('glPushMatrix', window.onScreen)
                Screen('glTranslate', window.onScreen, posX, posY)
                Screen('glRotate', window.onScreen, targAngle, 0, 0);
                Screen('glTranslate', window.onScreen, -posX, -posY)
                Screen('FillRect', window.onScreen, prefs.testGray, CenterRectOnPoint(baseRect, posX, posY));
                Screen('glPopMatrix', window.onScreen)
                
            end
            
            % Flip to the screen
            Screen('Flip', window.onScreen);
            WaitSecs(1); % wait for one frame... then carry on
            
        end
        
        
        
        
        
        
        
        
        % remove stimulus for retention interval
        returnToFixation(window, window.centerX, window.centerY, prefs.fixationSize, prefs);
        WaitSecs(prefs.retentionInterval);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                            PROBE LOGIC                          %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % prefs.ProbeType: How will we probe the subject to report the stimulus?
        %                  0 = show position of probed target
        %                  1 = show color of probed target (in center screen)
        %                  2 = show orientation of probed target (in center screen)
        
        % pick a target to test... get the color and angle of that target.
        %targetToTest(trialIndex) = Ranint(prefs.nItems);
        %probeOrder=randperm(prefs.nItems,prefs.nItems);
        probeOrder=randperm(prefs.nItems);
        
        for i=1:prefs.nItems
            targetToTest(trialIndex,i) = probeOrder(i);
            presentedColor(trialIndex,i) = colorsInDegrees(trialIndex, targetToTest(trialIndex,i));
            presentedAngle(trialIndex,i) = anglesInDegrees(trialIndex, targetToTest(trialIndex,i));
        end
        
        
        
        % prefs.taskType: Pick one of several task types
        
        %                  0 = probe is location, estimation = color
        %                  1 = probe is location, estimation = orientation
        %                  2 = probe is location, estimation = both (50/50 mix)
        
        
        %                  3 = probe is color, estimate = orientation
        %                  4 = probe is orientation, estimate = color
        %                  5 = mixture (50/50) of types 3 and 4.
        
        % each trial might have a different type of probe:
        
        
        
        
        % Determine which dimension will be used for estimation and which
        % dimension will be used for a probe.
        
        if prefs.taskType <=2
            probeType=2;
            if prefs.taskType==0|(prefs.taskType==0&rand>.5)
                estType  =0; %estimation of color
            else
                estType  =1; % estimation of orientation
            end
        elseif prefs.taskType ==3 | (prefs.taskType == 5 & rand > .5)
            probeType=1;
            estType  =1;  % estimation will be done on orientation
        elseif prefs.taskType == 10 | (prefs.taskType == 11)
            probeType = 2;
            estType = 0;
        else
            probeType=2;
            estType  =0 % estimation will be done on color
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                              SHOW PROBE                               %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        allProbeType(trialIndex)=probeType;
        allEstType(trialIndex)=estType;
        
        for i=1:prefs.nItems
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%         Probe            %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if prefs.doET && i==1
                io64(ioObject,LPT1address,prefs.triggerNum.cue1on);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
            end
            
            if prefs.doET && i==2
                io64(ioObject,LPT1address,prefs.triggerNum.cue2on);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
            end
            
            switch probeType
                case 0
                    
                    % show probed target in gray.
                    % choose a circle to test, then display response screen
                    colorsOfTest = repmat([10 10 10], prefs.nItems, 1);
                    if prefs.useCieLab
                        colorsOfTest(targetToTest(trialIndex), :) = [200 200 200];
                    else
                        colorsOfTest(targetToTest(trialIndex), :) = [145 145 145];
                    end
                    drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
                    Screen('FillRect', window.onScreen, colorsOfTest', rects);
                    drawColorWheel(window, prefs);
                    SetMouse(window.centerX, window.centerY);
                    ShowCursor('Arrow');
                    
                    % NEED TO STORE PROBE CENTER!!!
                    
                case 1
                    % show probed colored square in center of screen
                    defBox=[0 0 40 40];
                    % this should do it... though its not clear what shape we
                    % should show...
                    
                    probeColor=colorsToDisplay(:, targetToTest(trialIndex))
                    
                    Screen('FillRect', window.onScreen,  probeColor', CenterRectOnPoint(defBox, window.centerX, window.centerY));
                    probeCent=[window.centerX, window.centerY]
                    
                    
                case 2
                    % show probed orientation rectangle in gray at center of screen
                    k=targetToTest(trialIndex,i)
                    posX = coor(1, k);
                    posY = coor(2, k);
                    display(posX);
                    display(posY);
                    targAngle = anglesInDegrees(trialIndex, k);
                    % Translate, rotate, re-tranlate and then draw our square
                    Screen('glPushMatrix', window.onScreen)
                    Screen('glTranslate', window.onScreen, window.centerX, window.centerY)
                    %                 Screen('glRotate', window.onScreen, targAngle, 0, 0);
                    Screen('glTranslate', window.onScreen, -window.centerX, -window.centerY)
                    Screen('FillRect', window.onScreen, [128, 128, 128], CenterRectOnPoint(baseRect,posX, posY));
                    %                 Screen('FillRect', window.onScreen, prefs.testGray, CenterRectOnPoint(baseRect,window.centerX, window.centerY));
                    Screen('glPopMatrix', window.onScreen)
                    probeCent=[window.centerX, window.centerY]
                    
                    
            end
            
            
            Screen('Flip', window.onScreen);
            
            
            %% OK, subject should now see the probe.
            %  Next step involves the subject estimating either orientation or
            %  color based on the cue...
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%                   ALLOW SUBJECT TO MAKE ESTIMATE                      %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             Screen('FrameRect', window,[127],[200 200 400 400])
            %             WaitSecs(prefs.preEstPause);
            
            
            
            
            
            [x,inY,buttons] = GetMouse(window.onScreen);
            inY=inY./prefs.mouseScale;
            
            
            
            
            while any(buttons) % if already down, wait for release
                [x,y,buttons] = GetMouse(window.onScreen);
                y=y./prefs.mouseScale;
            end
            
            
            startPoint=round(rand.*360);
            %initRT(trialIndex)=nan;
            t1=GetSecs;
            estRT(trialIndex,i)=nan;
            MaxWait = 10;  %  maximum RT of  seconds
            StopTime = GetSecs + MaxWait;  % This computes the time at which you want to stop (max RT)
            mousePos = [];

            while ~any(buttons)  && (GetSecs<StopTime) % wait for press
                
                [x,y,buttons] = GetMouse(window.onScreen);
                mousePos = [mousePos,y];
                y=y./prefs.mouseScale;

                %                 if ~isfinite(initRT(trialIndex))&&abs(y-inY)>prefs.initRTthresh
                %                     initRT(trialIndex)= GetSecs-t1
                %                 end
                
                if ~isfinite(estRT(trialIndex,i))&&abs(y-inY)<prefs.responseRTthresh&&abs(y-inY)>prefs.initRTthresh
                    estRT(trialIndex,i)= GetSecs-t1;
                end
                
                if estType==0     % if we're estimating color, allow color graphics to be changed
                    nValues=360;
                    
                    chosenVal=mod(y-inY+startPoint, nValues);
                    %chosenVal=y-inY+startPoint;
                    
                    % DISPLAY "estimate"
                    if prefs.useCieLab
                        if isfield(prefs.monitorData, 'constantB')
                            
                            currentColor=round(hueTo_calRGB(chosenVal./360, prefs.monitorData))';
                        else
                            currentColor = round(255.*hueTo_calRGB(chosenVal./360, prefs.monitorData));
                        end
                    else
                        disp('HSV colors NOT enabled!!!!') % degrees on color wheel
                    end
                    
                    
                    
                    k=targetToTest(trialIndex,i); % set orientation to cue value
                    currentOrientation = anglesInDegrees(trialIndex, k);
                    
                elseif estType==1
                    nValues=180;
                    chosenVal=mod((y./2)+startPoint, nValues)
                    
                    
                    i=targetToTest(trialIndex);
                    currentColor=probeColor; % set color to cue value
                    currentOrientation = chosenVal;
                    
                    
                end
                
                
                % Recompute the angle every time through the loop? This seems
                % ridiculous...
                %  SUBMITTING ESTIMATE ---------------------------------------------------
                % show probed orientation rectangle in gray at center of screen
                
                % Translate, rotate, re-tranlate and then draw our square
                k=targetToTest(trialIndex,i)
                posX = coor(1, k);
                posY = coor(2, k);
                Screen('glPushMatrix', window.onScreen);
                Screen('glTranslate', window.onScreen, window.centerX, window.centerY);
                %             Screen('glRotate', window.onScreen, currentOrientation, 0, 0);
                Screen('glTranslate', window.onScreen, -window.centerX, -window.centerY);
                drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
                Screen('FillRect', window.onScreen, currentColor', CenterRectOnPoint(baseRect,posX, posY));
                Screen('glPopMatrix', window.onScreen);
                if trialIndex==1 && block==1 && i==1
                    instructionText='Now, you can move your cursor up and down to \n find the color you think you saw at this stimulus location.\n\n Once you find the color, you can select by left click the mouse.\n\n You will only have 10s to make this response!';
                    DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                elseif trialIndex==1 && block==1 && i==2
                    instructionText='Again, find the color you think you saw at this stimulus location \n by moving your cursor up and down.\n\n Once you find the color, you can select by left click the mouse.\n\n Again, you will only have 10s to make this response!';
                    DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                    
                else
                    instructionText='Reproduce the color of this square on this trial.';
                    DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                end
                %Screen('DrawText', window.onScreen, 'Estimate the color of this square on this trial.', 'center', window.centerY*0.2, 255);
                %Screen('FrameRect', window.onScreen,255,[100 150 700 450])
                Screen('Flip', window.onScreen);
                
                
            end
            
            allStartPoint(trialIndex,i) = startPoint;
            reportedValue(trialIndex,i) = chosenVal;
            allMousePos{trialIndex,i} = mousePos;

            while any(buttons) % wait for releasessssss
                [x,y,buttons] = GetMouse(window.onScreen);
                y=y./prefs.mouseScale;
            end
            
            
            
            if prefs.doET && i==1
                io64(ioObject,LPT1address,prefs.triggerNum.responseMade1);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
                %WaitSecs(0.5)
            end
            
            if prefs.doET && i==2
                io64(ioObject,LPT1address,prefs.triggerNum.responseMade2);
                WaitSecs(0.01)
                io64(ioObject,LPT1address,0);
            end
            % return to fixation
            returnToFixation_estimate(window, window.centerX, window.centerY, prefs.fixationSize, prefs);
            WaitSecs(prefs.ITI);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%                       Optional Feedback                         %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %if prefs.doFeedback
            
            if estType==0
                trialErr=abs(circ_dist(deg2rad(reportedValue(trialIndex,i)), deg2rad(presentedColor(trialIndex,i))));
                if trialErr<prefs.fdbk_cThresh
                    correct=1;
                else
                    correct=-1;
                end
            else
                
                trialErr=abs(circ_dist(2.*deg2rad(reportedValue(trialIndex,i)), 2.*deg2rad(anglesInDegrees(trialIndex, targetToTest(trialIndex,i)))));
                if trialErr<prefs.fdbk_oThresh
                    correct=1;
                else
                    correct=-1;
                end
            end
            
            
            
            
            % If required, flip some feedback at random.
            if rand<prefs.fdbkProp
                posFdbk=correct.*-1;
            else
                posFdbk=correct;
            end
            
            allTrialErr(trialIndex,i) = trialErr;
            allCorr(trialIndex,i)     = correct;
            allPosFdbk(trialIndex,i)  = posFdbk;
            %score for reprduction
            allScoreEst(trialIndex,i)    = correct;
            
            
            if prefs.doFeedback
                
                if posFdbk==1
                    if trialIndex==1
                        instructionText='Good job, you just correctly reproduced \n the color you saw for this square on this trial.\n\n Left click to continue.';
                        DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                        Screen('DrawTexture', window.onScreen, fdbk_textures{2}, [], CenterRectOnPoint([0 0 200 200],window.centerX, window.centerY));
                        
                        
                    else
                        % Show positive feedback
                        Screen('DrawTexture', window.onScreen, fdbk_textures{2}, [], CenterRectOnPoint([0 0 200 200],window.centerX, window.centerY));
                    end
                elseif posFdbk==-1
                    if trialIndex==1
                        instructionText='Sorry, the color you chose is not close to what you saw. \n\n Remember to only reproduce the color you saw at the prompted location. \n\n Left click to continue. ';
                        DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                        Screen('DrawTexture', window.onScreen, fdbk_textures{1}, [], CenterRectOnPoint([0 0 200 200],window.centerX, window.centerY));
                    else
                        % show negative feedback
                        Screen('DrawTexture', window.onScreen, fdbk_textures{1}, [], CenterRectOnPoint([0 0 200 200],window.centerX, window.centerY));
                    end
                end
                
                if trialIndex==1
                    Screen(window.onScreen, 'Flip');
                    pause(1)
                    [x,y,buttons]=GetMouse;
                    while ~any(buttons)
                        [x,y,buttons]=GetMouse
                        WaitSecs(.1)
                    end
                else
                    Screen('Flip', window.onScreen);
                    WaitSecs(prefs.fdbk_duration);
                end
                
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                       Optional PDW                              %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prefs.pdw_prop:      Proportion of trials that will have a post decision wager
        % prefs.pdw_vals:      Possible bet values [array of two integers]
        % prefs
        
        
        
        if rand<prefs.pdw_prop
            pdw(trialIndex)=true;
            
            while any(buttons) % if already down, wait for release
                [x,y,buttons] = GetMouse(window.onScreen);
                y=y./prefs.mouseScale;
                
            end
            
            
            SetMouse(window.centerX, window.centerY, window.onScreen)
            while ~any(buttons)|~isfinite(bet) % wait for decision and press.
                
                [x, y, buttons] = GetMouse(window.onScreen);
                locVar=x;
                
                if locVar>window.centerX+10
                    msg='Bet!';
                    bet=max(prefs.pdw_vals);
                elseif locVar<window.centerX-10
                    msg='Pass';
                    bet=min(prefs.pdw_vals);
                else
                    msg='Pass or Bet?';
                    bet=nan;
                end
                
                Screen('TextSize',window.onScreen,30);
                Screen('TextFont',window.onScreen,'Times');
                Screen('TextStyle',window.onScreen,0);
                Screen('TextColor',window.onScreen, window.white);
                DrawFormattedText(window.onScreen,msg,'center','center')
                Screen(window.onScreen, 'Flip');
                
            end
            
            
            while any(buttons) % wait for releasessssss
                [~,~,buttons] = GetMouse(window.onScreen);
            end
        else
            bet = 1;
            pdw(trialIndex)=false;
        end
        
        
        allBet(trialIndex)=bet;
        
        
        
        %% **************************************************************%%
        %                        PREDICTION PART                          %
        % ***************************************************************%%
        
        allPredErr = nan(trialIndex,prefs.nItems);
        allPredCorr= nan(trialIndex,prefs.nItems);
        
        if prefs.taskType~=0
            
            
            if prefs.doPredict && trialIndex < prefs.nTrials
                
                WaitSecs(prefs.preEstPause);
                
                predictOrder=randperm(prefs.nItems,prefs.nItems);
                
                for i=1:prefs.nItems
                    targetToPredict(trialIndex,i) = predictOrder(i);
                    toPredictColor(trialIndex,i) = colorsInDegrees(trialIndex+1, targetToPredict(trialIndex,i));
                    toPredictAngle(trialIndex,i) = anglesInDegrees(trialIndex+1, targetToPredict(trialIndex,i));
                end
                
                %                 [x,inY,buttons] = GetMouse(window.onScreen);
                %                 inY=inY./prefs.mouseScale;
                %
                %                 while any(buttons) % if already down, wait for release
                %                     [x,y,buttons] = GetMouse(window.onScreen);
                %                     y=y./prefs.mouseScale;
                %                 end
                %                 while ~any(buttons) % wait for press
                %
                %                     [x,y,buttons] = GetMouse(window.onScreen);
                %                     y=y./prefs.mouseScale;
                %
                %                     if ~isfinite(initRT(trialIndex,i))&&abs(y-inY)>prefs.initRTthresh
                %                         initRT(trialIndex,i)= GetSecs-t1
                %                     end
                %
                %                     Screen('DrawText', window.onScreen, 'Predict the color of this square on the next trial.', 100, 100, 255);
                %                     Screen('FrameRect', window.onScreen,255,[100 150 700 450])
                %                 end
                %                 Screen('DrawText', window.onScreen, 'Predict the color of this square on the next trial.', 100, 100, 255);
                %                 Screen('FrameRect', window.onScreen,255,[100 150 700 450])
                %                 dontClear=1;
                %                 Screen('Flip', window.onScreen,dontClear);
                for i=1:prefs.nItems
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%         Predict         %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if prefs.doET && i==1
                        io64(ioObject,LPT1address,prefs.triggerNum.predCue1on);
                        WaitSecs(0.01)
                        io64(ioObject,LPT1address,0);
                    end
                    
                    if prefs.doET && i==2
                        io64(ioObject,LPT1address,prefs.triggerNum.predCue2on);
                        WaitSecs(0.01)
                        io64(ioObject,LPT1address,0);
                    end
                    
                    %Screen('FrameRect', window.onScreen,127,[100 100 200 200])
                    %WaitSecs(prefs.preEstPause);
                    [x,inY,buttons] = GetMouse(window.onScreen);
                    inY=inY./prefs.mouseScale;
                    
                    
                    
                    
                    while any(buttons) % if already down, wait for release
                        [x,y,buttons] = GetMouse(window.onScreen);
                        y=y./prefs.mouseScale;
                    end
                    
                    
                    startPoint=round(rand.*360);
                    %initRT(trialIndex,i)=nan;
                    t1=GetSecs;
                    predRT(trialIndex,i)=nan;
                    MaxWait = 10;  %  maximum RT of  seconds
                    StopTime = GetSecs + MaxWait;  % This computes the time at which you want to stop (max RT)
                    mousePosPred = [];

                    while ~any(buttons)  && (GetSecs<StopTime) % wait for press
                        
                        
                        [x,y,buttons] = GetMouse(window.onScreen);
                        y=y./prefs.mouseScale;
                        mousePosPred = [mousePosPred,y];
                        
                        if ~isfinite(predRT(trialIndex,i))&&abs(y-inY)>prefs.initRTthresh&&abs(y-inY)<prefs.responseRTthresh
                            predRT(trialIndex,i)= GetSecs-t1
                        end
                        
                        if estType==0     % if we're estimating color, allow color graphics to be changed
                            nValues=360;
                            
                            predictedVal=mod(y-inY+startPoint, nValues);
                            
                            
                            % DISPLAY "estimate"
                            if prefs.useCieLab
                                if isfield(prefs.monitorData, 'constantB')
                                    
                                    currentColor=round(hueTo_calRGB(predictedVal./360, prefs.monitorData))';
                                else
                                    currentColor = round(255.*hueTo_calRGB(predictedVal./360, prefs.monitorData));
                                end
                            else
                                disp('HSV colors NOT enabled!!!!') % degrees on color wheel
                            end
                            
                            
                            
                            k=targetToPredict(trialIndex,i); % set orientation to cue value
                            currentOrientation = anglesInDegrees(trialIndex+1, k);
                            
                        elseif estType==1
                            nValues=180;
                            predictedVal=mod((y./2)+startPoint, nValues)
                            
                            
                            i=targetToPredict(trialIndex);
                            currentColor=probeColor; % set color to cue value
                            currentOrientation = predictedVal;
                            
                            
                        end
                        
                        
                        % Recompute the angle every time through the loop? This seems
                        % ridiculous...
                        %  SUBMITTING PREDCITION ---------------------------------------------------
                        % show probed orientation rectangle in gray at center of screen
                        
                        % Translate, rotate, re-tranlate and then draw our square
                        k=targetToPredict(trialIndex,i)
                        posX = coor(1, k);
                        posY = coor(2, k);
                        Screen('glPushMatrix', window.onScreen);
                        Screen('glTranslate', window.onScreen, window.centerX, window.centerY);
                        %             Screen('glRotate', window.onScreen, currentOrientation, 0, 0);
                        Screen('glTranslate', window.onScreen, -window.centerX, -window.centerY);
                        drawFixation(window, window.centerX, window.centerY, prefs.fixationSize,255);
                        Screen('FillRect', window.onScreen, currentColor', CenterRectOnPoint(baseRect,posX, posY));
                        Screen('glPopMatrix', window.onScreen);
                        instructionText='Predict the color of this square on the next trial.';
                        DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
                        %Screen('DrawText', window.onScreen, 'Predict the color of this square on the next trial.', 'center', window.centerY*0.2, 255);
                        
                        Screen('FrameRect', window.onScreen,255,[window.centerX*0.3 window.centerY*0.4 window.centerX*1.7 window.centerY*1.7])
                        %Screen('FrameRect', window.onScreen,255,[200 250 1700 800])
                        Screen('Flip', window.onScreen);
                        %
                    end
                    
                    
                    
                    
                    
                    allPredStartPoint(trialIndex,i) = startPoint;
                    predictedValue(trialIndex,i) = predictedVal;
                    allMousePosPred{trialIndex,i} = mousePosPred;
                    
                    while any(buttons) % wait for releasessssss
                        [x,y,buttons] = GetMouse(window.onScreen);
                        y=y./prefs.mouseScale;
                        
                    end
                    
                    if prefs.doET && i==1
                        io64(ioObject,LPT1address,prefs.triggerNum.predResp1);
                        WaitSecs(0.01)
                        io64(ioObject,LPT1address,0);
                    end
                    
                    if prefs.doET && i==2
                        io64(ioObject,LPT1address,prefs.triggerNum.predResp2);
                        WaitSecs(0.01)
                        io64(ioObject,LPT1address,0);
                    end
                    
                    % return to fixation
                    returnToFixation_wRect(window, window.centerX, window.centerY, prefs.fixationSize, prefs);
                    
                    
                    
                    
                    if estType==0
                        trialErr=abs(circ_dist(deg2rad(predictedValue(trialIndex,i)), deg2rad(toPredictColor(trialIndex,i))));
                        if trialErr<prefs.fdbk_cThresh
                            correctPredict=1;
                        else
                            correctPredict=-1;
                        end
                    else
                        
                        trialErr=abs(circ_dist(2.*deg2rad(predictedValue(trialIndex,i)), 2.*deg2rad(anglesInDegrees(trialIndex, targetToPredict(trialIndex)))));
                        if trialErr<prefs.fdbk_oThresh
                            correctPredict=1;
                        else
                            correctPredict=-1;
                        end
                    end
                    
                    display(trialErr)
                    display(correctPredict)
                    prefs.Correct = correctPredict
                    
                    allPredErr(trialIndex,i) = trialErr;
                    allPredCorr(trialIndex,i)     = correctPredict;
                    %score for prediction
                    allScorePred(trialIndex,i)    = correctPredict.*1;
            
                    
                    %         allPosFdbk(trialIndex)  = posFdbk;
                    
                    %         allScore(trialIndex)    = posFdbk.*bet;
                    returnToFixation_wRect(window, window.centerX, window.centerY, prefs.fixationSize, prefs);
                    WaitSecs(prefs.ITI);
                end
                %Screen('Flip', window.onScreen,);
                
            end
            
        end
        
        clear ans buttons x y minDistance minDistanceIndex
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                create and save data structure                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % data:          Structure containing the following fields of data:
    
    % prefs      :    structure of initial preferences.
    % colorArray :    matrix of all colors for each trial
    % orientArray:    matrix of all target orientations for each trial
    % isRand     :    were colors and orientations generated at random?
    % chosenTarg :    number of probed target
    % probeDim   :    Dimension that was used as a probe (0=location, 1=color, 2= orientation)
    % estimDim   :    1= color, 2 = orientation.
    % reportedValue:  in degrees
    % pdw        :    true indicates a PDW was made on this trial
    % bet        :    Value of bet (1 if there was no post decision wager)
    % error      :    absolute error (in radians).
    % correct    :    Was error smaller than relavant threshold?
    % fdbkType   :    pos = +1,  neg = -1.  (Can be different than correct if pseudofeedback is used)
    % score      :    bet*fdbkType
    % presentedColor: the color presented to participant based on order of stimulus
    
    % save mean of the prediction
    % and when changepoint/oddball happen (binary)
    
    data                = struct;
    data.prefs          = prefs;
    data.colorArray     = colorsInDegrees;
    data.orientArray    = anglesInDegrees;
    data.isRand         = isRand;
    data.chosenTargEst  = targetToTest;
    
    data.presentedColor = presentedColor;
    data.presentedAngle = presentedAngle;
    data.probeDim       = allProbeType;
    data.estimDim       = allEstType;
    data.pdw            = pdw;
    data.bet            = allBet;
    data.error          = allTrialErr;
    data.correct        = allCorr;
    data.predError      = allPredErr;
    data.predCorrect    = allPredCorr;
    data.fdbkType       = allPosFdbk;
    %data.score          = allScore;
    data.isRand         = isRand;
    data.locJit         = allLocJit;
    data.reportedValue  = reportedValue;
    data.condition      = condition;
    data.mousePos       = allMousePos;
    
    
    if prefs.doPredict
        data.chosenTargPredict=targetToPredict;
        data.predictedValue = predictedValue;
        data.toPredict      = toPredictColor;
        data.generativeMean = generativeMean;
        data.surprise       = surprise;
        data.predRT         = predRT;
        data.predStartPoint = allPredStartPoint;
        data.mousePosPred   = allMousePosPred;
    end
    
    
    
    data.startPoint     = allStartPoint;
    data.validStartPoint= true(size(allStartPoint)); % just a flag added
    
    % when i fixed the startpoint computation
    %data.initRT         = initRT;
    data.estRT         = estRT;
    
    
    allScoreEst(allScoreEst==-1)= prefs.negPointScale;
    %keyboard
    if prefs.doPredict
        allScorePred(allScorePred==-1)= prefs.negPointScale;
        allScore=cat(1,allScoreEst,allScorePred);
    else
        allScore=allScoreEst;
    end
    totScore=sum(nansum(allScore));
    
    data.allScore=allScore;
    data.totScore=totScore;
    data.sumScore=sumScore+totScore;
    
    
    
    save(saveFileName, 'data')
    
    if block~=1
        %Provide feedback:
        %######################## INSTRUCTIONS
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        feedbackText = sprintf('This block you earned: \n\n\n  %s points! \n\n\n click to continue', num2str(totScore));
        DrawFormattedText(window.onScreen,feedbackText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            
            WaitSecs(.1)
        end
    end
    
    if block==3 && prefs.doPredInstructions==0
        %Provide feedback:
        %######################## INSTRUCTIONS
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        feedbackText = sprintf('Congratulations on completing this block! \n\n The way the colors change over time will be different in the next block. \n\n click to continue');
        DrawFormattedText(window.onScreen,feedbackText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            
            WaitSecs(.1)
        end
    end
    
    if block==4 
        Screen('TextSize',window.onScreen,40);
        Screen('TextFont',window.onScreen,'Times');
        Screen('TextStyle',window.onScreen,0);
        Screen('TextColor',window.onScreen, 255);
        feedbackText = sprintf('You earned: %s points for %s bonus dollars!', num2str(data.sumScore),num2str(ceil(data.sumScore./prefs.pointsPerDollar)));
        %feedbackText = 'The end!';
        DrawFormattedText(window.onScreen,feedbackText,'center','center')
        Screen(window.onScreen, 'Flip');
        pause(1)
        [x,y,buttons]=GetMouse;
        while ~any(buttons)
            [x,y,buttons]=GetMouse
            
            WaitSecs(.1)
        end
    
    end
    
    
    
    
    if prefs.closeWindow
        if prefs.doET
            io64(ioObject,LPT1address,prefs.triggerNum.end);
            WaitSecs(0.01)
            io64(ioObject,LPT1address,0);
        end
        postpareEnvironment;
        if prefs.doET
            Eyelink('StopRecording'); %stop eyelink
            Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
            Eyelink('Command', 'clear_screen 0'); % Clear Host PC backdrop graphics at the end of the experiment
            WaitSecs(0.5); % Allow some time before closing and transferring file
            Eyelink('CloseFile'); 
            % download data file
            
            try
                fprintf('Receiving data file ''%s''\n', subStr );
                status=Eyelink('ReceiveFile'); %receive file from eyelink
                if status > 0
                    fprintf('ReceiveFile status %d\n', status);
                end
                if 2==exist(subStr, 'file')
                    fprintf('Data file ''%s'' can be found here: ''%s''\n', subStr, edfFile );
                end
            catch rdf
                fprintf('Problem receiving data file ''%s''\n', subStr );
                rdf;
            end
            [el] = EyelinkSetup(0,window.onScreen) %turn off eyelink
        end
    else
        prefs.window=window;
    end
    
    
catch
    postpareEnvironment;
    psychrethrow(psychlasterror);
end

function prepareEnvironment(prefs)
commandwindow;
a=clock;
rand('twister', a(6).*10000);
ListenChar(2);
Screen('Preference','SkipSyncTests', 0)





function postpareEnvironment

loadLinearGammaTable;
ShowCursor;
ListenChar(0);
Screen('CloseAll');

function instruct(window)
Screen('TextSize',window.onScreen,40);
Screen('TextFont',window.onScreen,'Times');
Screen('TextStyle',window.onScreen,0);
Screen('TextColor',window.onScreen, 255);
%Screen('TextSize', window.onScreen, window.fontsize);
instructionText='When you are ready, left click the mouse to begin.';
DrawFormattedText(window.onScreen,instructionText,'center','center');
%DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
%Screen('DrawText', window.onScreen, 'Take a break if you need one and then click to begin.', window.centerX*0.2, window.centerY*0.2, 255);
Screen('Flip', window.onScreen);
[clicks,x,y,whichButton] = GetClicks(window.onScreen);

function instructBetween(window)
Screen('TextSize',window.onScreen,40);
Screen('TextFont',window.onScreen,'Times');
Screen('TextStyle',window.onScreen,0);
Screen('TextColor',window.onScreen, 255);
%Screen('TextSize', window.onScreen, window.fontsize);
instructionText='Please take a break and close your eyes for a minute,\n\n then when you are ready, left click the mouse to begin.';
DrawFormattedText(window.onScreen,instructionText,'center','center');
%WaitSecs(30)
%DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
%Screen('DrawText', window.onScreen, 'Take a break if you need one and then click to begin.', window.centerX*0.2, window.centerY*0.2, 255);
Screen('Flip', window.onScreen);
%WaitSecs(30)
[clicks,x,y,whichButton] = GetClicks(window.onScreen);

function drawFixation(window, fixationX, fixationY, fixationSize,color)
Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, color);%255);


function offsets = circularArrayOffsets(n, centerX, centerY, radius, rotation)
degreeStep = 360/n;
offsets = [sind([0:degreeStep:(360-degreeStep)] + rotation)'.* radius, cosd([0:degreeStep:(360-degreeStep)] + rotation)'.* radius];

function [coor, jit] = circularArrayRects(rect, nItems, radius, centerX, centerY)
jit=round(rand.*360./nItems);
coor = (circularArrayOffsets(nItems, centerX, centerY, radius, jit) + repmat([centerX centerY], nItems, 1))';


function returnToFixation(window, fixationX, fixationY, fixationSize, prefs)
Screen('FillRect', window.onScreen, window.bcolor);
Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
Screen('Flip', window.onScreen);


function returnToFixation_wRect(window, fixationX, fixationY, fixationSize, prefs)
% Screen('TextSize',window.onScreen,prefs.fontSize);
% Screen('TextFont',window.onScreen,'Times');
% Screen('TextStyle',window.onScreen,0);
% Screen('TextColor',window.onScreen, 255);
Screen('FillRect', window.onScreen, window.bcolor);
Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
instructionText='Predict the color of this square on the next trial.';
DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
%Screen('DrawText', window.onScreen, 'Predict the color of this square on the next trial.', 200, 100, 255);
Screen('FrameRect', window.onScreen,255,[window.centerX*0.3 window.centerY*0.4 window.centerX*1.7 window.centerY*1.7])
%Screen('FrameRect', window.onScreen,255,[200 250 1700 800])
Screen('Flip', window.onScreen);

function returnToFixation_estimate(window, fixationX, fixationY, fixationSize, prefs)

Screen('FillRect', window.onScreen, window.bcolor);
Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
instructionText='Reproduce the color of this square on this trial.';
DrawFormattedText(window.onScreen,instructionText,'center',window.screenY*0.15);
%Screen('DrawText', window.onScreen, 'Estimate the color of this square on this trial.', 200, 100, 255);
Screen('Flip', window.onScreen);


function window = openWindow(prefs)
% get screen pointer:


window.screenNumber = max(Screen('Screens'));

% get relevant colors:
if prefs.useCieLab
    window.black=round(LabTo_calRGB([20, 0 0], prefs.monitorData));
    window.white=round(LabTo_calRGB([80, 0,0], prefs.monitorData));
    window.gray =round(LabTo_calRGB([50, 0, 0], prefs.monitorData));
end



% if we're in debug mode, open a small screen.
if prefs.debug
    screenRect =  [0,0,800,600]
else
    screenRect  = []; % screen rect
end


if isfield(window, 'gray')
    bg=window.gray;
else
    bg=[128 128 128];
end


if prefs.skipSyncTests
    Screen('Preference', 'SkipSyncTests', 1);
    disp('Careful: skipping sync tests!!!!')
end




% open screen:
[window.onScreen, rect] = Screen('OpenWindow', window.screenNumber, bg, screenRect,[],[],[]);
window.screenRect=rect;

if ~prefs.debug
    HideCursor(window.onScreen);
end



Screen('BlendFunction', window.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

if ~prefs.useCieLab
    window.black    = BlackIndex(window.onScreen);
    window.white    = WhiteIndex(window.onScreen);
    window.gray     = mean([window.black window.white]);
end
% get background color
window.fontsize = prefs.fontSize;
window.bcolor   = window.gray;


% Check some locations on the current screen:
[window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
window.centerX = window.screenX * 0.5; % center of screen in X direction
window.centerY = window.screenY * 0.5; % center of screen in Y direction
window.centerXL = floor(mean([0 window.centerX])); % center of left half of screen in X direction
window.centerXR = floor(mean([window.centerX window.screenX])); % center of right half of screen in X direction
%% if we have a gamma table, lets load it:
if isfield(prefs, 'monitorData')
    if isfield(prefs.monitorData, 'gammaTable') % if there is a gamma table load it
        Screen('LoadNormalizedGammaTable', window.onScreen, prefs.monitorData.gammaTable);
    end
end


function drawColorWheel(window, prefs, cwRotation)

colorWheelLocations = [cosd([1:360]+cwRotation).*prefs.colorWheelRadius + window.centerX; sind([1:360]+cwRotation).*prefs.colorWheelRadius + window.centerY];
colorWheelSizes = round(random('Uniform', 20, 20, 1, 360));
if prefs.useCieLab
    hue=[1:360]/360;
    RGB=hueTo_calRGB(hue, prefs.monitorData);
    Screen('DrawDots', window.onScreen, colorWheelLocations, colorWheelSizes, round(RGB.*255), [], 1);
else
    colorWheelColors = [[1:360]/360; ones(1, 360); ones(1, 360)];
    Screen('DrawDots', window.onScreen, colorWheelLocations, colorWheelSizes, round(hsv2rgb(colorWheelColors')*255)', [], 1);
end




