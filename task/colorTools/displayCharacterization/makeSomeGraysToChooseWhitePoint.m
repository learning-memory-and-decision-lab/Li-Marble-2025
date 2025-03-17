%% This code is for trying a few  possible neutral gray tones and chooosing one 
%  to use as a white point. 


probeVals=[0:17:255]; % We'll show stimuli at these values
nReps    = 1;               % We'll repeat each stimulus this many times


commandwindow;
window.screenNumber = max(Screen('Screens'));
[window.onScreen rect] = Screen('OpenWindow', window.screenNumber, [0 0 0],[],[],[],[]);
[window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
window.screenRect  = [0 0 window.screenX window.screenY]; % screen rect
window.centerX = window.screenX * 0.5; % center of screen in X direction
window.centerY = window.screenY * 0.5; % center of screen in Y direction
window.centerXL = floor(mean([0 window.centerX])); % center of left half of screen in X direction
window.centerXR = floor(mean([window.centerX window.screenX])); % center of right half of screen in X direction


nProbes  =length(probeVals);
rVals=repmat([1 0 0], nProbes,1).*repmat(probeVals', 1, 3);
gVals=repmat([0 1 0], nProbes,1).*repmat(probeVals', 1, 3);
bVals=repmat([0 0 1], nProbes,1).*repmat(probeVals', 1, 3);
allStim=[repmat(rVals, nReps, 1); repmat(gVals, nReps, 1); repmat(bVals, nReps, 1)]
randOrder=1:length(allStim);
reading=nan(size(randOrder));

% load identity color lookup table.
idClut=LoadIdentityClut(window.onScreen)
Screen('LoadNormalizedGammaTable', window.onScreen, idClut);

% create a box in the cente of the screen
boxRect = [0 0 400 400];
boxRect = CenterRect(boxRect,rect);



allGrays=[110, 100, 80; 100, 100, 100; 100, 80 , 100; 100, 80, 80]


% run through each test triplet... show a box of that color, take a
% reading, record it.
for i=1:length(allGrays)
    
    boxColor = allGrays(i,:)
   
    Screen('FillRect', window.onScreen, boxColor, boxRect);
    Screen('Flip', window.onScreen, 0, 1);
    
    
    %######################## INSTRUCTIONS
  %  Screen('TextSize',window.onScreen,30);
  %  Screen('TextFont',window.onScreen,'Times');
  %  Screen('TextStyle',window.onScreen,0);
  %  Screen('TextColor',window.onScreen, [255 255 255]);
  %  instructionText = num2str(allStim(i));
  %  DrawFormattedText(window.onScreen,instructionText,'center','center')
  %  Screen(window.onScreen, 'Flip');

    
    reading=nan;
    while ~isfinite(reading)
        reading = input('Enter photometer reading (cd/m2): ');
    end

end
Screen('CloseAll');

