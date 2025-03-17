% loadLinearGammaTable:

linGamma=[0:255; 0:255; 0:255]'./255;
window.screenNumber = max(Screen('Screens'));
[window.onScreen rect] = Screen('OpenWindow', window.screenNumber, [128 128 128],[],[],[],[]);
Screen('LoadNormalizedGammaTable', window.onScreen, linGamma);
Screen('CloseAll');
