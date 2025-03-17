% Color sandbox.

% most recent characterization data
% lab EEG testing room: labMonitorData_2016-03-17_315b.mat

cd('C:\Users\lncc\Dropbox\MattExperiments\vwm_task\code\calibrationInfo')
%load  labMonitorData_2016-03-17_315b.mat

addpath(genpath('../../'))

% OK, so lets mess around with showing some CieLAB colors on a
% characterized (but not gamma corrected) monitor. The information required
% will be the same information collected with Karen Schloss... you
% basically get the constant and slope fit to the log(relative voltage)
% versus log(luminance) for each color channel (ie RGB). You will also need
% the chromaticity of each channel (in x, y coordinate space). Together
% these values can be used to get relatively accurate reproduction of xyY
% colors using xyY2RGB.m ... however the exact method is still a bit of a
% mystery to me, so i'll probably spend some time digging around to better
% understand that at some point. 

% But for now the goal will be to stipulate a color in CIELAB space, which
% has fairly good perceptual variance properties, convert that color to xyY
% values, and then compute the RGB values that would allow the color to be
% reproduced on the characterized monitor. 

% NOTE: This is very different from the strategy that we've used previously
% in that it DOES NOT require gamma correcting the monitor.


% Choose some CIELAB colors: 

% Get these from CALIBRATION_FrankLab.xls
xR=.676      %monitorData.xR; % x chromaticity value for red channel
yR=.309      %monitorData.yR; % y chormaticity value for red channel
xG=.208      %monitorData.xG;
yG=.683       %monitorData.yG;
xB=.149     %monitorData.xB;
yB=.048       %monitorData.yB;

% slope and intercept from excel spreadsheet (get from graphs!)
constantR= 1.53                 %monitorData.constantR
slopeR=    .48                   %monitorData.slopeR
constantG= 1.33                   %monitorData.constantG
slopeG=    .482          %monitorData.slopeG
constantB=  1.83  %       monitorData.constantB
slopeB=     .476    %monitorData.slopeB

int=(2.*pi)./360;

% make a circle in CIELAB space
A_offset=0;
B_offset=0;
gain=55;    % If you get an error that says the circle is out of range, 
Lum =60;    % You will need to adjust the gain and/or luminance!


A=((sin(0:int:2.*pi) ).*gain) +A_offset;
B=(cos(0:int:2.*pi) ).*gain  +B_offset;
L=ones(size(A)).*Lum;
numMeas=length(L);

hold on
plot(0,0, '.r')
plot(A, B, '.')
ylabel('A*')
xlabel('B*')
close all

% project the circle in CIELAB into XYZ space.
whitePoint=[100, 100, 100]
[XYZ] = lab2xyz([L' A' B' ], 'user', whitePoint)


% Ugh. Looks like i have two separate functions for converting LAB to XYZ.
%  [XYZ_alt]=      LabToXYZ([L(1)', A(1)', B(1)']', whitePoint')

% The white point will be 0,0 in A*B* space... but don't take my word for
% it...
[lab] = xyz2lab(whitePoint, 'user', whitePoint)


% get from XYZ space to xyY notation (no transformation of color space,
% just specifying the location in a different way.
x=XYZ(:,1)./sum(XYZ, 2)
y=XYZ(:,2)./sum(XYZ, 2)
Y=XYZ(:,2);

% ok, lets get the same info but now for the white point. We'll also want
% a couple versions of white... one with luminance matched to our stimuli
% (ie the gray stimulus with no color info; stimulusLuminance) and one with
% a low luminance level that serves as a background (bgLuminance). 

stimulusLuminance=nanmean(Y);

bgLuminance=10;   % set the background luminance...
xyY_white= [whitePoint(1)./sum(whitePoint), whitePoint(2)./sum(whitePoint), whitePoint(2)];

[R, G, BB, GunPercentR, GunPercentG, GunPercentB] = ...
    xyY2RGB (x, y, Y, xR, yR, xG, yG, xB, yB, constantR, constantG, constantB, slopeR,slopeG, slopeB)

[Rgray, Ggray, Bgray] = ...
    xyY2RGB (xyY_white(1), xyY_white(2), stimulusLuminance, xR, yR, xG, yG, xB, yB, constantR, constantG, constantB, slopeR,slopeG, slopeB)


[Rw, Gw, Bw] = ...
    xyY2RGB (xyY_white(1), xyY_white(2), bgLuminance, xR, yR, xG, yG, xB, yB, constantR, constantG, constantB, slopeR,slopeG, slopeB)


if ~isreal(R) | ~isreal(G) | ~isreal(BB)
    disp('Problem: RGB values go out of bounds!!!')
    error('try making the circle smaller')
end



figure
hold on
plot([255, 255], [0 255], '--k')
plot([0 255],[255, 255],  '--k')
plot(R, BB, '.r')
plot(G, BB, '.g')
plot(Rgray,Bgray, '*r')
plot(Ggray,Bgray, '*g')
plot(Rw,Bw, '.m')
plot(Gw,Bw, '.c')
ylabel('blue')
xlabel('green/red')
%pause(5)
close all



% lets choose a background gray at the white point chromaticity but with
% low luminance. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Measure color reproduction                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

commandwindow;
window.screenNumber = max(Screen('Screens'));
[window.onScreen rect] = Screen('OpenWindow', window.screenNumber, [0 0 0],[],[],[],[]);
[window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
window.screenRect  = [0 0 window.screenX window.screenY]; % screen rect
window.centerX = window.screenX * 0.5; % center of screen in X direction
window.centerY = window.screenY * 0.5; % center of screen in Y direction
window.centerXL = floor(mean([0 window.centerX])); % center of left half of screen in X direction
window.centerXR = floor(mean([window.centerX window.screenX])); % center of right half of screen in X direction

sel=mod(1:numMeas, 36)==0;
rVals=R(sel);
gVals=G(sel);
bVals=BB(sel);


allStim=[[0 0 0];[rVals, gVals, bVals]; [Rgray, Ggray, Bgray]; [Rw, Gw, Bw]];
randOrder=1:length(allStim);
reading=nan(size(randOrder));

% load identity color lookup table.
idClut=LoadIdentityClut(window.onScreen)
Screen('LoadNormalizedGammaTable', window.onScreen, idClut);

% create a box in the cente of the screen
boxRect = [0 0 400 400];
boxRect = CenterRect(boxRect,rect);

% preallocate space for input:
Y_reading_raw=nan(length(allStim), 1);
x_reading_raw=nan(length(allStim), 1);
y_reading_raw=nan(length(allStim), 1);

% run through each test triplet... show a box of that color, take a
% reading, record it.
for i=1:length(randOrder)
    
    boxColor = allStim(randOrder(i),:);   % foreground color
    Screen('FillRect', window.onScreen, boxColor, boxRect);
    Screen('Flip', window.onScreen, 0, 1);
    
    while ~isfinite(Y_reading_raw(i))
        Y_reading_raw(i) = input('Enter photometer reading (cd/m2): ');
    end

    while ~isfinite(x_reading_raw(i))
        x_reading_raw(i) = input('Enter x reading: ');
    end

    while ~isfinite(y_reading_raw(i))
        y_reading_raw(i) = input('Enter y reading: ');
    end

    
end
Screen('CloseAll');




Y_reading=Y_reading_raw(2:end);
x_reading=x_reading_raw(2:end);
y_reading=y_reading_raw(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at intended versus actual colors in cieLAB space                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% measure coordinates in xyY space...
x=x_reading;
y=y_reading;
Y=Y_reading;

%  then convert to XYZ
X=(Y./y).*x;
Z=(Y./y).*(1-x-y);

% then convert to cieLAB
[lab] = xyz2lab([X Y Z], 'user', whitePoint)


%% OK, this is the big one. You should see a circle in blue. A few of the
%  points on the blue circle should be covered with cyan points indicating
%  where we've taken measurements. Hopefully the red points fall almost
%  directly on top of the cyan points, indicating that our measurements
%  match the colors that we intended to show.
close all
hold on
plot(A, B, '.')
plot(A(sel), B(sel), 'oc')

%plot(lab(:,2), lab(:,3), '.r')

plot(lab(1:end-2,2), lab(1:end-2,3), '.r')

plot(0,0, '*k')

plot(lab(end,2), lab(end, 3), '*m')
%plot(lab(end-1,2), lab(end-1, 3), '*g')

dateStr=date
ylabel('A*')
xlabel('B*')
saveas(gcf, sprintf('calibrationFigure_%s', dateStr));
input('hit enter to continue...')
close all




%% If the plot above looks good, save the monitor data.
monitorData.xR=xR;
monitorData.yR=yR;
monitorData.xG=xG;
monitorData.yG=yG;
monitorData.xB=xB;
monitorData.yB=yB;

monitorData.constantR=constantR;
monitorData.slopeR=slopeR;
monitorData.constantG=constantG;
monitorData.slopeG=slopeG;
monitorData.constantB=constantB;
monitorData.slopeB=slopeB;

% These are actually properties of the stimuli that we've checked... 
monitorData.scale     = gain;       % how big should the circle be (centered at A*=0, B*=0)
monitorData.whitePoint= whitePoint; % what is 0,0 in XYZ colors?
monitorData.lVal      = unique(L);  % luminance of stimuli.

Date=datestr(now, 'yyyy-mm-dd')
save(sprintf('labMonitorData_%s_315b.mat', Date), 'monitorData');


