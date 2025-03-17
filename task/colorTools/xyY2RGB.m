function [R, G, B, GunPercentR, GunPercentG, GunPercentB] = xyY2RGB (Oldx, Oldy, Y, XR, YR, XG, YG, XB, YB, constantR, constantG, constantB, slopeR,slopeG, slopeB)


%Oldx, Oldy, and Y are the CIE xyY coordinates, respectively.
%XR and YR are the cie x and y coordinates for the red gun
%XG and YG are the cie x and y coordinates for the green gun
%XB and YB are the cie x and y coordinates for the blue gun
%constantR, constantG, constantB, are constants from the linear equations
%that fit the functions of log voltage over log luminance for red, green
%and blue guns respectively
%slopeR, slopeG, and slopeB are the slopes of the fuctions described above.


%This code was adapted from a program written by Michael Webster,
%University of Nevada, Reno

step1 = (Oldx - XB) .* (YR - YB) - (Oldy - YB) .* (XR - XB);
step2 = (XG - XB) .* (YR - YB) - (YG - YB) .* (XR - XB);     
step3 = step1 ./ step2;
step4 = ((Oldx - XB) - step3 .* (XG - XB)) ./ (XR - XB);
step5 = 1.0 - step3 - step4;
step6 = step4 .* YR + step3 .* YG + step5 .* YB;


% OK... there is going to be some total amount of gun blasting... this
% solves for the proportiality across the guns? It seems that this 
GunPercentR = (step4 .* YR) ./ step6;
GunPercentG = step3 .* YG ./ step6;
GunPercentB = 1.0 - GunPercentR - GunPercentG;

% OK. this is pretty straight forward. 

R = round(10^constantR.*(GunPercentR .* Y).^slopeR);

G = round(10^constantG.*(GunPercentG .* Y).^slopeG);

B = round(10^constantB.*(GunPercentB .* Y).^slopeB);


