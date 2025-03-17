

% matlab way (uses projects into real/imaginary space):
x=deg2rad(359); % x (in radians)
y=deg2rad(1);   % y (in radians)
r = rad2deg(angle(exp(1i*x)./exp(1i*y)));

% alternate way: arc distance = (deltaX^2 + deltaY^2)^(1/2):
x=359
y=1;
cosDif=cosd(x)-cosd(y);
sinDif=sind(x)-sind(y);
arcDist= rad2deg(sqrt(cosDif^2 + sinDif^2))

% should be able to set a limit on the arc distance that should get around
% the 360 issue. 
