function RGB=hueTo_calRGB(hue, monitorData)
% Create circular CIELAB colors based on 0-1 "hue" value.
if isfield(monitorData, 'testParams')
    aStar=sin(hue.*pi.*2).*monitorData.testParams.scale ;
    bStar=cos(hue.*pi.*2).*monitorData.testParams.scale ;
    lStar=ones(size(hue)).*monitorData.testParams.lVal;

    Lab=[lStar' aStar' bStar'];
    
    % convert CIELAB colors to CIE XYZ colors:
    XYZ = LabToXYZ(Lab(:,:)',monitorData.testParams.whitePoint');
    
    % OTHERWISE do it the old way:
    % convert XYZ colors to RGB colors:
    RGB=nan(size(XYZ));
    for i = 1:size(XYZ, 2)
        RGB(:,i)=inv(monitorData.normMatrix)*XYZ(:,i);
    end
    
    % IF monitorData contains the slope and constants from the log/log fits
    % then use Karen Schloss' method for getting RGB values.
    
elseif isfield(monitorData, 'xR');
    aStar=sin(hue.*pi.*2).*monitorData.scale ;
    bStar=cos(hue.*pi.*2).*monitorData.scale ;
    lStar=ones(size(hue)).*monitorData.lVal;

    Lab=[lStar' aStar' bStar'];
    
    % convert CIELAB colors to CIE XYZ colors:
    XYZ = LabToXYZ(Lab(:,:)',monitorData.whitePoint')';

    
    x=XYZ(:,1)./sum(XYZ, 2);
    y=XYZ(:,2)./sum(XYZ, 2);
    Y=XYZ(:,2);

    [R, G, B] = ...
        xyY2RGB (x, y, Y, ...  % required x,y,Y values.
        monitorData.xR, monitorData.yR, monitorData.xG, monitorData.yG, monitorData.xB, monitorData.yB, ...%chromaticity measurements
        monitorData.constantR, monitorData.constantG, monitorData.constantB, ... %constant from log log fits
        monitorData.slopeR, monitorData.slopeG, monitorData.slopeB);  %slope from log log fits.
    RGB=[R,G,B];

end