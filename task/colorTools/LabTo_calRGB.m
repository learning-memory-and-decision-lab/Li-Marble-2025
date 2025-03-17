function RGB=LabTo_calRGB(Lab, monitorData);


if isfield(monitorData, 'testParams')
    
    % convert CIELAB colors to CIE XYZ colors:
    XYZ = LabToXYZ(Lab(:,:)',monitorData.testParams.whitePoint');
    
    % convert XYZ colors to RGB colors:
    RGB=nan(size(XYZ));
    for i = 1:size(XYZ, 2)
        RGB(:,i)=inv(monitorData.normMatrix)*XYZ(:,i);
    end
    
else
    XYZ = lab2xyz(Lab(:,:), 'user', monitorData.whitePoint)
    
    x=XYZ(:,1)./sum(XYZ, 2);
    y=XYZ(:,2)./sum(XYZ, 2);
    Y=XYZ(:,2);
    
    [R, G, B] = ...
        xyY2RGB (x, y, Y, ...  % required x,y,Y values.
        monitorData.xR, monitorData.yR, monitorData.xG, monitorData.yG, monitorData.xB, monitorData.yB, ...%chromaticity measurements
        monitorData.constantR, monitorData.constantG, monitorData.constantB, ... %constant from log log fits
        monitorData.slopeR, monitorData.slopeG, monitorData.slopeB)  %slope from log log fits.
    RGB=[R,G,B];
    
end


    
    
 