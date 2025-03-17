function [X, negLogLike, simplexP,uniformProb, genData, genData_altTargErrors]=fit_VWM_mixtureModelTrial(data, whichParams, startPoint, prior)

% The goal of this code is to fit the classic mixture of Gaussian and
% uniform model to data from a VWM estimation task with variable set size.

% INPUTS:

% data should be a structure including the following fields:

% data.signedError  --> signed errors (or just angles if not for VWM task)
% data.signedError_allTargs --> errors computed as if subject were
    % estimating all of the colors in the array (ie colorArray - subject
    % response). 
% data.doFit --> if true minimize neg log like, otherwise return likelihood
% of startpoint.

% data.zeroMeanBE = logical... if true binding errors are assumed to be
% centered on actual non-probed target colors. This is useful if you are
% trying to compute the bias in correct reports. 
    
    
% simplex order: gaussian, binding error, uniform



% OUTPUTS:
%
% params: maximum  likelihood model parameters, 
% 1= proportion gaussian,
% 2= concentration (ie 1./sigma^2) of von mises, 
% 3= mean of gaussian (should be zero... but who knows!)
% 4= proportion binding error (ie propr of gaussian that are evenly distributed across targets) 

% startPoint: initial parameter values
% whichParasms: which parameters do you want to fit?  others will remain at
% start point.

whichParams=logical(whichParams);

if nargin<3||isempty(startPoint)
    startPoint=[.9, 1, 0, 0];
end

if nargin<4||isempty(startPoint)
    wPrior=false;
else
    wPrior=true;
end

if ~(isfield(data, 'doFit'))
    data.doFit=true;
end

if ~(isfield(data, 'zeroMeanBE'))
    data.zeroMeanBE=false;
end



lb=[0 0 -pi, 0];
ub=[1 1000 pi, 1];

%X=fminsearch(@mixModel,startPoint(whichParams));
options = optimset('Algorithm','interior-point', 'MaxPCGIter', 5000, 'MaxProjCGIter', 5000, 'MaxSQPIter', 5000);

if data.doFit
[X] = fmincon(@mixModel, startPoint(whichParams), [], [], [], [], lb(whichParams), ub(whichParams), [], options);
else
X=startPoint(whichParams);
end
[negLogLike,uniformProb]=mixModel(X);

if nargout>3
% get categorical simplex:
uParam(whichParams)=X;
uParam(~whichParams)=startPoint(~whichParams);
% simplex order: gaussian, binding error, uniform
simplexP=[uParam(1)-uParam(1).*uParam(4), uParam(1).*uParam(4), 1-uParam(1)];
end


if nargout>4   
    %keyboard
    altTargLocs=circ_dist(repmat(data.signedError, 1, size(data.signedError_allTargs, 2)), data.signedError_allTargs);
    [genData, genData_altTargErrors]=genMixData(altTargLocs, [simplexP, uParam(2)]);
  %  keyboard;

end




    function [negLogLike,uniformProb]=mixModel(X)
        uParam(whichParams)=X;
        uParam(~whichParams)=startPoint(~whichParams);
        % compute probability on a von mises distribution for a given mean and
        % dispersion:
        
       % pVM1= circ_vmpdf(data.signedError,  uParam(3), uParam(2));

        pVM= circ_vmpdf_ov(data.signedError,  uParam(3), uParam(2));
       
        
        
        % compute probability on von mises distributions conditional on a
        % binding error:
        
        if data.zeroMeanBE==true
            be_mean=0;
        else
            be_mean=uParam(3);
        end
        
        
        
        % only compute binding error probability if there is the potential
        % for a binding error.
        if size( data.signedError_allTargs, 2)>0;
            pVM_BE=nan(size(data.signedError_allTargs));
            
            for i = 1:size(data.signedError_allTargs,2)
                pVM_BE(:,i)=circ_vmpdf_ov(data.signedError_allTargs(:,i), be_mean, uParam(2));
            end
            pVM_BE=nanmean(pVM_BE, 2);
            
            
            % weighted combination of gaussian likelihoods:
            pVM= pVM.*(1-uParam(4)) + pVM_BE.*uParam(4);
        else
            pVM= pVM.*(1-uParam(4));
        end
        
            
        
        if any(isinf(pVM))
            pVM(isinf(pVM))=10e10;
            disp('numeric overflow flag')
        end
        
        
        
        pU =    ones(size(pVM))./(2.*pi); % compute probability on a uniform        
        totProb= pVM.*uParam(1) + pU.*(1-uParam(1));
        
        uniformProb=(pU.*(1-uParam(1)))./totProb;
        negLogLike= -1.*nansum(log(totProb));
        
        if wPrior
            if ~isstruct(prior)
                % don't do anything?
            else
                % prior is a structure with fields:
                % .distr         = distribution for pdf
                % .p1            = parameter value 1
                % .p2            = parameter value 2
                % .whichParamNum = which parameter value should come from
                % this prior distribution?
                
                logPrior=0;
                for i = 1:length(prior.whichParamNum)
                    logPrior= logPrior + pdf(prior.distr{i},X(prior.whichParamNum(i)), prior.p1(i), prior.p2(i));
                end
                negLogLike=negLogLike-logPrior;
            end
        end 
        
        % Penalize for probabilities outside of the range [0 1];
            
        if ~isfinite(negLogLike)
          disp('got inf')
        end
    end
end






