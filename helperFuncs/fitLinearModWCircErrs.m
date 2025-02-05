function [params, negLogLike]=fitLinearModWCircErrs(data)

%% GOAL: Fit a regression-like model to data on a circle (Updates) with any
% set of predictors.

% data.Y              = ydata
% data.X              = xdata
% data.includeUniform = do you want to include a uniform mixture component?
% data.whichParams    = which parameters should we fit?
% data.startPoint     = where should we start parameter search
% data.lb             = lower bound
% data.ub             = upper bound
% data.priorWidth     = mean of gaussian parameter priors (should have one for each coefficient [NOT PRECISION OR MIXTURE parameters])
% data.priorWidth     = width of gaussian priors (same as above).
% data.nStart         = number of start points to use for optimizer


% Parameters =
% 1) concentration
% 2:(size(data.X, 2)+1)= coefficients
% end) if you include uniform -- then you get another parameter which is
%                                  the uniform probability

% vm concentration:


if ~isfield(data, 'nStart')
    nStart=10;
else
    nStart=data.nStart;
end


if ~isfield(data, 'startPoint')
    startPoint=[5, zeros(1, size(data.X, 2))];
    makeStartPoint=true;
else
    makeStartPoint=false;
    startPoint=data.startPoint;
end

if ~isfield(data, 'whichParams')
    data.whichParams=true(size(startPoint));
end


if isfield(data, 'lb')
    lb=data.lb;
    ub=data.ub;
else
    lb=[.0001, ones(1, size(data.X, 2)).*-5];
    ub=[100, ones(1, size(data.X, 2)).*5];
end


% Set startpoint and boundaries on uniform if you include it in mix.
if isfield(data, 'includeUniform')&data.includeUniform==1
    if makeStartPoint
        startPoint(end+1)=.5;
    end
    lb(end+1)=0;
    ub(end+1)=1;
elseif ~isfield(data, 'includeUniform')
    data.includeUniform=false;
end

% Setup optimization:
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon','objective',...
    @circLinPredMod,'x0',startPoint(data.whichParams),'lb',lb(data.whichParams),'ub',ub(data.whichParams),'options',opts);
ms=MultiStart;
[params,f]=run(ms,problem,nStart);
negLogLike=circLinPredMod(params);

% function to optimize:
    function [negLogLike]=circLinPredMod(params)
        % Get "params" from optimizer for parameters that you are fitting: 
        gParams(data.whichParams)=params(data.whichParams);
        gParams(~data.whichParams)=startPoint(~data.whichParams); % otherwise, use the "start point"
        
        % parameters 2:end are coefficients... unless you are including a
        % uniform mixture -- in which case they are 2:end-1:
        coeffs=gParams(2:end-data.includeUniform);
        %multiply coefficients by predictors to get yHats (ie. predictions for outcome)
        yHat=data.X*coeffs';
        
        
        
        if max(abs(yHat))> pi
            %keyboard
            %this is sometimes a useful keyboard to use for debugging local
            %minimum issues, but not currently doing anything. 
        end
        
        
        % get likelihood of data under von mises distribution with mean
        % yHat:
        allVM_likelihoods=circ_vmpdf(data.Y, yHat, gParams(1));
        
        % get overall likelihood which includes  uniform mixture component.
        if data.includeUniform==1
            allVM_likelihoods=allVM_likelihoods.*(1-gParams(end))+ (1./(2.*pi)).*(gParams(end));
        end
        
        % sum all trial likelihoods and multiply by negative 1 to get
        % negative Log Likelihood:
        negLogLike=-1.*sum(log(allVM_likelihoods));

        % if we are using a (normal) prior to regularize coefficients,
        % incorporate it here:
        if isfield(data, 'priorWidth')
            logPriorProb=sum(log(normpdf(coeffs(:), data.priorMean(:), data.priorWidth(:)))); % get log prior probability. 
            if logPriorProb <-1e300  % if log likelihood goes to negative infinity, replace it with something arbitrarily small. 
                logPriorProb=-1e300; % set some minimum prior probability (flat outside some range)
            end
            negLogLike=negLogLike-logPriorProb;
        end
        
        if ~isfinite(negLogLike)
            % keyboard
        end
        
    end
end