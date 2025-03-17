 function data = modelCP(vars)
%% Perceptual inference in a dynamic changing environment (changepoints)

% Notes : 

% We, as experimenters, do not have access to Y (subjects internal
% representation). We do have access to X -- as this is controlled in our
% task. So outside of this loop, we will create Y values by sampling normal
% random variables centered on X. Then we will repeat this procedure
% multiple times in order to get "average" behavior of the model (averaged
% across different exact patterns of internal noise). 

% overarching goal: create a generative environment to model perception in an
% environment that includes an unstable generative mean (occasional
% changepoints), variable sampling of actual stimuli from that
% distribution (with some stimulus variability around the true mean
% sigma_s), and an imperfect sensory representation of that stimulus (with
% perceptual variability equal to sigma_obs). 

% We're basically modeling the paradigm of krishnamurthy 2017. 

% Generative model:

% SAMPLE CHANGEPOINT:
% S_t    ~B(H)

% SAMPLE MU:
% if S_t = 0:
% mu_t   =mu_{t-1}
% elseif S_t = 1
% mu_t   ~ U(0-360)

% SAMPLE X:
% X_t    = N(mu_t, sigma_s)

% SAMPLE Y:
% Y_t    = N(X_t, sigma_obs)



%% STEP 1: Generate observations (Ys) from the process above.

%close all;
%This variable is to indicate whether the data will be obtained from an
%external variable called data or not. If it is set to true then we will
%just replace the values of Mu, X and Y with their respective variables in
%data. If false then we just do not.
replace_data = vars.replace;
replace_data = 1;
nTrials   = 2*vars.nTrials;   % number of trials
haz         = vars.H;   % Hazard rate
sigma_x   = vars.sigma_x;   % variance on distribution of stimuli/outcomes around mean
sigma_obs = vars.sigma_y;
X= ones(nTrials,1).*-1;
nDegrees = vars.nDegrees;
circSwitch =true;  % if true we use von mises distributions instead of gaussians

Mu = nan(nTrials,1);
S = nan(nTrials,1);
Y = nan(nTrials,1);


if replace_data == false
    
    for t = 1:nTrials
    
        if data.condition(1,1)==1
            Mu = [data.generativeMean(data.block==3,1) ; data.generativeMean(data.block==3,2)] ;                              %of data to use imported user data
            X = [data.colorArray(data.block==3,1) ;data.colorArray(data.block==3,2)]; %.* (pi/360);
        else
            Mu = [data.generativeMean(data.block==4,1) ; data.generativeMean(data.block==4,2)] ;                              %of data to use imported user data
            X = [data.colorArray(data.block==4,1) ;data.colorArray(data.block==4,2)]; %.* (pi/360);
        end

        % SAMPLE CHANGEPOINT:
        % S_t    ~B(H)
        
%         S(t)=rand<haz;
%     
%         % SAMPLE MU:
%         % if S_t = 0:
%         % mu_t   =mu_{t-1}
%         % elseif S_t = 1
%         % mu_t   ~ U(0-360)
%         
%         if S(t)==1 || t ==1
%             Mu(t)=rand*360; 
%         else
%             Mu(t)=Mu(t-1);
%         end
%         
%         % SAMPLE X:
%         % X_t    = N(mu_t, sigma_s
%         X(t) = normrnd(Mu(t), sigma_x);
%         while X(t)<0 || X(t)>360
%             X(t)=normrnd(Mu(t), sigma_x);
%         end
    
        % SAMPLE Y:
        % Y_t    = N(X_t, sigma_obs)
        
        Y(t)=normrnd(X(t), sigma_obs);
    

    end
    
else 

    S = zeros(1, 200); %Need to also save S in the colorOrient_estTask.m
    %zeros(1, 200); %Need to also save Mu in the colorOrient_estTask.m
    data = load(vars.path);%importdata('data.mat'); %Replace 'data.mat' with the file path 
    data = data.alldata;
    if data.condition(1,1)==1
        Mu = [data.generativeMean(data.block==3,1) ; data.generativeMean(data.block==3,2)] ;                              %of data to use imported user data
        X = [data.colorArray(data.block==3,1) ;data.colorArray(data.block==3,2)]; %.* (pi/360);
    else
        Mu = [data.generativeMean(data.block==4,1) ; data.generativeMean(data.block==4,2)] ;                              %of data to use imported user data
        X = [data.colorArray(data.block==4,1) ;data.colorArray(data.block==4,2)]; %.* (pi/360);
    end
    nTrials = size(X, 1);
    Y = zeros(nTrials, 1);   
    for i = 1:nTrials
        Y(i)=normrnd(X(i, 1), sigma_obs);
    end
end




%% Step 2: Do inference on generated Ys (eg. figure out mu and X on each trial)

% Now we get to see Y and need to infer X and Mu:


possMu = 1:360;
possX  = 1:360;

pMu=ones(length(possMu), 1)./length(possMu); %---> A bidimensinal array of size length(possMu)*1
                                             %     containing only the
                                             %     probability of getting a
                                             %     certain X from a uniform
                                             %     distribution containing
                                             %     360 X's.
transFunc=eye(length(possMu)); %---> Transition function: Identity square 
                               %     matrix whose size equals 
                               %     the number of possible mu's
transFunc(:,:,2)= repmat(pMu', length(possMu), 1);
                               %---> This is A three dimensional matrix
                               %whose 3rd dimension is just the value of
                               %getting a certain X from a uniform
                               %distribution with 100 possible X's

transFuncByS=cat(3, transFunc(:,:,1).*(1-haz), transFunc(:,:,2).*haz);


likelihoodX = []; %--> (1) = poss X
                 %--> (2) = poss Mu
 for mu  = possMu
        if circSwitch
        pXgivMu = circ_vmpdf(deg2rad(possX'), deg2rad(mu),   1./(deg2rad(sigma_x).^2))';      % 
        else
        pXgivMu = normpdf(possX, mu, sigma_x); %horizontal vector
        end
        pXgivMu = pXgivMu./sum(pXgivMu);
        likelihoodX = [likelihoodX; pXgivMu];
 end
 shannon = zeros(nTrials, 1); %--> Contains the shanon information content of each percieved outcome
             %    this is just how surprising each of the outcomes is.
 entropy = zeros(nTrials, 1);% -->  Contains the entropy of every trial.
 relative_entropy_Mu = zeros(nTrials, 1); %--> Contains the relative entropy between the 
                       %    probability distribution at time t and the one 
                       %    at time t + 1. This will measure how much the
                       %    model is updating from trial to trial
posterior = nan(nTrials, nDegrees);
posteriorPercept = nan(nTrials, nDegrees);

for t = 1:nTrials
    
    % Storing the entropy of trial t
    entropy(t) = sum(pMu .* log2(1 ./ pMu));

    % Step 1: multiply prior (t-1) by transition probability to get prior for (t):
    % Here we multiply P(Mu_t, S|Mu_t-1) .* P(Mu_t-1) to get P(Mu_t, S, Mu_t-1)
    pMu_t=transFuncByS.*repmat(pMu',length(pMu), 1,  2);
    
    
    % Here we will calculate the [E[P(Mu_t), S]
%     pMu_t = sum(pMu_t, 2);
    pMu_t = sum(pMu_t, 2);
    pMu_t=squeeze(pMu_t);    %---> (2) = S
                            %     (1) = MU_T
                            %      Each value is [E[P(Mu_t, S)]
    
    if circSwitch % do inference on curcular distribution:   
    % how likely would true color values have been to give rise to your
    % current color percept? 
    % possX = true color values that are possible
    % Y(t)  = what you saw on this trial
    likelihoodY = circ_vmpdf(deg2rad(Y(t)), deg2rad(possX), 1./(deg2rad(sigma_obs)^2)); %horizontal vector
    likelihoodY = likelihoodY./sum(likelihoodY);
    else
    likelihoodY = normpdf(Y(t), possX, sigma_obs); %horizontal vector
    likelihoodY = likelihoodY./sum(likelihoodY);
    end

    pXandYgivMu = likelihoodX .* likelihoodY;%--> (1) X
                                             %      each value is P(X, Y|Mu)                                          
    pXandYandMuandS(:,:,1) = pXandYgivMu .* pMu_t(:,1);
    pXandYandMuandS(:,:,2) = pXandYgivMu .* pMu_t(:,2) ;%P(X, y, Mu, S)
                                                       %(1): X
                                                       %(2): Mu
                                                       %(3): S
  
    surprise(t) = sum(sum(pXandYandMuandS(:,:,2),1),2)./  sum(pXandYandMuandS(:));
    
    pYandXandMu = sum(pXandYandMuandS, 3); % integrate over S... 
   
    shannon(t) = log2(1/sum(pYandXandMu(:)));

    
    % P(X,Mu|y):
    pXandMugivY=pYandXandMu./sum(pYandXandMu(:)); % Normalize to make probability distribution (this is like dividing by p(y))
    
  
    % Marginalize over X
    pMugivY = sum(pXandMugivY, 2); % integrate over X
    pXgivY  = sum(pXandMugivY, 1); % integrate over mu
    
    posterior(t,:) = pMugivY;
    posteriorPercept(t,:)=pXgivY;
    
    % Storing the sharon content of the perceived 
    % pXgivTrials = likelihoodX .* pMu; % p(X|Mu) * p(Mu|t, t-1...1) = p(X, Mu|t, t-1...1)
    % pXgivTrials = sum(pXgivTrials, 1); % Integrate over Mu = p(X|t, t-1...1)
                                       % --> 1 row
                                       % --> columns (poss X)
    
    % likelihood_Y = zeros(length(possX), length(possX)); % p(Y|X)
                                               %--> rows(poss Y)
                                               %--> columns (poss X)
    % for y = possX
    %     % MRN changed this to try to get a circular likelihood function:
    %     if circSwitch
    %     likelihood_Y(y,:) = circ_vmpdf(deg2rad(y), deg2rad(possX)',   1./(deg2rad(sigma_obs).^2));      % 
    %     else
    %     likelihood_Y(y,:) = normpdf(y, possX, sigma_obs);
    %     end
    % 
    %     likelihood_Y(y,:) = likelihood_Y(y,:) ./ sum(likelihood_Y(y,:));
    % end
    % %p(Y, X | Trials)
    % %--> rows(poss Y)
    % %--> columns(poss X)
    % pYandXgivTrials = likelihood_Y .* pXgivTrials;
    % 
    % %p(Y | Trials)
    % % one column
    % pYgivTrials = sum(pYandXgivTrials, 2);
    % 
    
        
    
    % Calculating the relative entropy between trial t-1 and trial t for X
    if t == 1
        temp = repmat((1/length(possMu)), length(possMu), 1);
        temp = log(temp' ./ posterior(t,:));
        relative_entropy_Mu(t,:) = sum((1/length(possMu)) .* temp);
    else
        temp =  log(posterior(t-1,:) ./  posterior(t,:));
        relative_entropy_Mu(t,:) =  sum(posterior(t-1,:) .* temp);
    end
%     if t == 1
%         temp = repmat((1/length(possMu)), length(possMu), 1);
%         temp = log(posterior(t,:) ./ temp');
%         relative_entropy_Mu(t,:) = sum(posterior(t,:) .* temp);
%     else
%         temp =  log(posterior(t,:) ./  posterior(t-1,:));
%         relative_entropy_Mu(t,:) =  sum(posterior(t,:) .* temp);
%     end
    
    pMu = pMugivY;
    
    
    
        
    
end 





%% 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % ANALYSIS OF DATA 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear maxLikePostMu maxLikePostX
maxLikePostX = nan(1, nTrials);
maxLikePostMu = nan(1, nTrials);
expectationX = nan(1, nTrials);
expectationMu = nan(1, nTrials);
for i = 1:size(posteriorPercept, 1)
   maxLikePostX(i)=possX(find(posteriorPercept(i,:)==max(posteriorPercept(i,:)), 1));
   maxLikePostMu(i)=possMu(find(posterior(i,:)==max(posterior(i,:)), 1));
   
   
   expectationX(i)=posteriorPercept(i,:)*possX';
   expectationMu(i) = posterior(i,:) * possX';
   
end


% Lets compare performance of our model to an unbiased max likelihood
% estimate (AKA: just say y). 

mean((Y'-X).^2);            % just reporting sensory experience directly
mean((maxLikePostX'-X).^2); % our (max a-posteriori -- has error due to grid approximation)
mean((expectationX'-X).^2); % our model (expectation over X)



% YAY!!! bayesian perception is better than max like perception!




%% STEP 3: MEASURING LEARNING AND PERCEPTUAL BIAS IN THIS MODEL

% One measure of how fast the model is "learning" would be how much mu is
% updated in response to a prediction error. 





maxLikePostMu= [mean(possMu), maxLikePostMu]; % Note = we are now indexing as if this is at the beginning of the trial!
muUpdate=circ_dist(deg2rad(maxLikePostMu(2:end)), deg2rad(maxLikePostMu(1:end-1)));
% [~, muUpdate, ~] = difference(deg2rad(X), deg2rad(maxLikePostMu(1:end-1)), 'polar');
% muUpdate = muUpdate';
predError=circ_dist(deg2rad(maxLikePostX), deg2rad(maxLikePostMu(1:end-1)));%subjective error

objPredError=circ_dist(deg2rad(X)', deg2rad(maxLikePostMu(1:end-1)));%objective error


% perceptual bias? 
% figure
% hold on
predictionErrorOnX=circ_dist(deg2rad(X(:,1)),deg2rad(maxLikePostMu(1:end-1)'));

perceptualErrorOnX=circ_dist( deg2rad(X(:,1)),deg2rad(maxLikePostX'));

noisyPrediction = normrnd(maxLikePostMu,sigma_obs);
noisyPerception = maxLikePostX;

noisyPredEOnX = circ_dist(deg2rad(X(:,1)),deg2rad(noisyPrediction(1:end-1)'));
noisySubPredEOnX = circ_dist(deg2rad(noisyPerception),deg2rad(noisyPrediction(1:end-1)))';
noisyPercEOnX = circ_dist(deg2rad(X(:,1)),deg2rad(noisyPerception'));

modelPredUpdate = circ_dist(deg2rad(maxLikePostMu(1:end-1)),deg2rad(maxLikePostMu(2:end)));
%%
%Storing data, all one dimensional data is stored in the form of horizontal
% vectors.
data.nTrialsCP = nTrials;
data.YCP = Y; 
data.X = X';
data.Mu = Mu;
data.SCP = S;
data.shannonCP = shannon';
data.klCP = relative_entropy_Mu';
data.entropyCP = entropy';
data.predictionErrorOnX = predictionErrorOnX';
data.perceptualErrorOnX = perceptualErrorOnX';
data.predErrorCP = predError;
data.objPredErrorCP=objPredError;
data.muUpdate = muUpdate;
data.posteriorPerceptCP = posteriorPercept;
data.posteriorCP = posterior;
data.expectationX = expectationX;
data.expectationMu = expectationMu;
data.maxLikePostMu = maxLikePostMu;
data.modelPredUpdateCP = modelPredUpdate;
data.maxLikePostX = maxLikePostX;
data.surpriseCP = surprise; 
data.noisyPredictionCP = noisyPrediction;
data.noisyPerceptionCP = noisyPerception;
data.noisyPredEOnX = noisyPredEOnX;
data.noisyPercEOnX = noisyPercEOnX;
data.noisySubPredEOnX = noisySubPredEOnX;

% save('../Desktop/neuro/changepoint_test_data', 'data');
data=straightStruct(data);

%%
% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % PLOTS
% % % % % % % % % % % % % % % % % % % % % % % % % 
% num        = 1;
% wid        = 17; % total width
% hts        = [8, 8 ]; % height of each row
% cols       = {1, 1}; % width of columns
% [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 4, 4, [], ''); % etc
% set(axs,'Units','normalized');
% set(gcf,'color','white')
% fig_.Position = [10 10 30 40];
% left_color = [0 0 255] / 256;
% right_color = [255 0 0] / 256;
% set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % draw in each panel, one at a time
% 
% 
% lw=1.5
% lw2=3
% exSub=12
% 
% threshold = 0.5;
% 
% 
% 
% for xx = 1:length(axs)
%     axes(axs(xx)); hold on; cla(gca)
%     if xx==1
%         small_error = abs(predError) < threshold;
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         small_b = regress(muUpdate(small_error)',[ones(size(predError(small_error)))' predError(small_error)']);
%         large_b = regress(muUpdate(~small_error)',[ones(size(predError(~small_error)))' predError(~small_error)']);
%         plot(predError(small_error) , muUpdate(small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'blue');
%         plot(predError(~small_error) , muUpdate(~small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'red');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%     elseif xx==2
%         small_error = abs(predictionErrorOnX) < threshold;
%         small_b = regress(perceptualErrorOnX(small_error),[ones(size(predictionErrorOnX(small_error))) predictionErrorOnX(small_error)]);
%         large_b = regress(perceptualErrorOnX(~small_error),[ones(size(predictionErrorOnX(~small_error))) predictionErrorOnX(~small_error)]);
%         plot(predictionErrorOnX(small_error), perceptualErrorOnX(small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'blue');
%         plot(predictionErrorOnX(~small_error), perceptualErrorOnX(~small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'red');
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%         ylabel('Perceptual Error')
%         xlabel('Prediction Error')
%     end
%     
%     
%     setPLOT_panelLabel(gca, xx);
% end
% 
%  kk=annotation('textbox')
%  set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
%  
% 
end











