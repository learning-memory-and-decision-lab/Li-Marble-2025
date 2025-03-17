%STEP 1: GENERATE DATA
% Generative process (Again within the range of 360 degrees):
function data = modelOB(vars)
% vars = shared_variables;
nTrials = 2*vars.nTrials;
replace_data = vars.replace;
replace_data = 1;
H = vars.H;
sigma_c = vars.sigma_c;
sigma_b = vars.sigma_x;
sigma_obs = vars.sigma_y;
circSwitch =true;

switchPoint = 50; % just for fitting piecewise regression model

B = nan(nTrials,1);
C = nan(nTrials,1);
Y = nan(nTrials,1);
S = nan(nTrials,1);

if replace_data == false
    for t = 1:nTrials
        %     If it is the first trial we need to create a random point for C to
        %     start
        if t == 1
            C(t) = rand*360;
        else
            C(t) = normrnd(C(t-1), sigma_c);
            while C(t)<0||C(t)>360
                C(t)=normrnd(C(t-1), sigma_c);
            end
        end
        
        %     Generate the S value.
        %        ==>  If S  = 1 with p=H, then the value of B is not
        %             going to be drawn from C but from a uniform distribution
        %        ==> .If S = 0, with probability 1-H then the vaue of B will be
        %        drawn from a normal distribution with mean B and variance of
        %        sigma_b
        if rand < H
            S(t) = 1;
            B(t) = rand*360;
        else
            S(t) = 0;
            B(t) = normrnd(C(t), sigma_b);
            while B(t)<0 || B(t)>360
                B(t)=normrnd(C(t), sigma_b);
            end
        end
        
        %     Generate the value of Y
        %       ==> It is drawn from a normal distributuon with mean of B and
        %       variance of sigma_obs
        Y(t) = normrnd(B(t), sigma_obs);
        
    end
else
    
    data = load(vars.path);
    data = data.alldata;
    if data.condition(1,1)==1
        C = [data.generativeMean(data.block==4,1); data.generativeMean(data.block==4,2)] ;
        B = [data.colorArray(data.block==4,1); data.colorArray(data.block==4,2)];
    else
        C = [data.generativeMean(data.block==3,1); data.generativeMean(data.block==3,2)] ;
        B = [data.colorArray(data.block==3,1); data.colorArray(data.block==3,2)];
    end
    
    
    nTrials = length(B);
    S = nan(nTrials,1);
    Y = nan(nTrials,1);
    for i = 1:nTrials
        Y(i)=normrnd(B(i), sigma_obs);
        
    end
    
end


%% STEP 2: MAKE INFERENCE ON B AND C GIVEN Y

possB = 1:360;
possC = 1:360;

pC = ones(length(possC), 1)./length(possC); %---> A bidimensinal array of size length(possC)*1
%     containing only the
%     probability of getting a
%     certain C from a uniform
%     distribution containing
%     360 C's. Note that this
%     vector will be all the
%     same only for the frist
%     trial, this will change
%     and be passed on to the
%     next trials with
%     different values.
pB = ones(length(possB), 1)./length(possB);
transFunc=eye(length(possB)); %---> Transition function: Identity square
%     matrix whose size equals
%     the number of possible B's
transFunc(:,:,2)= repmat(pB', length(possB), 1);
%---> This is A three dimensional matrix
%whose 3rd dimension is just the value of
%getting a certain X from a uniform
%distribution with 100 possible X's

% We need a 3-dimenassional matrix that will store the value of getting
% P(B | C,S) because we will use many times in the following iteration and we do not want
% to create it every time in every new iteration.
pBgivCandS = []; % ==> (1) = C
% ==> (2) = B
% ==> (3) = S
pCtgivC = []; % ==> (1) C_{t-1}
% ==> (2) C_{t}
for c = possC
    if circSwitch
        pBgivC = circ_vmpdf(deg2rad(possB), deg2rad(c), 1./(deg2rad(sigma_b).^2))';
    else
        pBgivC = normpdf(possB, c, sigma_b);
    end
    pBgivC = pBgivC./sum(pBgivC);
    
    pBgivCandS = [pBgivCandS; pBgivC];
    
    if circSwitch
        likelihoodCt = circ_vmpdf(deg2rad(possC), deg2rad(c), 1./(deg2rad(sigma_c).^2))';
    else
        likelihoodCt = normpdf(possC, c, sigma_c);
    end
    likelihoodCt = likelihoodCt./sum(likelihoodCt);
    pCtgivC = [pCtgivC;likelihoodCt]; % (1): C_{t-1}
    % (2): C_t
end

pBandSgivC = pBgivCandS * (1-H);
pBandSgivC(:,:,2) = repmat(pB'* H, length(possB), 1);



for t = 1:nTrials
    
    %storing entropy (uncertainty) of trial t
    entropy(t) = sum(pC .* log2(1 ./ pC));
    
    
    %Multiply P(C_t | C_{t-1}) * P(C_{t-1}) = P(C_t, C_{t-1})
    pCtandC = pCtgivC .* pC;
    
    pCt = sum(pCtandC, 1)';
    
    %multiply P(B,S | C) * P(C) = P(B, S, C)
    pBandCandS = pBandSgivC .* pCt;
    
    %   integrate over S (3rd dimension) in P(B, C, S)to get the value of P(B, C)
    pBandC = sum(pBandCandS, 3);
    
    %   Getting the likelihood of P(Y | B)
    if circSwitch
        pYgivB = circ_vmpdf(deg2rad(Y(t)), deg2rad(possB)', 1./(deg2rad(sigma_obs)^2));
        pYgivB = pYgivB./sum(pYgivB);
    else
        pYgivB = normpdf(Y(t), possB, sigma_obs);
        pYgivB = pYgivB./sum(pYgivB);
    end
    
    
    %   Calcuating the probabilit of P(B, C, Y = yt) by multiplying P(Y | B) *  P(B, C)
    pCandBandYandS = pBandCandS.* repmat(pYgivB', size(pBandC, 1), 1);
    pCandBandYandS = pCandBandYandS ./ sum(pCandBandYandS(:));
    surprise_t = sum(pCandBandYandS, 1);
    surprise_t = sum(surprise_t, 2);
    surprise_t = squeeze(surprise_t);
    surprise(t) = surprise_t(2);
    pCandBandY =   pBandC.* repmat(pYgivB', size(pBandC, 1), 1);
    pCandBgivY = pCandBandY./sum(pCandBandY(:));
    
    pCgivY = sum(pCandBgivY, 2);
    pC = pCgivY;
    
    pBgivY = sum(pCandBgivY, 1);
    
    
    
    posterior(t,:) = pCgivY;
    posteriorPercept(t,:)=pBgivY;
    %Calculating the shannon information content of each trial
    prob(t) = sum(pCandBandY(:));
    
    
    shannon(t) =log2(1/sum(pCandBandY(:)));
    
    
    
    % Calculating the relative entropy between trial t-1 and trial t for B
    if t == 1
        temp = repmat((1/length(possC)), length(possC), 1);
        temp = log(temp' ./ posterior(t,:));
        relative_entropy_C(t,:) = sum((1/length(possC)) .* temp);
    else
        temp =  log(posterior(t-1,:) ./  posterior(t,:));
        relative_entropy_C(t,:) =  sum(posterior(t-1,:) .* temp);
    end
    
    
end

%%
clear maxLikePostMu maxLikePostX
for i = 1:size(posteriorPercept, 1)
    maxLikePostB(i)=possB(find(posteriorPercept(i,:)==max(posteriorPercept(i,:)), 1));
    maxLikePostC(i)=possC(find(posterior(i,:)==max(posterior(i,:)), 1));
    
    
    expectationB(i)=posteriorPercept(i,:)* possB';
    expectationC(i) = posterior(i,:) * possB';
    
    
end


mean((Y-B).^2);            % just reporting sensory experience directly
mean((maxLikePostB'-B).^2); % our (max a-posteriori -- has error due to grid approximation)
mean((expectationB'-B).^2); % our model (expectation over B)



%% Measuring the models update and perceptual bias


% Measuring update
maxLikePostC= [mean(possC), maxLikePostC]; % Note = we are now indexing as if this is at the beginning of the trial!
cUpdate=circ_dist(deg2rad(maxLikePostC(2:end)'), deg2rad(maxLikePostC(1:end-1)'));
% [~, cUpdate, ~] = difference(deg2rad(B), deg2rad(maxLikePostC(1:end-1)), 'polar');
% cUpdate = cUpdate';
predError=circ_dist(deg2rad(maxLikePostB), deg2rad(maxLikePostC(1:end-1)));
 objPredError=circ_dist(deg2rad(B), deg2rad(maxLikePostC(1:end-1)'));
%objPredError=circ_dist(deg2rad(B), deg2rad(maxLikePostC(1:end-1)));


% Measuring perceptual bias
predictionErrorOnB=circ_dist(deg2rad(B),deg2rad(maxLikePostC(1:end-1)'));% for behave data
perceptualErrorOnB=circ_dist(deg2rad(B),deg2rad(maxLikePostB)'); % for behave data
% predictionErrorOnB=circ_dist( deg2rad(B),deg2rad(maxLikePostC(1:end-1)));
% perceptualErrorOnB=circ_dist(deg2rad(B),deg2rad(maxLikePostB) );

noisyPrediction = normrnd(maxLikePostC,sigma_obs);
noisyPerception = maxLikePostB;

noisyPredEOnB = circ_dist(deg2rad(B),deg2rad(noisyPrediction(1:end-1)'));
noisySubPredEOnB = circ_dist(deg2rad(noisyPerception),deg2rad(noisyPrediction(1:end-1)));
noisyPercEOnB = circ_dist(deg2rad(B),deg2rad(noisyPerception)');

%% Storing data, all data is stored in the form of horizontal vectors
data.nTrialsOB = nTrials;
data.YOB = Y;
data.B = B;
data.C = C;
data.SOB = S;
data.shannonOB = shannon;
data.klOB = relative_entropy_C';
data.entropyOB = entropy;
data.predictionErrorOnB = predictionErrorOnB';
data.perceptualErrorOnB = perceptualErrorOnB';
data.predErrorOB = predError;
data.objPredErrorOB=objPredError;
data.cUpdate = cUpdate;
data.posteriorPerceptOB = posteriorPercept;
data.posteriorOB = posterior;
data.expectationB = expectationB;
data.expectationC = expectationC;
data.maxLikePostC = maxLikePostC;
data.maxLikePostB = maxLikePostB;
data.noisyPredictionOB = noisyPrediction;
data.noisyPerceptionOB = noisyPerception;
data.noisyPredEOnB = noisyPredEOnB;
data.noisyPercEOnB = noisyPercEOnB;
data.noisySubPredEOnB = noisySubPredEOnB;
data.surpriseOB = surprise;
% save('../Desktop/neuro/odball_test_data', 'data');


data=straightStruct(data);


%% Plots:

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
%         small_b = regress(cUpdate(small_error)',[ones(size(predError(small_error)))' predError(small_error)']);
%         large_b = regress(cUpdate(~small_error)',[ones(size(predError(~small_error)))' predError(~small_error)']);
%         plot(predError(small_error) , cUpdate(small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'blue');
%         plot(predError(~small_error) , cUpdate(~small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'red');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%     elseif xx==2
%         small_error = abs(predictionErrorOnB) < threshold;
%         small_b = regress(perceptualErrorOnB(small_error),[ones(size(predictionErrorOnB(small_error))) predictionErrorOnB(small_error)]);
%         large_b = regress(perceptualErrorOnB(~small_error),[ones(size(predictionErrorOnB(~small_error))) predictionErrorOnB(~small_error)]);
%         plot(predictionErrorOnB(small_error), perceptualErrorOnB(small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'blue');
%         plot(predictionErrorOnB(~small_error), perceptualErrorOnB(~small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'red');
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%         ylabel('Perceptual Bias')
%         xlabel('prediction error')
%     end
%
%
%     setPLOT_panelLabel(gca, xx);
% end
%
%  kk=annotation('textbox')
%  set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
%

end








