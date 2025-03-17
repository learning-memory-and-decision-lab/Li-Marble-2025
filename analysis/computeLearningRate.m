function [LR, UP, PE]=computeLearningRate(outcomes, predictions, newBlock, coordSys)
% Compute update (UP), prediction error (PE), and learning rate (LR=UP/PE)
% for either cartesian or polar coordinates
% Last modified by Ryan Thorpe on 2020.01.24

% Define default coordinate system
if nargin<4
    coordSys = 'cartesian';
end

% find the last trial of each block
a=find(newBlock)-1;
a=a(a>1);

UP=nan(length(outcomes)-1,1);   %nans
UP_options=nan(length(outcomes)-1,3);

switch coordSys
    case 'cartesian'
        UP(1:end) = predictions(2:end)-predictions(1:end-1);
        PE = outcomes-predictions;
    case 'polar' %%% predictions [nsamples, 2]
        PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations 
        UP_options(:,1) = circ_dist(predictions(2:end),predictions(1:end-1)); % Distance one way around the circle
        UP_options(:,2) = UP_options(:,1)- 2*pi; % Distance the other way around the circle
        UP_options(:,3) = UP_options(:,1)+ 2*pi; % Distance the other way around the circle
        
        [~,opt_ind] = min(abs(UP_options-PE(1:end-1)),[],2); % Index
        
        % Select the update that travels in the direction of PE
        for ti=1:length(opt_ind)
            UP(ti) = UP_options(ti,1);
        end

    case 'polarHalfCorrect'
        PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations 
        UP_options(:,1) = circ_dist(predictions(2:end),predictions(1:end-1)); % Distance one way around the circle
        UP_options(:,2) = UP_options(:,1)- 2*pi; % Distance the other way around the circle
        UP_options(:,3) = UP_options(:,1)+ 2*pi; % Distance the other way around the circle
        
        [~,opt_ind] = min(abs(UP_options-(PE(1:end-1))),[],2); % Index
        kk = [PE(1:end-1),UP_options(:,1),opt_ind];
        %keyboard
        % Select the update that travels in the direction of PE
%         for ti=1:length(opt_ind)
%             UP(ti) = UP_options(ti,opt_ind(ti));
%         end
        for ti=1:length(opt_ind)
            if abs(UP_options(ti,1)-PE(ti)) > 3*pi/2
                UP(ti) = UP_options(ti,opt_ind(ti));
            else
                UP(ti) = UP_options(ti,1);
            end
        end
%         
    case 'polarNoCorrect'
         PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations 
         UP(1:end) = circ_dist(predictions(2:end),predictions(1:end-1)); % Distance one way around the circle
        
    otherwise
        error("computeLR error: didn't recognize coordinate system!! Use coordSys='cartesian' or 'polar'.")
end

% Learing rate
LR=UP./PE(1:end-1);

% Undefined for the first trial of each block
UP(a)=nan;
PE(a)=nan;
LR(a)=nan;
end
