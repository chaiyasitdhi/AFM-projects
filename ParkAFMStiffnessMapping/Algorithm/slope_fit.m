function [slopeF, rsq] = slope_fit(ZDetector,Force)
% choose only an approaching curve
zDist = ZDetector(1: length(ZDetector)/2); %Z distance
force = Force(1: length(Force)/2); %force

% use 2SD of the first 20% to define the threshold 
thermalDeflLength = round( length(force)*(0.2) );
thermalDefl = diff( force(1:thermalDeflLength) );
threshold = 2*( std(thermalDefl + mean(thermalDefl)) );

% find the initial and final position of the data we want to fit with
% a linear regression model
overThresholdPos = find(diff(force) >= threshold);
ThresholdPos = overThresholdPos(1); % initial position
FinalPos = overThresholdPos(end);   % final position

% extract data by using initial and final position
nForce = force(ThresholdPos:FinalPos);
nzDist = zDist(ThresholdPos:FinalPos);

% fitting with a linear regression model
fit_params = polyfit(nzDist,nForce,1);
slopeF = fit_params(1);
predictedForce = polyval(fit_params, nzDist);

% calculate R-square
yresid = nForce-predictedForce;
SSresid = sum(yresid.^2);
SStotal = (length(nForce)-1)*var(nForce);
rsq = 1 - SSresid/SStotal;


