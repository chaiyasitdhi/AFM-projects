function [E R2 ub] = CPD(z,defl,k,G,nu,model_geo,model_para)
 
%%_______________________________________________________________________%%
% CPD function is a MATLAB function based on the Contact Point Detection 
% (CPD) algorithm for automated AFM force curve analysis. 
% Required Inputs: 
    % z    [=] m       z is Z sensor data from AFM. 
    % defl [=] m       defl is Deflection data from AFM. 
    % k    [=] pN/nm   k is the AFM cantilever spring constant.
    % G    [=] 1       G is the threshold level defined as 5 in CPD 
    %                   algorithm.
    % nu   [=] 1       nu is Poisson's ratio. nu is 0.5 for incompressible
    %                   materials. This may be changed depending on
    %                   analysis.
    % model_geo        model_geo represents the Hertz model geometry. For
    %                   conical tip geometry, set model_geo to 'cone'. For 
    %                   spherical tip geometry, set model_geo to 'sphere'. 
    % model_para       model_para represents the tip geometry parameter. 
    %                   For conical tips, this is the half-angle opening in
    %                   degrees. For sphereical tips, this is the sphere
    %                   radius in micrometers.
% Outputs: 
    % figure(1)     Raw data plot 
    % figure(2)     Approach and retract data
    % figure(3)     Zeroing deflection
    % figure(4)     Zeroed deflection plot
    % figure(5)     Force versus separation
    % figure(6)     Force versus indentation with Hertz model
%%_______________________________________________________________________%%
 
% Start of code
if nargin == 0
    disp('Error: No inputs')
    close all
    return
end
  
if all(z) == 0 && all(defl) == 0
    E = 0;
    R2 = 0;
    ub = 0;
end

% % Convert data to nanometers and correct for constant compliance
z = z.*10^9;                    % [=] nm
defl = defl.*10^9;              % [=] nm
 
% Split the data into approach and retract
% Find the max defl and split the data into approach and retract
% max_loc = find(defl == max(defl));
% z_r = z(max_loc+1:end);         % z retract [=] nm
% defl_r = defl(max_loc+1:end);   % defl retract [=] nm
% z = z(1:max_loc);               % z approach [=] nm
% defl = defl(1:max_loc);         % defl approach [=] nm
 
% Find conversion between data points and nanometers
% Do this by taking the slope of z position versus index
% This is calculated for convenience and used in other parts of the code
points = [1:length(z)]';                
fit_points = polyfit(points,z,1);
slope = abs(fit_points(1));             % slope [=] nm/point
 
% Zero the deflection 
% First define a range of data (first 40%)
zero_defl_range = [z(1:round(0.4*length(z))) defl(1:round(0.4*length(z)))];
% Fit a line to this data
zero_defl_fit = polyfit(zero_defl_range(:,1),zero_defl_range(:,2),1);
 
% Calculate the line for all data and subtracts the line from data
zero_defl_eval = zero_defl_fit(1).*z + zero_defl_fit(2);
zero_defl = defl - zero_defl_eval;
 
% Zero separation
% Find the min and max deflection values
maxDefl = z(zero_defl == max(zero_defl));
minDefl = z(zero_defl == min(zero_defl));

if size(maxDefl,1) == 0
   E = 0;
   R2 = 0;
   ub = 0;
   return
end

% Take the average of the min and max and define separation
avg_MinMax = mean([maxDefl(end) minDefl(end)]);
sep = -(zero_defl - z);                     % separation [=] nm
zero_sep = sep - avg_MinMax;                % zeroed separation [=] nm
 
% Calculate force and plot force versus separation
force = zero_defl.*k./1000;                 % force [=] nN 
 
% Find the approximate contact point by detecting threshold level
% Take average and standard deviation of initial portion of data
avg = mean(force(1:round(0.4*length(zero_sep))));
stdev = std(force(1:round(0.4*length(zero_sep))));
for i = 20:length(zero_sep)
    if force(i) < avg + G*stdev
        P = i;
    end
end
if exist('P','var') == 0
   E = 0;
   R2 = 0;
   ub = 0;
   return
end

% Define the Hertz model 
if strcmp(model_geo,'sphere') == 1
        
    R = model_para./1000;        % Define the radius of sphere [=] nm
    f = @(x,s) (4/3).*x(1).*R.^(0.5).*((s-x(2)).^(3/2))./(1-nu.^2);
    f_1 = @(x,I) (4/3).*x(1).*R.^(0.5).*(I.^(3/2))./(1-nu.^2);
        
elseif strcmp(model_geo,'cone') == 1
        
    alpha = model_para.*pi/180;  % Define the half angle of the tip [=] rad
    f = @(x,s) (2/pi).*tan(alpha).*(x(1)./(1-nu.^2)).*((s-x(2)).^2);
    f_1 = @(x,I) (2/pi).*tan(alpha).*(x(1)./(1-nu.^2)).*(I.^2);
        
else
    disp('Error: input a valid model. Accepted models: "sphere" or "cone"')
    close all
    clear all
    clc
    return
        
end
 
% Provide an initial guess for the nlinfit function
% We found that the initial guess does not affect the fitting process
Ini_guess = 50/10^6;                % Initial modulus guess 50 kPa
 
% Check the slope from approximate contact point to max force
slope_fit = polyfit(zero_sep(P:end),force(P:end),1);
slope_a = slope_fit(1);                       % [=] nN/nm
 
% Define the indentation length based on a measure of how stiff the sample
% is
if slope_a > -0.1 && slope_a < -0.01
    ub = 100;            % ub is the indentation length(upper bound) [=] nm 
elseif slope_a < -0.1
    ub = 50;
elseif slope_a < -0.001 && slope_a > -0.01
    ub = 200;
else
    ub = 100;
end

% Check how much data is after the approximate contact point
% This check is to ensure that the indentation length defined does not 
% exceed the actual indentation length in the data
points_remaining = length(zero_sep) - P;
nanometers_remaining = points_remaining*slope;      % [=] nm
if nanometers_remaining < ub
    ub = floor(nanometers_remaining)-5;             % [=] nm
end
 
% Identify an initial fit region using the approximate contact point P
X = zero_sep(P:P+floor(ub/slope));  % X and Y are the data set for the fit
Y = force(P:P+floor(ub/slope));
guess = [Ini_guess zero_sep(P)];

    S = nlinfit(X,Y,f,guess);       % Perform non-linear fit for CP
    E_fit_1 = S(1)*10^6;            % This value of E is discarded
    CP_fit_1 = S(2);                % This is the fitted CP [=] nm

if isnan(CP_fit_1)
   E = 0;
   R2 = 0;
   ub = 0;
   return
elseif or(CP_fit_1 >= max(zero_sep), CP_fit_1 < min(zero_sep))
   E = 0;
   R2 = 0;
   ub = 0;
   return
end

% Define indentation
I = zero_sep - CP_fit_1;    % I = indentation [=] nm
for i = 1:length(I)-1
    if abs(zero_sep(i) - CP_fit_1) < abs(zero_sep(i+1) - CP_fit_1)
        loc = i;
    break
    end
end
force_zero = force - force(loc);

% Define the fit region as the same indentation length used previously
X_1 = I(loc:loc+floor(ub/slope));
Y_1 = force_zero(loc:loc+floor(ub/slope));
guess_1 = Ini_guess;
Q = nlinfit(X_1,Y_1,f_1,guess_1);   % perform the non-linear fit
E_fit_2 = Q(1)*10^6;
E = E_fit_2;    % E is the final reported elastic modulus [=] kPa
 
%Plot force versus indentation and the Hertz model fit
% figure(6);
% plot(I,force_zero,'k-')
% title('Force vs Indentation Contact Point Fitted ')
% xlabel('Indentation [nm]')
% ylabel('Force [nN]')
% hold on
% Plot the contact point
%plot(I(loc),force_zero(loc),'rx','MarkerSize',10)
%plot(X_1,Y_1,'bo');
% Plot the Hertz model for the indentation range
fit_val_1 = feval(f_1,E_fit_2/10^6,X_1);
% % Calculate an R^2 value for the fit
SSres = sum((Y_1 - fit_val_1).^2);
SStot = (length(Y_1)-1)*var(Y_1);
R2 = 1 - SSres/SStot;

% plot(X_1,fit_val_1,'r-')
% axis([-ub-50 10 -1 max(force)]);
% legend('Data','Fitted CP','Fit Region','Hertz Model Fit')
% %Display E, R^2, indentation, threshold, and nu
% annotation('textbox',[0.5 0.6 0.35 0.1], 'String', {[num2str(E_fit_2)...
%     ' kPa' ' R^2: ' num2str(R2)];[num2str(ub) 'nm ' 'G' num2str(G) ...
%     ' \nu:' num2str(nu)]});
%  
% End of code
end

