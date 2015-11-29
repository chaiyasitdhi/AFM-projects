function [best_elastic R2 ub]= elastic_subtraction(ZDetector,Force)
% ################# Elastic Modulus Extraction ############################
% ################# Atitheb Chaiyasitdhi, UPDATED: 17-12-2014 #############
% Here is the algorithm for extracting elastic modulus of data obtained 
% from the AFM force spectroscopy in a liquid phase ONLY. This algorithm 
% is most suitable for the model XE-120 from PARK Systems.

% ## General Paramaters ##
n_data = numel(ZDetector);                  % # of data points
force_unit = 10^-9;                         % nanoNewton (unit of force)
dist_unit = 10^-6;                          % microMeter (unit of distance)

%% ### Data Extraction ###

global v alpha
% *** For Fitting the Hertz Model to single data, SPECIFY the "v" and
% "alpha" below and UNCOMMENT the 3 lines.

% ########################################################################
v = 0.5; % poission ratio ~ 0.5, for incompressible materials(i.e. rubber)
face_angle =25; % face angle of the AFM tip in degree
alpha = deg2rad(face_angle); % convert the face angle into radian unit
spring_con = 600;
% ########################################################################

% Approaching and retracting curves should be the same. We decide to use
% only the approaching curve. 
rawZpost = dist_unit.*ZDetector(1:round(n_data)/2); % Approaching Z-position
rawforce = force_unit.*Force(1:round(n_data/2));   % Approaching Forces

% create an array of Z-distance 
Zdist = abs(rawZpost - rawZpost(1).*ones(size(rawZpost))); 

% ## Force Subtraction ## 
% We subtract the resulant forces by a median of the first 20% of the data.
int_def = rawforce(1:round(n_data*0.2),:);      % the first 20% of the data
subtF = rawforce - [median(int_def).*ones(n_data/2,1)];  %subtract the data
 
% ## Find the Initial Contact Point ##
initial_contact = contact_estimate(Zdist,subtF);    % initial contact point
 
% ## Find the Local Minima of the Elastic Modulus ##
% # Boundary Conditions for Hertz Model #
% Contact Parameters

if isempty(initial_contact)
    best_elastic = 0;
    R2 = 0;
    ub = 0;
    return
end

lowerb_contact = initial_contact(1);

[best_elastic R2 ub] = CPD(Zdist,subtF,spring_con,5,v,'cone',face_angle);
%[best_elastic R2 ub] = CPD(rawZpost,rawforce,spring_con,5,v,'cone',15);

%upperb_contact = in_depth*10^-9; % Hertz model only reliable within 10% of 
                            % the sample thickness

if lowerb_contact == 0
    lowerb_force = find(subtF > 0);
    contact_data = Zdist(subtF > 0);
else
    contact_data = [Zdist(Zdist >= lowerb_contact)];
    lowerb_force = find(Zdist >= contact_data(1));
end

end


