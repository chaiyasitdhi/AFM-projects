% ################### Park AFM Stiffness Mapping 1.1 #####################
% ####################### by Atitheb Chaiyasitdhi ########################
% ########### King's Mongkut University of Technology Thonburi ###########
% ################### create : 21-09-2015 ################################
% ################### update : 22-09-2015 ################################

clear
clc

% ############################# READ ME ! ###############################
% #######################################################################

% 1.) Export all point data and image reference into text files by using 
% XEI software (Park System) 
% *** Check a folder "HowToExport" for help *** 

% 2.) Specify a directory of point data 
% CAUTION! this directory must have only point data files, info file and
% image reference file!
DataLocation = ['/home/theb/Documents/ubeTestData/S4/'];

% 3.) Specify a directory of MATLAB scripts
AlgorithmLocation = ['/home/theb/Documents/ParkAFMStiffnessMapping/Algorithm/'];  

% 4.) Specify units for conversion
forceUnit = 10e-9;      % force unit in nanoNewton (10e-9)
deflUnit = 10e-6;       % delfection unit in micron (10e-6)

% 5.) Specify the option if you want to show reference image, but you have
% to export image file into the data directory
show_image = 0; % 0 = not show, 1 = show

% 6.) Colorbar range setting
init_colorbar = 1; % N/m
final_colorbar = 4; % N/m

% 7.) Fitting range setting (in percent)
init_data_range = 10; % in percent
final_data_range = 90; % in percent 

% 8.) Run this MATLAB script (Press F5)

% #######################################################################
% #######################################################################
% #######################################################################


%% ####################### File Import Setting ##########################
% #######################################################################

cd(DataLocation)    % move to the data directory

InfoFile = dir('*info.txt');    % searching for the info file
if isempty(InfoFile)
   fprintf('Error: Cannot find the *info.txt\n');
   break
end

% Import the "info" file into MATLAB workspace. The info file contains 
% coordinates of all force-distance curves. We use this information to 
% construct an elasticity map.
infoObject = importdata(InfoFile.name); % import data from the info file
infoData = infoObject.data;

Point = infoData(:,1);      %Number of Data Points 
PositionX = infoData(:,2);  %X coordinate (in um)
PositionY = infoData(:,3);  %Y coordinate (in um)

%N of point data
numPoint = numel(Point);

%resolution in pixels
sizeX = numel( find( PositionX == PositionX(1) ) );
sizeY = numel( find( PositionY == PositionY(1) ) );

%dimension in um
% find the dimension of the 2D map 
% (sometimes data points are not collected with a square frame
dimX = abs( PositionX(sizeX)-PositionX(1) );
dimY = abs( PositionY(end)-PositionY(1) );

File = dir('*.txt');    % acquire filenames with *.txt
index_count = 1;

for itr = 1:length({File.name})
    if ~isempty( regexp(File(itr).name,'^^S\d{1}\.tiff_info\.txt','match') )
        image_ref_filename = File(itr).name;
    else
        Filename{itr} = File(itr).name;
    end
end

% check if the REGEX pattern correcly matches point data files
if length(Filename) ~= numPoint
    fprintf('Error: Point data files mismatching\nCarefully check that the REGEX pattern is correct\n');
    break
end

%create empty arrays for Stiffness and R-square
stiffness = zeros(1,numPoint);
rsq = zeros(1,numPoint);

% #######################################################################
% #######################################################################
% #######################################################################
   
%% ##################### Stiffness Extraction #####################
% #######################################################################

for itr = 1:numPoint
LoadFile = Filename{itr};       % current filename 
cd(DataLocation)                % move to the data directory

% read data from the specified file
[idx ZScan ZDetector Force] = textread(LoadFile,...
    '%u %f %f %f', 'headerlines', 3, 'whitespace', '\t');
% idx = data index
% ZScan = Z-position from the extension of the piezo-electric scanner 
% ZDetector = Z-position from the Position-Sensitive-Photodiode(PSPD)
% 'headerlines' = 3, assgin the first 3 lines as header lines
% 'whitespace' = '\t', use the horizontal tab as the delimiter 

cd(AlgorithmLocation);      % move to the algorithm directory
% fitting by using linear regression with a slope_fit function
[stiffness(1,itr), rsq(1,itr)] = slope_fit(ZDetector.*deflUnit, Force.*forceUnit,...
    init_data_range, final_data_range);

% print out "current_index/total_index"
fprintf(' data No#: %i /%i \n',itr,numPoint); 
end

% create reference image
if show_image == 1
    cd(DataLocation)                % move to the data directory
    % import image reference metadata
    [X xUnit Y yUnit S sUnit] = textread(image_ref_filename,...
    '%f %s %f %s %f %s', 'headerlines', 7, 'whitespace', '\t');
    % X = position X
    % xUnit = unit of X (um)
    % Y = position Y
    % yUnit = unit of Y (um)
    % S = position stiffness (XEI's results)
    % xUnit = unit of stiffness (N/m)
    % 'headerlines' = 7, assgin the first 7 lines as header lines
    % 'whitespace' = '\t', use the horizontal tab as the delimiter 

    %image size
    im_sizeX = numel( find( X == X(1) ) );
    im_sizeY = numel( find( Y == Y(1) ) );

    image_filename = dir('*.png');
    image_ref_file = imread(image_filename.name);
    
    figure()
    imshow( imrotate(image_ref_file,-90) ); 
    title('Reference Image');
    xlabel('Pixel');
    ylabel('Pixel');
end

% create a stiffness map
stiffness_map = reshape(abs(stiffness),sizeX,sizeY);    
% create a r-square map
rsquare_map = reshape(rsq,sizeX,sizeY);

% #######################################################################
% #######################################################################
% #######################################################################

%% ########################## Visualization #############################
% #######################################################################

figure()
imagesc(stiffness_map); title('Stiffness Map'); 
xlabel('Pixel');
ylabel('Pixel');
c1 = colorbar; ylabel(c1, 'N/m');
set(gca, 'CLim', [init_colorbar, final_colorbar]); 
% setting a range of colorbar
daspect([1 1 1]);

figure()
imagesc(rsquare_map); title('RSquare Map'); 
xlabel('Pixel');
ylabel('Pixel');
c2 = colorbar(); caxis([0 1]); ylabel(c2, 'A.U.');

% #######################################################################
% #######################################################################
% #######################################################################

