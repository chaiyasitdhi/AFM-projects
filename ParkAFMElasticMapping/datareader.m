% ############################ Data Reader ################################
% ################# Atitheb Chaiyasitdhi, UPDATED: 17-12-2014 #############

% Here is the algorithm for elasticity mapping of data obtained from the 
% AFM force spectroscopy in a liquid phase ONLY. This algorithm is most 
% suitable for the model XE-120 from PARK Systems.
warning('off','all')

%% ### File Import Setting ###

% input the data directory
DataLocation = ['/home/theb/Documents/AFMresults/Spirulina_liquidphase_FDspectroscopy_29062015_centralLabVsC005/test4'];
% input the algorithm direcctory
AlgorithmLocation = ['/home/theb/Documents/ElasticNewAl'];                

%% ### Data Acquisition ###

% First, import the "info" file into MATLAB. The info file contains 
% coordinates of the data collected. We use information from the info 
% file to automatically construct an elasticity map.

% Point = Number of Data Points 
% PositionX = X coordinate 
% PositionY = Y coordinate

cd(DataLocation);  
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
    % check if filename matches default format with regular expression
    if ~isempty( regexp(File(itr).name,'^.+Spectroscopy-\d{3}_\d{4}\.txt','match') )
        % store matched filenames in a specified folder
        Filename(index_count) = {File(itr).name}; 
        index_count = index_count + 1;
    elseif ~isempty( regexp(File(itr).name,'^.+Spectroscopy-\d{3}_glass_\d{4}\.txt','match') )
        glass_ref_filename = File(itr).name;
    else
        if isempty( regexp(File(itr).name,'^.+Spectroscopy.+info\.txt','match') )
            image_ref_filename = File(itr).name;
        end
    end
end

%% ### Parameter Setting ###
% ZDetector = raw Z-distance (X)
% Force = raw force data (Y)
  
%% ### Elastic Modulus Extraction ###
% for in_depth = 100:10:300
% fprintf('fit: %i \n',in_depth);    
cd(DataLocation);

% read data from the specified file
[gidx gZScan gZDetector gForce] = textread(glass_ref_filename,...
    '%u %f %f %f', 'headerlines', 3, 'whitespace', '\t');

cd(AlgorithmLocation);              % move to the algorithm directory
[g_elastic_modulus gR2 g_init_con] = elastic_subtraction(gZDetector, gForce);

for itr = 1:numPoint
cd(DataLocation);                   % move to the data directory
LoadFile = Filename{itr};           % current filename 

% read data from the specified file
[idx ZScan ZDetector Force] = textread(LoadFile,...
    '%u %f %f %f', 'headerlines', 3, 'whitespace', '\t');

% idx = data index
% ZScan = Z-position from extension of the piezo-electric scanner 
% ZDetector = Z-position from the Position-Sensitive-Photodiode(PSPD)
% 'headerlines' = 3, assgin the first 3 lines as header lines
% 'whitespace' = '\t', use the horizontal tab as the delimiter 

cd(AlgorithmLocation);              % move to the algorithm directory

% extract elastic modulus of the specified data 
[elastic_modulus R2 init_con] = elastic_subtraction(ZDetector, Force);
%elastic_modulus = CPD(Zdist,subtF,0.06,lowerb_contact,v,'cone',alpha);
%clear ZDetector Force 

%%[Py Px] = find(mapPosition == itr); % locate the position of data point
%[EPy EPx] = find(mapPosition == itr);
%elastic_map(Px,Py) = elastic_modulus; % fill the elastic modulus to a map
elastic_map(itr,1) = elastic_modulus;
error_map(itr,1) = R2;
contact_map(itr,1) = init_con-g_init_con;

fprintf(' data No#: %i /%i \n',itr,numPoint);

end

%% ### Plot the Results ###
fimage = flipimage(elastic_map);% rescale to log(E)
eimage = flipimage(error_map);
iimage = flipimage(contact_map);


figure()
imagesc(fimage);
surf_elastic = imcrop;
close all

thresh_surf = min(min(surf_elastic));
fnimage = fimage;
%fnimage(eimage < 0.5) = NaN;
fnimage(fimage > thresh_surf) = NaN;

% fstack(:,:, (in_depth-90)/10 ) = fimage;
% estack(:,:,(in_depth-90)/10 ) = eimage;

figure()
subplot(1,3,1); im1 = pcolor(fnimage); 
daspect([1 1 1]); set(im1, 'edgecolor', 'none')    
subplot(1,3,2); imagesc(eimage); daspect([1 1 1]);
subplot(1,3,3); imagesc(iimage); daspect([1 1 1]);

%set(gca, 'YTickLabel', num2str(coor_y));    % real scale
%set(gca, 'XTickLabel', num2str(coor_x));    % real scale
% 
% title(['Elasticity Map ( log(Pa) ) ' num2str(idepth)]);
% xlabel('Scan Size (Pixel)');
% ylabel('Scan Size (Pixel)');
%colorbar();                                 % show colorbar

%saveas(h,['emap_' num2str(in_depth) 'um.png'])
%save(['edat_' num2str(in_depth) '.mat'],'fimage');
%close
% end




