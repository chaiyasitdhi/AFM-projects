clear
clc
close all
cd('~/Documents/PersistentLengthSpirulina/');
% c005_handle = dir('c005*.jpg');
% cent_handle = dir('cent*.jpg');
% 
% c005_filenames = {c005_handle.name};
% cent_filenames = {cent_handle.name};
% all_filenames = {[c005_filenames, cent_filenames]};

spatial_resolution_ii = 1.700; %pixel/um
spatial_resolution_iii = 1.788; %pixel/um only for iii experiment

%experiment II
cent_ii_handle = dir('cent_*_ii.jpg');
c005_ii_handle = dir('c005_*_ii.jpg');

%experiment III
cent_iii_handle = dir('cent_*_iii.jpg');
c005_iii_handle = dir('c005_*_iii.png');

%experiment IIII
cent_iiii_handle = dir('cent_*_iiii.png');
c005_iiii_handle = dir('c005_*_iiii.png');

experiment = 'str';
handle_to_analysis = cent_iiii_handle;
spatial_resolution = spatial_resolution_ii;

if strcmp(experiment, 'str') % str for straight
    [cPx, cPy, gauss_eudist_image] = getRealPosition(handle_to_analysis, spatial_resolution);
else
    [cPx, cPy, gauss_eudist_image] = fit_helix(handle_to_analysis, spatial_resolution);
end
    
[bodyWidth, contourLength, relative_distance_initial] = getWidthLength(cPx, cPy, gauss_eudist_image, spatial_resolution);
[tangent_corr, deltaT, fit_data, tanVec] = getTanCorrelationFunction(cPx, cPy, contourLength, spatial_resolution);

figure()
for file_i = 1:length(tangent_corr)
    if file_i < 23
        plot(deltaT{file_i}, tangent_corr{file_i}, 'r--o'); hold on;
    else
        plot(deltaT{file_i}, tangent_corr{file_i}, 'b--o'); hold on;
    end
    
    curvature(file_i) = norm(gradient(tanVec{file_i}));
    
end

reliable_Lp = fit_data.Lp(fit_data.rsq >= 0.50);
reliable_bodyWidth = bodyWidth(fit_data.rsq >= 0.50);

morethanusual_id = find(reliable_bodyWidth <= 15 & reliable_bodyWidth >= 6);
reliable_bodyWidth = reliable_bodyWidth(morethanusual_id);
reliable_Lp = reliable_Lp(morethanusual_id)';

meanLp = mean(reliable_Lp)
stdLp = std(reliable_Lp)
seLp = stdLp./sqrt(numel(reliable_Lp));
numel(reliable_Lp)

l2ree_ratio = contourLength(fit_data.rsq>0.5)./fit_data.Lp(fit_data.rsq>0.5);
% for file_i = 1:length(all_filenames{1})
%     
%     imRead = imread(all_filenames{1}{file_i});
%     bw_image = im2bw(imadjust(rgb2gray(imRead)));
%     eudist_image = bwdist(bw_image);
%     gauss_eudist_image = imgaussian(eudist_image, 2);
%     bw_gauss_image = im2bw(gauss_eudist_image);
%     
%     [cPx, cPy, n_image] = extractCenterLine(gauss_eudist_image);
%     ncPx = [cPx(1):cPx(end)];
%     ncPy = spline(cPx, cPy, ncPx);
%         
%     cPx_s{file_i} = ncPx./spatial_resolution;
%     cPy_s{file_i} = ncPy./spatial_resolution;
%     
%     %calculate contour length 
%     relative_dist{file_i} = sqrt(gradient(ncPx).^2 + gradient(ncPy).^2);
%     relative_dist_initial{file_i} = cumsum([0, relative_dist{file_i}]);
%     contour_length(file_i) = sum(relative_dist{file_i});
%    
%     %calculate body width
%     bodyWidth = widthCalc2(gauss_eudist_image')./spatial_resolution;
%     bodyWidth = bodyWidth(~isnan(bodyWidth));
%     bodyWidth_s{file_i} = bodyWidth(round(50*length(bodyWidth)/100):...
%        round(75*length(bodyWidth)/100));
%     avrBodyWith(file_i) = median(bodyWidth_s{file_i});   
%    clear n_image imRead bw_image eudist_image...
%        gauss_eudist_image bw_gauss_image ncPx ncPy bodyWidth
% 
% end
% 
%    %calculate tangent correlation function
%    minimal_displacement = min(contour_length);
%    minimal_dis_pixel = minimal_displacement*spatial_resolution;
%    
%    %fittype
%    myfittype = fittype('a*exp(-x/b)',...
%     'dependent', {'y'}, 'independent', {'x'},...
%     'coefficients', {'a', 'b'});
% 
% for file_i = 1:length(all_filenames{1})
%     
%     display(file_i);
%     data = [cPx_s{file_i}; cPy_s{file_i}]';
%     nData = size(data,1); %# number of data points
%     numberOfDeltaT = floor(nData/3); %# for MSD, dt should be up to 1/4 of number of data points
% 
%     tangent_corr = zeros(numberOfDeltaT,1); %# We'll store [mean, std, n]
% 
%     %# calculate msd for all deltaT's
%     Px = data(:,1);
%     Py = data(:,2);
%     dydx = gradient(Py)./gradient(Px);
%     theta = atan(dydx);
%     r = [cos(theta), sin(theta)];
%     
%     for dt = 1:numberOfDeltaT
%        deltaT(dt) = dt;
%        rr = dot(r(1+dt:end,:),r(1:end-dt,:),2);    
%        tangent_corr(dt) = mean(rr); %# average
%        
%        clear rr
%     end
%     
%      deltaT = deltaT./spatial_resolution;
% %     figure();
% %     plot(deltaT, tangent_corr, '-o');
% %     filename = all_filenames{1}{file_i};
% %     title(filename(1:end-4), 'interpreter', 'none');
% %     ylabel({'C_{t}<ds>', '(U.A.)'});
% %     xlabel('ds (\mum)');
% %     % Export Figure
% %     set(gcf, 'PaperPositionMode', 'auto');
% %     print('-dpng', '-r300', ['tangent_corr_' filename(1:end-4) '.png']);
% %    close
%     tangent_corr_s{file_i} = tangent_corr;
%     deltaT_s{file_i} = deltaT;
%   
%     [myfit{file_i} gof{file_i}] = fit(deltaT', tangent_corr, myfittype,...
%         'StartPoint', [1, 500]);
%     b(file_i) = myfit{file_i}.b;
%     rsq(file_i) = gof{file_i}.rsquare;
%     sef(file_i) = gof{file_i}.rmse;
%     clear deltaT tangent_corr Px Py
% end
% 
% reliable_b = b(rsq >= 0.85);
% 
% figure()
% for file_i = 1:length(all_filenames{1})
%     if file_i < 23
%         plot(deltaT_s{file_i}, tangent_corr_s{file_i}, 'r-.'); hold on;
%     else
%         plot(deltaT_s{file_i}, tangent_corr_s{file_i}, 'b-.'); hold on;
%     end
%     
%     filename = all_filenames{1}{file_i};
%     %title(filename(1:end-4), 'interpreter', 'none');
%     ylabel({'C_{t}<ds>', '(U.A.)'});
%     xlabel('ds (\mum)');
% end
%     