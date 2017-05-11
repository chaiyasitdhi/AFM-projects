function [n_bound_x, n_bound_fit, gauss_eudist_image, normResidual] = fit_helix(fileHandle, spatial_resolution)

 fileNames = {fileHandle.name};
    
    for file_i = 1:length(fileNames)
        imRead = imread(fileNames{file_i});
        
        if numel(size(imRead)) == 3
            bw_image = im2bw(imadjust(rgb2gray(imRead)));
        else
            bw_image = im2bw(imadjust(imRead));
        end

        eudist_image = bwdist(bw_image);
        gauss_eudist_image{file_i} = imgaussian(eudist_image, 2);
        bw_gauss_image = im2bw(gauss_eudist_image{file_i});

        [cPx, cPy, n_image] = extractCenterLine(gauss_eudist_image{file_i});
        ncPx = [cPx(1):cPx(end)];
        ncPy = spline(cPx, cPy, ncPx);

        cPx_s = ncPx./spatial_resolution;
        cPy_s = ncPy./spatial_resolution;

        %rotate
        params = polyfit(cPx_s, cPy_s, 1);
        yfit_original = polyval(params, cPx_s);
        yfit = polyval([0 params(2)], cPx_s);

        diffY = yfit - yfit_original;
        cPy_s = cPy_s + diffY;
% 
%         cPy_s = cPy;
%         cPx_s = cPx;
        
        %find peaks
        diffY = gradient(cPy_s)./gradient(cPx_s);
        diffY = imgaussian(diffY,5);
        [pks,locs] = findpeaks(diffY);
        [pks_neg start] = findpeaks(-diffY);

        [Y_peaks pks_post] = find(abs(diffY) < 0.05);
        bound_x = [cPx_s(locs), cPx_s(start), cPx_s(pks_post)];
        bound_y = [cPy_s(locs), cPy_s(start), cPy_s(pks_post)];
        [bound_fit_params S] = polyfit(bound_x, bound_y, 3);
        
        normResidual(file_i) = S.normr;
        
        %n_bound_x{file_i} = [cPx_s(20):cPx_s(end-20)];
        %n_bound_fit{file_i} = polyval(bound_fit_params, n_bound_x{file_i});
        n_bound_x{file_i} = cPx_s;
        n_bound_fit{file_i} = imgaussian(cPy_s, 60, length(cPx_s));
        
         figure()
          %plot(cPx_s, cPy_s, 'b'); daspect([1 1 1]);
%          imshow(~bw_gauss_image);
          hold on; plot(n_bound_x{file_i}, n_bound_fit{file_i}, 'r--', 'LineWidth', 2);
%           plot(cPx_s(locs), cPy_s(locs), 'ro', 'MarkerFaceColor', 'red');
%           plot(cPx_s(start), cPy_s(start), 'ro', 'MarkerFaceColor', 'red');
%           plot(cPx_s(pks_post), cPy_s(pks_post), 'go', 'MarkerFaceColor', 'green');
    end
            %mean(normResidual)
end





