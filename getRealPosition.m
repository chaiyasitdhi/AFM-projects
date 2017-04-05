function [cPx_s, cPy_s, gauss_eudist_image] = getRealPosition(fileHandle, spatial_resolution)

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

        cPx_s{file_i} = ncPx./spatial_resolution;
        cPy_s{file_i} = ncPy./spatial_resolution;
        
        clear n_image imRead bw_image eudist_image...
        bw_gauss_image ncPx ncPy bodyWidth
    end
    
end