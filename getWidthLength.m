function [avrBodyWidth, contour_length, relative_dist_initial] = getWidthLength(cPx, cPy, gauss_eudist_image, spatial_resolution)

    for file_i = 1:length(cPx)
        %calculate contour length 
        relative_dist{file_i} = sqrt(gradient(cPx{file_i}).^2 + gradient(cPy{file_i}).^2);
        relative_dist_initial{file_i} = cumsum([0, relative_dist{file_i}]);
        contour_length(file_i) = sum(relative_dist{file_i});
    
        
        %calculate body width
        bodyWidth = widthCalc2(gauss_eudist_image{file_i}')./spatial_resolution;
        bodyWidth = bodyWidth(~isnan(bodyWidth));
        bodyWidth_s{file_i} = bodyWidth(round(50*length(bodyWidth)/100):...
        round(75*length(bodyWidth)/100));
        avrBodyWidth(file_i) = median(bodyWidth_s{file_i});
    end
    
end