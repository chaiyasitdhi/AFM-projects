function [ntangent_corr, ndeltaT, fit_data, tangentVectors] = getTanCorrelationFunction(cPx, cPy, contour_length, spatial_resolution)
    
   %calculate tangent correlation function
   minimal_displacement = min(contour_length);
   minimal_dis_pixel = floor(minimal_displacement*spatial_resolution);
   
   for file_i = 1:length(contour_length)
       
       %fittype
       myfittype = fittype('a*exp(-x/b)',...
        'dependent', {'y'}, 'independent', {'x'},...
        'coefficients', {'a', 'b'});

        display(file_i);
        data = [cPx{file_i}; cPy{file_i}]';
        
        nData = size(data,1); %# number of data points
        numberOfDeltaT = floor(nData/3); %# for MSD, dt should be up to 1/4 of number of data points

        tangent_corr = zeros(numberOfDeltaT,1); %# We'll store [mean, std, n]

        %# calculate msd for all deltaT's
        Px = data(:,1);
        Py = data(:,2);
        dydx = gradient(Py)./gradient(Px);
        theta = atan(dydx);
        r = [cos(theta), sin(theta)];
        tangentVectors{file_i} = r;
            for dt = 1:numberOfDeltaT
               deltaT(dt) = dt;
               rr = dot(r(1+dt:end,:),r(1:end-dt,:),2);    
               tangent_corr(dt) = mean(rr); %# average

               clear rr
            end

        deltaT = deltaT./spatial_resolution;
%         deltaT 
%         minimal_displacement
%          deltaT = deltaT(deltaT < minimal_displacement-250);
%          tangent_corr = tangent_corr(1:length(deltaT));

        tangent_corr_s{file_i} = tangent_corr;
        deltaT_s{file_i} = deltaT;
        sizeDeltaT(file_i) = length(deltaT);
        
        clear deltaT tangent_corr
   end
   
   minDeltaT = min(sizeDeltaT);
   
   for file_i = 1:length(contour_length)
       ntangent_corr{file_i} = tangent_corr_s{file_i}(1:minDeltaT);
       ndeltaT{file_i} = deltaT_s{file_i}(1:minDeltaT);
   end
   
   for file_i = 1:length(contour_length)
        [myfit{file_i} gof{file_i}] = fit(ndeltaT{file_i}', ntangent_corr{file_i}, myfittype,...
            'StartPoint', [1, 500]);

        Lp(file_i) = myfit{file_i}.b;
        rsq(file_i) = gof{file_i}.rsquare;
        sef(file_i) = gof{file_i}.rmse;
    
   end
   
   fit_data.Lp = Lp;
   fit_data.rsq = rsq;
   fit_data.sef = sef;
   
end