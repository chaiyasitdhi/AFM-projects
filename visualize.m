% for layer = 1:size(elastic_3d, 3)
%     figure()
%     for i = 1:size(elastic_3d,1) 
%         for k = 1:size(elastic_3d,2)
%             if elastic_3d(i,k,layer) > 3000;
%                 plot3(i,k,NaN);
%             else
%                 plot3(i,k,elastic_3d(i,k,layer), 'o', 'MarkerSize', 5);
%                 hold on;
%             end
%         end
%     end
%     xlabel('Scan Size (Pixel)');
%     ylabel('Scan Size (Pixel)');
%     zlabel('Elasticity (kPa)');
% end
% 

% for layer = 1:size(error_3d, 3)
%     data = error_3d(:,:,layer);
%     data_linear = reshape(data,1,32*32);
%     data_linear = data_linear(data_linear<3000);
%     [Y(:,layer) x(:,layer)] = hist(data_linear,100);
%     meanDat(layer) = mean(data_linear);
%     stdDat(layer) = std(data_linear);
% end

% for i = 1:size(tt,1)
%     for k = 1:size(tt,2)
%         depth = i*5;
%         plot3(k,i, tt(k,i));
%     end
% end
