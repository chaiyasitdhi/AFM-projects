function map3d = tomograph(map2D, sizeX, sizeY)

for i = 1:length(map2D)
    maxm = length(map2D{1});
    if length(map2D{i}) > maxm
        maxm = length(map2D{i});
    end 
end

for itr = 1:maxm
    for i = 1:length(map2D)
        if length(map2D{i}) < itr
            Layer(i,1) = NaN;
        else
            Layer(i,1) = map2D{i}(itr);
        end
    end
 
    imfloor = reshape(Layer, sizeX, sizeY);
    imfloor(imfloor < 0 ) = NaN;
    map3d(:,:,itr) = imfloor;
end


end