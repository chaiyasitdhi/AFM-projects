clear concave_dat convex_dat

orient = 1; % 1 for horizontal, 2 for vertical.
concave = 2; % 1 for left ), 2 for right (.

if size(elastic_3d, 1) == size(elastic_3d, 2)
    lenDat = size(elastic_3d, 1);
else
    lenDat = 3;
end

concave_dat(1) = NaN;
convex_dat(1) = NaN;

for i = 1:lenDat
    if orient == 1
        height = Rconmap(:,i);
        edat = elastic_3d(:,i);
    else
        height = Rconmap(i,:);
        edat = elastic_3d(i,:);
    end
    
    coord = find(height > 0.5*10^-6 + min(height));
    top_coord = find(height == max(height));
    
    if concave == 1
        concave_dat(end+1:end+length(edat(coord(1):top_coord))) = edat(coord(1):top_coord);
        convex_dat(end+1:end+length(edat(top_coord+1:end))) = edat(top_coord+1 : end);
    else
        convex_dat(end+1:end+length(edat(top_coord+1:end))) = edat(top_coord+1 : end);
        concave_dat(end+1:end+length(edat(coord(1):top_coord))) = edat(coord(1):top_coord);
    end
    
end 
    concave_dat = concave_dat(2:end);
    convex_dat = convex_dat(2:end);
   
    
    