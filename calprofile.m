imagesc(fnimage); 
clear h ih
for i = 1:30
    h(i) = {improfile(10)};
    ih(:,i) = h{i};
end 