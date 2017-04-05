clear result
clear resultN

num = 1;
dimen = 1; % 1 horizontal 2 vertical
for i = 1:3
%for i = num:size(Rconmap,2)-num
    
    if dimen == 1
        data = Rconmap(:,i);
        max_data = find(data == max(data));
        max_data = max_data(1);
    
        edat = elastic_3d(:,i);
                
    else
        data = Rconmap(i,:);
        max_data = find(data == max(data));
        max_data = max_data(1);
    
        edat = elastic_3d(i,:);
       
    end
    
    if max_data == 1 
       result(i-num+1,:) = [NaN, edat(max_data:max_data+1)'];
    elseif max_data == size(Rconmap,2)
       result(i-num+1,:) = [edat(max_data-1:max_data)', NaN];
    else
       result(i-num+1,:) = edat(max_data-1:max_data+1);
    end

end

resultN = reshape(result, 1, size(result,1)*size(result,2));

    