function [ flipimg ] = flipimage( img )

imgHeight = sqrt(size(img,1));

for i = 1:imgHeight
 flipimg(i,:) = img( ( imgHeight*i-(imgHeight-1) ) : (imgHeight*i) );
    
%     if logical(mod(i,2))
%         img(:,i) = fliplr(img(:,i));
%     end
%         
end
      
end

