function img_new_connected = load_image(photos_number, foldernameTortuosity )

filePattern = fullfile(foldernameTortuosity, '*.bmp'); 
theFiles = dir(filePattern);

for n=1:photos_number    
    
    baseFileName = theFiles(n).name;
    fullFileName = fullfile(foldernameTortuosity, baseFileName)   ;     % Changes the number of the picture into string variable; 
                                                        % Pictures end with 0000,0001.. indexes
    img(:,:,n) = double(imread(fullFileName))/255;            % Loads all the pictues into arrays (and divides them by 255)

end


img_new = img;
img_new(1,:,:) = 1;
img_new(end,:,:) = 1;
clear img

img_new_connected_to_z0 = bwselect3(img_new == 1, 1, 1, 1);
img_new_connected_to_zend = bwselect3(img_new == 1, 280, 280, 280);
%clear img_new
img_new_connected = img_new_connected_to_z0.*img_new_connected_to_zend;


save load_image.mat