function img = load_image(photos_number, foldernameTortuosity )

filePattern = fullfile(foldernameTortuosity, '*.bmp'); 
theFiles = dir(filePattern);

%%%%%%%% dodac jako mic.mat, i wyciagnac z tej f.
for n=1:photos_number    
    
    baseFileName = theFiles(n).name;
    fullFileName = fullfile(foldernameTortuosity, baseFileName)   ;     % Changes the number of the picture into string variable; 
                                                        % Pictures end with 0000,0001.. indexes
    img(:,:,n) = double(imread(fullFileName))/255;            % Loads all the pictues into arrays (and divides them by 255)

end
%%%%%%%%
save load_image.mat