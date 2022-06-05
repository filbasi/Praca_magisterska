function [Y] = LoadImages(filePath, fileDataType, nx, ny, nz, phaseValue, direction, B)


%% Images Measurements 
%nx=250; ny=250; numberOfFiles=250;  
%nz=numberOfFiles;


%B = readRaw('./Microstructure.raw',nx, ny, nz,'int');  


%% Load Raw

% loaded images transform in matrix of values 0 or 1
% where: 0- phase under consideration;  1- other phases
%%%% w tym miejscu instrukcje warunkowe są całkowicie niepotrzebne. Proszę się zastanowić jak zamienić poniższe 60 linijek na 20
  if phaseValue == 0   
        %%%% Najpierw tworzone jest Y0 tylko po to, żeby być przepisane do Y, nie ma takiej potrzeby
        Y0 = double(B==0);       %Change value to select the phase to be consideration! 
                                %for example: ==0 Black/Pore; ==128 Grey/YSZ; ==255 White/Ni;
            Y = Y0;            
            %%%% mkdir Photos jest wykonywane w każdym przypadku, więc nie powinno znajdować się wewnątrz instrukcji warunkowych
            if direction == 'x'
                mkdir Photos
                for i = 1:nz
                imwrite(Y(:,:,i),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'y'
                mkdir Photos
                for i = 1:ny
                imwrite(reshape(Y(:,i,:), nx, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'z'
                mkdir Photos
                for i = 1:nx
                imwrite(reshape(Y(i,:,:), ny, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            end
  elseif phaseValue == 127
        Y127 = double(B==127);   %Change value to select the phase to be consideration!
                                %for example: ==0 Black/Pore; ==128 Grey/YSZ; ==255 White/Ni;
            Y = Y127;
            if direction == 'x'
                mkdir Photos
                for i = 1:nz
                imwrite(Y(:,:,i),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'y'
                mkdir Photos
                for i = 1:ny
                imwrite(reshape(Y(:,i,:), nx, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'z'
                mkdir Photos
                for i = 1:nx
                imwrite(reshape(Y(i,:,:), ny, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            end                               
 elseif phaseValue == 255 
    Y255 = double(B==255);       %Change value to select the phase to be consideration!     
                                %for example: ==0 Black/Pore; ==128 Grey/YSZ; ==255 White/Ni;
            Y = Y255;                    
            if direction == 'x'
                mkdir Photos
                for i = 1:nz
                imwrite(Y(:,:,i),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'y'
                mkdir Photos
                for i = 1:ny
                imwrite(reshape(Y(:,i,:), nx, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            elseif direction == 'z'
                mkdir Photos
                for i = 1:nx
                imwrite(reshape(Y(i,:,:), ny, nz ),['./Photos/Photo',num2str(i, '%04.f'),'.bmp'])
                end
            end
  end 
                         
save Mic.mat
end