%function [result2,result3,result4,result5,result6,result7,result11,result12] = Tortuosity (photos_number,walkers_number,steps,foldernameTortuosity)
function [result11] = Tortuosity (photos_number,walkers_number,steps,foldernameTortuosity)
%#################################################################################################################%
%######################################                          #################################################%
%######################################     RANDOM WALKER 3D     #################################################%
%######################################                          #################################################%
%#################################################################################################################%
%##########################    made by Filip Lyczek/Adam Zajaczewski         #####################################%
%#################################################################################################################%
%##########################                  ver. 1.3.2                      #####################################%
%#################################################################################################################%
%{
clc                                                                                          %####################%
clear all                                                                                    %####################%
close all                                                                                    %####################%
                                                                                             %####################%
% Parameters                                                                                 %####################%
photos_number = 292;            % Number of used photos, aka the 3rd dimension               %####################%
walkers_number = 10000;        % Number of walkers released in one test                     %####################%
steps = 10000;                 % Number of steps of each and every walker                   %####################%
display = 2;                    % Option for displaying images;                              %####################%
                                % '0'-maximized windows on top of each other;                %####################%
                                % '1'-photos arranged in order;                              %####################%
                                % '2'-no display;                                            %####################%
boundary = 3;                   % Walking area by colour:                                    %####################%
                                % '1' - white                                                %####################%
                                % '2' - black                                                %####################%
                                % '3' - no boundary                                          %####################%
                                % This program is made for only two phases - black or white  %####################%
                                                                                             %####################%
%#################################################################################################################%
%#################################################################################################################%

%}
% Boundries setting array 
usr = memory
display = 2;
boundary = 3;

low_boundaries      = [0.96, 0.00, 0.00]; % Low and high boundries are values for chosing starting point
high_boundaries     = [1.00, 0.04, 1.00];
neighbor_boundaries = [0.01, 0.01, 1.00]; % Neighbour boundary is for possibility to change color by set value
l_bound = low_boundaries(boundary);
h_bound = high_boundaries(boundary);
global bound;
bound = neighbor_boundaries(boundary);

%t = cputime;
StartTimer = tic; %%% Pommiar czasu
StartLoadImages = tic; %%% Pommiar czasu
% Loads all the photo files

load load_image.mat;
% filePattern = fullfile(foldernameTortuosity, '*.bmp'); 
% theFiles = dir(filePattern);
% img = zeros(280,280);
% %%%%%%%% dodac jako mic.mat, i wyciagnac z tej f.
% for n=1:photos_number    
%     
%     baseFileName = theFiles(n).name;
%     fullFileName = fullfile(foldernameTortuosity, baseFileName)   ;     % Changes the number of the picture into string variable; 
%                                                         % Pictures end with 0000,0001.. indexes
%     img(:,:,n) = double(imread(fullFileName))/255;            % Loads all the pictues into arrays (and divides them by 255)
% 
% end

%%%%%%%%
StopLoadImages = toc(StartLoadImages) %%% Pommiar czasu
StartTimer2 = tic; %%% Pommiar czasu
%[Y,X,Z]=size(img);
[Y,X,Z]=size(img_new_connected);
img = img_new_connected;

% for n=1:Z
% 
%     if display ~= 2
%         figure (n)
%         if display == 0
%             % Command to maximize each photo and highlight only the one in use
%             set(n,'units','normalized','outerposition',[0 0 1 1])
%         elseif display == 1
%             % Mini script to arrange photos on the screen in order
%             autoArrangeFigures();
%         else 
%             disp('Choose 0,1 or 2 for display parameter')
%             break;
%         end
% 
%         imagesc(img(:,:,n));
%         colormap('gray')
%         axis equal
%         axis tight
%         axis off
%         hold on
%         plot(1,1,'y.','MarkerSize',9) % Example of starting point (yellow dot)
%         plot(X,Y,'r.','MarkerSize',9) % Example of finishing point (red dot)
%         name2 = sprintf('File: %s',name); 
%         title(name2)
%     end
%     
% end

StopTimer2 = toc(StartTimer2) %%% Pommiar czasu
StartPlaceWalkers = tic; %%% Pommiar czasu
% Picks random starting places for walkers (random in middle area of extended picture)

rx = zeros(1,walkers_number);
ry = zeros(1,walkers_number);
rz = zeros(1,walkers_number);
rx2 = zeros(1,walkers_number);
ry2 = zeros(1,walkers_number);
rz2 = zeros(1,walkers_number);
xp = zeros(1,walkers_number);
yp = zeros(1,walkers_number);
zp = zeros(1,walkers_number);

for k=1:walkers_number

   x=randi([1,X]);                                   % Or you can set parameters to '1' so walkers will start from the starting green point
   y=randi([1,Y]);                                   % x=1; y=1;
   z=randi([1,Z]);                                   % z=1; % Starts from the first picture
   
   while img(y,x,z)<l_bound || img(y,x,z)>h_bound    % Walkers must be placed (and move) within examined area
        x=randi([1,X]); 
        y=randi([1,Y]); 
        z=randi([1,Z]); 
   end
                                                     % Chosen starting points for walkers

     rx(k) = x; 
     ry(k) = y; 
     rz(k) = z; 
     rx2(k) = x; 
     ry2(k) = y; 
     rz2(k) = z; 
     xp(k) = x;                                      % Starting points - saved for further calculations
     yp(k) = y;
     zp(k) = z;
end

Lx=zeros(walkers_number,1);                         % Creates indicators to track current location in virtual mirrored structure
Ly=zeros(walkers_number,1);
Lz=zeros(walkers_number,1);

% Walking - green for the trace, blue for the current position
StopPlaceWalkers = toc(StartPlaceWalkers) %%% Pommiar czasu
StartRS3d = tic; %%% Pommiar czasu

%xd = [];
%xd = zeros(1,9*steps*walkers_number);
% xd1 = zeros(steps, 9*walkers_number);
% xd1 = zeros(1,3*walkers_number);

% if display == 0 || display == 1
%     
%     % Loop 1 - normal random and drawing
%     
%     for i=1:steps
%         for j=1:walkers_number
%             figure (rz(j))
%             plot(rx(j),ry(j),'g.','MarkerSize',9)
%             [Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)] = randomstep3d(Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),img,rx(j),ry(j),rz(j));
%             figure (rz(j))
%             plot(rx(j),ry(j),'b.','MarkerSize',13)        
%         end
%         current_step = i; %shows progres of random walking
%     end
    
if  display == 2
    
    % Loop 2 - normal random without drawing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      for i=1:steps
%           for j=1:walkers_number
%             [Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)] = randomstep3d(Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),img,rx(j),ry(j),rz(j));
%             xd = [xd,Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)];
% %      for i=1:steps
% %           for j=1:walkers_number   
% %             [Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)] = randomstep3d(Lx(i),Ly(i),Lz(i),rx2(i),ry2(i),rz2(i),img,rx(i),ry(i),rz(i));
%           end
%         current_step = i %shows progres of random walking
%      end
%%%%%%%%%%%%%%%%%%%%%TEST PARALLELIZACJI%%%%%%%%%%%%%%%%%%%%%%%%
steps_MSD = zeros(steps,1);
    for i=1:steps
        steps_ave_MSD = 0;
          for j=1:walkers_number
            [Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)] = randomstep3d(Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),img,rx(j),ry(j),rz(j));
%            xd = [xd,Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)];
%             xd1(i,9*j-8:9*j) = [Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),rx(j),ry(j),rz(j)];
%            [rx2(j),ry2(j),rz2(j)] = randomstep3d(Lx(j),Ly(j),Lz(j),rx2(j),ry2(j),rz2(j),img,rx(j),ry(j),rz(j));
%            xd1(1,3*j-2:3*j) = [rx2(j),ry2(j),rz2(j)];

            steps_ave_MSD = steps_ave_MSD + (rx2(j)-xp(j)).^2+(ry2(j)-yp(j)).^2+(rz2(j)-zp(j)).^2;

          end
          steps_MSD(i,1) = steps_ave_MSD/walkers_number;
        current_step = i; %shows progres of random walking
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kontrukscja 3d tablicy
%xd2 = reshape(xd, [9,round(walkers_number),round(steps)]);
% xdd = reshape(xd1.', [1, 9*walkers_number*steps]);
% xd2 = reshape(xdd, [9,round(walkers_number),round(steps)]);
% xdd = reshape(xd1.', [1, 3*walkers_number*steps]);
% xd2 = reshape(xd1, [3,round(walkers_number),1]); % CHECK THIS IF SOMETHING NOT WORKING

else 
    disp('Something went wrong')      
end
StopRS3d = toc(StartRS3d) %%% Pommiar czasu
StartMSD = tic; %%% Pommiar czasu
% Starting points marking

% if display == 0 || display == 1
%     for j=1:walkers_number
%         figure (zp(j))
%         plot(xp(j),yp(j),'m.','MarkerSize',13) % Pink dots
%     end
% end

fprintf('\n\n');

% Mean square displacement MSD
% suma=0; sumax=0; sumay=0; sumaz=0; sumadt = 0; suma1=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rt = [];
% for nu=1:steps
%     rt2 = [];
%     for n=1:walkers_number 
%         rt2 = [rt2, (xd2(4,n,nu)-xp(n)).^2+(xd2(5,n,nu)-yp(n)).^2+(xd2(6,n,nu)-zp(n)).^2];
%     end
%     rt = [rt, sum(rt2) / walkers_number];
% 
% end
% wykres = rt;

%%%%%%%%%%%%%%%PARALLELIZATION TESTING%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%% wersja testowa oraz wszystko poniżej do usuniecia %%%%%%%%%%%%%
% rt = [];
% 
% xd24 = xd2(4,:,:);
% xd25 = xd2(5,:,:);
% xd26 = xd2(6,:,:);
% 
% 
% for n=1:walkers_number
%     rt2 = [];
%     for nu=1:steps 
%         rt2 = [rt2, (xd24(:,n,nu)-xp(n)).^2+(xd25(:,n,nu)-yp(n)).^2+(xd26(:,n,nu)-zp(n)).^2];
%     end
%     rt = [rt, rt2];
% 
% end
% 
% rt3 = reshape(rt, [steps, walkers_number]);
% rt4 = sum(rt3,2)/walkers_number;
% rt4 = reshape(rt4,[1,steps]);
% wykres = rt4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  PARALLELIZATION TESTING v2   %+%%%%%%%%%%%%%%%%%%%%%%

%rt = zeros(1,steps*walkers_number);

% rt = zeros(walkers_number,1);
% 
% % xd24 = xd2(4,:,:);
% % xd25 = xd2(5,:,:);
% % xd26 = xd2(6,:,:);
% xd24 = xd2(1,:,:);
% xd25 = xd2(2,:,:);
% xd26 = xd2(3,:,:);
% 
% for n=1:walkers_number
%     rt2 = zeros(1,steps);
%     for nu=1:steps 
%         rt2(nu) = (xd24(:,n,1)-xp(n)).^2+(xd25(:,n,1)-yp(n)).^2+(xd26(:,n,1)-zp(n)).^2;
%     end
% %     rt = [rt, sum(rt2) / walkers_number];
% %     rt(n) = rt2
%     %rt(1,1+(n-1)*steps:n*steps) = rt2  ;
%     rt(n,:) = rt2  ;
% 
% end
% 
% clear xdd xd2
% clear xd24 xd25 xd26
% 
% rt1d = reshape(rt.',1,[]);
% clear rt
% rt3 = reshape(rt1d, [steps, walkers_number]);
% clear rt1d
% rt4 = sum(rt3,2)/walkers_number;
% clear rt3
% wykres = reshape(rt4,[1,steps]);
StopMSD = toc(StartMSD) %%% Pomiar czasu
StopTimer = toc(StartTimer) %% Pomiar czasu

%%%%% Create and save data into new txt file %%%%%%% NEW %%%%%%%

% FileName = fullfile('..', 'dmr-quant-micros','ParallelData');
% %FileName = 'C:\studia\Praca_magisterska\dmr-quant-micros\ParallelData';
% %[fPath, fName, fExt] = fileparts(FileName);
% fDir = dir(fullfile(FileName, '*.txt'));
% fDirr = length(fDir);
% fDirrr = int2str(fDirr);
% 
% if exist(fullfile(FileName, 'Parallel.txt'), 'file')
%     fid = fopen(fullfile(FileName, strcat('Parallel',fDirrr, '.txt')), 'at');
%     fprintf(fid, '%f ', rt4);
%     fclose(fid);
% else
%     fid = fopen(fullfile(FileName, "Parallel.txt"), 'at');
%     fprintf(fid, '%f ', rt4);
%     fclose(fid); 
% end

% clear rt4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Saveing data into txt file %% OLD %%%

% %rt44 = rt4.';
% % fid = fopen("ParallelData.txt", "at");
% % fprintf(fid, '%f \n' , rt4);
% % fclose(fid);
% name = fullfile("ParallelData/", "ParallelData");
% name1 = nextname(name,'_1','.txt', true)
% %name = nextname('ParallelData','_1','.txt')
% fid = fopen(name1,'wt');
% fprintf(fid, '%f ' , rt4);
% fclose(fid);

%%%%%%%%%%Convert data to array

% projectdir = fullfile(".\ParallelData/");
% dinfo = dir( fullfile(projectdir, '*.txt'));
% nfiles = length(dinfo);
% assert(nfiles == 10, 'expected 10 files');
% filenames = fullfile(projectdir, {dinfo.name});
% data = zeros(nfiles, 10000);
% for K = 1 : nfiles
%     thisfile = filenames{K};
%     thisdata = load(thisfile);
%     data(K, :) = thisdata;
% end
% %columnMeans = mean(data, 1);
% wykres = mean(data, 1);


%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j=1:walkers_number
%     r2(j) = (rx2(j)-xp(j)).^2 + (ry2(j)-yp(j)).^2 + (rz2(j)-zp(j)).^2; % msd step1
%     
%     r2x(j) = (rx2(j).^2)-(xp(j).^2);% msd step1,  x direction
%     r2y(j) = (ry2(j).^2)-(yp(j).^2);
%     r2z(j) = (rz2(j).^2)-(zp(j).^2);
% 
%     r2x(j)= (rx2(j)-xp(j)).^2;
%     r2y(j) = (ry2(j)-yp(j)).^2;% msd step1,  y direction   % directional mean square displacements
%     r2z(j) = (rz2(j)-zp(j)).^2;% msd step1,  z direction
%     
%     r2(j) = r2x(j)+r2y(j)+r2z(j);
%     
%     rtrue(j)= r2(j)/walkers_number;
%     dt = diff(rtrue);
%     
%     result = ['Walker no. ',num2str(j),' scored: ',num2str(r2(j))];
%     disp(result)
%     
%     suma=suma+r2(j);    % summed msd for all walkers
%     suma1=suma1+(r2(j)/walkers_number);
%     wykres(j)=suma1;
%     sumax=sumax+r2x(j); % summed msd in single dimension; for each dimension: x,y,z
%     sumay=sumay+r2y(j);
%     sumaz=sumaz+r2z(j); 
%    
% end

%wykres = rt;
%z=1:walkers_number;
z=1:steps;
plot(z,steps_MSD,'r-','LineWidth',2);
%ylim([0 inf]);
ylabel('\it MSD \rm[-]');
xlabel('\it \vartheta \rm[-]');
aprox = polyfit(z,steps_MSD,1);
% Evaluate fit equation using polyval
y_est = polyval(aprox,z);
% Add trend line to plot
hold on 
plot(z,y_est,'black--','LineWidth',2)
hold off
legend('Funkcja rzeczywista', 'Prosta aproksymujaca')
txt1 = ['y = ' num2str(aprox(1)) 'x + ' num2str(aprox(2))];
text(x, y, txt1);

% msd = suma/walkers_number;    % average; msd 
%                  (j = 'walkers_number') at that point
% msdx = sumax/walkers_number;  
% msdy = sumay/walkers_number;  % average msds in x,y and z dimension
% msdz = sumaz/walkers_number;
% 
% 
%D = msd/(6*steps);     % Mean square displacement to time; D parameter!
%                        %Steps represent time                       
% Dx = msdx/(2*steps); 
% Dy = msdy/(2*steps);   % msd for x, y and z dimension
% Dz = msdz/(2*steps);
% 
% a2 = 1; % woksel

% sumadt = sum(dt);
% Tautrue=a2/sumadt;

Taugraph=1/aprox(1);
%Taugraphxyz=Taugraph/3;
% Results
% result2 = ['Score represents square distance between starting and finishing point. '];
% result3 = ['Average Walker score (MSD) is ', num2str(msd),'. Overall number of steps was ', num2str(steps),'. '];
% result4 = ['Mean Square Displacement to Time (D) is ', num2str(D),'. '];

% result5 = ['Mean Square Displacement (Dx) for x-axis is ', num2str(Dx),'. '];
% result6 = ['Mean Square Displacement (Dy) for y-axis is ', num2str(Dy),'. '];
% result7 = ['Mean Square Displacement (Dz) for z-axis is ', num2str(Dz),'. '];

% result8 = ['Tortuosity factor is ', num2str(Tautrue),'. '];
% result9 = ['Tortuosity factor in x,y,z direction is ', num2str(Tauxyz),'. '];
% result10 = ['Tortuosity factor without abs is ', num2str(Tautrue),'. '];

result11 = ['Tortuosity factor calculated from graph is ', num2str(Taugraph),'. '];
%result12 = ['Tortuosity factor calculated from graph in x,y,z direction is ', num2str(Taugraphxyz),'. '];

% disp(result2)
% fprintf('\n\n\n');
% disp('~~Results~~')
% disp(result3)
% disp(result4)
% fprintf('\n\n');
% disp(result5)
% disp(result6)
% disp(result7)

% disp(result8)
% disp(result9)
% disp(result10)

disp(['Trend line equation is y = ' num2str(aprox(1)) '*x + ' num2str(aprox(2))])
disp(result11)
%disp(result12)

%e = cputime-t
usr = memory
whos