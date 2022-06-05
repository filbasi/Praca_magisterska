function [PercolationResults]  = Percolation(nx,ny,nz,phaseValue,B)
%% Main Loop

% The code doesn't compare first and second images because percolated is
% cluster in contact with the external surface of the microstructure. 

% The code creates matrices that represent clusters beginning from the each
% external surface of the microstructure, so started from first and last
% images. After that compares both matrices to get a percolated cluster.
load Mic.mat
perkolated = true;
%%%% poniższe linijki (do końca while) powtarzają się później z lekkimi zmianami, 
%%%% można zrobić z tego osobną funkcję, która wyznaczy Z w zalezności od parametrów pfun oraz pfun2
Zz=zeros(nx,ny,nz);
Z2=zeros(nx,ny,nz);
C = 1;
while sum(C)~=0
    pfun=false;
    pfun2=true;
    [Z2]=fun(Z2,Y,nx,ny,nz,pfun);   
    [Z2]=fun2(Z2,Y,nx,ny,nz,pfun2);
    
    %%%% w tych dwóch linijkach nie ma potrzeby dawać (:)
    C = Z2(:)-Zz(:);
    Zz(:)=Z2(:); 
    %%%% wypisywanie sum(C) można tutaj usunąć
    sum(C)
end

Zz=zeros(nx,ny,nz);
Z1=zeros(nx,ny,nz);
C = 1;

while sum(C)~=0
    pfun=true;
    pfun2=false;
    [Z1]=fun2(Z1,Y,nx,ny,nz,pfun2);
    [Z1]=fun(Z1,Y,nx,ny,nz,pfun);
    C = Z1(:)-Zz(:);
    Zz(:)=Z1(:); 
    sum(C)
end

for i=1:1:nz
Z(:,:,i)=((Z1(:,:,i)==1) & (Z2(:,:,i)==1));   
end

%% Creation Images
Matrixq = CreationImages(nx,ny,nz,perkolated,phaseValue,Y,Z);

%% Calculation
PercolationResults = Calculation(nx,ny,nz,perkolated,Y,phaseValue,B,Matrixq);
end