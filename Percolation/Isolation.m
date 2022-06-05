function [IsolationResults]  = Isolation(nx,ny,nz,phaseValue)
% Form matrix Z:

% Z - Matrix of values 0/1, marks non-isolated clusters beginning from the first image
% where  1- non-isolated phase;  0- other phases;
load Mic.mat

Z=zeros(nx,ny,nz);
Z(:,:,2)=((Y(:,:,1)==1) & (Y(:,:,2)==1)); % Non-isolated clusters between first and second images




%% Main Loop
% The code compares matrix 'Z' with the next image (matrix 'Y+1'). When in theese two
% images, woxels with the same coordinates are equal '1' in both arrays, code
% marks this woxel in next matrix 'Z' (Z+1) by '1'. In next step code checks
% if the cluster is growing: in the matrix 'Y+1' code checks cells in all
% directions if the cell is '1', the code includes it to matrix 'Z+1',
% In that case the process repeats for this cell.

    
% A - Array informs us about action that have been taken in woxels
    % 0 - don't consider this cell
    % 1 - include this cell to matrix 'Z'
    % 2 - check growth of this woxel in every direction.
            
% The code looks for value of '1' cells in array 'A' , and checks adjecent cells 
% in all directions, if they also equal '1' : 
% includes this adjecent cells to matrix 'Z+1', 
% changes values of adjecent cells for '1' or '0' in array A 
% changes value cell under consideration to '2'

perkolated = false;
for i=2:1:nz-1 
            
    A=zeros(nx,ny);  
        
        Z(:,:,i+1)=((Z(:,:,i)==1) & (Y(:,:,i+1)==1));  
        A(:,:)=Z(:,:,i+1);          
        
    while true                              
        [rows,cols] = find(A==1);      % Repeat procces until all cells are equal 2 or 0
        if isempty(rows)
            break
        end
        
        for k=1:length(rows)            
                                 
              Z(rows(k),cols(k),i+1)=1;
             
                if (cols(k)~=1) &&  A(rows(k),cols(k)-1)~=2 && (Y(rows(k),cols(k)-1,i+1)==1) 
                    Z(rows(k),cols(k)-1,i+1)=1;
                    A(rows(k),cols(k)-1)=1;
                end   

                if (rows(k)~=nx) && A(rows(k)+1,cols(k))~=2 &&  (Y(rows(k)+1,cols(k),i+1)==1) 
                    Z(rows(k)+1,cols(k),i+1)=1;
                    A(rows(k)+1,cols(k))  =1;
                end              
                
                if (rows(k)~=1) && A(rows(k)-1,cols(k))~=2 &&  (Y(rows(k)-1,cols(k),i+1)==1) 
                    Z(rows(k)-1,cols(k),i+1)=1;
                    A(rows(k)-1,cols(k))=1; 
                end
                
                if (cols(k)~=ny) && A(rows(k),cols(k)+1)~=2 && (Y(rows(k),cols(k)+1,i+1)==1) 
                    Z(rows(k),cols(k)+1,i+1)=1;
                    A(rows(k),cols(k)+1)=1;
                end
                
                A(rows(k),cols(k))=2; 
          end
          
     end
 

end

% Single initialization is not enough!
% Repeats the whole process taking into concideration the current matrix 'Z'
% beginning from the opposite side and repeating all actions as long as the
% matrix 'Z' changes.

Zz=zeros(nx,ny,nz);
C = Z(:)-Zz(:);

while sum(C)~=0
    C = Z(:)-Zz(:);
    Zz(:)=Z(:);  
sum(C)
    pfun=true;
    pfun2=true;
    [Z]=fun2(Z,Y,nx,ny,nz,pfun2);
    [Z]=fun(Z,Y,nx,ny,nz,pfun);
end

%% Creation of Images
Matrixq = CreationImages(nx,ny,nz,perkolated,phaseValue,Y,Z);

%% Calculation 
IsolationResults = Calculation(nx,ny,nz,perkolated,Y,phaseValue,B,Matrixq);
end