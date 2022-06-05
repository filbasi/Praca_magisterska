function [Lxnew,Lynew,Lznew,dxnew,dynew,dznew,pxnew,pynew,pznew]=randomstep3d(Lx,Ly,Lz,dx,dy,dz,img,px,py,pz)
    global bound;
    
    % Taking value of current location
    
    [Y,X,Z]=size(img);
    curr=img(py,px,pz);
    
    % Checking of surrounding pixels
    
    if py>1 
        u=img(py-1,px,pz);
    end 
    if py==1
        u=curr;
    end          % Gives freedom to cross the boundary of sample
    if py<Y
        d=img(py+1,px,pz);
    end
    if py==Y
        d=curr;
    end          % Gives freedom to cross the boundary of sample
    if px>1
        l=img(py,px-1,pz);
    end
    if px==1
        l=curr;
    end          % Gives freedom to cross the boundary of sample
    if px<X
        r=img(py,px+1,pz);
    end
    if px==X
        r=curr;
    end           % Gives freedom to cross the boundary of sample
    if pz>1
        f=img(py,px,pz-1);
    end
    if pz==1
        f=curr;
    end           % Gives freedom to cross the boundary of sample
    if pz<Z
        b=img(py,px,pz+1);
    end
    if pz==Z
        b=curr;
    end

    % Comparing current value with surrounding pixels and estimating
    % probability to jump into them
    
    if abs(curr-u)>bound 
        prob_u=0;
    else
        prob_u=1;
    end

    if abs(curr-d)>bound
        prob_d=0;
    else
        prob_d=1;
    end
    
    if abs(curr-l)>bound
        prob_l=0;
    else
        prob_l=1;
    end
    
    if abs(curr-r)>bound
         prob_r=0;
    else 
        prob_r=1;
    end
    
    if abs(curr-f)>bound
        prob_f=0;
    else 
        prob_f=1;
    end
    
    if abs(curr-b)>bound
        prob_b=0;
    else 
        prob_b=1;
    end
    
    s = prob_u + prob_d + prob_l + prob_r + prob_b + prob_f;
    if s~=0
    prob_u=prob_u/s;
    prob_d=prob_d/s;
    prob_r=prob_r/s;
    prob_l=prob_l/s;
    prob_f=prob_f/s;
    prob_b=prob_b/s;
    
    % Creating array of probability and taking one random of possible ones
    
    prob=[prob_u prob_d prob_l prob_r prob_f prob_b];
    prob2=find(prob > 0);
    x1=randi(length(prob2));
    x=prob2(x1);
    
 if x==1 || x==2                  % Checking direction of drawn step
   
  if Ly==0                        % Checking current sector of virtual mirrored structure
    if py==Y && x==2              % Checking if current positions is on the edge of the sector 
        Ly=Ly+1;
        %py=py;
        dy=dy+1;
    elseif py==1 && x==1
        Ly=Ly-1;
        %py=py;
        dy=dy-1;
    elseif x==2
        
        py=py+1;
        dy=dy+1;
    elseif x==1
        
        py=py-1;
        dy=dy-1;
    end
 
  elseif Ly>0 && abs(rem(Ly,2))==1
    if py==Y && x==2
        Ly=Ly-1;
        %py=py;
        dy=dy-1;
    elseif py==1 && x==1
        Ly=Ly+1;
        %py=py;
        dy=dy+1;
    elseif x==2
        
        py=py+1;
        dy=dy-1;
    elseif x==1
        
        py=py-1;
        dy=dy+1;
    end
    
  elseif Ly<0 && abs(rem(Ly,2))==1
    if py==Y && x==2
        Ly=Ly-1;
        %py=py;
        dy=dy-1;
    elseif py==1 && x==1
        Ly=Ly+1;
        %py=py;
        dy=dy+1;
    elseif x==2
        
        py=py+1;
        dy=dy-1;
    elseif x==1
        
        py=py-1;
        dy=dy+1;
    end 

  elseif Ly>0 && abs(rem(Ly,2))==0
    if py==Y && x==2 
        Ly=Ly+1;
        %py=py;
        dy=dy+1;
    elseif py==1 && x==1
        Ly=Ly-1;
        %py=py;
        dy=dy-1;
    elseif x==2
        
        py=py+1;
        dy=dy+1;
    elseif x==1
        
        py=py-1;
        dy=dy-1;
    end
    
  elseif Ly<0 && abs(rem(Ly,2))==0
    if py==Y && x==2
        Ly=Ly+1;
        %py=py;
        dy=dy+1;
    elseif py==1 && x==1
        Ly=Ly-1;
        %py=py;
        dy=dy-1;
    elseif x==2
        
        py=py+1;
        dy=dy+1;
    elseif x==1
        
        py=py-1;
        dy=dy-1;
    end
  end  
end
 
 if x==3 || x==4
     
  if Lx==0
    if px==X && x==4
        Lx=Lx+1;
        %px=px;
        dx=dx+1;
    elseif px==1 && x==3
        Lx=Lx-1;
        %px=px;
        dx=dx-1;
    elseif x==4
        %Lx=Lx;
        px=px+1;
        dx=dx+1;
    elseif x==3
        %Lx=Lx;
        px=px-1;
        dx=dx-1;
    end

  elseif Lx>0 && abs(rem(Lx,2))==1
    if px==X && x==4
        Lx=Lx-1;
        %px=px;
        dx=dx-1;
    elseif px==1 && x==3
        Lx=Lx+1;
        %px=px;
        dx=dx+1;
    elseif x==4
        px=px+1;
        dx=dx-1;
    elseif x==3
        px=px-1;
        dx=dx+1;
    end

  elseif Lx<0 && abs(rem(Lx,2))==1
    if px==X && x==4
        Lx=Lx-1;
        %px=px;
        dx=dx-1;
    elseif px==1 && x==3
        Lx=Lx+1;
        %px=px;
        dx=dx+1;
    elseif x==4
        px=px+1;
        dx=dx-1;
    elseif x==3
        px=px-1;
        dx=dx+1;
    end 

  elseif Lx>0 && abs(rem(Lx,2))==0
    if px==X && x==4 
        Lx=Lx+1;
        %px=px;
        dx=dx+1;
    elseif px==1 && x==3
        Lx=Lx-1;
        %px=px;
        dx=dx-1;
    elseif x==4
        px=px+1;
        dx=dx+1;
    elseif x==3
        px=px-1;
        dx=dx-1;
    end

  elseif Lx<0 && abs(rem(Lx,2))==0
    if px==X && x==4
        Lx=Lx+1;
        %px=px;
        dx=dx+1;
    elseif px==1 && x==3
        Lx=Lx-1;
        %px=px;
        dx=dx-1;
    elseif x==4
        px=px+1;
        dx=dx+1;
    elseif x==3
        px=px-1;
        dx=dx-1;
    end 
  end
end
       
 if x==5 || x==6
        
   if Lz==0
    if pz==Z && x==6
        Lz=Lz+1;
        %pz=pz;
        dz=dz+1;
    elseif pz==1 && x==5
        Lz=Lz-1;
        %pz=pz;
        dz=dz-1;
    elseif x==6
        
        pz=pz+1;
        dz=dz+1;
    elseif x==5
        
        pz=pz-1;
        dz=dz-1;
    end
 
   elseif Lz>0 && abs(rem(Lz,2))==1
    if pz==Z && x==6
        Lz=Lz-1;
        %pz=pz;
        dz=dz-1;
    elseif pz==1 && x==5
        Lz=Lz+1;
        %pz=pz;
        dz=dz+1;
    elseif x==6
        
        pz=pz+1;
        dz=dz-1;
    elseif x==5
        
        pz=pz-1;
        dz=dz+1;
    end

   elseif Lz<0 && abs(rem(Lz,2))==1
    if pz==Z && x==6
        Lz=Lz-1;
        %pz=pz;
        dz=dz-1;
    elseif pz==1 && x==5
        Lz=Lz+1;
        %pz=pz;
        dz=dz+1;
    elseif x==6
        
        pz=pz+1;
        dz=dz-1;
    elseif x==5
        
        pz=pz-1;
        dz=dz+1;
    end 
    
   elseif Lz>0 && abs(rem(Lz,2))==0
    if pz==Z && x==6 
        Lz=Lz+1;
        %pz=pz;
        dz=dz+1;
    elseif pz==1 && x==5
        Lz=Lz-1;
        %pz=pz;
        dz=dz-1;
    elseif x==6
        
        pz=pz+1;
        dz=dz+1;
    elseif x==5
        
        pz=pz-1;
        dz=dz-1;
    end
    
   elseif Lz<0 && abs(rem(Lz,2))==0
    if pz==Z && x==6
        Lz=Lz+1;
        %pz=pz;
        dz=dz+1;
    elseif pz==1 && x==5
        Lz=Lz-1;
        %pz=pz;
        dz=dz-1;
    elseif x==6
        
        pz=pz+1;
        dz=dz+1;
    elseif x==5
        
        pz=pz-1;
        dz=dz-1;
    end  
  end
 end

    pynew=py;
    pxnew=px;
    pznew=pz;
    dxnew=dx;
    dynew=dy;
    dznew=dz;
    Lxnew=Lx;
    Lynew=Ly;
    Lznew=Lz;
    else
        
    pynew=py;
    pxnew=px;
    pznew=pz;
    dxnew=dx;
    dynew=dy;
    dznew=dz;
    Lxnew=Lx;
    Lynew=Ly;
    Lznew=Lz;
    end
        
end