function [ResultCalculation,Ni,YSZ,Pore]  = Calculation(nx,ny,nz,perkolated,Y,phaseValue,B,Matrixq)

%% Calculation
%%%% nie ma potrzeby dodatkowej zmiennej q, skoro istnieje Matrixq, zbędna linijka
q = Matrixq
per=0;
%%%% all jest funkcją wbudowaną matlaba, zmienić nazwę
all=0;
iso=0;

%%%% percolated 
%%%% wykorzystanie instrukcji warunkowej nie ma zbytnio sensu, w obu przypadkach zliczanie jest 
%%%% zrobione w prawie ten sam sposób. Nie do końca rozumiem wykorzystywanie zmiennej percolated, skoro w trakcie wykonywania funkcji Percolation 
%%%% wpływa ona tylko na miejsce zapisu wyznaczonych tablic. W każdym razie zliczanie raczej można zrobić poza if
if perkolated==true
    for i=1:1:nz
        for k=1:1:nx
            for w=1:1:ny
                if q{i}(w,k)==255
                per=per+1;
                end
            end
        end
    end
    
    all=sum(Y(:) == 1);

    %%%% instrukcje warunkowe całkowicie niepotrzebne
    % percentage of percolated phase 
    if phaseValue == 255
        ResultCalculation=(all-per)/all
    elseif phaseValue == 0
        ResultCalculation=(all-per)/all
    elseif phaseValue == 127
        ResultCalculation=(all-per)/all
    end

elseif perkolated==false
    for i=2:1:nz-1
        for k=1:1:nx
            for w=1:1:ny

                if q{i}(w,k)==255
                iso=iso+1;
                end
            end
        end
    end
  all=sum(Y(:) == 1);

% percentage of isolated phase 
    if phaseValue == 255
        ResultCalculation=(iso)/all
    elseif phaseValue == 0
        ResultCalculation=(iso)/all
    elseif phaseValue == 127
        ResultCalculation=(iso)/all
    end  
 
end
    

c=sum(B(:) == 127);
d=sum(B(:) == 0);
e=sum(B(:) == 255);


Ni= e/(c+d+e) % percentage of nickel phase 
YSZ=c/(c+d+e) % percentage of YSZ phase 
Pore=d/(c+d+e) % percentage of pore 
end