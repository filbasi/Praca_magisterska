function [q] = CreationImages(nx,ny,nz,perkolated,phaseValue,Y,Z)
%% Creation Images
 
for i=1:1:nz
%%%% brak wcięcia wewnątrz pętli for ogranicza czytelność
y{i}=Y(:,:,i);
z{i}=Z(:,:,i);

%%%% q może być zwykłą tablicą o 3 wymiarach zamiast tablicą komórek, będzie
%%%% szybciej i czytelniej tzn. na końcu tylko zapisywanie zdjęć do plików będzie w pętli 
%%%% trzbea też będzie zmienić funkcję Calculation.m która traktuje q jako tablicę komórek
q{i}=zeros(nx,ny);
q{i}(y{i}==1)=255;   
q{i}(z{i}==1) = 133;

   %%%% nazewnictwo plików jest trochę niekonsekwentne
   %%%% zamiast w każdej instr. warunkowej wykonywać funkcję imwrite, można najpierw wyznaczyć ścierzkę i przypisać do zmiennej wewnątrz instrukcji 
   %%%% warunkowych a dopiero po nich wykonać imwrite z wyznaczonym argumentem  
   if phaseValue==0 && perkolated==true
       %%%% tworzenie folderu w każdym przebiegu pętli nie ma sensu, taka niepozorna funkcja może znacznie zwalniać kod, ze względu na komunikację z systemem
       %%%% Należy przenieść tworzenie folderów przed pętlę for 
       mkdir Percolation_Pore;
       imwrite(uint8(q{i}),sprintf("./Percolation_Pore/Black%04d.png", i));
   elseif phaseValue==0 && perkolated==false
       mkdir Isolated_Pore;
       imwrite(uint8(q{i}),sprintf("./Isolated_Pore/PorIso%04d.png", i)); % sprintf("path to save images", i)
   end;

   if phaseValue==127 && perkolated==true
        mkdir Percolation_YSZ;
        imwrite(uint8(q{i}),sprintf("./Percolation_YSZ/Grey%04d.png", i))
   elseif phaseValue==127 && perkolated==false
        mkdir Isolated_YSZ;
        imwrite(uint8(q{i}),sprintf("./Isolated_YSZ/Szary%04d.png", i))
   end;

    if phaseValue==255 && perkolated==true
        mkdir Percolation_Ni;
        imwrite(uint8(q{i}),sprintf("./Percolation_Ni/White%04d.png", i));
    elseif phaseValue==255 && perkolated==false 
        mkdir Isolated_Ni;
        imwrite(uint8(q{i}),sprintf("./Isolated_Ni/Bialy%04d.png", i));
    end

end

end