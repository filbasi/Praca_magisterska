function [GrainSizeResult] = AverageGrainSize(foldername,direction)
%% AGI main script
% Author: Piotr Zurek
% Date of cerate: 08.12.2019
% Up-date: 08.12.209
%% Comand window
%% Inputs
% voxel size nm
%%%% Ta wartośc na pewno nie jest dobra, dodatkowo powinna być parametrem, więcej informacji w general_comments_SB.txt
v = 10;
% Folder path
% e.g. foldername = 'E:\Programy\Matlab2019b\Pore_1.bmp';
%foldername = 'C:\Users\Piotr\Desktop\SOFC\Percolation_Ni\White0001.png';
%%%% ta linijka wykonywana jest na potrzeby komend ok. 15 linijek niżej, więc tampowinna się znaleźć
%[filepath,name,ext] = fileparts(foldername);
% AGI
    [aveAGI,dev] = secantmethode(foldername,v,direction);
%% Results
R = [aveAGI,dev];
%%%% Niezby czytelny format wyników, lepiej chyba byłoby zwracać dwie tablice, jedną z wartościami rozmiarów ziaren R, drugą z odchyleniami devR 
GrainSizeResult = table(R(1:end,1), R(1:end,2),...
    'VariableNames', { 'GrainSize','AGIStandardDeviation'});
%%%% zapisywanie wyników nie jest funkcjonalnością funkcji AverageGrainSize, więc nie powinno się tutaj znajdować
%writetable(T, num2str(['AGI from ',name,'.txt']))
