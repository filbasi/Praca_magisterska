clear;
clc;
tic;
%%%% Przyałyby się komentarze oraz możliwe wartości (np dla direction i phase value)
filePath = 'Microstructure.raw'
nx = 280;
ny = 280;
nz = 280;
fileDataType='int32'
phaseValue = 127
direction = 'z'
%%%% Tutaj zamiast zmiennych filePath oraz fileDataType jest wpisany string
B = readRaw('./Microstructure.raw',nx, ny, nz,'int');
%%%% parametry dla różnych programów można oddzielić, np.
%%%% %% Tortuosity algorithm parameters

%%%% Dodatkowo warto zastanaowić się, czy te parametry (walkers_number, steps) na pewno powinny być łatwo modyfikowalne?
%%%% Jeżeli tak, to proszę dać komentarz jakie są odpowiednie zakresy i ewentualnie co oznaczają (np. choose from [a, b], higher value results in slower, but more precise computation)
walkers_number = 1000;
steps = 10000;
%%%% Nie jestem pewien co do tej funkcjonalności, wprowadza ona spore zamieszanie w kodzie do kretosci i służy bardziej do debugowania
%%%% Myślę że gdyby decyodwać się na jej zostawienie, lepiej byłoby zamienić to na zapis przebiegu algorytmu do plików 
%%%% Dodatkowo display jest już funkcją zdefiniowaną w Matlabie, i należy unikać takich przypadków

photos_number = nz;
%%%% tutaj tak samo, funkcja boundary jest już zdefiniowana w Matlab


microstructure = LoadImages(filePath, fileDataType, nx, ny, nz, phaseValue,direction, B);

%%%% Zamiast Results lepiej użyć więcej mówiacych nazw, np PercolatedMicrostructure, IsolatedMicrostructure, results mogłoby się równie dobrze odnościć do wartości perkolacji itp.
% PercolationResults = Percolation(nx, ny, nz, phaseValue, B);
%%%% Dlaczego Isolation nie przyjmuje mikrostruktury?
% IsoResults = Isolation(nx,ny,nz,phaseValue);
% PartiallyPercolated = 1 - PercolationResults - IsoResults;

%%%% Zawsze unikac ścieżek bezwzględnych. Do łączenia ścieżek lepiej wykorzystać uniwersalną funkcję fullfile, np.:
%%%% foldername = fullfile('..', 'SOFC', 'Photos.bmp')
%%%% wtedy kod staje się uniwersalny pomiędzy systemami 
%foldername = fullfile('..', 'dmr-quant-micros', 'Photos.bmp');

foldernameTortuosity = fullfile('..', 'dmr-quant-micros','Photos');
%foldernameTortuosity = fullfile('..', 'tortuosity','Photo*s');

%%%% URUCHOMIC W PRZYPADKU ZMIANY MIKROSTRUKTURY %%%%%%%%%%%%%%
%%%% Wczytywanie zdjęć i tworzenie pliku .mat %%%%%%%%%%%%%%%%%
zdj = load_image(photos_number, foldernameTortuosity);

%AGIRESULTS = AverageGrainSize(foldername,direction);
TortuosityRESULTS = Tortuosity(photos_number,walkers_number,steps,foldernameTortuosity);
toc;