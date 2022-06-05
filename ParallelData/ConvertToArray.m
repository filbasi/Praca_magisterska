% clc; clear all ;
% fid = fopen('ParallelData_1.txt') ;  % open the text file
% S = textscan(fid,'%s');   % text scan the data
% fclose(fid) ;      % close the file
% S = S{1} ;
% N = cellfun(@(x)str2double(x), S);  % convert the cell array to double 
% N(isnan(N))=[]; % Remove NaN's which were strings earlier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc; clear all ;
% 
% % dirName = fullfile(".\ParallelData/");   % folder path (please edit this as your own)
% % files = dir( fullfile(dirName,'*.txt') );   % list of all *.txt files
% % %files = {files.name}';                   % file names
% files = dir( fullfile(fullfile(".\ParallelData/"),'*.txt'))
% pat = (".txt");
% str = dir( fullfile(fullfile(".\ParallelData/"),'*.txt'));
% TF = endsWith(pat,str);
% 
% for filename = 1: length(files)
%     if TF
%         fid = fopen(dirName + '\\' + filename); % open the text file
%         S = textscan(fid,'%s');   % text scan the data
%         fclose(fid) ;      % close the file
%         S = S{1} ;
%         N = cellfun(@(x)str2double(x), S);  % convert the cell array to double 
%         N(isnan(N))=[]; % Remove NaN's which were strings earlier
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectdir = fullfile(".\ParallelData/");
dinfo = dir( fullfile(projectdir, '*.txt'));
nfiles = length(dinfo);
assert(nfiles == 10, 'expected 10 files');
filenames = fullfile(projectdir, {dinfo.name});
data = zeros(nfiles, 10000000);  %Insert steps_number
steps = 10000000;
for K = 1 : nfiles
    thisfile = filenames{K};
    thisdata = load(thisfile);
    data(K, :) = thisdata;
end

wykres = mean(data, 1);
%columnMeans = max(data);


%wykres = rt;
%z=1:walkers_number;
z=1:steps;
plot(z,wykres,'r-','LineWidth',2);
%ylim([0 inf]);
ylabel('\it MSD \rm[-]');
xlabel('\it \vartheta \rm[-]');
aprox = polyfit(z,wykres,1);
% Evaluate fit equation using polyval
y_est = polyval(aprox,z);
% Add trend line to plot
hold on 
plot(z,y_est,'black--','LineWidth',2)
hold off
legend('Funkcja rzeczywista', 'Prosta aproksymujaca')
txt1 = ['y = ' num2str(aprox(1)) 'x + ' num2str(aprox(2))];
%text(x, y, txt1);
Taugraph=1/aprox(1);
result11 = ['Tortuosity factor calculated from graph is ', num2str(Taugraph),'. '];
disp(['Trend line equation is y = ' num2str(aprox(1)) '*x + ' num2str(aprox(2))])
disp(result11)