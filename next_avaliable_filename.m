%%%%%%%%%%%%%
% Require "Next Available Filename" Toolbox by Stephen
%%%%%%%%%%%%%

if exist(filename,'file')
    [fPath, fName, fExt] = fileparts(filename);
    fDir = dir(fullfile(fPath, [fName,' (*)', fExt]));
    if isempty(fDir)
        filename = fullfile(fPath, [fName,' (1)', fExt]);
    else
        pattern = "(" + digitsPattern + ")" + fExt;
        hasThePattern = endsWith(extractfield(fDir,'name'),pattern);
        Extracted = extract(extractfield(fDir(hasThePattern),'name'),pattern);
        num = max(cell2mat(cellfun(@(C) textscan(C,'(%d)') , Extracted,'UniformOutput',true)));
        num = num+1;
        filename = fullfile(fPath, [fName,' (',num2str(num),')', fExt]);
    end
end