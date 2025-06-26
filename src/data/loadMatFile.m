function [data,info] = loadMatFile(fileName,rangeRows,modelVariables)
% loadMatFile load measurement data from files
%
% [data,info] = loadMatFile('path/fileName',[])
% loads all data
%
% [data,info] = loadMatFile('path/fileName',[startRow endRow])
% loads data from specified rows 
%
% [data,info] = loadMatFile('path/fileName',[],{'variable1';'variable2';...})
% loads all data and compares variable names in the data file with supplied list of variable names 
%
% [data,info] = loadMatFile('path/fileName',[startRow endRow],{'variable1';'variable2';...})
% loads data from specified rows and compares variable names in the data file with supplied list of variable names 


fromFile = load(fileName,'data','info');
info = fromFile.info;
info.selectedRange = rangeRows;
if isempty(rangeRows)
  data = fromFile.data;
else
  data.time=fromFile.data.time(rangeRows(1):rangeRows(2));
  data.data=fromFile.data.data(rangeRows(1):rangeRows(2),:);
end

if nargin>2

% check whether inputs from data file are the same as model inputs !!! also check units?
if ~all(string(modelVariables)==string(info.variableNames))
    disp(' ');
    disp('Warning: Variable names in data file are not identical to model variables.');
    disp(' ');
    disp('   data file:');   
    disp(['     ',fileName]); 
    disp(' ');
    disp('   data file variables:');
    disp([string(info.variableNames)]);
    disp('   model variables:');
    disp([string(modelVariables)]);
    pause();
end

end

clearvars fromFile; % remove large variable from memory
