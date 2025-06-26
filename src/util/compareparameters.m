function compareparameters(Parameters1,Parameters2)
%COMPAREPARAMETERS The function looks for differences in the model 
% parameters used for simulation and outputs them in the Command Window.
% 
% Syntax: 
%   COMPAREPARAMETERS(Parameters1,Parameters2)
%
% Inputs:
%   - used model parameters defined as a struct



% Get the field names of both structs
parameternames1 = fieldnames(Parameters1);
parameternames2 = fieldnames(Parameters2);

% Find fields that are missing in each struct
missingParameters1 = setdiff(parameternames2, parameternames1);
missingParameters2 = setdiff(parameternames1, parameternames2);

% Display missing parameters
if ~isempty(missingParameters1)
    disp('Parameters missing in first parameter set:');
    disp(missingParameters1);
end

if ~isempty(missingParameters2)
    disp('Parameters missing in second parameter set:');
    disp(missingParameters2);
end

% Compare the values of common parameters
commonParameters = intersect(parameternames1, parameternames2);

for iParameter = 1:length(commonParameters)
    parameter = commonParameters{iParameter};
    value1 = Parameters1.(parameter);
    value2 = Parameters2.(parameter);
    
    % Check if the values are functions (anonymous functions)
    if isa(value1, 'function_handle') && isa(value2, 'function_handle')
            % Convert functions to strings for comparison
            str1 = func2str(value1);
            str2 = func2str(value2);
            
            if ~isequal(str1, str2)
                fprintf('Difference in parameter "%s" (anonymous function):\n', parameter);
                fprintf('Parameter set 1: %s\n', str1);
                fprintf('Parameter set 2: %s\n', str2);
                fprintf('\n');
            end
    elseif ~isequal(value1, value2)
        fprintf('Difference in parameter "%s":\n', parameter);
        fprintf('Parameter set 1: ');
        disp(value1);
        fprintf('Parameter set 2: ');
        disp(value2);
    end
end

end