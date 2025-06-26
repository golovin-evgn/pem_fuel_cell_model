function stop = optimsaveresults(x,optimValues,state,filename)
%OPTIMSAVERESULTS is a custom output function for the optimization function, 
% that saves the optimization results at each iteration.
% The optimization function passes the values to the input arguments of the 
% output function at each iteration.
%
%   Syntax:
%   stop = OPTIMSAVERESULTS(x,optimValues,state,filename)
% 
%   Inputs:
%       x           - point computed by the optimization algorithm at the current iteration
%       optimValues - structure containing data from the current iteration
%       state       - current state of the optimization algorithm
%       filename    - Name of file, specified as a string scalar or character vector
% 
%   Outputs:
%       stop - a flag that is true or false depending on whether the 
%              optimization routine should stop (true) or continue (false)
% 
%   Example:
%   Create an options structure that will use OPTIMSAVERESULTS as the
%   output function
% 
%     outputfcn = @(x, optimValues, state) ...
%     optimsaveresults(x,optimValues,state,'optim_results.mat');
%
%     options = optimoptions('fminunc','OutputFcn',outputfcn);
%
%   Pass the options into an optimization problem to save the results of
%   each iteration in the specified file
% 
%     fminunc(@(x) 3*sin(x(1))+exp(x(2)),x0,options)

stop = false;

% Check if the file exists or create a new structure
if isfile(filename)
    load(filename,'results');
else
    results.optimizedParameters = [];    
    results.firstorderopt = [];
    results.funccount = [];
    results.fvalues = [];
    results.gradient = [];
    results.nIteration = [];
    results.stepsize = [];
end

% Append or overwrite current iteration results to the structure
results.optimizedParameters = x;
results.firstorderopt = [results.firstorderopt; optimValues.firstorderopt];
results.funccount = [results.funccount; optimValues.funccount];
results.fvalues = [results.fvalues; optimValues.fval];
results.gradient = optimValues.gradient;
results.nIteration = optimValues.iteration;
results.stepsize = [results.stepsize; optimValues.stepsize];

% Save the updated structure to the specified .mat file
save(filename, 'results');

end

