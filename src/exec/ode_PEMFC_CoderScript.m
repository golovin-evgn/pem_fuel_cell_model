function ode_PEMFC_CoderScript()
% ODE_PEMFC_CODERSCRIPT   Generate multisignature MEX-function ode_PEMFC 
% from ode_PEMFC.
% ode_PEMFC_CoderScript() generates mex file for nargs = 3 arguments
% (simulation, default) and for nargs = 4 arguments (parameter identification)
% 
%   basic script generated from project 'ode_PEMFC.prj' on 18-Jan-2023.
% 
%   See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.)

%% load parameters
p = mod_param_PEMFC();
p_ident = ident_param();

nstates = sum(p.states.variableDimensions);
ninputs = sum(p.inputs.variableDimensions);
npar = p_ident.npar;

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = true;

%% Define argument types for entry-point 'ode_PEMFC'.

ARGS = cell(2,1);

% simulation, default signature: nargs = 3
ARGS{1} = cell(3,1);                   % a total of three arguments
ARGS{1}{1} = coder.typeof(0);              % argument 1: time
ARGS{1}{2} = coder.typeof(1i,[nstates 1]); % argument 2: state vector, complex for complex step method in generated code
ARGS{1}{3} = coder.typeof(0,[ninputs 1]);  % argument 3: input vector

% parameter identification signature: nargs = 4
% default first three arguments 
ARGS{2} = [ARGS{1}; {coder.typeof(1i,[npar   1])}];  % argument 4: adaptable input parameter vector

%% Invoke MATLAB Coder.
% source files in model folder
% temporary files in dedicated codegen folder
% resulting mex file in binaries folder
codegen -config cfg -I ./src/model ode_PEMFC -args ARGS{1} -args ARGS{2} -d ./codegen/ode_PEMFC/ -o ./bin/ode_PEMFC

disp('Calls to ode_PEMFC() are now directed to the mex-file /bin/ode_PEMFC.mexXXX.');
