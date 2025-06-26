%%%
%%% example of parameter identification of pem fuel cell model based on the
%%% data generated from the model
%%%

%% initialization
clear all; close all; clc;

% use matlab coder for faster computation by generating mex-function of ode PEMFC
force_codegen = 0; % 0 - use presaved mex-file; 1 - new mex-file needed

if force_codegen 
    ode_PEMFC_CoderScript();
end

dequant_IO = 0; % 1 - use de-quantization filter for input trajectories; 0 - don't use it

%% Parameters initialization
% load model parameters
p = mod_param_PEMFC();

% initialize default model parameters 
p_ident = ident_param();
n_par = p_ident.npar;

% initial guess for identified parameters
dev_par = 30; % in percent
rng(1);
rand_init = (1-dev_par/100 + (1+dev_par/100 - 1+dev_par/100)*rand(n_par,1));
initialParameters = 0.5*ones(n_par,1) .* rand_init; % random initialization

% descale from normalized values to model values
initialParameters_descaled = p_ident.scaled2phys(initialParameters);
p_phys = initialParameters_descaled;
% from vector to struct
p_scal_var = p_ident.par2struct(p_phys);

% update parameters
p = param_update( mod_param_PEMFC(), p_scal_var );

%% Data Selection and Preprocessing
% load data of model inputs and outputs
% set needed initial variables
% perform a normalization of model outputs for optimization

inputFileName = 'dataset_15_OP_250622_v03';
outputFileName_in = [inputFileName,'_inputs'];
outputFileName_out = [inputFileName,'_outputs'];
rangeRows = []; % determined manually

% input trajectories data
[u_traj,u_traj_info] = loadMatFile(outputFileName_in,rangeRows,p.testbench.variableNames);

% test bench output data
[y_traj,y_traj_info] = loadMatFile(outputFileName_out,rangeRows);

% min-max-scaling
y_traj.y_min = min(y_traj.data);
y_traj.y_max = max(y_traj.data);
y_traj.data_scaled = p.output2scaled(y_traj.data,y_traj.y_min,y_traj.y_max);

%% pre-processing and interpolation for data
info('Pre-processing started.');

% de-quantize input trajectories
if dequant_IO
    for i1=1:length(u_traj_info.variableNames)
        [u_traj.data(:,i1),~] = deQuant(u_traj.data(:,i1));
    end
end

% convert testbench inputs to model inputs
u_traj_pp = testbenchTraj2inputTraj(u_traj, p);

% create interpolant for model inputs
uInterpolant_pp = griddedInterpolant(u_traj_pp.time, u_traj_pp.data,'pchip','nearest'); % spline, pchip, linear

% initial point from measurement file 
initialInput = p.testbench2struct(u_traj.data(1,:).');

info('Pre-processing complete.');

%% DAE options
% set options for ode15s-solver including tolerances, a mass matrix and a
% jacobian Pattern

options_ode.Mass = p.M();
% adjusted tolerances for higher gradient accuracy
options_ode.AbsTol = 1e-6; % default 1e-6
options_ode.RelTol = 1e-6; % default 1e-3, RelTol >= 1e-9
options_ode.MStateDependence = "none";


% JPattern
force_rec_jac = 1; % 0 - use presaved JPattern; 1 - new JPattern needed

JPattern = myJPattern();

rec_jac = force_rec_jac || p.nstates ~= size(JPattern,1);

if rec_jac == 0; options_ode.JPattern = JPattern; end

% generate JPattern if it is needed
% based on initial guess parameters and initialization step
if rec_jac == 1

    disp("JPattern calculation started.");

        initialStates = steady_state_PEMFC(p,initialInput,options_ode, ...
                                           initialParameters);

        % calculate model input vector from testbench input struct
        u_vec = p.inputs2vec(testbench2model(initialInput,p));

        % generate and save new Jpattern
        genJpattern(@(x) ode_PEMFC(0,x,u_vec,initialParameters),initialStates);
        options_ode.JPattern=myJPattern(); 

    disp("JPattern calculation complete.");

end



%% Optimization and Training the Network
% optimize adaptable and neural network parameters

startOptimization = true;   % true - start optimization
                             % false - load optimization results

filename_results = "res/optimization/test_parameter_ident_01.mat"; % file extension needed                            

if startOptimization == true && ~isfile(filename_results) % prevent file overwriting

    options_opt = optimoptions("fminunc");
    options_opt.Algorithm = "quasi-newton";
    options_opt.SpecifyObjectiveGradient = true;
    options_opt.HessianApproximation = "lbfgs";
    options_opt.OptimalityTolerance = 1e-16;
    options_opt.FunctionTolerance = 1e-8;
    options_opt.StepTolerance = 1e-16;
    options_opt.Diagnostics = "on";
    options_opt.Display = "iter-detailed";
    options_opt.MaxIterations = 10;
    options_opt.OutputFcn = @(x,optimValues,state) optimsaveresults( ...
        x,optimValues,state,filename_results);
    
    % problem solving
    info("Parametrization started.");
    
    [par_opt,fval,exitflag,output] = fminunc( ...
        @(par) objective_PEMFC(par,u_traj,y_traj,options_ode), ...
        initialParameters, ...
        options_opt);
    
    info("Parametrization complete.");

    save(filename_results,"output","exitflag","-append");

    optimizedParameters = par_opt;
end


%% Validation
% solution with optimized parameters
if ~startOptimization 
    optimizedParameters = initialParameters;
end
optimizedAnnParameters_descaled = p_ident.scaled2phys(optimizedParameters);
p_phys = optimizedAnnParameters_descaled;

% separate into parameter structures for each parameter category (hard code)
p_scal_var = p_ident.par2struct(p_phys);

% update parameters and activated ANNs
p = param_update( mod_param_PEMFC(), p_scal_var);
 
tspan = u_traj.time;

u_traj_pp = testbenchTraj2inputTraj(u_traj, p);

uInterpolant_pp = griddedInterpolant(u_traj_pp.time, u_traj_pp.data,"pchip","nearest"); % spline, pchip, linear

u = @(t) uInterpolant_pp(t).';

initialStates = steady_state_PEMFC(p,initialInput,options_ode,optimizedParameters);

% tic
[t, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t),optimizedParameters), tspan, initialStates, options_ode);
% tsim = toc

xStruct = p.state2struct(x.');

y_model = sys_output_PEMFC(x,u_traj_pp,p);
yStruct = p.output2struct(y_model.');

%% plots
yData = p.output2struct(y_traj.data.');

h1 = figure(1);
set(gcf, 'unit', 'normalized', 'position', ...
   [0.4589 0.1704 0.4964 0.6361]);

subplot(3,1,1);
plot(t,yStruct.U_cell); grid on; hold on;
plot(u_traj.time,yData.U_cell);
xlabel('$t ~/~ \mathrm{s}$','fontsize',12, 'Interpreter','latex');
ylabel('$U_{\mathrm{cell}} ~/~ \mathrm{V}$','fontsize',12, 'Interpreter','latex');
set(gca,'FontSize',14);

subplot(3,1,2);
plot(t,yStruct.a_H2O_a_out); grid on; hold on;
plot(u_traj.time,yData.a_H2O_a_out);
xlabel('$t ~/~ \mathrm{s}$','fontsize',12, 'Interpreter','latex');
ylabel('$a_{\mathrm{H_2 O,a}} ~/~ \mathrm{1}$','fontsize',12, 'Interpreter','latex');
set(gca,'FontSize',14);

subplot(3,1,3);
plot(t,yStruct.a_H2O_c_out); grid on; hold on;
plot(u_traj.time,yData.a_H2O_c_out);
xlabel('$t ~/~ \mathrm{s}$','fontsize',12, 'Interpreter','latex');
ylabel('$a_{\mathrm{H_2 O,c}} ~/~ \mathrm{1}$','fontsize',12, 'Interpreter','latex');
set(gca,'FontSize',14);

% exportgraphics(h1,['ident_results_01.jpeg']);