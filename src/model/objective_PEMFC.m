function [out, gradout] = objective_PEMFC(par,u_traj,y_traj,options_ode)

% -------------------------------------------------------------------------
% objective function for parameter identification
% fnc inputs
    % par           - scaled parameter vector
    % u_traj        - input trajectories of dataset
    % y_data_scaled - scaled output of dataset
    % options_ode   - ODE options
% fnc outputs
    % out           - objective function (mse/mse with regularization)
    % gradout       - objective function gradient
% -------------------------------------------------------------------------

% get configuration of identification task 
p_ident = ident_param();

% descale from normalized values to model values
p_phys = p_ident.scaled2phys(par);

% separate into parameter structures for each parameter category (hard code)
p_scal_var = p_ident.par2struct(p_phys);

% update parameters and activated ANNs
p = param_update( mod_param_PEMFC(), p_scal_var );

% ode solver options
if ~coder.target('MATLAB')
    options_ode.JPattern = @myJPattern;
else
    options_ode.JPattern = myJPattern();
end

% -------------------------------------------------------------------------
% initialization for current set of parameters

% initial point from measurement file 
initialInput = p.testbench2struct(u_traj.data(1,:).');

% initial conditions for experiments
try
    x0_initialization = steady_state_PEMFC(p,initialInput,options_ode,par);
catch
    % return NaN if a failure occur during model solving
    % force the optimization solver to attempt a different iterative step
    out = NaN;
    gradout = NaN*ones(size(par.'));
    return 
end

% -------------------------------------------------------------------------
% simulation results current set of parameters

% use only cell voltage and molar flux densities of water
ySelect = p.selectedOutputs == 1;

% min-max scaling
y_min = y_traj.y_min(ySelect);
y_max = y_traj.y_max(ySelect);
y_data_scaled = y_traj.data_scaled(:,ySelect);

% convert testbench inputs to model inputs
u_traj_pp = testbenchTraj2inputTraj(u_traj, p);
% create interpolant for model inputs
uInterpolant_pp = griddedInterpolant(u_traj_pp.time, u_traj_pp.data,'pchip','nearest'); % spline, pchip, linear
% input values from interpolant
u = @(t) uInterpolant_pp(t).';

% time span
tspan = u_traj.time;

if nargout>1
    % with gradient
    mymse=complex(zeros(size(par.')));
    stepSize=1e-100;

    % variables redefining for reduced memory overhead in parfor
    % (not necessary for thread-based environment)

    output2scaled = @p.output2scaled;
    p_constant = parallel.pool.Constant(p);

    % regularization parameter
    lambda_reg = 0.02;

    % try objective function evaluation and gradient calculation
    try
        parfor i1=1:length(par)
            % make complex step
            parStep = complex(par);
            parStep(i1) = parStep(i1) + 1i*stepSize;

            % -------------------------------------------------------------------------
            % solution with ode15s-solver
            [~, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t),parStep), tspan, x0_initialization, options_ode);

            % calculate output
            ptemp = p_constant.Value;
            y_model = sys_output_PEMFC(x,u_traj_pp,ptemp);
            y_model_scaled = output2scaled(y_model(:,ySelect),y_min,y_max);
    
            % calculate mean squared error
            % with regularization
            mymse(1,i1) = mean((y_model_scaled - y_data_scaled).^2,'all') ...
                        + mean((lambda_reg.*(parStep-par)).^2);
        end
    catch
            % return NaN if a failure occur during model solving
            % force the optimization solver to attempt a different iterative step
            out = NaN;
            gradout = NaN*ones(size(par.'));
            return
    end

    % objective value
    out=real(mymse(1,1));
    % numeric gradient
    gradout=imag(mymse)/stepSize;
    
else
    % without gradient
    try
        % try to solve model
        % solution with ode15s-solver
        [~, x] = ode15s(@(t,x) ode_PEMFC(t,x,u(t),par), tspan, x0_initialization, options_ode);
        % output equation
        y_model = sys_output_PEMFC(x,u_traj_pp,p);
    catch
        % return NaN if a failure occur during model solving
        % force the optimization solver to attempt a different iterative step
        out = NaN;
        return
    end

    % min-max-scaling of model outputs
    y_model_scaled = p.output2scaled(y_model(:,ySelect),y_min,y_max);

    % without regularization
    out = mean((y_model_scaled - y_data_scaled).^2,'all');

end

end

