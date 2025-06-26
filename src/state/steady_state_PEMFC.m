function [out] = steady_state_PEMFC(p,initialInput,varargin)
     
% --------------------------------------------------------------------------------------   
% calculates initial conditions for experiment
% (steady state for initial input trajectory)
% fnc inputs
    % p             - parameter structure
    % u_traj_init   - initial input trajectory
    % varargin      - structure of options_ode (optional)
    %               - vector of identified parameters (optional)
% fnc outputs
    % out           - steady state x
    
% -------------------------------------------------------------------------------------- 
use_par = 0; % solve ODE: 0 - for default parameter values, 1 - for input vector "par"

% -------------------------------------------------------------------------------------- 
% check and define the optional inputs
% two optional inputs are possible: structure "options_ode" and vector "par"
num_inputs = numel(varargin);
if num_inputs > 0
    if num_inputs == 1 && isa(varargin{1}, 'struct')
       options_ode_init = varargin{1};
    else
       par = varargin{1};
       use_par = 1;
       % options for initialization
       options_ode_init.Mass = p.M();
       options_ode_init.RelTol=1e-4;
       options_ode_init.AbsTol=1e-4;
       options_ode_init.MStateDependence = 'none';
    end
    if num_inputs == 2 && isa(varargin{1}, 'struct')
       options_ode_init = varargin{1};
       par = varargin{2};
       use_par = 1;
    elseif num_inputs == 2 && isvector(varargin{1})
       par = varargin{1};
       use_par = 1;
       options_ode_init = varargin{2};
    end
else
    % options for initialization
    options_ode_init.Mass = p.M();
    options_ode_init.RelTol=1e-4;
    options_ode_init.AbsTol=1e-4;
    options_ode_init.MStateDependence = 'none';
end

% -------------------------------------------------------------------------------------- 
% model input trajectories for initialization
[u_traj_init, ~] = initializationInputStruct(initialInput);
% convert to model inputs
u_traj_pp_init = testbenchTraj2inputTraj(u_traj_init, p);
% create interpolant
uInterpolant_pp_init = griddedInterpolant(u_traj_pp_init.time, u_traj_pp_init.data,'pchip','nearest');

% -------------------------------------------------------------------------------------- 
% close-to-feasible initial conditions
x0 = sys_states_PEMFC_initial(p,testbench2model( p.testbench2struct( u_traj_init.data(1,:).' ) ,p));

% -------------------------------------------------------------------------------------- 
% solution with ode15s-solver
if use_par
    options_ode_init.AbsTol=1e-4;
    [~, x] = ode15s(@(t,x) ode_PEMFC(t,x,uInterpolant_pp_init(t).',par), u_traj_init.time([1,end]), x0, options_ode_init);
else
    [~, x] = ode15s(@(t,x) ode_PEMFC(t,x,uInterpolant_pp_init(t).'), u_traj_init.time([1,end]), x0, options_ode_init);
end

% initial conditions for experiment
out = x(end,:).';

end
