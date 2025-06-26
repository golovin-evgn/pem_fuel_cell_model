function rhs_vec = ode_PEMFC(t, x_vec, u_vec, p_var)
%%% system of semi-explicit differential algebraic equations
% fnc inputs
    % t      - time points
    % x      - state vector (in fnc sys_states_PEMFC(x,p))
    % u_vec  - input vector, 'model' quantities #probably switchable in a future version#
  % optional:         
    % p_var  - non-constant parameters
% fnc outputs
    % rhs - right-hand side function (for M*dx_dt = rhs(x))

% -------------------------------------------------------------------------------------- 
% load parameters

if nargin == 4

    % get configuration of identification task 
    p_ident = ident_param();

    % descale from normalized values to model values
    p_phys = p_ident.scaled2phys(p_var);

    % separate into parameter structures for each parameter category (hard code)
    p_scal_var = p_ident.par2struct(p_phys);
    
    % update parameters
    p = param_update( mod_param_PEMFC(), p_scal_var);

else
    % use default parameters
    p = mod_param_PEMFC();
end

% --------------------------------------------------------------------------------------
% copy state vector to state struct

x = p.state2struct(x_vec);

% --------------------------------------------------------------------------------------   
% initialize right-hand side (for Matlab coder)

if coder.target('MATLAB')
    rhs=[];
else
    rhs = coder.nullcopy(p.rhs_init);
end

% --------------------------------------------------------------------------------------
% initialize input struct with internal fluxes (for Matlab coder)

if ~coder.target('MATLAB')
    u = coder.nullcopy(p.fluxes_init);
end

% --------------------------------------------------------------------------------------
% get external inputs from input trajectories

% get time-interpolated external inputs at current time from input dataset
% convert from test bench inputs to model inputs

u_ext = p.inputs2struct( u_vec );

% copy external inputs to input structure
u_ext_fieldNames=fieldnames(u_ext);
for i1 = 1:numel(u_ext_fieldNames)
 u.(u_ext_fieldNames{i1}) = u_ext.(u_ext_fieldNames{i1});
end

% calculate controlled inputs u.controlledInput=f(t,x) here

%-no controlled inputs yet-
 
% --------------------------------------------------------------------------------------   
% membrane

[rhs, u] = membrane(rhs, x, u, p);

% --------------------------------------------------------------------------------------   
% catalyst layers

[rhs, u] = anode_cat(rhs, x, u, p);        
[rhs, u] = cathode_cat(rhs, x, u, p);

% --------------------------------------------------------------------------------------   
% GDLs

[rhs] = anode_GDL(rhs, x, u, p);
[rhs] = cathode_GDL(rhs, x, u, p);

% --------------------------------------------------------------------------------------   
% channels

[rhs] = anode_gas_channel(rhs, x, u, p);
[rhs] = cathode_gas_channel(rhs, x, u, p);

% --------------------------------------------------------------------------------------   
% solid part

[rhs] = solid_part(rhs, x, u, p);

% --------------------------------------------------------------------------------------
% electrode

[rhs] = electrode(rhs, x, u, p);

% --------------------------------------------------------------------------------------   
% copy right-hand-side function to vector (for M*dx_dt = f(x))

rhs_vec = p.struct2state(rhs);

end 


% --------------------------------------------------------------------------------------
% original vectors

% rhs_vec = [dc_H2_a_dt; dc_H2O_a_dt; dc_O2_c_dt; dc_H2O_c_dt; dc_N2_c_dt; d_lambda_m_dt;...
%         d_delta_Phi_a_dt; d_delta_Phi_c_dt; drho_u_a_dt; drho_u_c_dt; drho_e_s_dt;...
%         eq_xi_H2_ca; eq_xi_H2O_ca; eq_xi_O2_cc; eq_xi_N2_cc; eq_xi_H2O_cc; eq_delta_Phi_m;...
%         eq_i; eq_T_a; eq_T_c; eq_T_s; eq_p_a; eq_p_c; eq_i_u];

% x_vec = [c_H2_a, c_H2O_a, c_O2_c, c_H2O_c, c_N2_c, lambda_m,...
%           delta_Phi_a, delta_Phi_c, rho_u_a, rho_u_c, rho_e_s,...
%           xi_H2_ca, xi_H2O_ca, xi_O2_cc, xi_N2_cc, xi_H2O_cc,...
%           delta_Phi_m, i_m, T_a, T_c, T_s, p_a, p_c, U_cell]

% fluxes = [n_dot_H2O_am, n_dot_H2O_cm,n_dot_H2_a, ...
%           n_dot_H2O_a,n_dot_O2_c,n_dot_H2O_c,n_dot_N2_c]
