function [out] = sys_states_PEMFC_initial(p,u)
     
% --------------------------------------------------------------------------------------   
% feasible initial values
% state_vector = sys_states_PEMFC_initial(parameter_struct)
%
% calculate (close to) feasible initial system states and algebraic variables
% requires 
%   T_cool and T_u equal to initial temperaure T_init
%   p_u equal to initial pressure p_init
%   n_dot_*_*_in such that relative humidity is below 1
%   n_dot_*_*_in small, but not too small (pressure gradient vs. convergence rate)
%   zero overall current
%
%  initial file version, based on sys_states_PEMFC.m
%  uses initial value T_ref in order to simplify things 
%  tested with v06

% state vector x = [dynamic states; algebraic variables]
%
%%% (11*p.N) dynamic states x_d  
% x_d = [c_H2_a; c_H2O_a; c_O2_c; c_H2O_c; c_N2_c; lambda_m;
%        delta_Phi_a; delta_Phi_c; rho_u_a; rho_u_c; rho_e_s]


p_init       = p.p_u; % Pa
x0.p_init = p_init;
T_init       = p.T_ref; % K
x0.T_init = T_init;
x_H2_a_init  = u.n_dot_H2_a_in/(u.n_dot_H2_a_in+u.n_dot_H2O_a_in); % 1
x0.x_H2_a_init = x_H2_a_init;
x_H2O_a_init = u.n_dot_H2O_a_in/(u.n_dot_H2_a_in+u.n_dot_H2O_a_in); % 1
x0.x_H2O_a_init = x_H2O_a_init;
x_O2_c_init  = u.n_dot_O2_c_in/(u.n_dot_O2_c_in+u.n_dot_H2O_c_in+u.n_dot_N2_c_in); % 1
x0.x_O2_c_init  = x_O2_c_init;
x_H2O_c_init = u.n_dot_H2O_c_in/(u.n_dot_O2_c_in+u.n_dot_H2O_c_in+u.n_dot_N2_c_in); % 1
x0.x_H2O_c_init = x_H2O_c_init;
x_N2_c_init  = u.n_dot_N2_c_in/(u.n_dot_O2_c_in+u.n_dot_H2O_c_in+u.n_dot_N2_c_in); % 1
x0.x_N2_c_init = x_N2_c_init;

% concentration vectors
c_H2_a = x_H2_a_init*p_init/p.R/T_init*ones(p.N,1);
x0.c_H2_a =c_H2_a;
c_H2O_a = x_H2O_a_init*p_init/p.R/T_init*ones(p.N,1);
x0.c_H2O_a = c_H2O_a;
c_O2_c = x_O2_c_init*p_init/p.R/T_init*ones(p.N,1);
x0.c_O2_c = c_O2_c;
c_H2O_c = x_H2O_c_init*p_init/p.R/T_init*ones(p.N,1);
x0.c_H2O_c = c_H2O_c;
c_N2_c = x_N2_c_init*p_init/p.R/T_init*ones(p.N,1);
x0.c_N2_c = c_N2_c;

% water content inside the membrane
lambda_m = p.lambda_m_x(x_H2O_a_init*p_init./p.p_sat(T_init), T_init)*ones(p.N,1);
x0.lambda_m = lambda_m;
% anode potential
x0.delta_Phi_a = p.delta_Phi_a_ref*ones(p.N,1); % if deltaG temperature-dependent, then only at T_ref
% cathode potential
x0.delta_Phi_c = p.delta_Phi_c_ref*ones(p.N,1); % if deltaG temperature-dependent, then only at T_ref
% internal energy density of anode side
x0.rho_u_a = c_H2_a .* p.h_T(p.h_H2_ref, p.c_p_H2, T_init) + c_H2O_a .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, T_init) - p.R * T_init .* (c_H2_a + c_H2O_a);
% internal energy density of cathode side
x0.rho_u_c = c_O2_c .* p.h_T(p.h_O2_ref, p.c_p_O2, T_init) + c_H2O_c .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, T_init)...
        + c_N2_c .* p.h_T(p.h_N2_ref, p.c_p_N2, T_init) - p.R * T_init .* (c_O2_c + c_H2O_c + c_N2_c);
% total energy for solid
x0.rho_e_s = 1/p.delta_s * ( (p.delta_s - p.delta_m)*p.rho_c_p_s + p.delta_m*p.rho_c_p_m ...
            + p.delta_m*(lambda_m.*p.rho_m(lambda_m).*p.x_m(lambda_m))*p.c_p_H2O_l ) .* (T_init - p.T_ref_s)...
            +(1/p.delta_s * p.C_a*p.delta_a_c*(p.delta_Phi_a_ref.^2)/2 + 1/p.delta_s * p.C_c*p.delta_c_c*(p.delta_Phi_c_ref.^2)/2)*ones(p.N,1); 

%%% (12*p.N + 1) algebraic variables x_a
% x_a = [xi_H2_ca; xi_H2O_ca; xi_O2_cc; xi_N2_cc; xi_H2O_cc;
%        delta_Phi_m; i_m; T_a; T_c; T_s; p_a; p_c; U_cell]
% molar fractions in the catalyst layer
x0.xi_H2_ca = x_H2_a_init*ones(p.N,1);
x0.xi_H2O_ca = x_H2O_a_init*ones(p.N,1);
x0.xi_O2_cc = x_O2_c_init*ones(p.N,1);
x0.xi_H2O_cc = x_H2O_c_init*ones(p.N,1);
x0.xi_N2_cc = x_N2_c_init*ones(p.N,1);
% membrane potential
x0.delta_Phi_m = 0*ones(p.N,1);
% membrane current density
x0.i_m = 0*ones(p.N,1); % actually, there are small non-zero current densities if negative current densities are allowed (numerical inaccuracy? Or due to H2O flow due to pressure gradients?)
% anode gas temperature
x0.T_a = T_init*ones(p.N,1);
% cathode gas temperature
x0.T_c = T_init*ones(p.N,1);
% solid temperature
x0.T_s = T_init*ones(p.N,1);
% anode pressure
x0.p_a = p_init*ones(p.N,1); % approximately zero gradient only for small flows
% cathode pressure
x0.p_c = p_init*ones(p.N,1); % approximately zero gradient only for small flows
% cell voltage
x0.U_cell = (p.delta_Phi_c_ref-p.delta_Phi_a_ref); % if deltaG temperature-dependent, then only at T_ref


% out=[c_H2_a; c_H2O_a; c_O2_c; c_H2O_c; c_N2_c; lambda_m;...
%           delta_Phi_a; delta_Phi_c; rho_u_a; rho_u_c; rho_e_s;...
%           xi_H2_ca; xi_H2O_ca; xi_O2_cc; xi_N2_cc; xi_H2O_cc;...
%           delta_Phi_m; i_m; T_a; T_c; T_s; p_a; p_c; U_cell];

out = p.struct2state(x0);

end

