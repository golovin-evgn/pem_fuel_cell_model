function [rhs] = cathode_gas_channel(rhs, x, u, p)
%%% cathode gas channel
% fnc inputs
    % c_O2_c            - concentration of H2 in cathode gas channel
    % c_H2O_c           - concentration of H2O in cathode gas channel
    % c_N2_c            - concentration of N2 in cathode gas channel
    % p_c               - pressure in cathode channel
    % n_dot_O2_c        - molar flux density of O2
    % n_dot_H2O_c       - molar flux density of H2O
    % n_dot_N2_c        - molar flux density of N2
    % T_c               - cathode gas temperature
    % T_s               - solid temperature
    % rho_u_c           - internal energy density for cathode channel
    % t                 - time
    % u                 - input structure
    % p                 - parameter structure
% fnc outputs
    % dc_O2_c_dt        - time derivative of c_H2_c
    % dc_H2O_c_dt       - time derivative of c_H2O_c
    % dc_N2_c_dt        - time derivative of c_N2_c
    % drho_u_c_dt       - time derivative of rho_u_c
    % eq_T_c            - algebraic equation
    % eq_p_c            - algebraic equation


% --------------------------------------------------------------------------------------
% mass balances in cathode channel (2)

% flow velocity: discretized (4)
v_c = -p.K_c/p.delta_z * [x.p_c(2:end) - x.p_c(1:end-1); 2*(u.p_c_out - x.p_c(end))];

% time derivative of O2 concentration: discretized (2)
v_c_O2_c = v_c .* x.c_O2_c;  
rhs.c_O2_c = - 1/p.delta_z .* [(v_c_O2_c(1,1) - u.n_dot_O2_c_in); v_c_O2_c(2:end) - v_c_O2_c(1:end-1)] ...
           -  u.n_dot_O2_c /(p.delta_c);
       
% time derivative of H2O concentration: discretized (2)
v_c_H2O_c = v_c .* x.c_H2O_c;  
rhs.c_H2O_c = - 1/p.delta_z .* [(v_c_H2O_c(1,1) - u.n_dot_H2O_c_in); v_c_H2O_c(2:end) - v_c_H2O_c(1:end-1)] ...
           -  u.n_dot_H2O_c /(p.delta_c);
       
% time derivative of N2 concentration: discretized (2)
v_c_N2_c = v_c .* x.c_N2_c;  
rhs.c_N2_c = - 1/p.delta_z .* [(v_c_N2_c(1,1) - u.n_dot_N2_c_in); v_c_N2_c(2:end) - v_c_N2_c(1:end-1)] ;

% algebraic equation for anode pressure calculation (6)
rhs.p_c = p.R * x.T_c .* (x.c_O2_c + x.c_H2O_c + x.c_N2_c) - x.p_c;   

% --------------------------------------------------------------------------------------
% balance of the internal energy in cathode gas channel (modified)

% term rho*u*v
rho_u_v_c_in = u.n_dot_O2_c_in * p.h_T(p.h_O2_ref, p.c_p_O2, u.T_c_in) + u.n_dot_H2O_c_in * p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, u.T_c_in)...
             + u.n_dot_N2_c_in * p.h_T(p.h_N2_ref, p.c_p_N2, u.T_c_in) - p.R * u.T_c_in .* (u.n_dot_O2_c_in + u.n_dot_H2O_c_in + u.n_dot_N2_c_in)...
             + 2*p.lambda_c/p.delta_z * (x.T_c(1) - u.T_c_in);
rho_u_v_c = x.c_O2_c .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c) .* v_c + x.c_H2O_c .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c) .* v_c ...
          + x.c_N2_c .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c) .* v_c - x.p_c .* v_c;

% time derivative of internal energy in cathode channel
% for control volume (CV) i = 1   
% for CV i = 2...N-1 
% for CV i = N 
iN = 2:length(x.T_c)-1;
rhs.rho_u_c =  - 1/p.delta_z .* [(rho_u_v_c(1) - rho_u_v_c_in);(rho_u_v_c(iN) - rho_u_v_c(iN-1));(rho_u_v_c(p.N) - rho_u_v_c(p.N-1))]...
               + p.lambda_c/p.delta_z^2 .* [(x.T_c(2) - 3*x.T_c(1) + 2*u.T_c_in);(x.T_c(iN+1) - 2*x.T_c(iN) + x.T_c(iN-1));(-x.T_c(p.N) + x.T_c(p.N-1))] ...
               + p.alpha_1/p.delta_c .* (x.T_s - x.T_c) - u.n_dot_O2_c /(p.delta_c) .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c) ...
               - u.n_dot_H2O_c /(p.delta_c) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c) - u.n_dot_N2_c /(p.delta_c) .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c)...
               + p.R/(p.delta_c) * x.T_c .* (u.n_dot_O2_c + u.n_dot_H2O_c + u.n_dot_N2_c); 
 

% algebraic equation for cathode gas temperature (10)
rhs.T_c = 1e-5*(  x.c_O2_c .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c) + x.c_H2O_c .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c)...
                + x.c_N2_c .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c) - x.p_c - x.rho_u_c );


end

