 function [rhs] = anode_gas_channel(rhs, x, u, p)
       
%%% anode gas channel
% fnc inputs
    % c_H2_a            - concentration of H2 in anode gas channel
    % c_H2O_a           - concentration of H2O in anode gas channel
    % p_a               - pressure of the anode channel
    % n_dot_H2_a        - molar flux density of H2
    % n_dot_H2O_a       - molar flux density of H2O
    % T_a               - anode gas temperature
    % T_s               - solid temperature
    % rho_u_a           - internal energy density for anode channel
    % t                 - time
    % u                 - input structure
    % p                 - parameter structure
% fnc outputs
    % dc_H2_a_dt        - time derivative of c_H2_a
    % dc_H2O_a_dt       - time derivative of c_H2O_a
    % drho_u_a_dt       - time derivative of rho_u_a
    % eq_T_a            - algebraic equation
    % eq_p_a            - algebraic equation

if p.counter_flow
    % counter-flow feeding

    % --------------------------------------------------------------------------------------
    % mass balances in anode channel (2)
    % flow velocity: discretized (4)
    v_a = -p.K_a/p.delta_z * [2*(u.p_a_out - x.p_a(1,1)); x.p_a(1:end-1) - x.p_a(2:end)];

    % time derivative of H2 concentration: discretized (2)
    v_c_H2_a = v_a .* x.c_H2_a;
    rhs.c_H2_a = - 1/p.delta_z .* [v_c_H2_a(1:end-1) - v_c_H2_a(2:end); (v_c_H2_a(end) - u.n_dot_H2_a_in)] ...
           -  u.n_dot_H2_a /(p.delta_a);

    % time derivative of H2O concentration: discretized (2)
    v_c_H2O_a = v_a .* x.c_H2O_a;
    rhs.c_H2O_a = - 1/p.delta_z .* [v_c_H2O_a(1:end-1) - v_c_H2O_a(2:end); (v_c_H2O_a(end) - u.n_dot_H2O_a_in)] ...
            - u.n_dot_H2O_a /(p.delta_a);

    % algebraic equation for anode pressure calculation (6)
    rhs.p_a = p.R * x.T_a .* (x.c_H2_a + x.c_H2O_a) - x.p_a;

    % --------------------------------------------------------------------------------------   
    % balance of the internal energy in anode gas channel (modified)
    % term rho*u*v
    rho_u_v_a = x.c_H2_a .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a) .* v_a + x.c_H2O_a .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) .* v_a ...
            - x.p_a .* v_a;
    rho_u_v_a_in = u.n_dot_H2_a_in * p.h_T(p.h_H2_ref, p.c_p_H2, u.T_a_in) + u.n_dot_H2O_a_in * p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, u.T_a_in)...
                 - p.R * u.T_a_in .* (u.n_dot_H2_a_in + u.n_dot_H2O_a_in) + 2*p.lambda_a/p.delta_z * (x.T_a(end) - u.T_a_in);

    % time derivative of internal energy in anode channel
    % for CV i = N
    iN = 2:length(x.T_a)-1;
    rhs.rho_u_a =  - 1/p.delta_z .* [(rho_u_v_a(1) - rho_u_v_a(2));(rho_u_v_a(iN) - rho_u_v_a(iN+1));(rho_u_v_a(p.N) - rho_u_v_a_in)]...
                   + p.lambda_a/p.delta_z^2 .* [(-x.T_a(1) + x.T_a(2));(x.T_a(iN-1) - 2*x.T_a(iN) + x.T_a(iN+1));(x.T_a(p.N-1) -3*x.T_a(p.N) + 2*u.T_a_in)]...
                   + p.alpha_1/p.delta_a .* (x.T_s - x.T_a) - u.n_dot_H2_a /(p.delta_a) .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a)...
                   - u.n_dot_H2O_a /(p.delta_a) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) + p.R/(p.delta_a) * x.T_a .* (u.n_dot_H2_a + u.n_dot_H2O_a);

    % algebraic equation for anode gas temperature (10)
    rhs.T_a = 1e-5*( x.c_H2_a .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a) + x.c_H2O_a .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) - x.p_a -  x.rho_u_a );
else
    % co-flow feeding

    % --------------------------------------------------------------------------------------
    % mass balances in anode channel (2)
    % flow velocity: discretized (4)
    v_a = -p.K_a/p.delta_z * [x.p_a(2:end) - x.p_a(1:end-1); 2*(u.p_a_out - x.p_a(end))];

    % time derivative of H2 concentration: discretized (2)
    v_c_H2_a = v_a .* x.c_H2_a;
    rhs.c_H2_a = - 1/p.delta_z .* [(v_c_H2_a(1,1) - u.n_dot_H2_a_in); v_c_H2_a(2:end) - v_c_H2_a(1:end-1)] ...
           -  u.n_dot_H2_a /(p.delta_a);

    % time derivative of H2O concentration: discretized (2)
    v_c_H2O_a = v_a .* x.c_H2O_a;
    rhs.c_H2O_a = - 1/p.delta_z .* [(v_c_H2O_a(1,1) - u.n_dot_H2O_a_in); v_c_H2O_a(2:end) - v_c_H2O_a(1:end-1)] ...
            - u.n_dot_H2O_a /(p.delta_a);

    % algebraic equation for anode pressure calculation (6)
    rhs.p_a = p.R * x.T_a .* (x.c_H2_a + x.c_H2O_a) - x.p_a;

    % --------------------------------------------------------------------------------------   
    % balance of the internal energy in anode gas channel (modified)
    % term rho*u*v
    rho_u_v_a = x.c_H2_a .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a) .* v_a + x.c_H2O_a .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) .* v_a ...
            - x.p_a .* v_a;
    rho_u_v_a_in = u.n_dot_H2_a_in * p.h_T(p.h_H2_ref, p.c_p_H2, u.T_a_in) + u.n_dot_H2O_a_in * p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, u.T_a_in)...
                 - p.R * u.T_a_in .* (u.n_dot_H2_a_in + u.n_dot_H2O_a_in) + 2*p.lambda_a/p.delta_z * (x.T_a(1) - u.T_a_in);

    % time derivative of internal energy in anode channel
    % for CV i = N
    iN = 2:length(x.T_a)-1;
    rhs.rho_u_a =  - 1/p.delta_z .* [(rho_u_v_a(1) - rho_u_v_a_in);(rho_u_v_a(iN) - rho_u_v_a(iN-1));(rho_u_v_a(p.N) - rho_u_v_a(p.N-1))]...
                   + p.lambda_a/p.delta_z^2 .* [(x.T_a(2) - 3*x.T_a(1) + 2*u.T_a_in);(x.T_a(iN+1) - 2*x.T_a(iN) + x.T_a(iN-1));(-x.T_a(p.N) + x.T_a(p.N-1))]...
                   + p.alpha_1/p.delta_a .* (x.T_s - x.T_a) - u.n_dot_H2_a /(p.delta_a) .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a)...
                   - u.n_dot_H2O_a /(p.delta_a) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) + p.R/(p.delta_a) * x.T_a .* (u.n_dot_H2_a + u.n_dot_H2O_a);

    % algebraic equation for anode gas temperature (10)
    rhs.T_a = 1e-5*( x.c_H2_a .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a) + x.c_H2O_a .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a) - x.p_a -  x.rho_u_a );
end
       
end

