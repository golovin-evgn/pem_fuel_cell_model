function [rhs] = solid_part(rhs, x, u, p)

%%% cathode gas channel
% fnc inputs
    % n_dot_H2_a        - molar flux density of H2
    % n_dot_H2O_a       - molar flux density of H2O
    % n_dot_O2_c        - molar flux density of O2
    % n_dot_H2O_c       - molar flux density of H2O
    % n_dot_N2_c        - molar flux density of N2
    % T_c               - cathode gas temperature
    % T_a               - anode gas temperature
    % T_s               - solid temperature
    % rho_e_s           - total energy for solid
    % U_cell            - cell voltage
    % i_m               - membrane current density
    % delta_Phi_a       - potential difference of anode double layer
    % delta_Phi_c       - potential difference of cathode double layer
    % lambda_m          - water content of membrane
    % t                 - time
    % u                 - input structure
    % p                 - parameter structure
% fnc outputs
    % drho_e_s_dt       - time derivative of rho_e_s
    % eq_T_s            - algebraic equation


% --------------------------------------------------------------------------------------
% coolant temperature system input
T_cool = u.T_cool*ones(p.N,1);

% --------------------------------------------------------------------------------------
% energy balance for the solid parts of fuel cell
        
% time derivative of total energy density for solid
% for control volume (CV) i = 1
% for CV i = 2...N-1
% for CV i = N
iN = 2:length(x.T_s)-1;
rhs.rho_e_s = [ 1/p.delta_s .* ( u.n_dot_H2_a(1) .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a(1)) + u.n_dot_H2O_a(1) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a(1)) ...
                + u.n_dot_O2_c(1) .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c(1)) + u.n_dot_H2O_c(1) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c(1)) ...
                + u.n_dot_N2_c(1) .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c(1)) ) + p.alpha_1/p.delta_s .* (x.T_a(1) - x.T_s(1)) + p.alpha_1/p.delta_s .* (x.T_c(1) - x.T_s(1))...
                + p.alpha_2/p.delta_s .* (T_cool(1) - x.T_s(1)) + p.lambda_s/p.delta_z^2 .* (x.T_s(2) -  x.T_s(1)) ...
                - 1/p.delta_s * x.U_cell * x.i_m(1);...
                1/p.delta_s .* ( u.n_dot_H2_a(iN) .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a(iN)) + u.n_dot_H2O_a(iN) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a(iN)) ...
                + u.n_dot_O2_c(iN) .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c(iN)) + u.n_dot_H2O_c(iN) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c(iN)) ...
                + u.n_dot_N2_c(iN) .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c(iN)) ) + p.alpha_1/p.delta_s .* (x.T_a(iN) - x.T_s(iN)) + p.alpha_1/p.delta_s .* (x.T_c(iN) - x.T_s(iN))...
                + p.alpha_2/p.delta_s .* (T_cool(iN) - x.T_s(iN)) + p.lambda_s/p.delta_z^2 .* (x.T_s(iN+1) - 2*x.T_s(iN) + x.T_s(iN-1)) ...
                - 1/p.delta_s * x.U_cell * x.i_m(iN);...
                1/p.delta_s .* ( u.n_dot_H2_a(p.N) .* p.h_T(p.h_H2_ref, p.c_p_H2, x.T_a(p.N)) + u.n_dot_H2O_a(p.N) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_a(p.N)) ...
                + u.n_dot_O2_c(p.N) .* p.h_T(p.h_O2_ref, p.c_p_O2, x.T_c(p.N)) + u.n_dot_H2O_c(p.N) .* p.h_T(p.h_H2O_g_ref, p.c_p_H2O_g, x.T_c(p.N)) ...
                + u.n_dot_N2_c(p.N) .* p.h_T(p.h_N2_ref, p.c_p_N2, x.T_c(p.N)) ) + p.alpha_1/p.delta_s .* (x.T_a(p.N) - x.T_s(p.N)) + p.alpha_1/p.delta_s .* (x.T_c(p.N) - x.T_s(p.N))...
                + p.alpha_2/p.delta_s .* (T_cool(p.N) - x.T_s(p.N)) + p.lambda_s/p.delta_z^2 .* (-x.T_s(p.N) + x.T_s(p.N-1)) ...
                - 1/p.delta_s * x.U_cell * x.i_m(p.N)];
            
% algebraic equation for solid temperature (39) -> (38)
        
rhs.T_s = 1e-7*( 1/p.delta_s * ( (p.delta_s - p.delta_m)*p.rho_c_p_s + p.delta_m*p.rho_c_p_m ...
                + p.delta_m*(x.lambda_m.*p.rho_m(x.lambda_m).*p.x_m(x.lambda_m))*p.c_p_H2O_l ) .* (x.T_s - p.T_ref_s) ...
                + 1/p.delta_s * p.C_a*p.delta_a_c*(x.delta_Phi_a.^2)/2 + 1/p.delta_s * p.C_c*p.delta_c_c*(x.delta_Phi_c.^2)/2 - x.rho_e_s );

end

