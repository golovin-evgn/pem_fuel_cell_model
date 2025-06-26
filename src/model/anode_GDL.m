function [rhs] = anode_GDL( rhs, x, u, p)
%%% anode GDL
% fnc inputs
    % n_dot_H2_a        - molar flux density of H2
    % n_dot_H2O_a       - molar flux density of H2O
    % c_H2_a            - concentration of H2 in anode gas channel
    % c_H2O_a           - concentration of H2O in anode gas channel
    % p_a               - pressure of anode channel
    % T_s               - solid temperature
    % xi_H2_ca          - molar fraction of H2 in catalyst layer
    % xi_H2O_ca         - molar fraction of H2O in catalyst layer
    % p                 - parameter structure
% fnc outputs
    % eq_xi_H2_ca       - algebraic equation
    % eq_xi_H2O_ca      - algebraic equation

% --------------------------------------------------------------------------------------   
% mole fractions of species in the anode channel
xi_H2_a = x.c_H2_a ./ (x.c_H2_a + x.c_H2O_a);
xi_H2O_a = x.c_H2O_a ./ (x.c_H2_a + x.c_H2O_a);

% total gas concentration in the anode GDL (14)
c_ca = x.p_a ./ (p.R * x.T_s);

% algebraic equations for mole fractions of species in anode GDL (11)
rhs.xi_H2_ca = (x.xi_H2_ca - xi_H2_a)/p.delta_a_g + (0.5*(x.xi_H2O_ca + xi_H2O_a) .* u.n_dot_H2_a ...
    - 0.5*(x.xi_H2_ca + xi_H2_a) .* u.n_dot_H2O_a) ./ (c_ca*p.D_eff_H2_H2O_g);

rhs.xi_H2O_ca = 1 - x.xi_H2_ca - x.xi_H2O_ca;

end

