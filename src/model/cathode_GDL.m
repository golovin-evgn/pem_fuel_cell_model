function [rhs] = cathode_GDL(rhs, x, u, p)
%%% cathode GDL
% fnc inputs
    % n_dot_O2_c        - molar flux density of O2
    % n_dot_N2_c        - molar flux density of N2
    % n_dot_H2O_c       - molar flux density of H2O
    % c_O2_c            - concentrtion of O2 in cathode gas channel
    % c_N2_c            - concentrtion of N2 in cathode gas channel
    % c_H2O_c           - concentrtion of H2O in cathode gas channel
    % p_c               - pressure of cathode channel
    % T_s               - solid temperature
    % xi_O2_cc          - molar fraction of O2 in catalyst layer
    % xi_N2_cc          - molar fraction of N2 in catalyst layer
    % xi_H2O_cc         - molar fraction of H2O in catalyst layer
    % p                 - parameter structure
% fnc outputs
    % eq_xi_O2_ca       - algebraic equation
    % eq_xi_N2_ca       - algebraic equation
    % eq_xi_H2O_ca      - algebraic equation

% --------------------------------------------------------------------------------------
% mole fractions of species in the cathode channel
xi_O2_c = x.c_O2_c ./ (x.c_O2_c + x.c_H2O_c + x.c_N2_c);
xi_N2_c = x.c_N2_c ./ (x.c_O2_c + x.c_H2O_c + x.c_N2_c);
xi_H2O_c = x.c_H2O_c ./ (x.c_O2_c + x.c_H2O_c + x.c_N2_c);

% total gas concentration in the cathode GDL (14)
c_cc = x.p_c ./ (p.R .* x.T_s);

% algebraic equations for mole fractions of species in cathode GDL (11)
rhs.xi_O2_cc = (x.xi_O2_cc - xi_O2_c)/p.delta_c_g + (0.5*(x.xi_H2O_cc + xi_H2O_c).*u.n_dot_O2_c ...
    - 0.5*(x.xi_O2_cc + xi_O2_c).*u.n_dot_H2O_c) ./ (c_cc*p.D_eff_O2_H2O_g) ...
    + (0.5*(x.xi_N2_cc + xi_N2_c).*u.n_dot_O2_c ...
    - 0.5*(x.xi_O2_cc + xi_O2_c).*u.n_dot_N2_c) ./ (c_cc*p.D_eff_O2_N2);

rhs.xi_N2_cc = (x.xi_N2_cc - xi_N2_c)/p.delta_c_g + (0.5*(x.xi_O2_cc + xi_O2_c).*u.n_dot_N2_c ...
    - 0.5*(x.xi_N2_cc + xi_N2_c).*u.n_dot_O2_c) ./ (c_cc*p.D_eff_O2_N2) ...
    + (0.5*(x.xi_H2O_cc + xi_H2O_c).*u.n_dot_N2_c ...
    - 0.5*(x.xi_N2_cc + xi_N2_c).*u.n_dot_H2O_c) ./ (c_cc*p.D_eff_N2_H2O_g);

rhs.xi_H2O_cc = 1 - x.xi_O2_cc - x.xi_N2_cc - x.xi_H2O_cc;

end

