function [rhs, u] = cathode_cat(rhs, x, u, p)
%%% cathode catalyst layer
% inputs
    % delta_Phi_c       - potential difference of cathode double layer
    % p_c               - pressure of cathode channel
    % xi_H2_cc          - molar fraction of O2 in catalyst layer
    % i_m               - membrane current density
    % n_dot_H2O_cm      - water flux between membrane and cathode gas bulk
    % T_s               - solid temperature
    % p                 - parameter structure
% outputs
    % d_delta_Phi_c_dt  - time derivative of delta_Phi_c
    % n_dot_O2_c        - O2 flux to catalyst layer
    % n_dot_H2O_c       - H2O flux to catalyst layer
    % n_dot_N2_c        - N2 flux to catalyst layer

% --------------------------------------------------------------------------------------   
% cathodic reaction rate with Butler-Volmer kinetics (21)

delta_Phi_c_ref = p.delta_Phi_c_ref;

r_c = p.f_v * p.i_0_ref_c / (2*p.F) * (x.p_c .* x.xi_O2_cc./p.p_O2_ref) .* exp(-p.E_cr ./ (p.R*x.T_s).*(1 - x.T_s/p.T_ref)) .*...
     ( (exp(-(p.alpha_c*2*p.F)./(p.R*x.T_s) .* ( x.delta_Phi_c - delta_Phi_c_ref ) ) ) ...
      - ( exp(((1-p.alpha_c)*2*p.F)./(p.R*x.T_s) .* ( x.delta_Phi_c - delta_Phi_c_ref) ) ) ) ;


% --------------------------------------------------------------------------------------   
% charge balance at the cathide double layer (41)

rhs.delta_Phi_c = 1/(p.C_c*p.delta_c_c) * (-x.i_m + 2*p.F*r_c);

% --------------------------------------------------------------------------------------   
% catalyst layer mass balances (17), (18), (19)

u.n_dot_O2_c = 1/2 * r_c;
u.n_dot_H2O_c = u.n_dot_H2O_cm - r_c;
u.n_dot_N2_c = zeros(size(r_c));

end

