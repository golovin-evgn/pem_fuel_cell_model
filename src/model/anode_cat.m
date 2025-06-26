function [rhs, u] = anode_cat(rhs, x, u, p)
%%% anode catalyst layer
% fnc inputs 
    % delta_Phi_a       - potential difference of anode double layer
    % p_a               - pressure of the anode channel
    % xi_H2_ca          - molar fraction of H2 in catalyst layer
    % i_m               - membrane current density
    % n_dot_H2O_am      - water flux between membrane and anode gas bulk
    % T_s               - solid temperature
    % p                 - parameter structure
% fnc outputs
    % d_delta_Phi_a_dt  - time derivative of delta_Phi_a
    % n_dot_H2_a        - H2 flux to catalyst layer
    % n_dot_H2O_a       - H2O flux to catalyst layer

% --------------------------------------------------------------------------------------   
% anodic reaction rate with Butler-Volmer kinetics (20)

r_a = p.f_v * p.i_0_ref_a / (2*p.F) * ( (exp((p.alpha_a*2*p.F)./(p.R*x.T_s) .* ( x.delta_Phi_a - p.delta_Phi_a_ref )) )...
     .* (x.p_a .* x.xi_H2_ca./p.p_H2_ref) - ( exp((-(1-p.alpha_a)*2*p.F)./(p.R*x.T_s) .* ( x.delta_Phi_a - p.delta_Phi_a_ref )) ) ) ;

% --------------------------------------------------------------------------------------  
% charge balance at the anode double layer (40)

rhs.delta_Phi_a = 1/(p.C_a*p.delta_a_c) * (x.i_m - 2*p.F*r_a);

% --------------------------------------------------------------------------------------  
% catalyst layer mass balances (15), (16)

u.n_dot_H2_a = r_a;
u.n_dot_H2O_a = u.n_dot_H2O_am;

end

