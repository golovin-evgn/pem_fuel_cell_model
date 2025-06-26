function [rhs] = electrode(rhs, x, u, p)

% integral charge balance over the length of the membrane (42)
rhs.U_cell = p.L_x * p.delta_z * sum(x.i_m) - u.I_cell;

% membrane potential (43)
rhs.delta_Phi_m = x.delta_Phi_c - x.delta_Phi_a - x.delta_Phi_m - x.U_cell * ones(p.N,1);

end