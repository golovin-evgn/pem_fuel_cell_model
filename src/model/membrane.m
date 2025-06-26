function [rhs, u] = membrane(rhs, x, u, p)
%%% membrane
% fnc inputs
    % lambda_m          - water content in membrane
    % xi_H2O_ca         - molar fraction of H2O in catalyst layer
    % xi_H2O_cc         - molar fraction of H2O in catalyst layer
    % T_s               - solid temperature
    % delta_Phi_m       - potential difference of membrane
    % i_m               - membrane current density
    % p_a               - anode gas pressure
    % p_c               - cathode gas pressure
    % p                 - parameter structure
% fnc outputs
    % d_lambda_m_dt     - time derivative of water content in membrane
    % eq_i              - algebraic equation
    % n_dot_H2O_am      - water flux between membrane and anode gas bulk
    % n_dot_H2O_cm      - water flux between membrane and cathode gas bulk


% --------------------------------------------------------------------------------------   
% concentrations of water on sides of membrane
% c_H2O_ca = x.xi_H2O_ca .* x.p_a ./ (p.R * x.T_s);
% c_H2O_cc = x.xi_H2O_cc .* x.p_c ./ (p.R * x.T_s);

% relative humidities on the sides of membrane
a_H2O_ca = x.xi_H2O_ca .* x.p_a ./ p.p_sat(x.T_s);
a_H2O_cc = x.xi_H2O_cc .* x.p_c ./ p.p_sat(x.T_s); 

% xi_H2O_a = x.c_H2O_a ./ (x.c_H2_a + x.c_H2O_a);
% xi_H2O_c = x.c_H2O_c ./ (x.c_O2_c + x.c_H2O_c + x.c_N2_c);
% 
% % relative humidities on the sides of membrane
% a_H2O_ca = xi_H2O_a .* x.p_a ./ p.p_sat(x.T_s);
% a_H2O_cc = xi_H2O_c .* x.p_c ./ p.p_sat(x.T_s); 

% % correction for more simulation robustness in invalid model region
% % (model invalid fo condensing water vapor)
% % -> changes model behavior for a_H2O_a/c>0.95
scaling = 300; % scaling parameter smoothness vs. model error 300
a_H2O_ca = -smoothmax(-a_H2O_ca,-0.97,scaling); %!!!
a_H2O_cc = -smoothmax(-a_H2O_cc,-0.97,scaling); %!!!


% mean water content of the membrane
lambda_m = x.lambda_m;

lambda_am = p.lambda_m_x(a_H2O_ca, x.T_s).*p.k_lambda;

lambda_cm = p.lambda_m_x(a_H2O_cc, x.T_s).*p.k_lambda;


% -----------------------------------------------------------------------------------
% help variables
xi_H2O_m  = p.xi_H2O(lambda_m);
xi_H2O_am = p.xi_H2O(lambda_am);
xi_H2O_cm = p.xi_H2O(lambda_cm);
xi_H_am   = p.xi_H(lambda_am);
xi_H_cm   = p.xi_H(lambda_cm);


% -------------------------------------------------------------------------------------- 
% gradients of the chemical potentials (29), (30)

% overall, anodic, and cathodic for water
grad_mu_H2O_m  = (p.R * x.T_s)./(0.5*xi_H2O_cm + 0.5*xi_H2O_am) .* (xi_H2O_cm - xi_H2O_am)/(p.delta_m);
grad_mu_H2O_am = (p.R * x.T_s)./(xi_H2O_m)                      .* (xi_H2O_m - xi_H2O_am) /(0.5*p.delta_m);
grad_mu_H2O_cm = (p.R * x.T_s)./(xi_H2O_m)                      .* (xi_H2O_m - xi_H2O_cm) /(0.5*p.delta_m);
% overall, anodic, and cathodic for H+
grad_mu_H_m  =   (p.R * x.T_s)./(0.5*xi_H_cm + 0.5*xi_H_am)     .* (xi_H_cm - xi_H_am)/(p.delta_m) + p.F * (-x.delta_Phi_m/p.delta_m);
grad_mu_H_am =   grad_mu_H_m;
grad_mu_H_cm =  -grad_mu_H_m;

% -------------------------------------------------------------------------------------- 
% water fluxes through the membrane (28)

% boundary fraction of total thickness
eps_m = 0.01; 

% k=1/resistance
% c_H2O = lambda * p.rho_m(lambda) * p.x_m(lambda)
k12_m  = (p.k_t.*p.t_w(lambda_m)  .*p.k_kappa.*p.kappa(lambda_m, x.T_s))  / p.F^2;
k12_am = (p.k_t.*p.t_w(lambda_am) .*p.k_kappa.*p.kappa(lambda_am, x.T_s)) / p.F^2;
k12_cm = (p.k_t.*p.t_w(lambda_cm) .*p.k_kappa.*p.kappa(lambda_cm, x.T_s)) / p.F^2;
k11_m  =  p.k_D_w.*p.D_w(lambda_m, x.T_s)  .*lambda_m   .* p.rho_m(lambda_m)  .* p.x_m(lambda_m)  ./ (p.R * x.T_s);
k11_am =  p.k_D_w.*p.D_w(lambda_am, x.T_s) .*lambda_am  .* p.rho_m(lambda_am) .* p.x_m(lambda_am) ./ (p.R * x.T_s);
k11_cm =  p.k_D_w.*p.D_w(lambda_cm, x.T_s) .*lambda_cm  .* p.rho_m(lambda_cm) .* p.x_m(lambda_cm) ./ (p.R * x.T_s);
k22_m  =  p.k_kappa.*p.kappa(lambda_m, x.T_s)  / p.F^2;
k22_am =  p.k_kappa.*p.kappa(lambda_am, x.T_s) / p.F^2;
k22_cm =  p.k_kappa.*p.kappa(lambda_cm, x.T_s) / p.F^2;

% k = 1/sum_i(eps_i * resistance_i) with resistance_i = 1/k_i 
% anodic flux
u.n_dot_H2O_am = -1./( (1-eps_m)./k12_m + eps_m./k12_am ) .* grad_mu_H_am ...
                 -1./( (1-eps_m)./k11_m + eps_m./k11_am ) .* grad_mu_H2O_am;

% cathodic flux
u.n_dot_H2O_cm = -1./( (1-eps_m)./k12_m + eps_m./k12_cm ) .* grad_mu_H_cm ...
                 -1./( (1-eps_m)./k11_m + eps_m./k11_cm ) .* grad_mu_H2O_cm;

% -------------------------------------------------------------------------------------- 
% proton flux through membrane (35)

n_dot_H = -1./( (1-2*eps_m)./k22_m + eps_m./k22_am + eps_m./k22_cm ) .* grad_mu_H_m ...
          -1./( (1-2*eps_m)./k12_m + eps_m./k12_am + eps_m./k12_cm ) .* grad_mu_H2O_m;

% -------------------------------------------------------------------------------------- 
% electrical current through the membrane (34)

rhs.i_m = x.i_m / p.F - n_dot_H;

% -------------------------------------------------------------------------------------- 
% time derivative of water content in the membrane (24), (25)
% d(c_H2O)/dt = sum(ndot_H2O)/delta_m    (for constant delta_m)
% d(lambda*rho_m*x_m)/dt = sum(ndot_H2O)/delta_m
% rho_m*x_m*d(lambda)/dt + lambda*d(rho_m*x_m)/dt = sum(ndot_H2O)/delta_m
% rho_m*x_m*d(lambda)/dt + lambda*d(rho_m*x_m)/d(lambda)*d(lambda)/dt = sum(ndot_H2O)/delta_m
% => d(lambda)/dt = ( (rho_m*x_m + lambda*d(rho_m*x_m)/d(lambda)) )^-1 * sum(ndot_H2O)/delta_m

rhs.lambda_m = (p.rho_H2O_l + p.M_H2O*lambda_m*p.rho_m_dry*p.x_m_dry).^2./(p.rho_H2O_l^2*p.rho_m_dry*p.x_m_dry) .* (u.n_dot_H2O_am + u.n_dot_H2O_cm) ./ p.delta_m ;

end
