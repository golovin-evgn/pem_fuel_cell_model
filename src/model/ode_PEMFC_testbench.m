function rhs_vec = ode_PEMFC_testbench(t, x_vec, u_vec_testbench, p_var)
%%% wrapper for testbench inputs of PEMFC model
% 
% test bench input quantities:
%  p_So_C       bar
%  FN_Si_Air_C  Nl/min      dry, 21% O2, 79% N2, 101325 Pa, 273.15 K
%  DPT_Si_C     °C
%  T_Si_C       °C
%  p_So_A       bar
%  FN_Si_H2_A   Nl/min      dry, 101325 Pa, 273.15 K, no nitrogen
%  DPT_Si_A     °C
%  T_Si_A       °C
%  T_Si_CL      °C                 
%  FN_Si_CL     liter/min   water  
%  I_S          A

% -------------------------------------------------------------------------------------- 
% load parameters

% use default parameters as long as adaptable parameters are not required
p = mod_param_PEMFC();

% --------------------------------------------------------------------------------------
% get input pressures from state vector

x = p.state2struct(x_vec);

if p.counter_flow
    p_a_in = x.p_a(end);
else
    p_a_in = x.p_a(1); 
end

p_c_in = x.p_c(1);

% -------------------------------------------------------------------------------------- 
% include input pressures and convert to model inputs 

u_struct_testbench = p.testbench2struct([u_vec_testbench; p_a_in*1e-5; p_c_in*1e-5]); % pressures from Pa to bar

u_struct = testbench2model(u_struct_testbench,p);

u_vec = p.inputs2vec(u_struct);

% -------------------------------------------------------------------------------------- 
% call PEMFC model with model inputs 

if nargin == 4
    rhs_vec = ode_PEMFC(t, x_vec, u_vec, p_var);
else
    rhs_vec = ode_PEMFC(t, x_vec, u_vec);
end

end 

