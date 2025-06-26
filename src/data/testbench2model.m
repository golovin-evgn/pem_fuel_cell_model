function out = testbench2model(in,p)
% testbench2model(in,p) Converts test bench input variable struct 'in' to model
% input struct 'out'. The model parameter struct 'p' is required for the saturation
% vapor pressure 'p_sat(T)' and geometry.
% 
% test bench input quantities (as measured, '_set' for actual set points):
%  p_So_C       bar
%  FN_Si_Air_C  Nl/min      dry, 21% O2, 79% N2, 101325 Pa, 273.15 K
%  DPT_Si_C     °C
%  T_Si_C       °C
%  p_So_A       bar
%  FN_Si_H2_A   Nl/min      dry, 101325 Pa, 273.15 K, no nitrogen
%  DPT_Si_A     °C
%  T_Si_A       °C
%  T_Si_CL      °C                 *** (T_Si_CL °C)
%  FN_Si_CL     liter/min   water  *** (T_So_CL °C)
%  I_S          A
% 
% *** current model version assumes constant coolant temperature along channel 
%     -> flow rate FN_Si_CL not used in current model  
%     option 1: T_cool(z) = T_cool_inlet <-
%     option 2: T_cool(z) = 0.5*(T_cool_inlet+T_cool_outlet)
%     option 3: discretized dynamic model for T_cool(z), with T_cool(0)=T_cool_inlet
%               and coolant mass flow
% 
% additional test bench quantities: measured states
% 
%  p_Si_A           bar
%  p_Si_C           bar
% 
% model input quantities:
%  n_dot_H2_a_in    mol/s/m²    inlet molar fluxes: molar flow divided by inlet cross-section area
%  n_dot_H2O_a_in   mol/s/m²
%  n_dot_O2_c_in    mol/s/m²    ratio of N2 to O2 fixed by air composition
%  n_dot_H2O_c_in   mol/s/m²
%  n_dot_N2_c_in    mol/s/m²    ratio of N2 to O2 fixed by air composition
%  T_a_in           K
%  T_c_in           K
%  T_cool           K
%  p_a_out          Pa
%  p_c_out          Pa
%  I_cell           A

% test bench stack geometry
N_cells = p.N_cells; % number of cells in the stack
A_cell  = p.A_cell;  % active cell area in m^2

% number of channels for translation to one-channel model with the same
% current density
N_channels = A_cell/(p.L_x*p.L_z); % A_cell=A_channel*N_channels

% help variables for H2O molar flow:
%  humidifier in front of fuel cell inlets
psat_a_in = p.p_sat(in.DPT_Si_A + 273.15);
psat_c_in = p.p_sat(in.DPT_Si_C + 273.15); 
%  pressure in humidifier:
%   is actually pressure at inlet + 1m of pipe (10-50mbar),
%   -> the pipe is not included in the model
p_a = in.p_Si_A *1e5;
p_c = in.p_Si_C *1e5; 

% ndot*A_cross=p_norm/R/T_norm*(norm volume flow in m^3/s)
% gases are distributed equally between the cells in the stack and channels in the cells -> '/N_cells/N_channels'
helpVar = 1/1000/60 * 101325/p.R/273.15 /N_cells/N_channels;
out.n_dot_H2_a_in  =                               in.FN_Si_H2_A  *helpVar /(p.delta_a*p.L_x);
out.n_dot_H2O_a_in = psat_a_in./(p_a-psat_a_in) .* in.FN_Si_H2_A  *helpVar /(p.delta_a*p.L_x);
out.n_dot_O2_c_in  =                       0.21 .* in.FN_Si_Air_C *helpVar /(p.delta_c*p.L_x);
out.n_dot_H2O_c_in = psat_c_in./(p_c-psat_c_in) .* in.FN_Si_Air_C *helpVar /(p.delta_c*p.L_x);
out.n_dot_N2_c_in  =                       0.79 .* in.FN_Si_Air_C *helpVar /(p.delta_c*p.L_x);

out.T_a_in = in.T_Si_A + 273.15;      
out.T_c_in = in.T_Si_C + 273.15;  
% coolant temperature models
out.T_cool = in.T_Si_CL + 273.15;
% out.T_cool = 0.5*(in.T_Si_CL+in.T_So_CL) + 273.15; % mean of measured temperatures

out.p_a_out = in.p_So_A * 1e5;
out.p_c_out = in.p_So_C * 1e5;

% channels are electrically connected in parallel -> I_channel = I_cell/N_channels
% Matlab model has only one channel -> I_cell=I_channel
out.I_cell = in.I_S /N_channels; 

end

% {'p_So_C', 'FN_Si_Air_C', 'DPT_Si_C', 'T_Si_C', 'p_So_A', 'FN_Si_H2_A', 'DPT_Si_A', 'T_Si_A', 'T_Si_CL', 'FN_Si_CL', 'I_S', 'p_Si_A', 'p_Si_C'}
