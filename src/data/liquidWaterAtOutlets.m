function [n_dot_H2Ol_c_out,n_dot_H2Ol_a_out] = liquidWaterAtOutlets(in,varargin)
%liquidWaterAtOutlets calculates liquid water flows at the gas outlets
%based on test bench input struct 'in' using rough assumptions
% 
%   [n_dot_H2Ol_c_out,n_dot_H2Ol_a_out] = liquidWaterAtOutlets(in)
%      assumes that T_So_CL = T_Si_CL
% 
%   [n_dot_H2Ol_c_out,n_dot_H2Ol_a_out] = liquidWaterAtOutlets(in,U_cell)
%      calculates T_So_CL from a simplified energy balance using a constant
%      average cell voltage U_cell
% 
%   [n_dot_H2Ol_c_out,n_dot_H2Ol_a_out] = liquidWaterAtOutlets(in,@U_cell)
%      calculates T_So_CL from a simplified energy balance using an average
%      U_cell-I_cell curve U_cell(I_cell)
% 
% output quantities:
%   n_dot_H2Ol_c_out  mol s^-1 for entire stack
%   n_dot_H2Ol_a_out  mol s^-1 for entire stack
% 
% argument 'in':
% structure containing the
% test bench input quantities ('_set' for actual set points) as fields:
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

% load parameter struct, set additional parameters
p = mod_param_PEMFC();
p.h_H2O_l_ref = -285830;    % J mol^-1 enthaply of formation for liquid water at 298.15 K
p.rho_H2O_l = 997;          % kg m^-3 density of water at 298.15 K
p.deltah_vap_H2O = p.h_H2O_g_ref-p.h_H2O_l_ref; % in J mol^-1, enthalpy of vaporization at 298.15 K
p.p_ref  = 101325;          % 101325 Pa, reference pressure
p.T_ref  = p.T_ref;         % 298.15 K, reference temperature
p.p_norm = 101325;          % in Pa, norm pressure for norm volumes
p.T_norm = 273.15;          % in K, norm temperature for norm volumes

% stack geometry
N_cells = 20; % number of cells in the stack

% dew points at stack inlets
DPT_Si_A = in.DPT_Si_A + 273.15; % in K
DPT_Si_C = in.DPT_Si_C + 273.15; % in K

% vapor partial pressure at stack inlets
psat_a_in = p.p_sat(DPT_Si_A);
psat_c_in = p.p_sat(DPT_Si_C); 

% pressure in humidifier:
%   is actually pressure at stack inlet + 1m of pipe (10-50mbar?, 100-500Pa?),
%   -> the pipe is not included in the model
p_a_in = in.p_So_A *1e5 ; % in Pa, p_Si_A \approx p_So_A_set
p_c_in = in.p_So_C *1e5 ; % in Pa, p_Si_C \approx p_So_A_set

% approximation of coolant temperature at coolant outlet
if isempty(varargin)
% 1. 
T_So_CL = in.T_Si_CL;

elseif isnumeric(varargin{1})

% 2. heavily simplified model for energy balance with constant voltage(!):
%    direct calculation of produced heat with thermodynamic voltage and
%    measured voltage and current (condensation not included! (at reference
%    pressure and temperature, no effects of mechanical work)
U_cell = varargin{1};
Wth_voltage = (-p.h_H2O_g_ref/2/p.F*N_cells-N_cells*U_cell).*in.I_S;
T_So_CL = in.T_Si_CL + Wth_voltage ./ (in.FN_Si_CL*1/60/1000 *p.rho_H2O_l * p.c_p_H2O_l/p.M_H2O) ;

elseif isa(varargin{1}, 'function_handle')

% 3. heavily simplified model for energy balance with U_cell-I_cell curve:
%    direct calculation of produced heat with thermodynamic voltage and
%    measured voltage and current (condensation not included! (at reference
%    pressure and temperature, no effects of mechanical work)
U_cell = varargin{1};
Wth_voltage = (-p.h_H2O_g_ref/2/p.F*N_cells-N_cells*U_cell(in.I_S)).*in.I_S;
T_So_CL = in.T_Si_CL + Wth_voltage ./ (in.FN_Si_CL*1/60/1000 *p.rho_H2O_l * p.c_p_H2O_l/p.M_H2O) ;

end

% temperatures at stack inlets and outlets
T_So_C =    T_So_CL+273.15; % in K, T_So_C \approx T_So_CL
T_So_A = in.T_Si_CL+273.15; % in K, T_So_A \approx T_Si_CL

% molar flows at stack inlet
%   assumption: dry air/hydrogen before humidifier
n_dot_O2_c_in  =                          0.21 * p.p_norm/p.R/p.T_norm .* (in.FN_Si_Air_C /1000/60); % mol/s entire stack
n_dot_N2_c_in  =                          0.79 * p.p_norm/p.R/p.T_norm .* (in.FN_Si_Air_C /1000/60); % mol/s entire stack
n_dot_H2O_c_in = psat_c_in./(p_c_in-psat_c_in) * p.p_norm/p.R/p.T_norm .* (in.FN_Si_Air_C /1000/60); % mol/s entire stack
n_dot_H2_a_in  =                                 p.p_norm/p.R/p.T_norm .* (in.FN_Si_H2_A  /1000/60); % mol/s entire stack
n_dot_H2O_a_in = psat_a_in./(p_a_in-psat_a_in) * p.p_norm/p.R/p.T_norm .* (in.FN_Si_H2_A  /1000/60); % mol/s entire stack

% mass balances with reaction
% assumptions: no diffusion or drag through the membrane, no short circuit currents
%              no side reactions, steady-state 
n_dot_O2_c_out  = n_dot_O2_c_in  - in.I_S/4/p.F * N_cells; % in mol/s
n_dot_N2_c_out  = n_dot_N2_c_in; % in mol/s
n_dot_H2O_c_out = n_dot_H2O_c_in + in.I_S/2/p.F * N_cells; % in mol/s
n_dot_H2O_a_out = n_dot_H2O_a_in; % in mol/s
n_dot_H2_a_out  = n_dot_H2_a_in  - in.I_S/2/p.F * N_cells; % in mol/s

% maximum transport capacity of water vapor at stack outlets
n_dot_H2O_c_out_max = p.p_sat(T_So_C)./(in.p_So_C*1e5 - p.p_sat(T_So_C)) .* (n_dot_O2_c_out+n_dot_N2_c_out); % in mol/s
n_dot_H2O_a_out_max = p.p_sat(T_So_A)./(in.p_So_A*1e5 - p.p_sat(T_So_A)) .* (n_dot_H2_a_out); % in mol/s

% liquid water at stack outlets
n_dot_H2Ol_c_out = max(0, n_dot_H2O_c_out - n_dot_H2O_c_out_max); % in mol/s
n_dot_H2Ol_a_out = max(0, n_dot_H2O_a_out - n_dot_H2O_a_out_max); % in mol/s

end
