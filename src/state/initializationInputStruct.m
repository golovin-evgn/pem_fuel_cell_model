function [data, info] = initializationInputStruct(in)
% initializationInputStruct(initialInputs) 
% input trajectory for simulation to steady state for test bench inputs 'in'
%
% 'FN_Si_CL' coolant flow not used in this version


info.variableNames =      {'p_So_C'    ;'FN_Si_Air_C'    ;'DPT_Si_C'    ;'T_Si_C'    ;'p_So_A'    ;'FN_Si_H2_A'    ;'DPT_Si_A'    ;'T_Si_A'    ;'T_Si_CL'    ;'FN_Si_CL' ;'I_S'    ;'p_Si_A'    ;'p_Si_C'}; 
info.variableUnits =      {'bara'      ;'Nl/min'         ;'°CTP'        ;'°C'        ;'bara'      ;'Nl/min'        ;'°CTP'        ;'°C'        ;'°C'         ;'l/min'    ;'A'      ;'bara'      ;'bara'}; % units for pressure from measurement files not correct
info.variableDimensions = [1           ;1                ;1             ;1           ;1           ;1               ;1             ;1           ;1            ;1          ;1        ;1           ;1];
info.timeUnits = {'s'};
info.fileName = {'initializationInputStruct'};


data.time = [0 .1 10 800 20000 30000].';

data.data = [1.01325    1               15.8        25          1.01325     1              15.8        25          25          0 0         1.01325     1.01325 ;... % reference temperature 25°C, small flows, (close to) zero current, 1atm ->bar
             1.01325    1               15.8        in.T_Si_C   1.01325     1              15.8        in.T_Si_A   in.T_Si_CL  0 5e-3      1.01325     1.01325 ;... % heat up first to avoid flooding (drying out is okay in the model)
             1.01325    in.FN_Si_Air_C  15.8        in.T_Si_C   1.01325     in.FN_Si_H2_A  15.8        in.T_Si_A   in.T_Si_CL  0 5e-3      1.01325     1.01325 ;... % increase gas flows to setpoints
             in.p_So_C  in.FN_Si_Air_C  in.DPT_Si_C in.T_Si_C   in.p_So_A   in.FN_Si_H2_A  in.DPT_Si_A in.T_Si_A   in.T_Si_CL  0 5e-3      in.p_Si_A   in.p_Si_C ;... % increase pressures and dew points to setpoints
             in.p_So_C  in.FN_Si_Air_C  in.DPT_Si_C in.T_Si_C   in.p_So_A   in.FN_Si_H2_A  in.DPT_Si_A in.T_Si_A   in.T_Si_CL  0 in.I_S    in.p_Si_A   in.p_Si_C ;... % increase current to set point
             in.p_So_C  in.FN_Si_Air_C  in.DPT_Si_C in.T_Si_C   in.p_So_A   in.FN_Si_H2_A  in.DPT_Si_A in.T_Si_A   in.T_Si_CL  0 in.I_S    in.p_Si_A   in.p_Si_C ]; % keep inputs for a while to ensure steady state


