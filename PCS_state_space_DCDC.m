clc;
clear;

%% DC link

% state variables
syms i_bat v_dc i_bat_ref duty_cycle;

% other parameters
syms C_dc L_dc; % filter parameters
syms v_ocv R_bat; % battery parameters
syms kpv_dc kiv_dc kpi_dc kii_dc; % Controller parameters
syms v_dc_ref; % reference dc voltage
syms i_o; % output current

%% DC link equation

v_bat = v_ocv - R_bat * i_bat;

% System equation

% Inductance
di_bat = (v_bat - (1 - duty_cycle) * v_dc) / L_dc;
% Capacitor
dv_dc = ((1 - duty_cycle) * i_bat - i_o) / C_dc;

% Voltage PI controller : i_bat_ref = kpv_dc * (v_dc_ref - v_dc) + kiv_dc * integral(v_dc_ref - v_dc)
di_bat_ref = - kpv_dc * dv_dc + kiv_dc * (v_dc_ref - v_dc);
% % Current PI controller: duty_cycle = kpi_dc * (i_bat_ref - i_bat) + kii_dc * integral(i_bat_ref - i_bat) + 1  
dduty_cycle = kpi_dc * (di_bat_ref - di_bat) + kii_dc * (i_bat_ref - i_bat);


%% Caculate the state matrix
state_dc = [i_bat; v_dc; i_bat_ref; duty_cycle];
f_xu_dc = [di_bat; dv_dc; di_bat_ref; dduty_cycle];


Amat_dc = jacobian(f_xu_dc,state_dc);

%% Set numerical number of parameters
C_dc = 0.224;
L_dc = 0.000125;

v_ocv = 0.5989;
R_bat = 0.0625;

wv_dc = 50 * 2 * pi;
kpv_dc = C_dc*wv_dc;
kiv_dc = C_dc*wv_dc^2/4;

wi_dc = 500 * 2 * pi;
kpi_dc = L_dc*wi_dc;
kii_dc = L_dc*(wi_dc^2)/4;

v_dc_ref = 1;
i_o = 0.25;

%% Set steady state of state variables
i_bat = 0.4174;
v_dc = 1;
i_bat_ref = 0.4174;
duty_cycle = 0.4011;


%% Replace symbolic by numerical number
Amat_dc = subs(Amat_dc,'C_dc',C_dc);
Amat_dc = subs(Amat_dc,'L_dc',L_dc);

Amat_dc = subs(Amat_dc,'v_ocv',v_ocv);
Amat_dc = subs(Amat_dc,'R_bat',R_bat);

Amat_dc = subs(Amat_dc,'kpi_dc',kpi_dc);
Amat_dc = subs(Amat_dc,'kii_dc',kii_dc);
Amat_dc = subs(Amat_dc,'kpv_dc',kpv_dc);
Amat_dc = subs(Amat_dc,'kiv_dc',kiv_dc);

Amat_dc = subs(Amat_dc,'v_dc_ref',v_dc_ref);
Amat_dc = subs(Amat_dc,'i_o',i_o);

Amat_dc = subs(Amat_dc,'i_bat',i_bat);
Amat_dc = subs(Amat_dc,'v_dc',v_dc);
Amat_dc = subs(Amat_dc,'i_bat_ref',i_bat_ref);
Amat_dc = subs(Amat_dc,'duty_cycle',duty_cycle);

Amat_dc = double(Amat_dc);

%% plot
EigVec = eig(Amat_dc);
EigVecHz = EigVec/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,9999);
