clc;clear;

%% Base value
BaseValue();

%% DC link variables

% state variables
syms i_bat v_dc i_bat_ref duty_cycle;

% other parameters
syms C_dc L_dc; % filter parameters
syms v_ocv R_bat SOC; % battery parameters
syms kpv_dc kiv_dc kpi_dc kii_dc; % Controller parameters
syms v_dc_ref; % reference dc voltage


%% AC link variables

% state variables
syms idr iqr ed eq id iq vd vq igd igq w delta

% other variables
syms Lg Rg % Line impedance
syms Lf Cf Rf % LC filter
syms vgD vgQ wg % Grid parameters
syms wf Dw Pr W0 % Droop control 
syms vd_ref vq_ref % reference ac voltage
syms kpv_ac kiv_ac kpi_ac kii_ac % controller parameters


%% AC link equation

vgd = vgD*cos(delta) + vgQ*sin(delta); 
vgq = -vgD*sin(delta) + vgQ*cos(delta);

q = vq*igd - vd*igq;
p = vd*id + vq*iq;


% System equation

% LC filter
did = (ed - vd + w * Lf * iq - Rf * id) / Lf;
diq = (eq - vq - w * Lf * id - Rf * iq) / Lf;
dvd = (id - igd + w * Cf * vq) / Cf;
dvq = (iq - igq - w * Cf * vd) / Cf;
digd = (vd - vgd + w * Lg * igq - Rg * igd) / Lg;
digq = (vq - vgq - w * Lg * igd - Rg * igq) / Lg;

% Voltage PI controller
didr = - kpv_ac * dvd + kiv_ac * (vd_ref - vd);
diqr = - kpv_ac * dvq + kiv_ac * (vq_ref - vq);

% Current PI controller
ded = kpi_ac * (didr - did) + kii_ac * (idr - id);
deq = kpi_ac * (diqr - diq) + kii_ac * (iqr - iq);

% Droop control
dw = ((Pr-p)*Dw/SOC + W0 - w)*wf;
ddelta = w - wg;

%% DC link equation

v_bat = v_ocv - R_bat * i_bat;

% System equation

% Inductance
di_bat = (v_bat - (1 - duty_cycle) * v_dc) / L_dc;
% Capacitor
dv_dc = ((1 - duty_cycle) * i_bat - p / v_dc) / C_dc;

% Voltage PI controller
di_bat_ref = - kpv_dc * dv_dc + kiv_dc * (v_dc_ref - v_dc);

% % Current PI controller
dduty_cycle = kpi_dc * (di_bat_ref - di_bat) + kii_dc * (i_bat_ref - i_bat);

%% Caculate the state matrix
state_dc = [i_bat; v_dc; i_bat_ref; duty_cycle];
f_xu_dc = [di_bat; dv_dc; di_bat_ref; dduty_cycle];

state_ac = [idr; iqr; ed; eq; id; iq; vd; vq; igd; igq; w; delta];
f_xu_ac = [didr; diqr; ded; deq; did; diq; dvd; dvq; digd; digq; dw; ddelta];

state = [state_dc; state_ac];
f_xu = [f_xu_dc; f_xu_ac];

Amat = jacobian(f_xu,state);

%% Set numerical number
C_dc = 0.224;
L_dc = 0.000125;
v_ocv = 0.625;
R_bat = 0.0625;
SOC = 0.5;
wv_dc = 50 * 2 * pi;
kpv_dc = C_dc*wv_dc;
kiv_dc = C_dc*wv_dc^2/4;
wi_dc = 500 * 2 * pi;
kpi_dc = L_dc*wi_dc;
kii_dc = L_dc*(wi_dc^2)/4;
v_dc_ref = 1;
i_o = 0.25;

Lg = 0.25/100/pi*2;
Rg = 0.25/5*2;
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;
vgD = 1;
vgQ = 0;
wg = 0.95 * Wbase;
wf = 2*pi*10;
Dw = 0.08*Wbase/Sbase;
Pr = 0;
W0 = Wbase;
wv_ac = 250*2*pi;
kpv_ac = Cf*wv_ac;
kiv_ac = Cf*wv_ac^2/4*100;
wi_ac = 1000*2*pi;
kpi_ac = Lf*wi_ac;
kii_ac = Lf*(wi_ac^2)/4;
vd_ref = 1;
vq_ref = 0;


%% Set steady state of state variables
eqn = [di_bat == 0, dv_dc == 0, di_bat_ref == 0, dduty_cycle == 0, didr == 0, diqr == 0, ded == 0, deq == 0, did == 0, diq == 0, dvd == 0, dvq == 0, digd == 0, digq == 0, dw == 0, ddelta == 0];
eqn = subs(eqn,'C_dc',C_dc);
eqn = subs(eqn,'L_dc',L_dc);
eqn = subs(eqn,'v_ocv',v_ocv);
eqn = subs(eqn,'R_bat',R_bat);
eqn = subs(eqn,'SOC',SOC);
eqn = subs(eqn,'kpi_dc',kpi_dc);
eqn = subs(eqn,'kii_dc',kii_dc);
eqn = subs(eqn,'kpv_dc',kpv_dc);
eqn = subs(eqn,'kiv_dc',kiv_dc);
eqn = subs(eqn,'v_dc_ref',v_dc_ref);
eqn = subs(eqn,'Rg',Rg);
eqn = subs(eqn,'Lg',Lg);
eqn = subs(eqn,'Cf',Cf);
eqn = subs(eqn,'Lf',Lf);
eqn = subs(eqn,'Rf',Rf);
eqn = subs(eqn,'vgD',vgD);
eqn = subs(eqn,'vgQ',vgQ);
eqn = subs(eqn,'wg',wg);
eqn = subs(eqn,'Dw',Dw);
eqn = subs(eqn,'wf',wf);
eqn = subs(eqn,'Pr',Pr);
eqn = subs(eqn,'W0',W0);
eqn = subs(eqn,'kpv_ac',kpv_ac);
eqn = subs(eqn,'kiv_ac',kiv_ac);
eqn = subs(eqn,'kpi_ac',kpi_ac);
eqn = subs(eqn,'kii_ac',kii_ac);
eqn = subs(eqn,'vd_ref',vd_ref);
eqn = subs(eqn,'vq_ref',vq_ref);

init_value = [0;1;0;0;0;0;1;0;0;0;1;0;0;0;Wbase;0];

[solve_i_bat,solve_v_dc, solve_i_bat_ref, solve_duty_cycle, solve_idr,...
    solve_iqr, solve_ed, solve_eq, solve_id, solve_iq, solve_vd, solve_vq, solve_igd,...
    solve_igq, solve_w, solve_delta] = vpasolve(eqn,[i_bat v_dc i_bat_ref duty_cycle idr iqr ed eq id iq vd vq igd igq w delta],init_value);

i_bat = double(solve_i_bat);
v_dc = double(solve_v_dc);
i_bat_ref = double(solve_i_bat_ref);
duty_cycle = double(solve_duty_cycle);
idr = double(solve_idr);
iqr = double(solve_iqr);
ed = double(solve_ed);
eq = double(solve_eq);
id = double(solve_id);
iq = double(solve_iq);
vd = double(solve_vd);
vq = double(solve_vq);
igd = double(solve_igd);
igq = double(solve_igq);
w = double(solve_w);
delta = double(solve_delta);

% i_bat = 0.5297;
% v_dc = 1;
% i_bat_ref = 0.5297;
% duty_cycle = 0.4081;
% 
% 
% idr = 0.6234;
% iqr = 0.06925;
% ed = 1.001;
% eq = 0.02689;
% id = 0.3108;
% iq = 0.0678;
% vd = 0.9996;
% vq = 0.02705;
% igd = 0.3113;
% igq = 0.04949;
% w = 0.95 * Wbase;
% delta = 0.1547;

%% Replace symbolic by numerical number
Amat = subs(Amat,'C_dc',C_dc);
Amat = subs(Amat,'L_dc',L_dc);

Amat = subs(Amat,'v_ocv',v_ocv);
Amat = subs(Amat,'R_bat',R_bat);
Amat = subs(Amat,'SOC',SOC);

Amat = subs(Amat,'kpi_dc',kpi_dc);
Amat = subs(Amat,'kii_dc',kii_dc);
Amat = subs(Amat,'kpv_dc',kpv_dc);
Amat = subs(Amat,'kiv_dc',kiv_dc);
Amat = subs(Amat,'v_dc_ref',v_dc_ref);

Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);
Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);

Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);
Amat = subs(Amat,'wg',wg);

Amat = subs(Amat,'Dw',Dw);
Amat = subs(Amat,'wf',wf);

Amat = subs(Amat,'kpv_ac',kpv_ac);
Amat = subs(Amat,'kiv_ac',kiv_ac);
Amat = subs(Amat,'kpi_ac',kpi_ac);
Amat = subs(Amat,'kii_ac',kii_ac);

Amat = subs(Amat,'vd_ref',vd_ref);
Amat = subs(Amat,'vq_ref',vq_ref);


Amat = subs(Amat,'i_bat',i_bat);
Amat = subs(Amat,'v_dc',v_dc);
Amat = subs(Amat,'i_bat_ref',i_bat_ref);
Amat = subs(Amat,'duty_cycle',duty_cycle);

Amat = subs(Amat,'idr',idr);
Amat = subs(Amat,'iqr',iqr);
Amat = subs(Amat,'ed',ed);
Amat = subs(Amat,'eq',eq);
Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);
Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);
Amat = subs(Amat,'w',w);
Amat = subs(Amat,'delta',delta);



%% Sweep parameters
Amat = double(Amat);

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,9999);