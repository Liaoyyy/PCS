clear all;
clc;

%% Base Value
Fbase = 50;
Wbase = 2*pi*Fbase;
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;
Vbase_dc = 1;
Sbase_dc = 1;
Ibase_dc = Sbase_dc / Vbase_dc;
Zbase_dc = Vbase_dc / Ibase_dc;
Ybase_dc = 1 / Zbase_dc;

%% AC link parameters
% Grid voltage
syms vgD vgQ wg;

% LCL filter
syms Lf Rf Cf Lg Rg

% System states
syms w theta;
syms id iq vd vq igd igq;

% Droop Controller
syms w0 m P0 wf; % wf - the bandwidth of LPF

% Voltage loop controller
syms vd_ref vq_ref; % reference voltage
syms kpv_ac kiv_ac;
syms vd_err vq_err; % The integral of the error between reference and current value

% Current loop controller
syms kpi_ac kii_ac;
syms id_err iq_err;

%% DC link parameters
% Battery parameters
syms SOC v_OCV R_bat MaxQ; %v_OCV - the open circuit voltage of the battery, MaxQ - the rated capacity (Ah)

% LC filter
syms L_dc C_dc;  

% System states
syms i_bat v_dc;

% Voltage loop controller
syms v_dc_ref v_err;
syms kpv_dc kiv_dc;

% Current loop controller
syms i_err;
syms kpi_dc kii_dc;


%% AC link Equation

% Inverse frame transformation
vgd = vgD*cos(theta) + vgQ*sin(theta); 
vgq = -vgD*sin(theta) + vgQ*cos(theta);

% Power calculation
p = vd * igd + vq * igq;
% p = vd * id + vq * iq;
q = -vd * igq + vq * igd;

% Droop control
dw = ((w0 - w) + m / SOC *(P0-p)) * wf; % adaptive droop control
% dw = ((w0 - w) + m  *(P0-p)) / Ts;


% Angle difference
dtheta = (w - wg)*Wbase;
% dtheta = (w - wg);

% Voltage controller
dvd_err = vd_ref - vd;
dvq_err = vq_ref - vq;
id_ref = kpv_ac * dvd_err + kiv_ac * vd_err;
iq_ref = kpv_ac * dvq_err + kiv_ac * vq_err;

% Current controller
did_err = id_ref - id;
diq_err = iq_ref - iq;
vid = kpi_ac * did_err + kii_ac * id_err; % Inverter output voltage
viq = kpi_ac * diq_err + kii_ac * iq_err; 

% Inverter-side inductor
did = (vid - vd - id * Rf + w * iq * Lf) / Lf;
diq = (viq - vq - iq * Rf - w * id * Lf) / Lf;

% Capacitor
dvd = (id - igd + w * Cf * vq) / Cf;
dvq = (iq - igq - w * Cf * vd) / Cf;

% Grid-side inductor
digd = (vd - vgd + w * Lg * igq - Rg * igd) / Lg;
digq = (vq - vgq - w * Lg * igd - Rg * igq) / Lg;


%% DC link Equation

% Battery
v_bat = v_OCV - R_bat * i_bat;

% Voltage Controller
dv_err = v_dc_ref - v_dc;
i_bat_ref = kpv_dc * dv_err + kiv_dc * v_err;

% Current Controller
di_err = i_bat_ref - i_bat;
duty_cycle = kpi_dc * di_err + kii_dc * i_err + 1;

% inductance
di_bat = (v_bat - (1 - duty_cycle) * v_dc) / L_dc;

% Capacitor
i_o = (vid * id + viq * iq) / v_dc;
% i_o = p / v_dc;
dv_dc = ((1 - duty_cycle) * i_bat - i_o) / C_dc;



%% Calculate the state matrix
state_ac = [vd_err; vq_err; id_err; iq_err; id; iq; vd; vq; igd; igq; w; theta];
f_xu_ac = [dvd_err; dvq_err; did_err; diq_err; did; diq; dvd; dvq; digd; digq; dw; dtheta];

state_dc = [v_err; i_err; i_bat; v_dc];
f_xu_dc = [dv_err; di_err; di_bat; dv_dc];

state = [state_ac; state_dc];
f_xu = [f_xu_ac; f_xu_dc];

Amat_ac = jacobian(f_xu_ac, state_ac);
Amat = jacobian(f_xu, state);

%% Set numerical number

% AC link
% Cf = 0.02/Wbase;
% Lf = 0.05/Wbase;
% Rf = 0.01;
% Xg = 0.3;
% Lg = Xg/Wbase;
% Rg = Xg/5;
% wg = 0.95;
Lf = 0.05 / Wbase;
Rf = 0.01;
Cf = 0.02 / Wbase;
Lg = 2*0.00079577;
Rg = 0.1;
wg = 0.95;

% Droop
wf = 2 * pi * 10;
m = 0.08;

% Controller
vd_ref = 1;
vq_ref = 0;
wv = 250*2*pi;
kpv_ac = Cf*wv;
kiv_ac = Cf*wv^2/4 * 100;

wi_ac = 1000*2*pi;
kpi_ac = Lf*wi_ac;
kii_ac = Lf*(wi_ac^2)/4;

% Circuit parameter
vgD = 1;
vgQ = 0;

% steady-state value
vd_err = 7.9603e-5;
vq_err = 1.4764e-5;
id_err = 6.3598e-4;
iq_err = 2.4847e-5;
id = 0.3108;
iq = 0.0678;
vd = 0.9996;
vq = 0.02705;
igd = 0.3113;
igq = 0.04949;
w = 0.95;
% theta = 0.1547;
theta = 0.14;

% DC link
v_OCV = 0.625 * Vbase_dc;
R_bat = 0.0625 * Zbase_dc;
SOC = 50 / 100; 

C_dc = 0.224;
L_dc = 0.000125;

v_dc_ref = 1;

wv_dc = 50 * 2 * pi;
kpv_dc = C_dc*wv_dc;
kiv_dc = C_dc*wv_dc^2/4;

wi_dc = 500 * 2 * pi;
kpi_dc = L_dc*wi_dc;
kii_dc = L_dc*(wi_dc^2)/4;

% steady-state value
v_err = 9.5839e-5;
i_err = -0.0019;
i_bat = 0.5297;
v_dc = 1;

%% Replace symbolic by numerical number
Amat = subs(Amat,'kpi_ac',kpi_ac);
Amat = subs(Amat,'kii_ac',kii_ac);
Amat = subs(Amat,'kpv_ac',kpv_ac);
Amat = subs(Amat,'kiv_ac',kiv_ac);
Amat = subs(Amat,'kpi_dc',kpi_dc);
Amat = subs(Amat,'kii_dc',kii_dc);
Amat = subs(Amat,'kpv_dc',kpv_dc);
Amat = subs(Amat,'kiv_dc',kiv_dc);

Amat = subs(Amat,'m',m);
Amat = subs(Amat,'wf',wf);

Amat = subs(Amat,'vd_ref',vd_ref);
Amat = subs(Amat,'vq_ref',vq_ref);
Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);

Amat = subs(Amat,'v_dc_ref',v_dc_ref);

Amat = subs(Amat,'v_OCV',v_OCV);
Amat = subs(Amat,'SOC',SOC);
Amat = subs(Amat,'R_bat',R_bat);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);
Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);
Amat = subs(Amat,'C_dc',C_dc);
Amat = subs(Amat,'L_dc',L_dc);

Amat = subs(Amat,'vd_err',vd_err);
Amat = subs(Amat,'vq_err',vq_err);
Amat = subs(Amat,'id_err',id_err);
Amat = subs(Amat,'iq_err',iq_err);
Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);
Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);
Amat = subs(Amat,'w',w);
Amat = subs(Amat,'wg',wg);
Amat = subs(Amat,'theta',theta);

Amat = subs(Amat,'v_err',v_err);
Amat = subs(Amat,'i_err',i_err);
Amat = subs(Amat,'i_bat',i_bat);
Amat = subs(Amat,'v_dc',v_dc);

Amat = double(Amat);

%% plot
EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,9999);
