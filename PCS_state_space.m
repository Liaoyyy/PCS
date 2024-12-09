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
syms vgd vgq;

% LCL filter
syms Lf Rf Cf Lg Rg

% System states
syms w theta;
syms id iq vd vq igd igq;

% Droop Controller
syms w0 m P0 Ts; % Ts represents the time constant of LPF

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

% Decoupling gain



%% AC link Equation

% Inverse frame transformation
vgd = vgD*cos(theta) + vgQ*sin(theta); 
vgq = -vgD*sin(theta) + vgQ*cos(theta);

% Power calculation
% p = vd * igd + vq * igq;
p = vd * id + vq * iq;
q = -vd * igq + vq * igd;

% Droop control
% dw = ((w0 - w) + m / SOC *(P0-p))/Ts; % adaptive droop control
dw = ((w0 - w) + m  *(P0-p))/Ts;

% Angle difference
dtheta = w - wg;

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

% input power
pi = vid * id + viq * iq;


%% DC link Equation

% Battery
v_bat = v_OCV - R_bat * i_bat;
dSOC = - i_bat / 3600;

% Voltage Controller
dv_err = v_dc_ref - v_dc;
i_bat_ref = kpv_dc * dv_err + kiv_dc * v_err;

% Current Controller
di_err = i_bat_ref - i_bat;
duty_cycle = kpi_dc * di_err + kii_dc * i_err;

% inductance
di_bat = (v_bat - (1 - duty_cycle) * v_dc) / L_dc;

% Capacitor
i_o = pi / v_dc;
dv_dc = ((1 - duty_cycle) * i_bat - i_o) / C_dc;



%% Calculate the state matrix
state_ac = [vd_err; vq_err; id_err; iq_err; id; iq; vd; vq; igd; igq; w; theta];
f_xu_ac = [dvd_err; dvq_err; did_err; diq_err; did; diq; dvd; dvq; digd; digq; dw; dtheta];

state_dc = [SOC; v_err; i_err; i_bat; v_dc];
f_xu_dc = [dSOC; dv_err; di_err; di_bat; dv_dc];

state = [state_ac; state_dc];
f_xu = [f_xu_ac; f_xu_dc];

Amat_ac = jacobian(f_xu_ac, state_ac);
Amat = jacobian(f_xu, state);

%% Set numerical number

% AC link
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;
Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

% Droop
Ts = 1 / (2 * pi * 10);
m = 0.05*Wbase/Sbase;

% Controller
wv_ac = 250*2*pi;
kpv_ac = Cf*wv_ac;
kiv_ac = Cf*wv_ac^2/4*50;

wi_ac = 1000*2*pi;
kpi_ac = Lf*wi_ac;
kii_ac = Lf*(wi_ac^2)/4;

% Circuit parameter
vgD = 1;
vgQ = 0;

vd = 1;
vq = 0;
P = 0.5;
Q = 0.2;
igd = P/vd;
igq = -Q/vd;
id = igd;
iq = igq;
vd_ref = vd;
vq_ref = vq;

theta = 0/180*pi;


% DC link
v_OCV = 0.625 * Vbase_dc;
R_bat = 0.0625 * Zbase_dc;
SOC = 50 / 100; 

C_dc = 0.224;
L_dc = 0.000125;

v_dc_ref = 1;

wv_dc = 50 * 2 * pi;
kpv_dc = C_dc*wv_dc;
kiv_dc = C_dc*wv_dc^2/4*50;

wi_dc = 500 * 2 * pi;
kpi_dc = L_dc*wi_dc;
kii_dc = L_dc*(wi_dc^2)/4;


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
Amat = subs(Amat,'Ts',Ts);

Amat = subs(Amat,'vd_ref',vd_ref);
Amat = subs(Amat,'vq_ref',vq_ref);
Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);
Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);
Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);

Amat = subs(Amat,'v_dc_ref',v_dc_ref);

Amat = subs(Amat,'w',Wbase);
Amat = subs(Amat,'wg',Wbase);
Amat = subs(Amat,'theta',theta);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);
Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

