%clc;
clear;

%% Base value
BaseValue();

%% AC link

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
dw = ((Pr-p)*Dw + W0 - w)*wf;
ddelta = w - wg;


%% Caculate the state matrix
state = [idr; iqr; ed; eq; id; iq; vd; vq; igd; igq; w; delta];
f_xu = [didr; diqr; ded; deq; did; diq; dvd; dvq; digd; digq; dw; ddelta];

Amat = jacobian(f_xu,state);

%% Set numerical number
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
W0 = 0;

wv_ac = 250*2*pi;
kpv_ac = Cf*wv_ac;
kiv_ac = Cf*wv_ac^2/4*100;

wi_ac = 1000*2*pi;
kpi_ac = Lf*wi_ac;
kii_ac = Lf*(wi_ac^2)/4;

vd_ref = 1;
vq_ref = 0;


%% Set steady state of state variables
% idr = 0.6234;
% iqr = 0.06925;
% ed = 1.004;
% eq = 0.03765;
% id = 0.6234;
% iq = 0.06925;
% vd = 0.9996;
% vq = 0.02709;
% igd = 0.6239;
% igq = 0.0511;
% w = 0.95 * Wbase;
% delta = 0.3066;
idr = -0.0005;
iqr = 0.01944;
ed = 0.9992;
eq = 0.0165;
id = -0.0005;
iq = 0.01944;
vd = 0.9996;
vq = 0.02843;
igd = 0.0000;
igq = 0.0000;
w = 1 * Wbase;
delta = 0.001658;

%% Replace symbolic by numerical number
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