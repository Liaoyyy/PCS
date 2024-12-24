% This code proves the conclusion: 
% Nonlinear transformation between two state variables will change the eigenvalues of a nonlinear system!!
% Please refer to the nonlinear.md for theoretical proof

% The following is an example of a Grid Following converter

clc;clear;

% parameters
syms Lg Rg Lf Cf Rf
syms kp_pll ki_pll
syms kpi kii idr iqr
syms wg
syms vgD vgQ

% state variables 1 -- based on grid side
syms idi iqi iD iQ vD vQ igD igQ vqi delta

% f1
id = iD * cos(delta) + iQ * sin(delta);
iq = - iD * sin(delta) + iQ * cos(delta);
vq = -vD*sin(delta) + vQ*cos(delta);

dvqi = vq;
ddelta = (kp_pll*vq + ki_pll*vqi) - wg;

didi = idr - id;
diqi = iqr - iq;
ed = kpi*didi + kii*idi;
eq = kpi*diqi + kii*iqi;

eD = ed*cos(delta) - eq*sin(delta);
eQ = ed*sin(delta) + eq*cos(delta);

diD = (eD - vD + wg*Lf*iQ - Rf*iD)/Lf;
diQ = (eQ - vQ - wg*Lf*iD - Rf*iQ)/Lf;
dvD = (iD-igD + wg*Cf*vQ)/Cf;
dvQ = (iQ-igQ - wg*Cf*vD)/Cf;
digD = (vD - vgD + wg*Lg*igQ - Rg*igD)/Lg;
digQ = (vQ - vgQ - wg*Lg*igD - Rg*igQ)/Lg;


x1 = [idi; iqi; iD; iQ; vD; vQ; igD; igQ; vqi; delta];
f1 = [didi; diqi; diD; diQ; dvD; dvQ; digD; digQ; dvqi; ddelta];

% state variables 2 -- based on controller side
syms id iq vd vq igd igq;

% f2
vgd = vgD*cos(delta) + vgQ*sin(delta);
vgq = -vgD*sin(delta) + vgQ*cos(delta);

dvqi = vq;
w = (kp_pll*vq + ki_pll*vqi);
ddelta = w - wg;

didi = idr - id;
diqi = iqr - iq;
ed = kpi*didi + kii*idi;
eq = kpi*diqi + kii*iqi;

did = (ed - vd + w*Lf*iq - Rf*id)/Lf;
diq = (eq - vq - w*Lf*id - Rf*iq)/Lf;
dvd = (id-igd + w*Cf*vq)/Cf;
dvq = (iq-igq - w*Cf*vd)/Cf;
digd = (vd - vgd + w*Lg*igq - Rg*igd)/Lg;
digq = (vq - vgq - w*Lg*igd - Rg*igq)/Lg;

x2 = [idi; iqi; id; iq; vd; vq; igd; igq; vqi; delta];
f2 = [didi; diqi; did; diq; dvd; dvq; digd; digq; dvqi; ddelta];

% x2 = h(x1)
id = iD * cos(delta) + iQ * sin(delta);
iq = - iD * sin(delta) + iQ * cos(delta);
vd = vD * cos(delta) + vQ * sin(delta);
vq = - vD * sin(delta) + vQ * cos(delta);
igd = igD * cos(delta) + igQ * sin(delta);
igq = - igD * sin(delta) + igQ * cos(delta);
h = [idi; iqi; id; iq; vd; vq; igd; igq; vqi; delta];

% x1 = h^{-1}(x2)
syms id iq vd vq igd igq
iD = id*cos(delta) - iq*sin(delta);
iQ = id*sin(delta) + iq*cos(delta);
vD = vd*cos(delta) - vq*sin(delta);
vQ = vd*sin(delta) + vq*cos(delta);
igD = igd*cos(delta) - igq*sin(delta);
igQ = igd*sin(delta) + igq*cos(delta);
h_inverse = [idi; iqi; iD; iQ; vD; vQ; igD; igQ; vqi; delta];


%% step 1 : Prove the transformation relationship between f1 and f2
f1_sub = f1;
f1_sub = subs(f1_sub,'iD', h_inverse(3));
f1_sub = subs(f1_sub,'iQ', h_inverse(4));
f1_sub = subs(f1_sub,'vD', h_inverse(5));
f1_sub = subs(f1_sub,'vQ', h_inverse(6));
f1_sub = subs(f1_sub,'igD', h_inverse(7));
f1_sub = subs(f1_sub,'igQ', h_inverse(8));

n = length(x1);
h_diff = [];
h_inverse_diff = [];
% get dx2/dx1 = h'(x1) dx1/dx2 = (h^{-1}(x2))'
for i = 1:n
    h_diff = [h_diff, diff(h,x1(i))];
    h_inverse_diff = [h_inverse_diff, diff(h_inverse, x2(i))];
end
h_diff_sub = h_diff;
h_diff_sub = subs(h_diff_sub,'iD', h_inverse(3));
h_diff_sub = subs(h_diff_sub,'iQ', h_inverse(4));
h_diff_sub = subs(h_diff_sub,'vD', h_inverse(5));
h_diff_sub = subs(h_diff_sub,'vQ', h_inverse(6));
h_diff_sub = subs(h_diff_sub,'igD', h_inverse(7));
h_diff_sub = subs(h_diff_sub,'igQ', h_inverse(8));
h_diff_sub = simplify(h_diff_sub);
h_inverse_diff = simplify(h_inverse_diff); % h_diff_sub * h_inverse_diff = Identity Matrix

% f2_cal = f2
f2_cal = simplify(h_diff_sub * f1_sub);
f2 = simplify(f2);

%% step 2 : Prove that the nonlinear transformation results in the slight difference between two A matrixs

% case 1 : use state variables 1
Amat1 = jacobian(f1,x1);

% case 2 : use state variables 2
Amat2 = jacobian(f2,x2);

% calculate Amat2 throught nonlinear tranformation
% Amat2_cal = d(h'(h^{-1}(x2))*f1(h^{-1}(x2),u))/dx2
% define 
%       g1(x2) = h'(h^{-1}(x2)) 
%       g2(x2) = f1(h^{-1}(x2),u)
% then
%       Amat2_cal = d(g1(x2)*g2(x2))/dx2

g1 = h_diff_sub; % size: 10 x 10
g2 = f1_sub; % size: 10 x 1
g1_diff = [];
g2_diff = [];
for i=1:n
    g1_diff = [g1_diff, diff(g1,x2(i))];
    g2_diff = [g2_diff, diff(g2,x2(i))];
end

% Amat2_cal = g1_diff * g2 + g1 * g2_diff;
Amat2_cal = simplify(g1 * g2_diff);
Amat1_sub = Amat1;
Amat1_sub = subs(Amat1_sub,'iD', h_inverse(3));
Amat1_sub = subs(Amat1_sub,'iQ', h_inverse(4));
Amat1_sub = subs(Amat1_sub,'vD', h_inverse(5));
Amat1_sub = subs(Amat1_sub,'vQ', h_inverse(6));
Amat1_sub = subs(Amat1_sub,'igD', h_inverse(7));
Amat1_sub = subs(Amat1_sub,'igQ', h_inverse(8));
Amat1_sub = simplify(Amat1_sub);
Amat2_cal2 = simplify(g1 * Amat1_sub * h_inverse_diff);

% set numerical number to verify
Wbase = 2 * 50 * pi;
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;
Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

wpll = 10*2*pi;
kp_pll = wpll;
ki_pll = wpll^2/4;

wi = 500*2*pi;
kpi = Lf*wi;
kii = Lf*(wi^2)/4;

P = 0.5;
Q = 0.2;
idr = P / 1.0817;
iqr = -Q / 1.0817;
wg = Wbase;

vgD = 1;
vgQ = 0;

% steady state
% idi; 
% iqi; 
id = idr; 
iq = iqr; 
vd = 1.0817; 
vq = 0;  
igd = idr; 
igq = iqr; 
vqi = (wg - kp_pll * vq)/ki_pll; 
delta = 3.6268/180*pi;
% idi = (vd + Rf * id - wg * Lf * iq)/kii; 
% iqi = (vq + Rf * iq + wg * Lf * id)/kii;
% Attention！ The steady-state values used in the GFL_Steady are different, as follows 
idi = (vd)/kii;
iqi = (vq)/kii;

% Replace symbolic by numerical number
Amat2 = subs(Amat2,'Cf',Cf);
Amat2 = subs(Amat2,'Lf',Lf);
Amat2 = subs(Amat2,'Rf',Rf);
Amat2 = subs(Amat2,'Lg',Lg);
Amat2 = subs(Amat2,'Rg',Rg);

Amat2 = subs(Amat2,'kp_pll',kp_pll);
Amat2 = subs(Amat2,'ki_pll',ki_pll);
Amat2 = subs(Amat2,'kpi',kpi);
Amat2 = subs(Amat2,'kii',kii);

Amat2 = subs(Amat2,'idr',idr);
Amat2 = subs(Amat2,'iqr',iqr);

Amat2 = subs(Amat2,'wg',wg);

Amat2 = subs(Amat2,'vgD',vgD);
Amat2 = subs(Amat2,'vgQ',vgQ);

Amat2 = subs(Amat2,'id',id);
Amat2 = subs(Amat2,'iq',iq);
Amat2 = subs(Amat2,'vd',vd);
Amat2 = subs(Amat2,'vq',vq);
Amat2 = subs(Amat2,'igd',igd);
Amat2 = subs(Amat2,'igq',igq);
Amat2 = subs(Amat2,'vqi',vqi);
Amat2 = subs(Amat2,'delta',delta);

Amat2 = double(Amat2);

EigVec1 = eig(Amat2);
EigVecHz1 = EigVec1/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz1,ZoomInAxis,100);
% The eigenvalues should be conincide with those in GFL_Swing.h

Amat2_cal = subs(Amat2_cal,'Cf',Cf);
Amat2_cal = subs(Amat2_cal,'Lf',Lf);
Amat2_cal = subs(Amat2_cal,'Rf',Rf);
Amat2_cal = subs(Amat2_cal,'Lg',Lg);
Amat2_cal = subs(Amat2_cal,'Rg',Rg);
Amat2_cal = subs(Amat2_cal,'kp_pll',kp_pll);
Amat2_cal = subs(Amat2_cal,'ki_pll',ki_pll);
Amat2_cal = subs(Amat2_cal,'kpi',kpi);
Amat2_cal = subs(Amat2_cal,'kii',kii);
Amat2_cal = subs(Amat2_cal,'idr',idr);
Amat2_cal = subs(Amat2_cal,'iqr',iqr);
Amat2_cal = subs(Amat2_cal,'wg',wg);
Amat2_cal = subs(Amat2_cal,'vgD',vgD);
Amat2_cal = subs(Amat2_cal,'vgQ',vgQ);

Amat2_cal = subs(Amat2_cal,'idi',idi);
Amat2_cal = subs(Amat2_cal,'iqi',iqi);
Amat2_cal = subs(Amat2_cal,'id',id);
Amat2_cal = subs(Amat2_cal,'iq',iq);
Amat2_cal = subs(Amat2_cal,'vd',vd);
Amat2_cal = subs(Amat2_cal,'vq',vq);
Amat2_cal = subs(Amat2_cal,'igd',igd);
Amat2_cal = subs(Amat2_cal,'igq',igq);
Amat2_cal = subs(Amat2_cal,'vqi',vqi);
Amat2_cal = subs(Amat2_cal,'delta',delta);

Amat2_cal = double(Amat2_cal);

EigVec2 = eig(Amat2_cal);
EigVecHz2 = EigVec2/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz2,ZoomInAxis,100);
% The eigenvalues should be conincide with those in GFL_Steady.h
