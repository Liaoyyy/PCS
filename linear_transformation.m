% This code proves the conclusion: 
% Linear transformation between two state variables won't change the eigenvalues of a nonlinear system!!
% The theoretical proof is as follows:
% assume that 
%           x1^{dot} = f1(x1,u) ---equ(1)
%           x2^{dot} = f2(x2,u)
%           x2 = T * x1 + B => x1 = inv(T) * (x2 - B) ( B refers to a constant matrix)
% then
%           left hand side of equ(1) = inv(T) * x2
%           right hand side of equ(2) = f1(inv(T) * (x2 - B), u)
% so
%           x2^{dot} = T * f1(inv(T) * (x2 - B), u)
%                    = f2(x2, u)
% obviously
%           Amat2 = T * Amat * inv(T)


% The following is an example of a boost DC/DC converter

clc;clear;

% state variable 1
syms i_bat v_dc v_err i_err;

% state variable 2
syms i_bat_ref duty

% parameter
syms L C;
syms kpv kpi kiv kii;

% input parameter
syms i_o v_bat v_dc_ref

% x1^{dot} = f1(x1,u)
dv_err = v_dc_ref - v_dc;
out_v = kpv * dv_err + kiv * v_err;

di_err = out_v - i_bat;
out_i = kpi * di_err + kii * i_err + 1;

di_bat = (v_bat - (1-out_i) * v_dc) / L;
dv_dc = ((1-out_i)*i_bat - i_o) / C;

% transformation x2 = T * x1 + B(constant)
i_bat_ref = kpv * dv_err + kiv * v_err;
duty = kpi * di_err + kii * i_err + 1;

x1 = [i_bat; v_dc; v_err; i_err];
x2 = [i_bat; v_dc; i_bat_ref; duty];
f1 = [di_bat;dv_dc;dv_err;di_err];
f1_trans = [di_bat;dv_dc;dv_err;di_err];

T = jacobian(x2,x1);
B = expand(x2 - T * x1);

% x1 = inv(T) * (x2 -B) redefine x2
syms i_bat_ref duty;
x2 = [i_bat; v_dc; i_bat_ref; duty];
x1_trans = inv(T) * (x2 - B);

% get the x2 state function through transformation
f1_trans = subs(f1_trans, 'v_err', x1_trans(3));
f1_trans = subs(f1_trans, 'i_err', x1_trans(4));
f2_calculate = expand(T * f1_trans);

Amat2_cal = jacobian(f2_calculate,x2);


% f2
di_bat = (v_bat - (1-duty) * v_dc) / L;
dv_dc = ((1-duty)*i_bat - i_o) / C;
di_bat_ref = - kpv * dv_dc + kiv * (v_dc_ref - v_dc);
dduty = kpi * (di_bat_ref - di_bat) + kii * (i_bat_ref - i_bat);
f2 = [di_bat;dv_dc;di_bat_ref;dduty];
Amat2 = expand(jacobian(f2,x2));

% Amat2 == Amat2_cal!
disp(Amat2)
disp(Amat2_cal)



