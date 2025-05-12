clc;
clear;

%% Base value
BaseValue();
enable_state_participation = 0;
enable_para_participation = 0;

%% DC link variables

% state variables
syms i_bat v_dc i_bat_ref duty_cycle;

% other parameters
syms C_dc L_dc; % filter parameters
syms v_ocv R_bat SOC; % battery parameters
% syms kpv_dc kiv_dc kpi_dc kii_dc; % Controller parameters
syms wv_dc wi_dc; % bandwidth of controller
syms v_dc_ref; % reference dc voltage


%% AC link variables

% state variables
syms idr iqr ed eq id iq vd vq igd igq w delta;

% other variables
syms Lg Rg % Line impedance
syms Lf Cf Rf % LC filter
syms vgD vgQ wg % Grid parameters
syms wf Dw Pr W0 % Droop control 
syms vd_ref vq_ref % reference ac voltage
syms wv_ac wi_ac; % bandwidth of controller
% syms kpv_ac kiv_ac kpi_ac kii_ac % controller parameters


%% AC link equation

vgd = vgD*cos(delta) + vgQ*sin(delta); 
vgq = -vgD*sin(delta) + vgQ*cos(delta);

q = vq*igd - vd*igq;
p = vd*id + vq*iq;

kpv_ac = Cf*wv_ac;
kiv_ac = Cf*wv_ac^2/4*100;
kpi_ac = Lf*wi_ac;
kii_ac = Lf*(wi_ac^2)/4;

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
kpv_dc = C_dc*wv_dc;
kiv_dc = C_dc*wv_dc^2/4;
kpi_dc = L_dc*wi_dc;
kii_dc = L_dc*(wi_dc^2)/4;

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

Amat_syms = jacobian(f_xu,state);

%% Set numerical number
C_dc = 0.224;
L_dc = 0.000125;
v_ocv = 0.625;
R_bat = 0.0625;
SOC = 1;
wv_dc = 50 * 2 * pi;
wi_dc = 500 * 2 * pi;
v_dc_ref = 1;

Lg = 0.5/Wbase;
% Rg = 0.1;
Rg_test = linspace(0.001,0.1,20);
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;
vgD = 1;
vgQ = 0;
wg = 1 * Wbase;
wf = 2*pi*10;
Dw = 0.08*Wbase/Sbase;
Pr = 0;
W0 = Wbase;
wv_ac = 250*2*pi;
wi_ac = 1000*2*pi;
vd_ref = 1;
vq_ref = 0;

eqn_syms = [di_bat == 0, dv_dc == 0, di_bat_ref == 0, dduty_cycle == 0, didr == 0, diqr == 0, ded == 0, deq == 0, did == 0, diq == 0, dvd == 0, dvq == 0, digd == 0, digq == 0, dw == 0, ddelta == 0];

for i=1:20

    Rg=Rg_test(i);

    %% Caculate steady state of state variables
    eqn=eqn_syms;
    syms i_bat v_dc i_bat_ref duty_cycle;
    syms idr iqr ed eq id iq vd vq igd igq w delta
    eqn = subs(eqn,'C_dc',C_dc);
    eqn = subs(eqn,'L_dc',L_dc);
    eqn = subs(eqn,'v_ocv',v_ocv);
    eqn = subs(eqn,'R_bat',R_bat);
    eqn = subs(eqn,'SOC',SOC);
    eqn = subs(eqn,'wv_dc',wv_dc);
    eqn = subs(eqn,'wi_dc',wi_dc);
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
    eqn = subs(eqn,'wv_ac',wv_ac);
    eqn = subs(eqn,'wi_ac',wi_ac);
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

    %% Replace symbolic by numerical number
    Para = [C_dc,L_dc,v_ocv,R_bat,SOC,wv_dc,wi_dc,v_dc_ref,Rg,Lg,Cf,Lf,Rf,vgD,vgQ,wg,Dw,wf,wv_ac,wi_ac,vd_ref,vq_ref];
    steadyState = [i_bat,v_dc,i_bat_ref,duty_cycle,idr,iqr,ed,eq,id,iq,vd,vq,igd,igq,w,delta];
    Amat = substitution(Amat_syms,Para,steadyState);


    %% Sweep parameters
    [phi,~] = eig(Amat);
    EigVec = eig(Amat);
    EigVecHz = EigVec/(2*pi);
    box on;
    target_mode = EigVecHz([13,14]);
    target_ZoomInAxis = [-10,0,-50,50];
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.45]);
    if i == 20
        hold on; grid on;
        s=scatter(real(target_mode),imag(target_mode),'x','LineWidth',1.5); 
        s.MarkerFaceColor = 'k';
        s.MarkerEdgeColor = 'k';
    else
        hold on; grid on;
        s= scatter(real(target_mode),imag(target_mode),'x','LineWidth',1.5); 
        s.MarkerFaceColor = [0.8500 0.3250 0.0980];
        s.MarkerEdgeColor = [0.8500 0.3250 0.0980];
    end
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis(target_ZoomInAxis);

end

Rg = 0.1;
wv_ac_test=linspace(50*2*pi,250*2*pi,20);
for i=1:20

    wv_ac=wv_ac_test(i);

    %% Caculate steady state of state variables
    eqn=eqn_syms;
    syms i_bat v_dc i_bat_ref duty_cycle;
    syms idr iqr ed eq id iq vd vq igd igq w delta
    eqn = subs(eqn,'C_dc',C_dc);
    eqn = subs(eqn,'L_dc',L_dc);
    eqn = subs(eqn,'v_ocv',v_ocv);
    eqn = subs(eqn,'R_bat',R_bat);
    eqn = subs(eqn,'SOC',SOC);
    eqn = subs(eqn,'wv_dc',wv_dc);
    eqn = subs(eqn,'wi_dc',wi_dc);
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
    eqn = subs(eqn,'wv_ac',wv_ac);
    eqn = subs(eqn,'wi_ac',wi_ac);
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

    %% Replace symbolic by numerical number
    Para = [C_dc,L_dc,v_ocv,R_bat,SOC,wv_dc,wi_dc,v_dc_ref,Rg,Lg,Cf,Lf,Rf,vgD,vgQ,wg,Dw,wf,wv_ac,wi_ac,vd_ref,vq_ref];
    steadyState = [i_bat,v_dc,i_bat_ref,duty_cycle,idr,iqr,ed,eq,id,iq,vd,vq,igd,igq,w,delta];
    Amat = substitution(Amat_syms,Para,steadyState);


    %% Sweep parameters
    [phi,~] = eig(Amat);
    EigVec = eig(Amat);
    EigVecHz = EigVec/(2*pi);
    box on;
    target_mode = EigVecHz([13,14]);
    target_ZoomInAxis = [-10,0,-50,50];
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.45]);
    if i == 1
        % hold on; grid on;
        % s=scatter(real(target_mode),imag(target_mode),'x','LineWidth',1.5); 
        % s.MarkerFaceColor = 'k';
        % s.MarkerEdgeColor = 'k';
    else
        hold on; grid on;
        s= scatter(real(target_mode),imag(target_mode),'x','LineWidth',1.5); 
        s.MarkerFaceColor = 'b';
        s.MarkerEdgeColor = 'b';
    end
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis(target_ZoomInAxis);

end





function Amat = substitution(Amat_syms,Para,steadyState)
    Amat = Amat_syms;

    Amat = subs(Amat,'C_dc',Para(1));
    Amat = subs(Amat,'L_dc',Para(2));

    Amat = subs(Amat,'v_ocv',Para(3));
    Amat = subs(Amat,'R_bat',Para(4));
    Amat = subs(Amat,'SOC',Para(5));

    Amat = subs(Amat,'wv_dc',Para(6));
    Amat = subs(Amat,'wi_dc',Para(7));
    Amat = subs(Amat,'v_dc_ref',Para(8));


    Amat = subs(Amat,'Rg',Para(9));
    Amat = subs(Amat,'Lg',Para(10));
    Amat = subs(Amat,'Cf',Para(11));
    Amat = subs(Amat,'Lf',Para(12));
    Amat = subs(Amat,'Rf',Para(13));

    Amat = subs(Amat,'vgD',Para(14));
    Amat = subs(Amat,'vgQ',Para(15));
    Amat = subs(Amat,'wg',Para(16));

    Amat = subs(Amat,'Dw',Para(17));
    Amat = subs(Amat,'wf',Para(18));

    Amat = subs(Amat,'wv_ac',Para(19));
    Amat = subs(Amat,'wi_ac',Para(20));

    Amat = subs(Amat,'vd_ref',Para(21));
    Amat = subs(Amat,'vq_ref',Para(22));

    Amat = subs(Amat,'i_bat',steadyState(1));
    Amat = subs(Amat,'v_dc',steadyState(2));
    Amat = subs(Amat,'i_bat_ref',steadyState(3));
    Amat = subs(Amat,'duty_cycle',steadyState(4));
    Amat = subs(Amat,'idr',steadyState(5));
    Amat = subs(Amat,'iqr',steadyState(6));
    Amat = subs(Amat,'ed',steadyState(7));
    Amat = subs(Amat,'eq',steadyState(8));
    Amat = subs(Amat,'id',steadyState(9));
    Amat = subs(Amat,'iq',steadyState(10));
    Amat = subs(Amat,'vd',steadyState(11));
    Amat = subs(Amat,'vq',steadyState(12));
    Amat = subs(Amat,'igd',steadyState(13));
    Amat = subs(Amat,'igq',steadyState(14));
    Amat = subs(Amat,'w',steadyState(15));
    Amat = subs(Amat,'delta',steadyState(16));

    Amat = double(Amat);
end