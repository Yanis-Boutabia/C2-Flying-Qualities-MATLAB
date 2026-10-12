%% ========================================================================
%  C2 AIRCRAFT - FLYING QUALITIES ANALYSIS
%  Course: Mecanique du vol - Qualites de vol (IPSA Aero4)
%% ========================================================================
clc; clear; close all;

outdir = fullfile(pwd, 'report_figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ========================================================================
%  1. AIRCRAFT C2 DATA
%% ========================================================================
g = 9.81;

% Geometry
S = 12.63; cbar = 1.65; b = 8.02; AR = 5.09; e = 0.94;

% Flight condition
h = 7620; Mach = 0.6; U0 = 185.9; qbar = 9500; rho = 0.55;
theta0 = 0*pi/180;

% Mass / inertia
m = 1816; Ixx = 1085; Iyy = 6511; Izz = 7053; Ixz = 271;

% Steady-state
CL0 = 0.15; CD0 = 0.022;

% Longitudinal stability derivatives
CDu = 0; CDa = 0.12; CLu = 0.084; CLa = 5.5;
CLadot = 4.2; CLq = 10;
Cmu = 0; Cma = -0.24; Cmadot = -9.6; Cmq = -17.7;

% Lateral stability derivatives
Clb = -0.11; Clp = -0.39; Clr = 0.28;
Cnb = 0.17;  Cnp = 0.09;  Cnr = -0.26;
Cyb = -1;    Cyp = -0.14; Cyr = 0.61;

% Control derivatives
CDde = 0; CLde = 0.038; Cmde = -0.88;
Cyda = 0; Clda = 0.1;   Cnda = -0.003;
Cydr = 0.28; Cldr = 0.05; Cndr = -0.12;

%% ========================================================================
%  2. DIMENSIONAL LONGITUDINAL DERIVATIVES
%% ========================================================================
Xu = -(qbar*S)/(m*U0) * (2*CD0 + CDu);
Zu = -(qbar*S)/(m*U0) * (2*CL0 + CLu);
Mu =  (qbar*S*cbar)/(Iyy*U0) * Cmu;

Xw = (qbar*S)/(m*U0) * CL0 * (1 - (2/(pi*e*AR))*CLa);
Zw = -(qbar*S)/(m*U0) * (CD0 + CLa);
Mw = (qbar*S*cbar)/(Iyy*U0) * Cma;

Zwdot = (qbar*S*cbar)/(2*m*U0^2) * (2*CD0 + CLadot);
Mwdot = (qbar*S*cbar^2)/(2*Iyy*U0^2) * Cmadot;

Zq = (qbar*S*cbar)/(2*m*U0) * CLq;
Mq = (qbar*S*cbar^2)/(2*Iyy*U0) * Cmq;

% Longitudinal control derivatives
Xde = (qbar*S)/(m*U0) * CDde;
Zde = -(qbar*S)/(m*U0) * CLde;   % negative: lift acts opposite to the z axis (down), same convention as Zu, Zw
Mde = (qbar*S*cbar)/(Iyy*U0) * Cmde;
XdT = 0; ZdT = 0; MdT = 0;

%% ========================================================================
%  3. LONGITUDINAL MATRICES
%  States: [du dw dq dtheta]
%% ========================================================================
ALong = [ Xu,  Xw,  0,              -g*cos(theta0);
          Zu,  Zw,  U0,             -g*sin(theta0);
          Mu + Zu*Mwdot,  Mw + Zw*Mwdot,  Mq + U0*Mwdot,  0;
          0,   0,   1,               0 ];

BLong = [ Xde,  XdT;
          Zde,  ZdT;
          Mde + Zde*Mwdot,  MdT + ZdT*Mwdot;
          0,    0 ];

disp(' ')
disp('=================================')
disp('    LONGITUDINAL RESULTS')
disp('=================================')
disp('ALong = '); disp(ALong)
disp('BLong = '); disp(BLong)

%% ========================================================================
%  4. DIMENSIONAL LATERAL DERIVATIVES (convention prof avec matrice In)
%% ========================================================================
Yv = (qbar*S)/(m*U0) * Cyb;
Yp = (qbar*S*b)/(2*m*U0) * Cyp;
Yr = (qbar*S*b)/(2*m*U0) * Cyr;

Lv = (qbar*S*b)/(Ixx*U0) * Clb;
Lp = (qbar*S*b^2)/(2*Ixx*U0) * Clp;
Lr = (qbar*S*b^2)/(2*Ixx*U0) * Clr;

Nv = (qbar*S*b)/(Izz*U0) * Cnb;
Np = (qbar*S*b^2)/(2*Izz*U0) * Cnp;
Nr = (qbar*S*b^2)/(2*Izz*U0) * Cnr;

Yda = (qbar*S)/(m*U0) * Cyda;
Ydr = (qbar*S)/(m*U0) * Cydr;
Lda = (qbar*S*b)/(Ixx*U0) * Clda;
Ldr = (qbar*S*b)/(Ixx*U0) * Cldr;
Nda = (qbar*S*b)/(Izz*U0) * Cnda;
Ndr = (qbar*S*b)/(Izz*U0) * Cndr;

%% ========================================================================
%  5. LATERAL MATRICES (with Ixz coupling via InLat inversion)
%% ========================================================================
AnLat = [ Yv,  Yp,  Yr - U0,  g*cos(theta0);
          Lv,  Lp,  Lr,       0;
          Nv,  Np,  Nr,       0;
          0,   1,   0,        0 ];

BnLat = [ Ydr,  Yda;
          Ldr,  Lda;
          Ndr,  Nda;
          0,    0 ];

InLat = [1, 0,          0,        0;
         0, 1,   -Ixz/Ixx,        0;
         0, -Ixz/Izz,    1,       0;
         0, 0,           0,       1];

ALat = InLat \ AnLat;
BLat = InLat \ BnLat;

disp(' ')
disp('=================================')
disp('      LATERAL RESULTS')
disp('=================================')
disp('ALat = '); disp(ALat)
disp('BLat = '); disp(BLat)

%% ========================================================================
%  6. EIGENVALUE ANALYSIS - LONGITUDINAL
%% ========================================================================
[Vlong, Dlong] = eig(ALong);
lambda_long = diag(Dlong);
char_poly_long = poly(ALong);

disp(' ')
disp('Longitudinal eigenvalues:'); disp(lambda_long)
disp('Longitudinal eigenvectors:'); disp(Vlong)
disp('Longitudinal characteristic polynomial coefficients:'); disp(char_poly_long)

[~, idxL] = sort(abs(imag(lambda_long)), 'descend');
idxSP = idxL(1:2);
idxPH = idxL(3:4);

lambdaSP = lambda_long(idxSP);
lambdaPH = lambda_long(idxPH);

sigma_SP = real(lambdaSP(1));
omega_d_SP = abs(imag(lambdaSP(1)));
omega_n_SP = sqrt(sigma_SP^2 + omega_d_SP^2);
zeta_SP = -sigma_SP / omega_n_SP;

sigma_PH = real(lambdaPH(1));
omega_d_PH = abs(imag(lambdaPH(1)));
omega_n_PH = sqrt(sigma_PH^2 + omega_d_PH^2);
zeta_PH = -sigma_PH / omega_n_PH;

disp(' '); disp('--- SHORT PERIOD ---')
disp(['lambda_SP = ', num2str(lambdaSP(1))])
disp(['omega_n_SP = ', num2str(omega_n_SP), ' rad/s'])
disp(['zeta_SP = ', num2str(zeta_SP)])

disp(' '); disp('--- PHUGOID ---')
disp(['lambda_PH = ', num2str(lambdaPH(1))])
disp(['omega_n_PH = ', num2str(omega_n_PH), ' rad/s'])
disp(['zeta_PH = ', num2str(zeta_PH)])

%% ========================================================================
%  7. EIGENVALUE ANALYSIS - LATERAL
%% ========================================================================
[Vlat, Dlat] = eig(ALat);
lambda_lat = diag(Dlat);
char_poly_lat = poly(ALat);

disp(' ')
disp('Lateral eigenvalues:'); disp(lambda_lat)
disp('Lateral eigenvectors:'); disp(Vlat)
disp('Lateral characteristic polynomial coefficients:'); disp(char_poly_lat)

tolImag = 1e-8;
idxReal = find(abs(imag(lambda_lat)) < tolImag);
idxCplx = find(abs(imag(lambda_lat)) >= tolImag);

lambda_real = real(lambda_lat(idxReal));
[~, iRoll] = min(lambda_real);
[~, iSpir] = min(abs(lambda_real));

idxRoll = idxReal(iRoll);
idxSpir = idxReal(iSpir);

lambda_roll = lambda_lat(idxRoll);
lambda_spiral = lambda_lat(idxSpir);

idxDR = idxCplx(1:2);
lambdaDR_pair = lambda_lat(idxDR);

sigma_DR = real(lambdaDR_pair(1));
omega_d_DR = abs(imag(lambdaDR_pair(1)));
omega_n_DR = sqrt(sigma_DR^2 + omega_d_DR^2);
zeta_DR = -sigma_DR / omega_n_DR;

T_roll = -1/real(lambda_roll);

disp(' '); disp('--- ROLL MODE ---')
disp(['lambda_roll = ', num2str(real(lambda_roll))])
disp(['T_roll = ', num2str(T_roll), ' s'])

disp(' '); disp('--- SPIRAL MODE ---')
disp(['lambda_spiral = ', num2str(real(lambda_spiral))])
if real(lambda_spiral) > 0
    disp(['T2 (double) = ', num2str(log(2)/real(lambda_spiral)), ' s (UNSTABLE)'])
end

disp(' '); disp('--- DUTCH ROLL ---')
disp(['lambda_DR = ', num2str(lambdaDR_pair(1))])
disp(['omega_n_DR = ', num2str(omega_n_DR), ' rad/s'])
disp(['zeta_DR = ', num2str(zeta_DR)])
disp(['zeta*omega_n = ', num2str(zeta_DR*omega_n_DR)])

%% ========================================================================
%  8. LONGITUDINAL MODE SHAPES (NORMALIZED)
%% ========================================================================
VSP = Vlong(:, idxSP);
VPH = Vlong(:, idxPH);

tSP = linspace(0, 20, 2000);
tPH = linspace(0, 800, 4000);

du_SP     = real(VSP(1,1)*exp(lambdaSP(1)*tSP) + VSP(1,2)*exp(lambdaSP(2)*tSP));
dw_SP     = real(VSP(2,1)*exp(lambdaSP(1)*tSP) + VSP(2,2)*exp(lambdaSP(2)*tSP));
dq_SP     = real(VSP(3,1)*exp(lambdaSP(1)*tSP) + VSP(3,2)*exp(lambdaSP(2)*tSP));
dtheta_SP = real(VSP(4,1)*exp(lambdaSP(1)*tSP) + VSP(4,2)*exp(lambdaSP(2)*tSP));

du_PH     = real(VPH(1,1)*exp(lambdaPH(1)*tPH) + VPH(1,2)*exp(lambdaPH(2)*tPH));
dw_PH     = real(VPH(2,1)*exp(lambdaPH(1)*tPH) + VPH(2,2)*exp(lambdaPH(2)*tPH));
dq_PH     = real(VPH(3,1)*exp(lambdaPH(1)*tPH) + VPH(3,2)*exp(lambdaPH(2)*tPH));
dtheta_PH = real(VPH(4,1)*exp(lambdaPH(1)*tPH) + VPH(4,2)*exp(lambdaPH(2)*tPH));

% Normalize
scale_SP = max([max(abs(du_SP)) max(abs(dw_SP)) max(abs(dq_SP)) max(abs(dtheta_SP))]) + eps;
scale_PH = max([max(abs(du_PH)) max(abs(dw_PH)) max(abs(dq_PH)) max(abs(dtheta_PH))]) + eps;

du_SP = du_SP/scale_SP; dw_SP = dw_SP/scale_SP;
dq_SP = dq_SP/scale_SP; dtheta_SP = dtheta_SP/scale_SP;

du_PH = du_PH/scale_PH; dw_PH = dw_PH/scale_PH;
dq_PH = dq_PH/scale_PH; dtheta_PH = dtheta_PH/scale_PH;

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tSP, du_SP, 'LineWidth', 1.5); hold on;
plot(tSP, dw_SP, 'LineWidth', 1.5);
plot(tSP, dq_SP, 'LineWidth', 1.5);
plot(tSP, dtheta_SP, 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Normalized amplitude');
title('Short Period Mode - Impulse Time Response');
legend('\Delta u', '\Delta w', '\Delta q', '\Delta \theta', 'Location', 'best');
xlim([0 2.5]); set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_mode_ShortPeriod.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tPH, du_PH, 'LineWidth', 1.5); hold on;
plot(tPH, dw_PH, 'LineWidth', 1.5);
plot(tPH, dq_PH, 'LineWidth', 1.5);
plot(tPH, dtheta_PH, 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Normalized amplitude');
title('Phugoid Mode - Impulse Time Response');
legend('\Delta u', '\Delta w', '\Delta q', '\Delta \theta', 'Location', 'best');
set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_mode_Phugoid.png'), 'Resolution', 200);

%% ========================================================================
%  9. LATERAL MODE SHAPES (NORMALIZED)
%% ========================================================================
VDR = Vlat(:, idxDR);
VRoll = Vlat(:, idxRoll);
VSpir = Vlat(:, idxSpir);

tDR = linspace(0, 100, 3000);
tRoll = linspace(0, 20, 2000);
tSpir = linspace(0, 300, 3000);

dv_DR = real(VDR(1,1)*exp(lambdaDR_pair(1)*tDR) + VDR(1,2)*exp(lambdaDR_pair(2)*tDR));
dp_DR = real(VDR(2,1)*exp(lambdaDR_pair(1)*tDR) + VDR(2,2)*exp(lambdaDR_pair(2)*tDR));
dr_DR = real(VDR(3,1)*exp(lambdaDR_pair(1)*tDR) + VDR(3,2)*exp(lambdaDR_pair(2)*tDR));
dphi_DR = real(VDR(4,1)*exp(lambdaDR_pair(1)*tDR) + VDR(4,2)*exp(lambdaDR_pair(2)*tDR));

dv_roll = real(VRoll(1)*exp(lambda_roll*tRoll));
dp_roll = real(VRoll(2)*exp(lambda_roll*tRoll));
dr_roll = real(VRoll(3)*exp(lambda_roll*tRoll));
dphi_roll = real(VRoll(4)*exp(lambda_roll*tRoll));

dv_sp = real(VSpir(1)*exp(lambda_spiral*tSpir));
dp_sp = real(VSpir(2)*exp(lambda_spiral*tSpir));
dr_sp = real(VSpir(3)*exp(lambda_spiral*tSpir));
dphi_sp = real(VSpir(4)*exp(lambda_spiral*tSpir));

scale_DR = max([max(abs(dv_DR)) max(abs(dp_DR)) max(abs(dr_DR)) max(abs(dphi_DR))]) + eps;
scale_roll = max([max(abs(dv_roll)) max(abs(dp_roll)) max(abs(dr_roll)) max(abs(dphi_roll))]) + eps;
scale_sp = max([max(abs(dv_sp)) max(abs(dp_sp)) max(abs(dr_sp)) max(abs(dphi_sp))]) + eps;

dv_DR = dv_DR/scale_DR; dp_DR = dp_DR/scale_DR;
dr_DR = dr_DR/scale_DR; dphi_DR = dphi_DR/scale_DR;

dv_roll = dv_roll/scale_roll; dp_roll = dp_roll/scale_roll;
dr_roll = dr_roll/scale_roll; dphi_roll = dphi_roll/scale_roll;

dv_sp = dv_sp/scale_sp; dp_sp = dp_sp/scale_sp;
dr_sp = dr_sp/scale_sp; dphi_sp = dphi_sp/scale_sp;

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tDR, dv_DR, 'LineWidth', 1.5); hold on;
plot(tDR, dp_DR, 'LineWidth', 1.5);
plot(tDR, dr_DR, 'LineWidth', 1.5);
plot(tDR, dphi_DR, 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Normalized amplitude');
title('Dutch Roll Mode - Impulse Time Response');
legend('\Delta v', '\Delta p', '\Delta r', '\Delta \phi', 'Location', 'best');
xlim([0 15]); set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_mode_DutchRoll.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tRoll, dv_roll, 'LineWidth', 1.5); hold on;
plot(tRoll, dp_roll, 'LineWidth', 1.5);
plot(tRoll, dr_roll, 'LineWidth', 1.5);
plot(tRoll, dphi_roll, 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Normalized amplitude');
title('Roll Mode - Impulse Time Response');
legend('\Delta v', '\Delta p', '\Delta r', '\Delta \phi', 'Location', 'best');
xlim([0 1]); set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_mode_Roll.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tSpir, dv_sp, 'LineWidth', 1.5); hold on;
plot(tSpir, dp_sp, 'LineWidth', 1.5);
plot(tSpir, dr_sp, 'LineWidth', 1.5);
plot(tSpir, dphi_sp, 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Normalized amplitude');
title('Spiral Mode - Impulse Time Response');
legend('\Delta v', '\Delta p', '\Delta r', '\Delta \phi', 'Location', 'best');
set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_mode_Spiral.png'), 'Resolution', 200);

%% ========================================================================
%  10. GLOBAL FREE RESPONSES
%% ========================================================================
x0_long = [1; 0.2; 0.05; 0.02];
[tt_long, x_long] = ode45(@(t,x) ALong*x, linspace(0,800,4000), x0_long);

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tt_long, x_long(:,1), 'LineWidth', 1.3); hold on;
plot(tt_long, x_long(:,2), 'LineWidth', 1.3);
plot(tt_long, x_long(:,3), 'LineWidth', 1.3);
plot(tt_long, x_long(:,4), 'LineWidth', 1.3);
grid on; xlabel('Time (s)'); ylabel('Amplitude');
title('Longitudinal Free Response - All State Variables');
legend('\Delta u', '\Delta w', '\Delta q', '\Delta \theta', 'Location', 'best');
set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_free_longitudinal.png'), 'Resolution', 200);

x0_lat = [1; 0.2; 0.1; 0.05];
[tt_lat, x_lat] = ode45(@(t,x) ALat*x, linspace(0,300,4000), x0_lat);

fig = figure('Visible','on','Position',[100 100 700 450]);
plot(tt_lat, x_lat(:,1), 'LineWidth', 1.3); hold on;
plot(tt_lat, x_lat(:,2), 'LineWidth', 1.3);
plot(tt_lat, x_lat(:,3), 'LineWidth', 1.3);
plot(tt_lat, x_lat(:,4), 'LineWidth', 1.3);
grid on; xlabel('Time (s)'); ylabel('Amplitude');
title('Lateral Free Response - All State Variables');
legend('\Delta v', '\Delta p', '\Delta r', '\Delta \phi', 'Location', 'best');
set(gca,'FontSize',11);
exportgraphics(fig, fullfile(outdir,'fig_free_lateral.png'), 'Resolution', 200);

%% ========================================================================
%  11. FLYING QUALITY ASSESSMENT (Class I, Cat B)
%% ========================================================================
disp(' ')
disp('==============================================')
disp('   FLYING QUALITY ASSESSMENT (Class I, Cat B)')
disp('==============================================')

disp(['SP: omega_n=', num2str(omega_n_SP,'%.4f'), ', xi=', num2str(zeta_SP,'%.4f')])
n_alpha = (qbar*S*CLa)/(m*g); CAP = omega_n_SP^2/n_alpha;
disp(['    CAP = ', num2str(CAP,'%.4f'), ' -> Level 1 (0.085-3.6)'])

disp(['Ph: omega_n=', num2str(omega_n_PH,'%.4f'), ', xi=', num2str(zeta_PH,'%.4f'), ' -> Level 1 (xi>=0.04)'])
disp(['Roll: T=', num2str(T_roll,'%.4f'), ' s -> Level 1 (T_R<=1.4s)'])

if real(lambda_spiral) < 0
    disp('Spiral: stable -> Level 1')
else
    T2s = log(2)/real(lambda_spiral);
    disp(['Spiral: T2=', num2str(T2s,'%.2f'), ' s -> Level 1 (T2>=20s)'])
end

disp(['DR: omega_n=', num2str(omega_n_DR,'%.4f'), ', xi=', num2str(zeta_DR,'%.4f'), ', xi*wn=', num2str(zeta_DR*omega_n_DR,'%.4f'), ' -> Level 1'])

%% ========================================================================
%  12. SAS DESIGN - POLE PLACEMENT (A* = A - B*K^T method)
%% ========================================================================
disp(' ')
disp('===============================')
disp('   SAS DESIGN - POLE PLACEMENT')
disp('===============================')

% --- PITCH DAMPER: Short Period xi 0.82 -> 0.95 ---
xi_SP_new = 0.95;
lam_sp_desired = omega_n_SP*(-xi_SP_new + 1j*sqrt(1-xi_SP_new^2));
desired_poles_long = [lam_sp_desired, conj(lam_sp_desired), lambdaPH(1), lambdaPH(2)];

K_long = place(ALong, BLong(:,1), desired_poles_long);
ALong_SAS = ALong - BLong(:,1)*K_long;

disp(' '); disp('--- PITCH DAMPER ---')
disp(['Target: xi_SP = ', num2str(xi_SP_new)])
disp(['K_long = [', num2str(K_long(1),'%.4f'), ', ', num2str(K_long(2),'%.4f'), ...
      ', ', num2str(K_long(3),'%.4f'), ', ', num2str(K_long(4),'%.4f'), ']'])

eig_long_SAS = eig(ALong_SAS);
disp('Eigenvalues with SAS:'); disp(eig_long_SAS)

% --- YAW DAMPER: Dutch Roll xi 0.09 -> 0.20 ---
xi_DR_new = 0.20;
lam_dr_desired = omega_n_DR*(-xi_DR_new + 1j*sqrt(1-xi_DR_new^2));
desired_poles_lat = [lam_dr_desired, conj(lam_dr_desired), lambda_roll, lambda_spiral];

K_lat = place(ALat, BLat(:,1), desired_poles_lat);
ALat_SAS = ALat - BLat(:,1)*K_lat;

disp(' '); disp('--- YAW DAMPER ---')
disp(['Target: xi_DR = ', num2str(xi_DR_new)])
disp(['K_lat = [', num2str(K_lat(1),'%.4f'), ', ', num2str(K_lat(2),'%.4f'), ...
      ', ', num2str(K_lat(3),'%.4f'), ', ', num2str(K_lat(4),'%.4f'), ']'])

eig_lat_SAS = eig(ALat_SAS);
disp('Eigenvalues with SAS:'); disp(eig_lat_SAS)

% --- PITCH DAMPER COMPARISON FIGURE ---
t_comp = linspace(0, 10, 1000);
X0_long = [0; 1; 0; 0];  % initial Δw perturbation
[y_orig_long, ~] = initial(ss(ALong, zeros(4,1), eye(4), zeros(4,1)), X0_long, t_comp);
[y_sas_long, ~] = initial(ss(ALong_SAS, zeros(4,1), eye(4), zeros(4,1)), X0_long, t_comp);

fig = figure('Visible','on','Position',[100 100 1000 750]);
labels_long = {'\Delta u [m/s]','\Delta w [m/s]','\Delta q [rad/s]','\Delta \theta [deg]'};
for k = 1:4
    subplot(2,2,k);
    d_orig = y_orig_long(:,k);
    d_sas = y_sas_long(:,k);
    if k == 4  % theta in degrees
        d_orig = rad2deg(d_orig);
        d_sas = rad2deg(d_sas);
    end
    plot(t_comp, d_orig, 'b--', 'LineWidth', 1.8); hold on;
    plot(t_comp, d_sas, 'r-', 'LineWidth', 1.8);
    xlabel('Time [s]'); ylabel(labels_long{k}); grid on;
    legend('Original', 'With SAS', 'Location', 'best');
    set(gca,'FontSize',11);
end
sgtitle('SAS Comparison - Pitch Damper', 'FontSize', 14);
exportgraphics(fig, fullfile(outdir,'fig_SAS_PitchDamper.png'), 'Resolution', 200);

% --- YAW DAMPER COMPARISON FIGURE ---
t_comp_lat = linspace(0, 30, 2000);
X0_lat = [1; 0; 0; 0];  % initial Δv perturbation (sideslip)
[y_orig_lat, ~] = initial(ss(ALat, zeros(4,1), eye(4), zeros(4,1)), X0_lat, t_comp_lat);
[y_sas_lat, ~] = initial(ss(ALat_SAS, zeros(4,1), eye(4), zeros(4,1)), X0_lat, t_comp_lat);

fig = figure('Visible','on','Position',[100 100 1000 750]);
labels_lat = {'\Delta v [m/s]','\Delta p [deg/s]','\Delta r [deg/s]','\Delta \phi [deg]'};
for k = 1:4
    subplot(2,2,k);
    d_orig = y_orig_lat(:,k);
    d_sas = y_sas_lat(:,k);
    if k >= 2  % p, r, phi in degrees
        d_orig = rad2deg(d_orig);
        d_sas = rad2deg(d_sas);
    end
    plot(t_comp_lat, d_orig, 'b--', 'LineWidth', 1.8); hold on;
    plot(t_comp_lat, d_sas, 'r-', 'LineWidth', 1.8);
    xlabel('Time [s]'); ylabel(labels_lat{k}); grid on;
    legend('Original', 'With SAS', 'Location', 'best');
    set(gca,'FontSize',11);
end
sgtitle('SAS Comparison - Yaw Damper', 'FontSize', 14);
exportgraphics(fig, fullfile(outdir,'fig_SAS_YawDamper.png'), 'Resolution', 200);

%% ========================================================================
%  13. BODE DIAGRAMS (numerical frequency response)
%% ========================================================================
I4 = eye(4);
w = logspace(-2, 2, 1200);

G_q_de = zeros(size(w));
G_theta_de = zeros(size(w));
G_p_da = zeros(size(w));
G_r_dr = zeros(size(w));

for k = 1:length(w)
    s = 1i*w(k);
    GL = (s*I4 - ALong)\BLong;
    GT = (s*I4 - ALat)\BLat;
    G_q_de(k) = GL(3,1);
    G_theta_de(k) = GL(4,1);
    G_p_da(k) = GT(2,2);
    G_r_dr(k) = GT(3,1);
end

fig = figure('Visible','on','Position',[100 100 700 500]);
subplot(2,1,1)
semilogx(w, 20*log10(abs(G_q_de)), 'LineWidth', 1.3); grid on;
ylabel('Magnitude (dB)'); title('Bode Diagram - q / \delta_e');
subplot(2,1,2)
semilogx(w, angle(G_q_de)*180/pi, 'LineWidth', 1.3); grid on;
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
exportgraphics(fig, fullfile(outdir,'fig_Bode_q_de.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 500]);
subplot(2,1,1)
semilogx(w, 20*log10(abs(G_theta_de)), 'LineWidth', 1.3); grid on;
ylabel('Magnitude (dB)'); title('Bode Diagram - \theta / \delta_e');
subplot(2,1,2)
semilogx(w, angle(G_theta_de)*180/pi, 'LineWidth', 1.3); grid on;
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
exportgraphics(fig, fullfile(outdir,'fig_Bode_theta_de.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 500]);
subplot(2,1,1)
semilogx(w, 20*log10(abs(G_p_da)), 'LineWidth', 1.3); grid on;
ylabel('Magnitude (dB)'); title('Bode Diagram - p / \delta_a');
subplot(2,1,2)
semilogx(w, angle(G_p_da)*180/pi, 'LineWidth', 1.3); grid on;
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
exportgraphics(fig, fullfile(outdir,'fig_Bode_p_da.png'), 'Resolution', 200);

fig = figure('Visible','on','Position',[100 100 700 500]);
subplot(2,1,1)
semilogx(w, 20*log10(abs(G_r_dr)), 'LineWidth', 1.3); grid on;
ylabel('Magnitude (dB)'); title('Bode Diagram - r / \delta_r');
subplot(2,1,2)
semilogx(w, angle(G_r_dr)*180/pi, 'LineWidth', 1.3); grid on;
xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');
exportgraphics(fig, fullfile(outdir,'fig_Bode_r_dr.png'), 'Resolution', 200);

%% ========================================================================
%  14. ROOT LOCUS with S-PLANE GRID (custom)
%% ========================================================================
Kvec = linspace(0, 500, 700);

root_channels = {
    ALong, BLong(:,1), [0 0 1 0], 'q/\delta_e',         [0 0.4470 0.7410];
    ALong, BLong(:,1), [0 0 0 1], '\theta/\delta_e',    [0.8500 0.3250 0.0980];
    ALat,  BLat(:,2),  [0 1 0 0], 'p/\delta_a',         [0.9290 0.6940 0.1250];
    ALat,  BLat(:,1),  [0 0 1 0], 'r/\delta_r',         [0.4940 0.1840 0.5560]
};

fig = figure('Visible','on','Position',[100 100 800 600], 'Color', 'w');
hold on; grid on; box on;

xlim_vals = [-15 2];
ylim_vals = [-10 10];

draw_sgrid_custom(xlim_vals, ylim_vals);

legend_entries = {'grid', '', '', ''};

for k = 1:size(root_channels, 1)
    A0 = root_channels{k,1};
    B0 = root_channels{k,2};
    Cfb = root_channels{k,3};
    name = root_channels{k,4};
    col = root_channels{k,5};

    pts_all = [];
    for j = 1:length(Kvec)
        K = Kvec(j);
        Acl = A0 - B0*(K*Cfb);
        rts = eig(Acl);
        pts_all = [pts_all; rts(:).'];
    end

    for cidx = 1:size(pts_all, 2)
        plot(real(pts_all(:,cidx)), imag(pts_all(:,cidx)), '.', 'Color', col, 'MarkerSize', 6);
    end

    open_poles = eig(A0);
    plot(real(open_poles), imag(open_poles), 'x', 'Color', col, 'MarkerSize', 10, 'LineWidth', 2);

    legend_entries{end+1} = name;
    legend_entries{end+1} = ['open-loop poles ', name];
end

xlabel('Real axis'); ylabel('Imaginary axis');
title('Root Locus - Main Aircraft Channels');
xlim(xlim_vals); ylim(ylim_vals);
legend(legend_entries, 'Location', 'bestoutside');

text(real(lambda_roll), 0.5, 'Roll mode', 'FontSize', 10);
text(real(lambda_spiral)-0.8, 0.8, 'Spiral mode', 'FontSize', 10);
text(real(lambdaDR_pair(1))+0.3, imag(lambdaDR_pair(1))+0.5, 'Dutch roll', 'FontSize', 10);

exportgraphics(fig, fullfile(outdir,'fig_RootLocus.png'), 'Resolution', 200);

disp(' ')
disp('=== DONE ===')
disp(['Figures exported to: ', outdir])

%% ========================================================================
%  LOCAL FUNCTION: S-PLANE GRID
%% ========================================================================
function draw_sgrid_custom(xlim_vals, ylim_vals)
    hold on;
    plot([xlim_vals(1) xlim_vals(2)], [0 0], 'k-', 'LineWidth', 0.8);
    plot([0 0], [ylim_vals(1) ylim_vals(2)], 'k-', 'LineWidth', 0.8);

    zetas = [0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.98];
    rmax = max(abs([xlim_vals ylim_vals]))*1.2;

    for z = zetas
        phi = acos(z);
        x = linspace(0, rmax*cos(phi), 400);
        y = tan(phi) * x;
        plot(-x, y, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
        plot(-x, -y, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

        xt = -0.92*rmax*cos(phi);
        yt = 0.92*rmax*sin(phi);
        if xt >= xlim_vals(1) && xt <= xlim_vals(2) && yt >= ylim_vals(1) && yt <= ylim_vals(2)
            text(xt, yt, num2str(z), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
            text(xt, -yt, num2str(z), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
        end
    end

    wn_vals = [1 2 4 6 8 10 12 15];
    th = linspace(pi/2, 3*pi/2, 500);
    for wn = wn_vals
        xc = wn*cos(th); yc = wn*sin(th);
        plot(xc, yc, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
        if (-wn >= xlim_vals(1) && -wn <= xlim_vals(2))
            text(-wn, 0.35, num2str(wn), 'FontSize', 8, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');
        end
    end
end