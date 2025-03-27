% MATLAB source codes of simulations in paper
% "A Fully Actuated System Approach to Interval Observer Design with
% Applications in Fault Detection”
% The 3rd Conference on Fully Actuated System Theory and Applications (FASTA2024)
% -------------------------------------------------------------------------
% Additional Function Needed: gse.m
% Additional Toolbox Needed:  None
% Additional Solver Needed:   None
% -------------------------------------------------------------------------
% Version:              1.0
% Author:               Weijie Ren
% Contact:              weijie.ren@outlook.com
% Initial modified:     Jan. 25, 2024
% Last modified:        
% -------------------------------------------------------------------------
% Function Description:
%       "gse.m": solve generalized Sylvester matrix equation (SVD approach)

clc,clear
close all

C = eye(3);

% Compute controller parametric matrices
s = [-1 -2 -3];         % expected pole
F = diag(s);            % A+BK~F
Z = [-1 2 -1];          % free matrix
V = [Z; Z*F; Z*F^2];    % transformation matrix
A = -Z * F^3 / V;       % parametric matrix in the high-order linear system
A0 = A(1);
A1 = A(2);
A2 = A(3);


% Compute observer 1's gain
Phi = [  0   1   0;
         0   0   1;
       -A0 -A1 -A2];    % construct state-space form
Gamma = [0; 0; 1];
s1 = [-15 -16 -17];     % expected poles of interval observer
F1 = diag(s1);
Z1 = [-1  0  1;
       0 -1  1;
       3  0 -1];
[L1,V1,~,~,~] = gse(Phi',C',F1,Z1);
L1 = -L1';  V1 = V1';
Phi - L1 * C        % Check if observer 1's error state matrix is Meztler
% eig(Phi - L1 * C)
% P = V1 * (Phi - L * C) / V1
% Phi - L1 * C


% Compute observer 2's gain
s2 = s1;
F2 = F1;
Z2 = [1 0 0;
      0 1 0;
      0 0 1];
[L2,V2,~,~,~] = gse(Phi',C',F2,Z2);
L2 = -L2';  V2 = V2';
% eig(Phi - L * C)
% P = V2 * (Phi - L * C) / V2
% Phi - L * C
tPhi = V2 * Phi / V2;
tL = V2 * L2;
tGamma = V2 * Gamma;


% Initial value
x0(1:3,1) = 1 * (2*rand(3,1)-1);  % initial system state
x0(4:6,1) = 2 * ones(3,1);        % initial upper estimation of observer 1
x0(7:9,1) = -x0(4:6);             % initial lower estimation of observer 1
V2p = max(0,V2);  V2m = V2p - V2;
x0(10:12,1) = V2p * x0(4:6,1) - V2m * x0(7:9,1);    % initial upper estimation of observer 2
x0(13:15,1) = V2p * x0(7:9,1) - V2m * x0(4:6,1);    % initial lower estimation of observer 2


% Solve ODE
tspan = [0 10];
[t,x_all] = ode45(@(t,x) ode_fcn(t,x,A,C,Phi,L1,Gamma,tPhi,tL,tGamma,V2),tspan,x0);
x_all = x_all';
x = x_all(1:3,:);       % system state
x1_u = x_all(4:6,:);    % upper estimation bound of observer 1
x1_l = x_all(7:9,:);    % lower estimation bound of observer 1
z_u = x_all(10:12,:);   % intermidiate upper bound of observer 2
z_l = x_all(13:15,:);   % intermidiate lower bound of observer 2
iV2 = inv(V2);
V2p = max(0,iV2); V2m = V2p - iV2;
x2_u = V2p * z_u - V2m * z_l;   % upper estimation bound of observer 2
x2_l = V2p * z_l - V2m * z_u;   % lower estimation bound of observer 2


% Plot interval estimation by observer 1
figure
subplot(311)
plot(t,x(1,:),'-',t,x1_u(1,:),'--',t,x1_l(1,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
h = legend('$x$','$\overline{x}$','$\underline{x}$');
set(h,'Interpreter','latex')
title('State estimation by interval observer 1')

subplot(312)
plot(t,x(2,:),'-',t,x1_u(2,:),'--',t,x1_l(2,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
h = legend('$\dot{x}$','$\dot{\overline{x}}$','$\dot{\underline{x}}$');
set(h,'Interpreter','latex')

subplot(313)
plot(t,x(3,:),'-',t,x1_u(3,:),'--',t,x1_l(3,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$\ddot{x}$','Interpreter','latex')
h = legend('$\ddot{x}$','$\ddot{\overline{x}}$','$\ddot{\underline{x}}$');
set(h,'Interpreter','latex')


% Plot interval estimation by observer 2
figure
subplot(311)
plot(t,x(1,:),'-',t,x2_u(1,:),'--',t,x2_l(1,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
h = legend('$x$','$\overline{x}$','$\underline{x}$');
set(h,'Interpreter','latex')
title('State estimation by interval observer 2')

subplot(312)
plot(t,x(2,:),'-',t,x2_u(2,:),'--',t,x2_l(2,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
h = legend('$\dot{x}$','$\dot{\overline{x}}$','$\dot{\underline{x}}$');
set(h,'Interpreter','latex')

subplot(313)
plot(t,x(3,:),'-',t,x2_u(3,:),'--',t,x2_l(3,:),'-.','LineWidth',1)
xlabel('$t$','Interpreter','latex')
ylabel('$\ddot{x}$','Interpreter','latex')
h = legend('$\ddot{x}$','$\ddot{\overline{x}}$','$\ddot{\underline{x}}$');
set(h,'Interpreter','latex')


% Compare interval estimation
% figure
% subplot(311)
% plot(t,x(1,:),t,x1_u(1,:),'--',t,x1_l(1,:),'--','LineWidth',1)
% hold on
% plot(t,x2_u(1,:),':',t,x2_l(1,:),':','LineWidth',1)
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$','Interpreter','latex')
% h = legend('$x$','upper bound by IO 1','lower bound by IO 1', ...
%            'upper bound by IO 2','lower bound by IO 2');
% set(h,'Interpreter','latex')
% title('Compare interval estimation')
% 
% subplot(312)
% plot(t,x(2,:),t,x1_u(2,:),'--',t,x1_l(2,:),'--','LineWidth',1)
% hold on
% plot(t,x2_u(2,:),':',t,x2_l(2,:),':','LineWidth',1)
% xlabel('$t$','Interpreter','latex')
% ylabel('$\dot{x}$','Interpreter','latex')
% 
% subplot(313)
% plot(t,x(3,:),t,x1_u(3,:),'--',t,x1_l(3,:),'--','LineWidth',1)
% hold on
% plot(t,x2_u(3,:),':',t,x2_l(3,:),':','LineWidth',1)
% xlabel('$t$','Interpreter','latex')
% ylabel('$\ddot{x}$','Interpreter','latex')


% Fault detection
y = C*x;

% Compute residual bounds of interval observer 1
y1_u = C*x1_u;
y1_l = C*x1_l;
r1_u = y1_u - y;
r1_l = y1_l - y;

% Compute residual bounds of interval observer 2
y2_u = C*x2_u;
y2_l = C*x2_l;
r2_u = y2_u - y;
r2_l = y2_l - y;

% Plot fault detection result by interval observer 1
figure
subplot(311)
plot(t,r1_u(1,:),'b--',t,r1_l(1,:),'r--','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k','LineWidth',1)
ylim([-1 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_1$','Interpreter','latex')
h = legend('$\overline{r}_1$','$\underline{r}_1$','$0$');
set(h,'Interpreter','latex')
title('Fault detection by interval observer 1')

subplot(312)
plot(t,r1_u(2,:),'b--',t,r1_l(2,:),'r--','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k','LineWidth',1)
ylim([-1 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_2$','Interpreter','latex')
h = legend('$\overline{r}_2$','$\underline{r}_2$','$0$');
set(h,'Interpreter','latex')

subplot(313)
plot(t,r1_u(3,:),'b-.',t,r1_l(3,:),'r-.','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k-','LineWidth',1)
ylim([-1.5 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_3$','Interpreter','latex')
h = legend('$\overline{r}_3$','$\underline{r}_3$','$0$');
set(h,'Interpreter','latex')


% Plot fault detection result by interval observer 2
figure
subplot(311)
plot(t,r2_u(1,:),'b--',t,r2_l(1,:),'r--','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k','LineWidth',1)
ylim([-1 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_1$','Interpreter','latex')
h = legend('$\overline{r}_1$','$\underline{r}_1$','$0$');
set(h,'Interpreter','latex')
title('Fault detection by interval observer 2')

subplot(312)
plot(t,r2_u(2,:),'b--',t,r2_l(2,:),'r--','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k','LineWidth',1)
ylim([-1 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_2$','Interpreter','latex')
h = legend('$\overline{r}_2$','$\underline{r}_2$','$0$');
set(h,'Interpreter','latex')

subplot(313)
plot(t,r2_u(3,:),'b-.',t,r2_l(3,:),'r-.','LineWidth',1)
hold on
plot([0 t(end)],[0 0],'k-','LineWidth',1)
ylim([-1.5 1])
xlabel('$t$','Interpreter','latex')
ylabel('$r_3$','Interpreter','latex')
h = legend('$\overline{r}_3$','$\underline{r}_3$','$0$');
set(h,'Interpreter','latex')


function dx = ode_fcn(t,x,A,C,Phi,L,Gamma,tPhi,tL,tGamma,V2)
% x(1~3): real system state vector
% x(4~6): upper estimation vector of interval observer 1
% x(7~9): lower estimation vector of interval observer 1
% x(10~12): upper estimation vector of interval observer 2
% x(13~15): lower estimation vector of interval observer 2

% Some signals and matrices
D = [1 -1 -1]';  % unknown input coefficient matrix, i.e., D = [D_1; ...; D_{m-2}; D] in the paper
wb = 4;          % magnitude of unknown inputs
w = wb * (2*rand(1)-1);
% w = wb * sin(t);
v = 20*sin(t)*1 + 20*1; % external signal

% Fully actuated system approach-based control                                                      。
u = -A*[x(1); x(2); x(3)] - 2*x(2)^2 - 2*x(1)*x(3) + v; % control input
if t >= 2
    f = 0;      % set fault here
else
    f = 0;
end
dx(1:3,1) = [x(2);
             x(3);
             2*(x(2)^2 + x(1)*x(3)) + u + f] + D*w;
y = C * x(1:3);   % measurement output

% Interval observer 1
w_u = wb; w_l = -wb;        % bounds of unknown inputs
D_p = max(0,D); D_m = D_p - D;
dx(4:6,1) = Phi * x(4:6) + L * (y - C * x(4:6)) + Gamma * v + D_p * w_u - D_m * w_l;
dx(7:9,1) = Phi * x(7:9) + L * (y - C * x(7:9)) + Gamma * v + D_p * w_l - D_m * w_u;

% Interval observer 2
VD_p = max(0,V2*D);  VD_m = VD_p - V2*D;
dx(10:12,1) = tPhi * x(10:12) + tL * (y - C / V2 * x(10:12)) + tGamma * v + VD_p * w_u - VD_m * w_l;
dx(13:15,1) = tPhi * x(13:15) + tL * (y - C / V2 * x(13:15)) + tGamma * v + VD_p * w_l - VD_m * w_u;
end

