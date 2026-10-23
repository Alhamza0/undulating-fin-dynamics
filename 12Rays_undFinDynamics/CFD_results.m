%this code consist of two parts partA: In order to solve and simulates the
%actuation forces and torques, an inverse dynamics approach is used to
%solve (8). partB;
%Parametric analysis to see what the effect of undulating parameters on
%hydrodynamic forces.
clc 
clear all
close all
%time
%------------------------------------------------------------------------
q= readmatrix('joints_angles.xlsx');
T_normal= readmatrix('Forces_phi50_6s_normal.xlsx');
T_coarse= readmatrix('phi50_6s_coarse.xlsx');
Cd=readmatrix('cd.xlsx');
figure(1)
subplot(2,1,1)
plot(q(:,1),q(:,2),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(q(:,1),q(:,3),'linewidth',1.5,'Color',[0 0 0]+0.6);
plot(q(:,1),q(:,4),'-.','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30','-20','-10','0','10','20','30',})
grid on 
ylabel('angle [\circ]')
%title('(a) motion profile of the rays');


subplot(2,1,2)
plot(T_normal(:,1),T_normal(:,2),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(T_normal(:,1),T_normal(:,3),'linewidth',1.5,'Color',[0 0 0]+0.6);
plot(T_normal(:,1),T_normal(:,4),'-.','Color',[0 0 0]);
hold off
legend('Sway','Heave','Surge')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
%title('(b) hydrodynamic force components');
xlabel('t [s]');

figure(2)
subplot(4,1,1)
plot(q(:,1),q(:,2),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(q(:,1),q(:,3),'linewidth',1.5,'Color',[0 0 0]+0.6);
plot(q(:,1),q(:,4),'-.','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ytixr = get(gca, 'YTick');
ytixlbl = regexp(sprintf('%.1f°\n',ytixr*180/pi), '\n', 'split');
set(gca, 'YTick', ytixr, 'YTickLabel',ytixlbl);
grid on 
ylabel('angle [\circ]')
title('(a) motion profile of the rays');

subplot(4,1,2)
plot(q(:,1),T_coarse(:,3),'Color',[0 0 0]+0.4);
hold on
plot(q(:,1),T_normal(:,3),'-.','Color',[0 0 0]);
hold off
legend('Coarse','Normal');
grid on
ylabel('Force [N]')
title('(b) Heave force');

subplot(4,1,3)
plot(q(:,1),T_coarse(:,2),'Color',[0 0 0]+0.4);
hold on
plot(q(:,1),T_normal(:,2),'-.','Color',[0 0 0]);
hold off
legend('Coarse','Normal')
grid on
ylabel('Force [N]')
title('(c) Sway force');

subplot(4,1,4)
plot(q(:,1),T_coarse(:,4),'Color',[0 0 0]+0.4);
hold on
plot(q(:,1),T_normal(:,4),'-.','Color',[0 0 0]);
hold off
legend('Coarse','Normal')
grid on
xlabel('t [s]');
ylabel('Force [N]')
title('(d) Surge force');

figure(3)
plot(Cd(:,1),Cd(:,2),'linewidth',1.5,'Color',[0 0 0]);
grid on
%title(' Drag Coefficient');
yticks([0 1 2 3 10 15 20 27])
ylabel('Cd')
xlabel('t [s]');


