%This code is used to compare CFD & analytical model results for Two_ray fin.
clear all
close all
%------------------------------------------------------------------------
global Cn rho Phi A f N X hmin
%------------------------------------------------------------------------
% Traveling wave parameters(sinusoidal path)
N = 3;%number of rays 
A = deg2rad(30);%Amplitude
f = 1;%frequency
Phi = deg2rad(30);%the inter-ray phase shift
X= deg2rad(0);%angular offset 
%------------------------------------------------------------------------
%time
ncycle=6; %Number of loading cycles
tMax = ncycle/f;
dt = 0.01/f;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1000;%[kg/m3]
hmin=0.034;
hmax=0.174;
D=0.096;
%------------------------------------------------------------------------
qdq= readmatrix('qdq.xlsx');%joints angles and velocities
q1_20=qdq(:,1);%angle of joint 1
q2_20=qdq(:,2);%angle of joint 2
dq1_20=qdq(:,3);%velcity of joint 1
dq2_20=qdq(:,4);%velcity of joint 2

q1_30=qdq(:,5);%angle of joint 1
q2_30=qdq(:,6);%angle of joint 2
dq1_30=qdq(:,7);%velcity of joint 1
dq2_30=qdq(:,8);%velcity of joint 2

q1_40=qdq(:,9);%angle of joint 1
q2_40=qdq(:,10);%angle of joint 2
dq1_40=qdq(:,11);%velcity of joint 1
dq2_40=qdq(:,12);%velcity of joint 2

q1_50=qdq(:,13);%angle of joint 1
q2_50=qdq(:,14);%angle of joint 2
dq1_50=qdq(:,15);%velcity of joint 1
dq2_50=qdq(:,16);%velcity of joint 2

q_2= readmatrix('joints_q_dq.xlsx');%joints angles and velocities
q3=q_2(:,2);%angle of joint 1
q4=q_2(:,3);%angle of joint 2
dq3=q_2(:,5);%velcity of joint 1
dq4=q_2(:,6);%velcity of joint 2
 %------------------------------------------------------------------------
 for (indT = 1:nT)
    
    qq1 = q1_20(indT);
 
    qq2 = q2_20(indT);
    
    dqq1=dq1_20(indT);
    
    dqq2=dq2_20(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq2]);
 FN_phi0=hydrodynamicforces_Func(D,hmax,[qq1 qq1],[dqq1 dqq1]);
 F_N_20(indT, :)=FN;
 end
 
  for (indT = 1:nT)
    
    qq1 = q1_30(indT);
 
    qq2 = q2_30(indT);
    
    dqq1=dq1_30(indT);
    
    dqq2=dq2_30(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq2]);
 F_N_30(indT, :)=FN;
  end
  for (indT = 1:nT)
    
    qq1 = q1_40(indT);
 
    qq2 = q2_40(indT);
    
    dqq1=dq1_40(indT);
    
    dqq2=dq2_40(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq2]);
 F_N_40(indT, :)=FN;
  end

  for (indT = 1:nT)
    
    qq1 = q1_50(indT);
 
    qq2 = q2_50(indT);
    
    dqq1=dq1_50(indT);
    
    dqq2=dq2_50(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq2]);
 FN_phi0=hydrodynamicforces_Func(D,hmax,[qq1 qq1],[dqq1 dqq1]);
 F_N_50(indT, :)=FN;
 end
 
%T_coarse= readmatrix('phi50_12s_coarse_0.5f');
T_normal= readmatrix('Forces_phi50_6s_normal.xlsx');% FN from CFD
T_normal_0= readmatrix('phi0_6s_normal.xlsx');% FN from CFD
ddq1=q3-q4;
%--------------------------------------------------------------------

figure(1)

subplot(2,1,1)
plot(timeVec,q_2(:,2),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,q_2(:,3),'Color',[0 0 0]+0.6);
plot(timeVec,ddq1,':','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30\circ','-20\circ','-10\circ','0','10\circ','20\circ','30\circ',})
grid on 
ylabel('angle [\circ]')
title('(a) motion profile of the rays');

subplot(2,1,2)
plot(timeVec,F_N_35(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,F_N_35(:,1),'-.','Color',[0 0 0]);
plot(timeVec,F_N_35(:,2),':','Color',[0 0 0]);
hold off
legend('F_{Surge}','F_{Heave}','F_{Sway}')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
xlabel('t [s]');
title('(c) hydrodynamic force components');

%--------------------------------------------------------------------
figure(2)
subplot(2,1,1)
plot(timeVec,q_2(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,q_2(:,3),'Color',[0 0 0]+0.4);
plot(timeVec,ddq1,':','Color',[0 0 0]);
%y1= linspace(min(q),max(ddq1),10);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30\circ','-20\circ','-10\circ','0','10\circ','20\circ','30\circ',})
grid on 
ylabel('angle [\circ]')
%xline(0.56,':',{'Max q_{1}-q_{2} '});
%title('(a) motion profile of the rays');

subplot(2,1,2)
plot(timeVec,F_N_50(:,2),'Color',[0 0 0]);
hold on
%plot(timeVec,T_coarse(:,3));
plot(timeVec,T_normal(:,3),'-.','Color',[0 0 0]);

hold off
legend('Analytical','CFD')%,'CFD-Normal')
grid on
ylabel('Heave[N]')
%title('(b) Sway force');
xlabel('t [s]');
%-------------------------------------------------------
figure(3)

subplot(2,1,1)
plot(timeVec,q_2(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,q_2(:,3),'Color',[0 0 0]+0.4);
plot(timeVec,ddq1,':','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30\circ','-20\circ','-10\circ','0','10\circ','20\circ','30\circ',})
grid on 
ylabel('angle [\circ]')
%title('(a) motion profile of the rays');

subplot(2,1,2)
plot(timeVec,F_N_50(:,1),'Color',[0 0 0]);
hold on
%plot(timeVec,T_coarse(:,2),);
plot(timeVec,T_normal(:,2),'-.','Color',[0 0 0]);
hold off
legend('Analytical','CFD')%,'CFD-Normal')
grid on
ylabel('Sway[N]')
xlabel('t [s]');
%title('(c) Heave force');
%-------------------------------------------------------
figure(4)

subplot(2,1,1)
plot(timeVec,q_2(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,q_2(:,3),'Color',[0 0 0]+0.4);
plot(timeVec,ddq1,':','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30\circ','-20\circ','-10\circ','0','10\circ','20\circ','30\circ',})
grid on 
ylabel('angle [\circ]')
%title('(a) motion profile of the rays');

subplot(2,1,2)
plot(timeVec,F_N_50(:,3),'Color',[0 0 0]);
hold on
%plot(timeVec,T_coarse(:,4));
plot(timeVec,T_normal(:,4),'-.','Color',[0 0 0]);

hold off
legend('Analytical','CFD')%,'CFD-Normal')
grid on
xlabel('t [s]');
ylabel('Surge[N]')
%title('(d) Surge force');
%-------------------------------------------------------
figure(5)
subplot(3,1,1)
plot(timeVec,F_N_50(:,2),'-.','Color',[0 0 0]);
hold on
plot(timeVec,F_N_35(:,2),':','Color',[0 0 0]);
hold on
%plot(timeVec,T_coarse(:,3));
plot(timeVec,T_normal(:,3),'Color',[0 0 0]);

hold off
legend('Analytical:\phi=50\circ','Analytical:\phi=35\circ','CFD:\phi=50\circ')%,'CFD-Normal')grid on
%xlim([4 6])
ylabel('Heave[N]')
%title('(b) Sway force');
%-------------------------------------------------------
subplot(3,1,2)
plot(timeVec,F_N_50(:,1),'-.','Color',[0 0 0]);
hold on
plot(timeVec,F_N_35(:,1),':','Color',[0 0 0]);
%plot(timeVec,T_coarse(:,2),);
plot(timeVec,T_normal(:,2),'Color',[0 0 0]);
hold off
legend('Analytical:\phi=50\circ','Analytical:\phi=35\circ','CFD:\phi=50\circ')%,'CFD-Normal')grid on
%xlim([4 6])
ylabel('Sway[N]')
%title('(c) Heave force');
%-------------------------------------------------------
subplot(3,1,3)
plot(timeVec,F_N_50(:,3),'-.','Color',[0 0 0]);
hold on
plot(timeVec,F_N_35(:,3),':','Color',[0 0 0]);
%plot(timeVec,T_coarse(:,4));
plot(timeVec,T_normal(:,4),'Color',[0 0 0]);
hold off
legend('Analytical:\phi=50\circ','Analytical:\phi=35\circ','CFD:\phi=50\circ')%,'CFD-Normal')
grid on
%xlim([4 6])
xlabel('t [s]');
ylabel('Surge[N]')
%title('(d) Surge force');
%-------------------------------------------------------