
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
ncycle=3; %Number of loading cycles
tMax = ncycle/f;
dt = 0.01/f;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.7;
rho=1000;%[kg/m3]
hmin=0.034;
hmax=0.174;
D=0.096;
%------------------------------------------------------------------------
qdq= readmatrix('qdq_step.xlsx');%joints angles and velocities
dq= readmatrix('dq.xlsx');%joints velocities
q1_20=qdq(:,1);%angle of joint 1
q2_20=qdq(:,2);%angle of joint 2
dq1_20=dq(:,1);%velcity of joint 1
dq2_20=dq(:,2);%velcity of joint 2

q1_30=qdq(:,5);%angle of joint 1
q2_30=qdq(:,6);%angle of joint 2
dq1_30=dq(:,3);%velcity of joint 1
dq2_30=dq(:,4);%velcity of joint 2

q1_40=qdq(:,9);%angle of joint 1
q2_40=qdq(:,10);%angle of joint 2
dq1_40=dq(:,5);%velcity of joint 1
dq2_40=dq(:,6);%velcity of joint 2

q1_50=qdq(:,13);%angle of joint 1
q2_50=qdq(:,14);%angle of joint 2
dq1_50=dq(:,7);%velcity of joint 1
dq2_50=dq(:,8);%velcity of joint 2

 %------------------------------------------------------------------------
 for (indT = 1:nT)
    
    qq1 = q1_20(indT);
 
    qq2 = q2_20(indT);
    
    dqq1=dq1_20(indT);
    
    dqq2=dq2_20(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq1]);
 FN_phi0=hydrodynamicforces_Func(D,hmax,[qq1 qq1],[dqq1 dqq1]);
 F20(indT, :)=FN;
 F0(indT, :)=FN_phi0;
 end
 
  for (indT = 1:nT)
    
    qq1 = q1_30(indT);
 
    qq2 = q2_30(indT);
    
    dqq1=dq1_30(indT);
    
    dqq2=dq2_30(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq1]);
 F30(indT, :)=FN;
  end
  for (indT = 1:nT)
    
    qq1 = q1_40(indT);
 
    qq2 = q2_40(indT);
    
    dqq1=dq1_40(indT);
    
    dqq2=dq2_40(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq1]);
 F40(indT, :)=FN;
  end

  for (indT = 1:nT)
    
    qq1 = q1_50(indT);
 
    qq2 = q2_50(indT);
    
    dqq1=dq1_50(indT);
    
    dqq2=dq2_50(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq1]);
 F50(indT, :)=FN;
  end
FCFD_50= readmatrix('F_phi50.xlsx');
FCFD_40= readmatrix('F_phi40.xlsx');
FCFD_30= readmatrix('F_phi30.xlsx');
FCFD_20= readmatrix('F_phi20.xlsx');
FCFD_0= readmatrix('F_phi0.xlsx');
%--------------------------------------------------------------------
figure(1)
%phi=50;
subplot(3,1,1)
%plot(timeVec,F50(:,2),'Color',[0 0 0]);
plot(timeVec,F50(:,2),'linewidth',3);
hold on
%plot(timeVec,FCFD_50(:,1),'-.','Color',[0 0 0]);
plot(timeVec,FCFD_50(:,1),'linewidth',3);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
%plot(timeVec,F50(:,1),'Color',[0 0 0]);
plot(timeVec,F50(:,1),'linewidth',3);
hold on
%plot(timeVec,FCFD_50(:,2),'-.','Color',[0 0 0]);
plot(timeVec,FCFD_50(:,2),'linewidth',3);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
%plot(timeVec,F50(:,3),'Color',[0 0 0]);
plot(timeVec,F50(:,3),'linewidth',3);
hold on
%plot(timeVec,FCFD_50(:,3),'-.','Color',[0 0 0]);
plot(timeVec,FCFD_50(:,3),'linewidth',3);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
figure(2)
%phi=40;
subplot(3,1,1)
plot(timeVec,F40(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_40(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F40(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_40(:,2),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F40(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_40(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
figure(3)
%phi=20;
subplot(3,1,1)
plot(timeVec,F30(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_30(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F30(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_30(:,2),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F30(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_30(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');


figure(4)
%phi=20;
subplot(3,1,1)
plot(timeVec,F20(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_20(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F20(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_20(:,2),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F20(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_20(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');

%--------------------------------------------------------------------
figure(5)
%phi=20;
subplot(3,1,1)
plot(timeVec,F0(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_0(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F0(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_0(:,2),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F0(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,FCFD_0(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.', 'CFD')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
