%this code consist of two parts partA: In order to solve and simulates the
%actuation forces and torques, an inverse dynamics approach is used to
%solve (8). partB;
%Parametric analysis to see what the effect of undulating parameters on
%hydrodynamic forces.
%% partA
clear all
close all
%------------------------------------------------------------------------
global Cn rho  A f N X hmin
%------------------------------------------------------------------------
% Traveling wave parameters(sinusoidal path)
N = 3;%number of rays 
A = deg2rad(30);%Amplitude
f = 1.02;%frequency
%Phi = deg2rad(20);%the inter-ray phase shift
X= deg2rad(0);%angular offset
%------------------------------------------------------------------------
%T_1= readmatrix('phi60_6s_normal.xlsx');
T_normal= readmatrix('forces_coarse -2s.xlsx');
%time
ncycle=3; %Number of loading cycles
tMax = ncycle;
dt = 1/(360*f);
timeVec = T_normal(:,1).';
nT = length(timeVec);
%------------------------------------------------------------------------
%phase shift
phiMax = deg2rad(50);
dphi = deg2rad(10);
phiVec = [deg2rad(20):dphi:phiMax].';
nphi = length(phiVec);
%------------------------------------------------------------------------
%frequency
fMax = 1.018;
df = 0.018;
fVec = [1:df:fMax].';
nf = length(fVec);
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1000;%[kg/m3]
hmin=0.034;
hmax=0.174;
D=0.096;
%------------------------------------------------------------------------
for (indphi = 1:nphi)
    Phi= phiVec(indphi);
for (indT = 1:nT)
    
    tt = timeVec(indT);
 %------------------------------------------------------------------------
    for (i = 1:N)
        %To define the generalized coordinates of the fin system
        % q=[q1 q2 q3]
        
        q_curr(i) = A*sin(2*pi*f*tt - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*tt);
        ddq_curr(i) = 4*A*f^2*pi^2*sin(Phi*(i - 1) - 2*pi*f*tt);
    end
  %------------------------------------------------------------------------
% the generalized coordinates of the system (joints position and velocity)
        %q=[q1 q2 qn]
        
    q(indT,:) = (q_curr);
    dq(indT,:) = (dq_curr);
    ddq(indT,:) = (ddq_curr);
 %------------------------------------------------------------------------
 %hydrodynamic torques on each joint
 hT= hydrodynamictorques_Func(D,hmax,q_curr,dq_curr);
    T_phi(indT, :) = hT;
 %------------------------------------------------------------------------
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2));
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F(indT, :)=FN;
end
F1(:,indphi)=F(:,1);%when f=1.018hz
F2(:,indphi)=F(:,2);
F3(:,indphi)=F(:,3);
end

for (indphi = 1:nphi)
    Phi= phiVec(indphi);
for (indT = 1:nT)
    f=1;
    tt = timeVec(indT);
 %------------------------------------------------------------------------
    for (i = 1:N)
        %To define the generalized coordinates of the fin system
        % q=[q1 q2 q3]
        
        q_curr(i) = A*sin(2*pi*f*tt - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*tt);
        ddq_curr(i) = 4*A*f^2*pi^2*sin(Phi*(i - 1) - 2*pi*f*tt);
    end
  %------------------------------------------------------------------------
% the generalized coordinates of the system (joints position and velocity)
        %q=[q1 q2 qn]
        
    q(indT,:) = (q_curr);
    dq(indT,:) = (dq_curr);
    ddq(indT,:) = (ddq_curr);
 %------------------------------------------------------------------------
 %hydrodynamic torques on each joint
 hT= hydrodynamictorques_Func(D,hmax,q_curr,dq_curr);
    T_phi(indT, :) = hT;
 %------------------------------------------------------------------------
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2));
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F(indT, :)=FN;
end
F11(:,indphi)=F(:,1);
F21(:,indphi)=F(:,2);
F31(:,indphi)=F(:,3);
end
FCFD_50= readmatrix('F_phi50_f1_02hz.xlsx');
FCFD_40= readmatrix('F_phi40.xlsx');
FCFD_30= readmatrix('F_phi30.xlsx');
FCFD_20= readmatrix('F_phi20.xlsx');
%--------------------------------------------------------------------
figure(1)
%phi=20;
subplot(3,1,1)
plot(timeVec,F2(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_20(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F1(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_20(:,2),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F3(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_20(:,3),'-.','Color',[0 0 0]);
%plot(timeVec,T_normal(:,9),':','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
figure(2)
%phi=30;
subplot(3,1,1)
%for diffrent f
plot(timeVec,F2(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,2),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F1(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,2),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,2),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F3(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,2),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
%-------------------------------------------------------
figure(3)
%phi=40;
subplot(3,1,1)
plot(timeVec,F2(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,3),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F1(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,3),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,2),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F3(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,3),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
%-------------------------------------------------------
figure(4)
%phi=50;
subplot(3,1,1)
%for diffrent f
plot(timeVec,F2(:,4),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,4),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_50(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
xlabel('t [s]');

subplot(3,1,2)
plot(timeVec,F1(:,2),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,4),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_50(:,2),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');

subplot(3,1,3)
plot(timeVec,F3(:,4),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,4),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_50(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.. (\it f = 1.02Hz)','An.. (\it f = 1Hz)',...
   'CFD (\it f = 1Hz)')
grid on
%xlim([1 6])
ylabel('Surge[N]')
xlabel('t [s]');
%-------------------------------------------------------
%-------------------------------------------------------
figure(5)
%ana (phi20) CFD (40)
subplot(3,1,1)
plot(timeVec,F2(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =20\circ ; \it f = 1.02Hz)','An..(\phi =20\circ ; \it f = 1Hz)',...
   'CFD (\phi =40\circ)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
%title('(b) Sway force');
xlabel('t [s]');
%-------------------------------------------------------
subplot(3,1,2)
plot(timeVec,F1(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,2),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =20\circ ; \it f = 1.02Hz)','An..(\phi =20\circ ; \it f = 1Hz)',...
   'CFD (\phi =40\circ)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');
%title('(c) Heave force');
%-------------------------------------------------------
subplot(3,1,3)
plot(timeVec,F3(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_40(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =20\circ ; \it f = 1.02Hz)','An..(\phi =20\circ ; \it f = 1Hz)',...
   'CFD (\phi =40\circ)')
grid on
%xlim([1 6])
xlabel('t [s]');
ylabel('Surge[N]')
%title('(d) Surge force');

%-------------------------------------------------------
%-------------------------------------------------------
figure(6)
%ana (phi15) CFD (30)
subplot(3,1,1)
plot(timeVec,F2(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F21(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,1),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =15\circ ; \it f = 1.02Hz)','An..(\phi =15\circ ; \it f = 1Hz)',...
   'CFD (\phi =30\circ ;\it f = 1.02Hz)')
grid on
%xlim([1 6])
ylabel('Heave[N]')
%title('(b) Sway force');
xlabel('t [s]');
%-------------------------------------------------------
subplot(3,1,2)
plot(timeVec,F1(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F11(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,2),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =15\circ ; \it f = 1.02Hz)','An..(\phi =15\circ ; \it f = 1Hz)',...
   'CFD (\phi =30\circ ;\it f = 1.02Hz)')
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');
%title('(c) Heave force');
%-------------------------------------------------------
subplot(3,1,3)
plot(timeVec,F3(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F31(:,1),'Color',[0 0 0]+0.6);
plot(timeVec,FCFD_30(:,3),'-.','Color',[0 0 0]);
hold off
legend('An.. (\phi =15\circ ; \it f = 1.02Hz)','An..(\phi =15\circ ; \it f = 1Hz)',...
   'CFD (\phi =30\circ ;\it f = 1.02Hz)')
grid on
%xlim([1 6])
xlabel('t [s]');
ylabel('Surge[N]')
%title('(d) Surge force');