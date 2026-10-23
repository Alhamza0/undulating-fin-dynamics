%this code consist of two parts partA:Solve the hydrodynamic force components
%part B Parametric analysis to see what the effect of undulating parameters on
%hydrodynamic forces.
%% partA
clear all
close all
%------------------------------------------------------------------------
global Cn rho hmin  Phi A f N X
%------------------------------------------------------------------------
% Traveling wave parameters(sinusoidal path)
N = 12;%number of rays 
A = deg2rad(30);%Amplitude
f = 1;%frequency
Phi = deg2rad(60);%the inter-ray phase shift 
X= deg2rad(0);%angular offset 
%------------------------------------------------------------------------
%time
tMax = 2;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%------------------------------------------------------------------------
%oscillation frequency,
fMax = 4;
df = 0.01;
freqVec = [0:df:fMax].';
nF = length(freqVec);
Famax = 2;
dFa = 0.5;
FaVec = [dFa:dFa:Famax].';
nFa = length(FaVec);
%-----------------------------------
%oscillation amplitude,
Amax = deg2rad(40);
dA = 0.01;
AmVec = [0:dA:Amax].';
nA = length(AmVec);
%oscillation amplitude for anlaysis 2,
Aamax = deg2rad(40);
dAa = deg2rad(10);
AaVec = [dAa:dAa:Amax].';
nAa = length(AaVec);
%-----------------------------------
%inter-ray phase shift
PhiMax = deg2rad(80);
dPhi = deg2rad(0.8);
PhiVec = [0:dPhi:PhiMax].';
nPhi = length(PhiVec);
%-----------------------------------
%Angular offset
XMax = deg2rad(90);
XMin=deg2rad(0);
dX = 0.01;
XVec = [XMin:dX:XMax].';
nX = length(XVec);
%Angular offset(figure6)
Xamax = deg2rad(40);
dXa = deg2rad(10);
XaVec = [dXa:dXa:Xamax].';
nXa = length(XaVec);
%-----------------------------------
%The membrane parameters
hmax = 0.0015; %[m]
Hmax = 0.005;
dH = 0.0005;%[m]
HVec = [dH:dH:Hmax].';
nH = length(HVec);
hmin = 0; %[m]
%-----------------------------------
D = 0.0015; %[m]
Dmax = 0.002;
dD = 0.0005;
DVec = [dD:dD:Dmax].';
nD = length(DVec);
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1000;%[kg/m3]
%------------------------------------------------------------------------
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
 FN=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))...
     +hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3))...
     +hydrodynamicforces_Func(D,hmax,q_curr(3:4),dq_curr(3:4))...
     +hydrodynamicforces_Func(D,hmax,q_curr(4:5),dq_curr(4:5))...
     +hydrodynamicforces_Func(D,hmax,q_curr(5:6),dq_curr(5:6))...
     +hydrodynamicforces_Func(D,hmax,q_curr(6:7),dq_curr(6:7))...
     +hydrodynamicforces_Func(D,hmax,q_curr(7:8),dq_curr(7:8))...
     +hydrodynamicforces_Func(D,hmax,q_curr(8:9),dq_curr(8:9))...
     +hydrodynamicforces_Func(D,hmax,q_curr(9:10),dq_curr(9:10))...
     +hydrodynamicforces_Func(D,hmax,q_curr(10:11),dq_curr(10:11))...
     +hydrodynamicforces_Func(D,hmax,q_curr(11:12),dq_curr(11:12));
 F_N(indT, :)=FN;
%------------------------------------------------------------------------
end
%------------------------------------------------------------------------
figure(1)
%plot(timeVec,F_N(:,1),'linewidth',1.5,'Color',[0 0 0]);
plot(timeVec,F_N(:,1),'linewidth',3);
hold on
%plot(timeVec,F_N(:,2),'linewidth',1.5,'Color',[0 0 0]+0.6);
plot(timeVec,F_N(:,2),'linewidth',3);
%plot(timeVec,F_N(:,3),':','Color',[0 0 0]);
plot(timeVec,F_N(:,3),'g','linewidth',3);
hold off
legend('F_{Sway}','F_{Heave}','F_{surge}')
grid on
xlabel('t [s]');
ylabel('Force [N]')
title('(c) hydrodynamic force components');
%------------------------------------------------------------------------
%% Part B
%-----------------------------------
PhiFN=PhiFN(nPhi,PhiVec,timeVec,nT,D,hmax);
%-------------------------------------------------------------------
FN_freq = freqFN(nF,freqVec,timeVec,nT,D,hmax);
%-------------------------------------------------------------------
FN_Amp = AmpFN(nA,AmVec,timeVec,nT,D,hmax);
%-------------------------------------------------------------------
XFN = XFN(nX,XVec,timeVec,nT,D,hmax);
%-------------------------------------------------------------------
Phi_A =Phi_A_FUN(nPhi,PhiVec,timeVec,nT,nAa,AaVec,D,hmax);
%-------------------------------------------------------------------
Phi_fr = Phi_F_FUN(nPhi,PhiVec,timeVec,nT,nFa,FaVec,D,hmax);
%-------------------------------------------------------------------
Phi_D = Phi_D_FUN(nPhi,PhiVec,timeVec,nT,nD,DVec,hmax);
%-------------------------------------------------------------------
Phi_H = Phi_H_FUN(nPhi,PhiVec,timeVec,nT,D,nH,HVec);
%-------------------------------------------------------------------
figure (2)
plot(freqVec,FN_freq(:,1),'Color',[0 0 0]);
hold on
plot(freqVec,FN_freq(:,2),':','Color',[0 0 0]);
plot(freqVec,FN_freq(:,3),'-.','Color',[0 0 0]);
hold off
legend('F_{Sway}','F_{Heave}','F_{surge}')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0.5 fMax])
grid on
ylabel('Force [N]')
xlabel('\it f [Hz]' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)
%-------------------------------------------------------------------
figure (2)
plot(AmVec,FN_Amp(:,1),'Color',[0 0 0]);
hold on
plot(AmVec,FN_Amp(:,2),':','Color',[0 0 0]);
plot(AmVec,FN_Amp(:,3),'-.','Color',[0 0 0]);
hold off
legend('F_{heave}','F_{sway}','F_{surge}')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 Amax])
xticks([0 deg2rad(10) deg2rad(20) deg2rad(30) deg2rad(40)])
xticklabels({'0','10','20','30','40'})
grid on
ylabel('Force [N]')
xlabel('\it A [\circ]')
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)
%-------------------------------------------------------------------
figure (3)
[maxy,idx] = max(PhiFN(:,3));
FX=PhiFN(:,3);
plot(PhiVec,PhiFN(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,PhiFN(:,2),':','Color',[0 0 0]);
plot(PhiVec,PhiFN(:,3),'-.','Color',[0 0 0]);
plot(PhiVec(idx),FX(idx),'*k');
hold off
legend('\it F_{heave}','\it F_{sway}','\it F_{surge}', sprintf(' \\Phi= %0.1f',rad2deg(PhiVec(idx))) );
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('Force [N]')
xlabel('\it \phi [\circ]')
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)
%-------------------------------------------------------------------
figure (4)
plot(XVec,XFN(:,1),'Color',[0 0 0]);
hold on
plot(XVec,XFN(:,2),':','Color',[0 0 0]);
plot(XVec,XFN(:,3),'-.','Color',[0 0 0]);
hold off
legend('F_{heave}','F_{sway}','F_{surge}')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([XMin XMax])
xticks([deg2rad(0) deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('Force [N]')
xlabel('\it X [\circ]')
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)



figure(5)
[maxy1,idx1] = max(Phi_A(:,1));
[maxy2,idx2] = max(Phi_A(:,2));
[maxy3,idx3] = max(Phi_A(:,3));
[maxy4,idx4] = max(Phi_A(:,4));
F1=Phi_A(:,1);
F2=Phi_A(:,2);
F3=Phi_A(:,3);
F4=Phi_A(:,4);
plot(PhiVec,Phi_A(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_A(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_A(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_A(:,4),'--','Color',[0 0 0]);
plot(PhiVec(idx1),F1(idx1),'*k');
plot(PhiVec(idx2),F2(idx2),'*k');
plot(PhiVec(idx3),F3(idx3),'*k');
plot(PhiVec(idx4),F4(idx4),'*k');
hold off
legend('\it A=10\circ','\it A=20\circ','\it A=30\circ','\it A=40\circ');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge} [N]')
xlabel('\it \phi [\circ]')
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl



figure (6)
[maxy1,idx1] = max(Phi_fr(:,1));
[maxy2,idx2] = max(Phi_fr(:,2));
[maxy3,idx3] = max(Phi_fr(:,3));
[maxy4,idx4] = max(Phi_fr(:,4));
F1=Phi_fr(:,1);
F2=Phi_fr(:,2);
F3=Phi_fr(:,3);
F4=Phi_fr(:,4);
plot(PhiVec,Phi_fr(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_fr(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_fr(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_fr(:,4),'--','Color',[0 0 0]);
plot(PhiVec(idx1),F1(idx1),'*k');
plot(PhiVec(idx2),F2(idx2),'*k');
plot(PhiVec(idx3),F3(idx3),'*k');
plot(PhiVec(idx4),F4(idx4),'*k');
hold off
legend('\it f=0.5Hz','\it f=1Hz','\it f=1.5Hz','\it f=2Hz');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')

figure (7)
[maxy1,idx1] = max(Phi_D(:,1));
[maxy2,idx2] = max(Phi_D(:,2));
[maxy3,idx3] = max(Phi_D(:,3));
[maxy4,idx4] = max(Phi_D(:,4));
F1=Phi_D(:,1);
F2=Phi_D(:,2);
F3=Phi_D(:,3);
F4=Phi_D(:,4);
plot(PhiVec,Phi_D(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_D(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_D(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_D(:,4),'--','Color',[0 0 0]);
plot(PhiVec(idx1),F1(idx1),'*k');
plot(PhiVec(idx2),F2(idx2),'*k');
plot(PhiVec(idx3),F3(idx3),'*k');
plot(PhiVec(idx4),F4(idx4),'*k');
hold off
legend('\it D=0.5mm','\it D=1mm','\it D=1.5mm','\it D=2mm');
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')

figure(8)
[maxy2,idx2] = max(Phi_H(:,2));
[maxy3,idx3] = max(Phi_H(:,3));
[maxy4,idx4] = max(Phi_H(:,4));
[maxy5,idx5] = max(Phi_H(:,5));
F2=Phi_H(:,2);
F3=Phi_H(:,3);
F4=Phi_H(:,4);
F5=Phi_H(:,5);
plot(PhiVec,Phi_H(:,2),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_H(:,3),':','Color',[0 0 0]);
plot(PhiVec,Phi_H(:,4),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_H(:,5),'--','Color',[0 0 0]);
plot(PhiVec(idx2),F2(idx2),'*k');
plot(PhiVec(idx3),F3(idx3),'*k');
plot(PhiVec(idx4),F4(idx4),'*k');
plot(PhiVec(idx5),F5(idx5),'*k');
hold off
legend('\it H=1mm','\it H=1.5mm','\it H=2mm','\it H=2.5mm');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')