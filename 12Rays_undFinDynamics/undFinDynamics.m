%this code consist of two parts partA:Solve the hydrodynamic force components
%part B Parametric analysis to see what the effect of undulating parameters on
%hydrodynamic forces.
%% partA
clear all
close all
%------------------------------------------------------------------------
global l v Cn rho hmin alpha Phi A f N X
%------------------------------------------------------------------------
% Traveling wave parameters(sinusoidal path)
N = 3;%number of rays 
A = deg2rad(30);%Amplitude
f = 2;%frequency
Phi = deg2rad(40);%the inter-ray phase shift 
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
 %number of rays,
NMax = 6;
dN = 2;
%NVec = [2:dN:NMax].';
NVec = [2;3;4;5;6].';
nN = length(NVec);
%-----------------------------------
%The membrane parameters
hmax = 0.175; %[m]
Hmax = 0.2;
dH = 0.05;
HVec = [dH:dH:Hmax].';
nH = length(HVec);
hmin = 0; %[m]
%-----------------------------------
D = 0.10; %[m]
Dmax = 0.2;
dD = 0.05;
DVec = [dD:dD:Dmax].';
nD = length(DVec);
%-----------------------------------
l = 0.08;
v = 0.5;
EK = 360; %[N/m]
alpha=1;
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1022;%[kg/m3]
%------------------------------------------------------------------------
%%Fin ray parameters
g = 9.81; %gravitational accelaration [m/s^2]
mi = 0.01;m0 = mi;m1 = mi;m2 = mi;m3=mi;%mass of the link
li = 0.12;l0 = li;l1 = li;l2 = li;l3=li;%length of the ray
Izi = 72e-6;Iz0 = Izi;Iz1 = Izi;Iz2 = Izi;Iz3 = Izi;%moment of inertia of the link 
bi = 0.02;b0 = bi;b1 = bi;b2 = bi;b3 = bi;%the viscous friction coefficient
kci = 0.001;kc0 = kci;kc1 = kci;kc2 = kci;kc3 = kci;%the Coulomb friction
g1 = (l1/2)*m1*g;g2 = (l2/2)*m2*g;g3 = (l3/2)*m3*g;
%------------------------------------------------------------------------
%the mass matrix of the fin system
M0 = m1 + m2+m3;
M1 = m1*(l1/2)^2 + Iz1;
M2 = m2*(l2/2)^2 + Iz2;
M3 = m3*(l3/2)^2 + Iz3;
M = diag([M1,M2,M3]);
%------------------------------------------------------------------------
B = diag([b1,b2,b3]);
Kc = diag([kc1,kc2,kc3]);
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
 %the generalized gravity vector 
    G = [g1*sin(q(indT,1));g2*sin(q(indT,2)) ;g3*sin(q(indT,3))];
 %------------------------------------------------------------------------
 %hydrodynamic torques on each joint
 hT= hydrodynamictorques_Func(D,hmax,q_curr,dq_curr);
    T_phi(indT, :) = hT;
 %------------------------------------------------------------------------
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3));
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F_N(indT, :)=FN;
%------------------------------------------------------------------------
%the force on each joint, due to gravity, inertia and fraction.
L(indT,:) = M*(ddq(indT,:).') + B*(dq(indT,:).') + Kc*sign(dq(indT,:).') + G ;
%------------------------------------------------------------------------
%the generalized force vector for the joints due to the elastic...
%deformation of the membrane 
    for (i = 1:N)
        if i < N
            j=1;
        tauG(:,i) = (-1)^(i+1) * sin(q(indT,j) - q(indT,j+1)) * EK * (DeltaFunc(D,hmax,q(indT,[j:j+1])) - DeltaFunc(D,hmin,q(indT,[j:j+1])));
        else
            j=2;
        tauG(:,i) = (-1)^(i+1) * sin(q(indT,j) - q(indT,j+1)) * EK * (DeltaFunc(D,hmax,q(indT,[j:j+1])) - DeltaFunc(D,hmin,q(indT,[j:j+1])));
        end
    end
    T_G(indT,:) = tauG;
%------------------------------------------------------------------------
%total torque on each joint
    T(indT,:) = L(indT,:) + tauG + hT';
end

ddq1=q(:,1)-q(:,2);
ddq2=q(:,2)-q(:,3);
ddq3=q(:,1)-q(:,3);
figure(1)
subplot(6,1,1)
plot(timeVec,q(:,1),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,q(:,2),'linewidth',1.5,'Color',[0 0 0]+0.4);
plot(timeVec,q(:,3),'linewidth',1.5,'Color',[0 0 0]+0.6);
hold off
legend('q_{1}','q_{2}','q_{3}')
ytixr = get(gca, 'YTick');
ytixlbl = regexp(sprintf('%.1f�\n',ytixr*180/pi), '\n', 'split');
set(gca, 'YTick', ytixr, 'YTickLabel',ytixlbl);
grid on 
ylabel('angle [\circ]')
title('(a) motion profile of the rays');

%subplot(6,1,2)
%plot(timeVec,ddq1,'-k',timeVec,ddq2,'-b',timeVec,ddq3,'-r','LineWidth',2.0);
%legend('q_1-q_2','q_2-q_3','q_1-q_3')
%grid on 
%ylabel('angle [rad]')
%title('(b) Angle separation');

subplot(6,1,2)
plot(timeVec,F_N(:,1),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,F_N(:,2),'linewidth',1.5,'Color',[0 0 0]+0.4);
plot(timeVec,F_N(:,3),'linewidth',1.5,'Color',[0 0 0]+0.6);
hold off
legend('F_{heave}','F_{sway}','F_{surge}')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
title('(c) hydrodynamic force components');

subplot(6,1,3)
plot(timeVec,T_phi(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,T_phi(:,2),'Color',[0 0 0]+0.4);
plot(timeVec,T_phi(:,3),'Color',[0 0 0]+0.6);
hold off
legend('\tau_{\Phi,1}','\tau_{\Phi,2}','\tau_{\Phi,3}')
grid on 
ylabel('Torque [Nm]')
title('(d) fin ray joint torques due to hydrodynamic loads');

subplot(6,1,4)
plot(timeVec,T_G(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,T_G(:,2),'Color',[0 0 0]+0.4);
plot(timeVec,T_G(:,3),'Color',[0 0 0]+0.6);
hold off
legend('\tau_{\Gamma,1}','\tau_{\Gamma,1}','\tau_{\Gamma,1}')
title('(e) fin ray joint torques due to the elastic membrane')
grid on
ylabel('Torque [Nm]')

subplot(6,1,5)
plot(timeVec,L(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,L(:,2),'Color',[0 0 0]+0.4);
plot(timeVec,L(:,3),'Color',[0 0 0]+0.6);
hold off
legend('\tau_{L,1}','\tau_{L,2}','\tau_{L,2}');
title('(e) fin ray joint torques due to inertia, gravity and friction');
grid on 
ylabel('Torque [Nm]')


subplot(6,1,6)
plot(timeVec,T(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,T(:,2),'Color',[0 0 0]+0.4);
plot(timeVec,T(:,3),'Color',[0 0 0]+0.6);
hold off
legend('\tau_{1}','\tau_{2}','\tau_{3}')
grid on
xlabel('t [s]'),ylabel('Torque [Nm]')
%axis([0 1.5 -0.5 0.5])
title('(f) total fin ray joint torques');

%sgtitle('h_{min}=0.116 : h_{max}=0.12 m : D=0.08 : A=30^{째} : f=2 Hz : 횠_{0}');
%figure(3)
%subplot(2,1,1)
%plot(timeVec,T(:,2),'-b');
%subplot(2,1,2)
%plot(timeVec,T(:,3),'-b');

figure(2)
subplot(4,1,1)
plot(timeVec,q(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,q(:,2),':','Color',[0 0 0]+0.4);
plot(timeVec,q(:,3),'--','Color',[0 0 0]+0.6);
hold off
legend('q_{1}','q_{2}','q_{3}')
ytixr = get(gca, 'YTick');
ytixlbl = regexp(sprintf('%.1f�\n',ytixr*180/pi), '\n', 'split');
set(gca, 'YTick', ytixr, 'YTickLabel',ytixlbl);
grid on 
ylabel('angle [\circ]')
title('(a) motion profile of the rays');

subplot(4,1,2)
plot(timeVec,ddq1,'Color',[0 0 0]);
hold on
plot(timeVec,ddq2,':','Color',[0 0 0]+0.4);
plot(timeVec,ddq3,'--','Color',[0 0 0]+0.4);
legend('q_{1}-q_{2}','q_{2}-q_{3}','q_{1}-q_{3}')
ytixr = get(gca, 'YTick');
ytixlbl = regexp(sprintf('%.1f�\n',ytixr*180/pi), '\n', 'split');
set(gca, 'YTick', ytixr, 'YTickLabel',ytixlbl);
grid on 
ylabel('angle [\circ]')
title('(b) Angle separation');

subplot(4,1,3)
plot(timeVec,F_N(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,F_N(:,2),'--','Color',[0 0 0]);
plot(timeVec,F_N(:,3),':','Color',[0 0 0]);
hold off
legend('F_{Sway}','F_{Heave}','F_{surge}')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
title('(c) hydrodynamic force components');

subplot(4,1,4)
plot(timeVec,T_phi(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,T_phi(:,2),'--','Color',[0 0 0]);
plot(timeVec,T_phi(:,3),':','Color',[0 0 0]);
hold off
legend('\tau_{\Phi,1}','\tau_{\Phi,2}','\tau_{\Phi,3}')
grid on 
ylabel('Torque [Nm]')
title('(d) fin ray joint torques due to hydrodynamic loads');
xlabel('t [s]'),ylabel('Torque [Nm]')

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
FN_N = NFN_func(NVec,nN,NMax,nT,timeVec,D,hmax);
%-------------------------------------------------------------------
Phi_A =Phi_A_FUN(nPhi,PhiVec,timeVec,nT,nAa,AaVec,D,hmax);
%-------------------------------------------------------------------
Phi_N =Phi_N_FUN(nPhi,PhiVec,timeVec,nT,NMax,D,hmax);
%-------------------------------------------------------------------
Phi_fr = Phi_F_FUN(nPhi,PhiVec,timeVec,nT,nFa,FaVec,D,hmax);
%-------------------------------------------------------------------
Phi_D = Phi_D_FUN(nPhi,PhiVec,timeVec,nT,nD,DVec,hmax);
%-------------------------------------------------------------------
Phi_H = Phi_H_FUN(nPhi,PhiVec,timeVec,nT,D,nH,HVec);
%-------------------------------------------------------------------
Phi_N_heave =Phi_N_heave(nPhi,PhiVec,timeVec,nT,NMax,D,hmax);
%-------------------------------------------------------------------
Phi_N_sway =Phi_N_Sway(nPhi,PhiVec,timeVec,nT,NMax,D,hmax);
%-------------------------------------------------------------------
figure (3)
subplot(1,2,1)
plot(freqVec,FN_freq(:,1),'Color',[0 0 0]);
hold on
plot(freqVec,FN_freq(:,2),':','Color',[0 0 0]);
plot(freqVec,FN_freq(:,3),'-.','Color',[0 0 0]);
hold off
legend('F_{Sway}','F_{Heave}','F_{surge}')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 fMax])
ylim([0 25])
grid on
ylabel('Force [N]')
xlabel('\it f [Hz]' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

subplot(1,2,2)
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
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

figure (4)
subplot(1,2,1)
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
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

subplot(1,2,2)
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
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

figure (5)
plot(NVec,FN_N(:,1),'Color',[0 0 0]);
hold on
plot(NVec,FN_N(:,2),':','Color',[0 0 0]);
plot(NVec,FN_N(:,3),'-.','Color',[0 0 0]);
hold off
legend('F_{heave}','F_{sway}','F_{surge}')
xticks([2 3 4 5 6])
xticklabels({'2','3','4','5','6'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it Number of Rays','FontName', 'Times' )

figure(6)
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
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl

figure(7)
plot(PhiVec,Phi_N(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_N(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_N(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_N(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Phi_N(:,5),'-o','Color',[0 0 0]);
legend('\it N=2','\it N=3','\it N=4','\it N=5','\it N=6');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f째\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl

figure (8)
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

figure (9)
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
legend('\it D=0.05m','\it D=0.1m','\it D=0.15m','\it D=0.2m');
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')

figure(10)
[maxy1,idx1] = max(Phi_H(:,1));
[maxy2,idx2] = max(Phi_H(:,2));
[maxy3,idx3] = max(Phi_H(:,3));
[maxy4,idx4] = max(Phi_H(:,4));
F1=Phi_H(:,1);
F2=Phi_H(:,2);
F3=Phi_H(:,3);
F4=Phi_H(:,4);
plot(PhiVec,Phi_H(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_H(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_H(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_H(:,4),'--','Color',[0 0 0]);
plot(PhiVec(idx1),F1(idx1),'*k');
plot(PhiVec(idx2),F2(idx2),'*k');
plot(PhiVec(idx3),F3(idx3),'*k');
plot(PhiVec(idx4),F4(idx4),'*k');
hold off
legend('\it H=0.05m','\it H=0.1m','\it H=0.15m','\it H=0.2m');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{surge}[N]')
xlabel('\it \phi [\circ]')

figure(11)
plot(PhiVec,Phi_N(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_N_heave(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_N_heave(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_N_heave(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Phi_N_heave(:,5),'-o','Color',[0 0 0]);
legend('\it N=2','\it N=3','\it N=4','\it N=5','\it N=6');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{heave}[N]')
xlabel('\it \phi [\circ]')

figure(12)
plot(PhiVec,Phi_N_sway(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Phi_N_sway(:,2),':','Color',[0 0 0]);
plot(PhiVec,Phi_N_sway(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Phi_N_sway(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Phi_N_sway(:,5),'-o','Color',[0 0 0]);
legend('\it N=2','\it N=3','\it N=4','\it N=5','\it N=6');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{heave}[N]')
xlabel('\it \phi [\circ]')