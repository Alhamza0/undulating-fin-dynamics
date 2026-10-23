%this code consist of two parts
%partA:
%In order to solve and simulates the actuation forces and torques, 
%an inverse dynamics approach is used to solve (8).
%partB;
%%Parametric analysis to see what the effect of 
%undulating parameters on hydrodynamic forces.

clear all
close all
%------------------------------------------------------------------------
global D l v Cn rho hmax hmin alpha Phi A f N X
%------------------------------------------------------------------------
%Traveling wave parameters(sinusoidal path)
N = 6;%number of rays 
A = deg2rad(30);%Amplitude
f = 2;%frequency
Phi = deg2rad(70);%the inter-ray phase shift 
X= deg2rad(0);%angular offset 

%------------------------------------------------------------------------
%time
tMax = 3;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%------------------------------------------------------------------------
%oscillation frequency,
fMax = 4;
df = 0.01;
freqVec = [0:df:fMax].';
nF = length(freqVec);

%-----------------------------------
%oscillation amplitude,
Amax = deg2rad(40);
dA = 0.01;
AmVec = [0:dA:Amax].';
nA = length(AmVec);
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
%-----------------------------------
 %number of rays,
NMax = 6;
dN = 2;
%NVec = [2:dN:NMax].';
NVec = [2;3;4;5;6].';
nN = length(NVec);
%-----------------------------------
%The membrane parameters
hmax = 0.16; %[m]
hmin = 0; %[m]
D = 0.08; %[m]
l = 0.08;
v = 0.5;
EK = 360; %[N/m]
alpha=1;
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.8;
rho=1022;%[kg/m3]
%------------------------------------------------------------------------
%%Fin ray parameters
g = 9.81; %gravitational accelaration [m/s^2]
mi = 0.01;m0 = mi;m1 = mi;m2 = mi;m3=mi;m4=mi;m5=mi;m6=mi;%mass of the link
li = 0.12;l0 = li;l1 = li;l2 = li;l3=li;l4 = li;l5 = li;l6=li;%length of the ray
Izi = 72e-6;Iz0 = Izi;Iz1 = Izi;Iz2 = Izi;Iz3 = Izi;Izi;Iz4 = Izi;Iz5 = Izi;Iz6 = Izi;%moment of inertia of the link 
bi = 0.02;b0 = bi;b1 = bi;b2 = bi;b3 = bi;b4 = bi;b5 = bi;b6 = bi;%the viscous friction coefficient
kci = 0.001;kc0 = kci;kc1 = kci;kc2 = kci;kc3 = kci;kc4 = kci;kc5 = kci;kc6 = kci;%the Coulomb friction
g1 = (l1/2)*m1*g;g2 = (l2/2)*m2*g;g3 = (l3/2)*m3*g;g4 = (l2/2)*m2*g;g6 = (l3/2)*m3*g;g5 = (l3/2)*m3*g;
%------------------------------------------------------------------------
%the mass matrix of the fin system
M0 = m1 + m2+m3;
M1 = m1*(l1/2)^2 + Iz1;
M2 = m2*(l2/2)^2 + Iz2;
M3 = m3*(l3/2)^2 + Iz3;
M4 = m4*(l1/2)^2 + Iz1;
M5 = m5*(l2/2)^2 + Iz2;
M6 = m6*(l3/2)^2 + Iz3;
M = diag([M1,M2,M3,M4,M5,M6]);
%------------------------------------------------------------------------
B = diag([b1,b2,b3,b4,b5,b6]);
Kc = diag([kc1,kc2,kc3,kc4,kc5,kc6]);
%------------------------------------------------------------------------

%%partA:

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
    G = [g1*sin(q(indT,1));g2*sin(q(indT,2)) ;g3*sin(q(indT,3)) ;g4*sin(q(indT,4));g5*sin(q(indT,5)) ;g6*sin(q(indT,6))];
 %------------------------------------------------------------------------
 %hydrodynamic torques on each joint
 hT= hydrodynamictorques_Func(q_curr,dq_curr);
    T_phi(indT, :) = hT;
 %------------------------------------------------------------------------
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3));
 FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F_N(indT, :)=FN;
%------------------------------------------------------------------------
%the force on each joint, due to gravity, inertia and fraction.
L(indT,:) = M*(ddq(indT,:).') + B*(dq(indT,:).') + Kc*sign(dq(indT,:).') + G ;
%------------------------------------------------------------------------
%the generalized force vector for the joints due to the elastic...

%------------------------------------------------------------------------
%total torque on each joint
    T(indT,:) = L(indT,:); %+ hT';
end

ddq1=q(:,1)-q(:,2);
ddq2=q(:,2)-q(:,3);
ddq3=q(:,1)-q(:,3);
figure(1)
subplot(2,1,1)
plot(timeVec,q(:,1),'-k',timeVec,q(:,2),'-b',timeVec,q(:,3),'-r','LineWidth',2.0);
legend('q_1','q_2','q_3')
%axis([0 1 -0.6 0.6])
%yticks([-0.5 0 0.5])
grid on 
ylabel('angle [rad]')
title('(a) motion profile of the rays');

%subplot(6,1,2)
%plot(timeVec,ddq1,'-k',timeVec,ddq2,'-b',timeVec,ddq3,'-r','LineWidth',2.0);
%legend('q_1-q_2','q_2-q_3','q_1-q_3')
%grid on 
%ylabel('angle [rad]')
%title('(b) Angle separation');



subplot(2,1,2)
plot(timeVec,F_N(:,1),'-k',timeVec,F_N(:,2),'-b',timeVec,F_N(:,3),'r','LineWidth',2.0);
legend('F_{heave}','F_{sway}','F_{surge}')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
title('(c) hydrodynamic force components');
xlabel('t [s]')

