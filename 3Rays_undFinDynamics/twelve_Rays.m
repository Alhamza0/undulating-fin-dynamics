%this code consist of two parts partA: In order to solve and simulates the
%actuation forces and torques, an inverse dynamics approach is used to
%solve (8). partB;
%Parametric analysis to see what the effect of undulating parameters on
%hydrodynamic forces.
%% partA
clear all
close all
%------------------------------------------------------------------------
global Cn rho Phi A f N X hmin
%------------------------------------------------------------------------
% Traveling wave parameters(sinusoidal path)
N = 12;%number of rays 
A = deg2rad(30);%Amplitude
f = 1;%frequency
Phi = deg2rad(40);%the inter-ray phase shift
X= deg2rad(0);%angular offset
%------------------------------------------------------------------------
%time
ncycle=3; %Number of loading cycles
tMax = ncycle/f;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
ntt=tMax/dt;
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1000;%[kg/m3]
hmin=0;
hmax=0.01;%m
D=0.002;
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
 %hT= hydrodynamictorques_Func(D,hmax,q_curr,dq_curr);
%    T_phi(indT, :) = hT;
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

end


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