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
N = 3;%number of rays 
A = deg2rad(30);%Amplitude
f = 1;%frequency
Phi = deg2rad(50);%the inter-ray phase shift
X= deg2rad(0);%angular offset
%------------------------------------------------------------------------
%T_1= readmatrix('phi60_6s_normal.xlsx');
T_normal= readmatrix('F_phi50_f2hz.xlsx');
%time
ncycle=3; %Number of loading cycles
tMax = ncycle/2;
dt = 0.005;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
ntt=tMax/dt;
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%hydrodynamic parameters
Cn=2.65;
rho=1000;%[kg/m3]
hmin=0.034;
hmax=0.174;
D=0.096;
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
 FN=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2));
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F_N(indT, :)=FN;

end
h = figure(1);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
ddq1=q(:,1)-q(:,2);
ddq2=q(:,2)-q(:,3);
figure(1)
subplot(3,1,1)
%plot(timeVec,q(:,1),'linewidth',1.5,'Color',[0 0 0]);
plot(timeVec,q(:,1),'linewidth',3);
hold on
%plot(timeVec,q(:,2),'linewidth',1.5,'Color',[0 0 0]+0.6);
plot(timeVec,q(:,2),'linewidth',3);
%plot(timeVec,ddq1,':','Color',[0 0 0]);
plot(timeVec,ddq1,':g','linewidth',3);
hold off
legend('q_{1}','q_{2}','q_{1}-q_{2}')
ylim([deg2rad(-30) deg2rad(30)])
yticks([deg2rad(-30) deg2rad(-20) deg2rad(-10) 0 deg2rad(10) deg2rad(20) deg2rad(30)])
yticklabels({'-30\circ','-20\circ','-10\circ','0','10\circ','20\circ','30\circ',})
grid on 
xlabel('t [s]');
ylabel('angle [\circ]')
%title('(a) motion profile of the rays');


subplot(3,1,2)
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
%title('(c) hydrodynamic force components');

subplot(3,1,3)
%plot(timeVec,T_phi(:,1),'Color',[0 0 0]);
plot(timeVec,T_phi(:,1),'linewidth',3);
hold on
%plot(timeVec,T_phi(:,2),'Color',[0 0 0]+0.6);
plot(timeVec,T_phi(:,2),'linewidth',3);
%plot(timeVec,T_phi(:,3),'Color',[0 0 0]+0.6);
hold off
legend('\tau_{\Phi,1}','\tau_{\Phi,2}','\tau_{\Phi,3}')
grid on
xlabel('t [s]');
ylabel('Torque [Nm]');
%title('(d) fin ray joint torques due to hydrodynamic loads');
drawnow
 % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
%--------------------------------------------------------------------
figure(3)
subplot(3,1,1)
plot(timeVec,F_N(:,2),'Color',[0 0 0]);
hold on
%plot(timeVec,T_coarse(:,3));
plot(timeVec,T_normal(:,3),'-.','Color',[0 0 0]);
hold off
legend('Analytical (\phi =50\circ)','CFD (\phi =50\circ)','CFD (\phi =40\circ)','CFD (\phi =30\circ)','CFD (\phi =20\circ)');
grid on
%xlim([1 6])
ylabel('Heave[N]')
%title('(b) Sway force');
xlabel('t [s]');
%-------------------------------------------------------
subplot(3,1,2)
plot(timeVec,F_N(:,1),'Color',[0 0 0]);
hold on
plot(timeVec,T_normal(:,2),'-.','Color',[0 0 0]);
hold off
legend('Analytical (\phi =50\circ)','CFD (\phi =50\circ)','CFD (\phi =40\circ)','CFD (\phi =30\circ)','CFD (\phi =20\circ)');
grid on
%xlim([1 6])
ylabel('Sway[N]')
xlabel('t [s]');
%title('(c) Heave force');
%-------------------------------------------------------
subplot(3,1,3)
plot(timeVec,F_N(:,3),'Color',[0 0 0]);
hold on
plot(timeVec,T_normal(:,4),'-.','Color',[0 0 0]);
%plot(timeVec,T_1(:,4),'Color',[0 0 0]+0.6);
hold off
legend('Analytical (\phi =50\circ)','CFD (\phi =50\circ)','CFD (\phi =40\circ)','CFD (\phi =30\circ)','CFD (\phi =20\circ)');
grid on
%xlim([1 6])
xlabel('t [s]');
ylabel('Surge[N]')
%title('(d) Surge force');