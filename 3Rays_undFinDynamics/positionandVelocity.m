clear all; close all; clc;
%this code is for calculate the posistion , velocity and accelration of the
%robot body 

global hmin hmax D Cn rho
%Traveling wave parameters(sinusoidal path)
N = 3;%number of rays 
Ar = deg2rad(30);%Amplitude
Al = deg2rad(30);%Amplitude
fr = 1;%frequency
fl = 1;%frequency
Phir = deg2rad(-20);%the inter-ray phase shift 
Phil = deg2rad(20);%the inter-ray phase shift 
Xr= deg2rad(0);%angular offset 
Xl= deg2rad(0);%angular offset 

%The membrane parameters
hmax = 0.175; %[m]
hmin = 0; %[m]
D = 0.1; %[m]
Cn=2.65;
rho=1022;
Cd=1.1;%for sphere
Cd2=0.70; %for the cylinder
%physicsal parameters 
R=0.1; 
L=0.22;
m=2.5; %kg mass of the robot body 
max=2/3*rho*pi*R^3;
may=rho*pi*R^2*L;
fbz=-rho*9.81*pi*R^2*L;
mx=m+max;
my=m+may;
mz=m+may;
x0=0; %robot initial posistion(m)
vx0=0; %robot initial velocity (m/s)
y0=0; %robot initial posistion(m)
vy0=0; %robot initial velocity (m/s)
z0=0; %robot initial posistion(m)
vz0=0; %robot initial velocity (m/s)

%time step
tMax = 10;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
ntt=tMax/dt;

F_Nl=zeros(nT,3);
F_Nr=zeros(nT,3);
for (indT = 1:nT)
    
    tt = timeVec(indT);
 %------------------------------------------------------------------------
    for (i = 1:N)
        %To define the generalized coordinates of the fin system
        % q=[q1 q2 q3]
        
        qr(i) = Ar*sin(2*pi*fr*tt - (i-1)*Phir)+Xr;
        %send current desired position (q(indT,i) for motor i
        
        dqr(i) = 2*Ar*fr*pi*cos(Phir*(i - 1) - 2*pi*fr*tt);   
    end
    Tr= hydrodynamictorques_Func(D,hmax,qr,dqr);
    T_r(indT, :) = Tr;
    
    FNr=hydrodynamicforces_Func(D,hmax,qr(1:2),dqr(1:2))+hydrodynamicforces_Func(D,hmax,qr(2:3),dqr(2:3));
    %FNr_6=hydrodynamicforces_Func(D,hmax,qr(1:2),dqr(1:2))+hydrodynamicforces_Func(D,hmax,qr(2:3),dqr(2:3))+hydrodynamicforces_Func(D,hmax,qr(3:4),dqr(3:4))+hydrodynamicforces_Func(D,hmax,qr(4:5),dqr(4:5))+hydrodynamicforces_Func(D,hmax,qr(5:6),dqr(5:6));
    %FNr_8=hydrodynamicforces_Func(D,hmax,qr(1:2),dqr(1:2))+hydrodynamicforces_Func(D,hmax,qr(2:3),dqr(2:3))+hydrodynamicforces_Func(D,hmax,qr(3:4),dqr(3:4))+hydrodynamicforces_Func(D,hmax,qr(4:5),dqr(4:5))+hydrodynamicforces_Func(D,hmax,qr(5:6),dqr(5:6))+hydrodynamicforces_Func(D,hmax,qr(6:7),dqr(6:7))+hydrodynamicforces_Func(D,hmax,qr(7:8),dqr(7:8));

    F_Nr(indT, :)=FNr;
    %F_Nr(indT, :)=FNr_8;
end

for (indT = 1:nT)
    
    tt = timeVec(indT);
 %------------------------------------------------------------------------
    for (i = 1:N)
        %To define the generalized coordinates of the fin system
        % q=[q1 q2 q3]
        
        ql(i) = Al*sin(2*pi*fl*tt - (i-1)*Phil)+Xl;
        %send current desired position (q(indT,i) for motor i
        
        dql(i) = 2*Al*fl*pi*cos(Phil*(i - 1) - 2*pi*fl*tt);   
    end
    Tl= hydrodynamictorques_Func(D,hmax,ql,dql);
    T_l(indT, :) = Tl;
    FNl=hydrodynamicforces_Func(D,hmax,ql(1:2),dql(1:2))+hydrodynamicforces_Func(D,hmax,ql(2:3),dql(2:3));
    %FNl_6=hydrodynamicforces_Func(D,hmax,ql(1:2),dql(1:2))+hydrodynamicforces_Func(D,hmax,ql(2:3),dql(2:3))+hydrodynamicforces_Func(D,hmax,ql(3:4),dql(3:4))+hydrodynamicforces_Func(D,hmax,ql(4:5),dql(4:5))+hydrodynamicforces_Func(D,hmax,ql(5:6),dql(5:6));
    %FNl_8=hydrodynamicforces_Func(D,hmax,ql(1:2),dql(1:2))+hydrodynamicforces_Func(D,hmax,ql(2:3),dql(2:3))+hydrodynamicforces_Func(D,hmax,ql(3:4),dql(3:4))+hydrodynamicforces_Func(D,hmax,ql(4:5),dql(4:5))+hydrodynamicforces_Func(D,hmax,ql(5:6),dql(5:6))+hydrodynamicforces_Func(D,hmax,ql(6:7),dql(6:7))+hydrodynamicforces_Func(D,hmax,ql(7:8),dql(7:8));

    F_Nl(indT, :)=FNl;
    %F_Nl(indT, :)=FNl_8;
end
TL=T_l(:,1)+T_l(:,2)+T_l(:,3);
TR=T_r(:,1)+T_r(:,2)+T_r(:,3);
T=TL-TR;
FX=F_Nr(:,3)+F_Nl(:,3);
FZ=F_Nr(:,1)-F_Nl(:,1);
FY=F_Nr(:,2)+F_Nl(:,2);
x=zeros(nT,1);
x(1,:)=x0;
vx=zeros(nT,1);
vx(1,:)=vx0;
fdx=zeros(nT,1);
fdx(1,:)=0.5*rho*vx0^2*Cd*pi*R^2;
ax=zeros(nT,1);
ax(1,:)=(F_Nr(1,3)+F_Nl(1,3)-fdx(1))/mx;
y=zeros(nT,1);
y(1,:)=y0;
vy=zeros(nT,1);
vy(1,:)=vy0;
ay=zeros(nT,1);
fdy=zeros(nT,1);
fdy(1,:)=0.5*rho*vy0^2*Cd2*2*pi*R*L;
ay(1,:)=(F_Nr(1,1)-F_Nl(1,1)-fdy(1))/my;

z=zeros(nT,1);
z(1,:)=z0;
vz=zeros(nT,1);
vz(1,:)=vz0;
fdz=zeros(nT,1);
fdz(1,:)=0.5*rho*vz0^2*Cd2*2*pi*R*L;
az=zeros(nT,1);
az(1,:)=(F_Nr(1,2)+F_Nl(1,2))/mz;




for n=2:(ntt+1)
    x(n,:)=x(n-1,:)+vx(n-1,:)*dt; % robot position
    vx(n,:)=vx(n-1,:)+ax(n-1,:)*dt;
    fdx(n,:)=0.5*rho*(vx(n-1,:))^2*Cd*pi*R^2;
    ax(n,:)=(F_Nr(n,3)+F_Nl(n,3)-fdx(n,:))/mx;
    
    y(n,:)=y(n-1,:)+vy(n-1,:)*dt; % robot position
    vy(n,:)=vy(n-1,:)+ay(n-1,:)*dt;
    fdy(n,:)=0.5*rho*(vy(n-1,:))^2*Cd2*2*pi*R*L;
    ay(n,:)=(F_Nr(n,1)-F_Nl(n,1)-fdy(n,:))/my;
    
    z(n,:)=z(n-1,:)+vz(n-1,:)*dt; % robot position
    vz(n,:)=vz(n-1,:)+az(n-1,:)*dt;
    fdz(n,:)=0.5*rho*(vz(n-1,:))^2*Cd2*2*pi*R*L;
    az(n,:)=(F_Nr(n,2)+F_Nl(n,2))/mz;
end

X=sqrt(x.^2+y.^2+z.^2);
V=sqrt(vx.^2+vy.^2);
A=sqrt(ax.^2+ay.^2+az.^2);

figure(2)
subplot (2,1,1)
plot(timeVec,FX,'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,FY,'-.','Color',[0 0 0]);
plot(timeVec,FZ,':','Color',[0 0 0]);
legend('F_{Surge}','F_{Heave}','F_{Sway}')
grid on
ylabel('Force [N]')
xlabel('time[s]')

subplot (2,1,2)
plot(timeVec,T(:,1),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,TR,'-.','Color',[0 0 0]);
plot(timeVec,TL,':','Color',[0 0 0]);
hold off
legend('\tau_{\Phi,T}','\tau_{\Phi,R}','\tau_{\Phi,L}')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Torque [Nm]')
xlabel('time[s]')

figure(3)
subplot(3,1,1)
plot(timeVec,x,'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,y,'-.','Color',[0 0 0]);
plot(timeVec,z,':','Color',[0 0 0]);
hold off
legend('x','y','z')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Position [m]')
title('(a) Robot Position in Surge(X), Sway(Y), Heave(Z)');

subplot(3,1,2)
plot(timeVec,vx,'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,vy,'-.','Color',[0 0 0]);
plot(timeVec,vz,':','Color',[0 0 0]);
hold off
legend('vx','vy','vz')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('velocity [m/s]')
title('(b) Robot velocity in Surge (X), Sway (Y), Heave (Z)');
xlabel('time[s]')
subplot(3,1,3)

plot(timeVec,ax,'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,ay,'-.','Color',[0 0 0]);
plot(timeVec,az,':','Color',[0 0 0]);
hold off
legend('ax','ay','az')
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('accleartion [m/s^{2}]')
title('(b) Robot accleartion in Surge (X), Sway (Y), Heave (Z)');
xlabel('time[s]')

