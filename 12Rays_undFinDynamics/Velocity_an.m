clear all
close all
%------------------------------------------------------------------------
global Cn rho A X hmin
%------------------------------------------------------------------------
%The membrane parameters
hmax = 0.175; %[m]
hmin = 0; %[m]
D = 0.1; %[m]
Cn=2.6;
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
tMax = 6;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
ntt=tMax/dt;
%------------------------------------------------------------------------
q= readmatrix('joints_q_dq.xlsx');
q1=q(:,2);
q2=q(:,3);
q3=q(:,4);
dq1=q(:,5);
dq2=q(:,6);
dq3=q(:,7);
 %------------------------------------------------------------------------
 for (indT = 1:nT)
    
    qq1 = q1(indT);
 
    qq2 = q2(indT);
    qq3 = q3(indT);

    dqq1=dq1(indT);
    
    dqq2=dq2(indT);
    dqq3=dq3(indT);
%overall hydrodynamic forces Fn=[Fx Fy Fz]
 FN=hydrodynamicforces_Func(D,hmax,[qq1 qq2],[dqq1 dqq2])+hydrodynamicforces_Func(D,hmax,[qq2 qq3],[dqq2 dqq3]);
 %FN=hydrodynamicforces_Func(q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(q_curr(5:6),dq_curr(5:6));
 F_N(indT, :)=FN;
 end
FX=F_N(:,3)+F_N(:,3);
FZ=F_N(:,1)-F_N(:,1);
FY=F_N(:,2)+F_N(:,2);
x=zeros(nT,1);
x(1,:)=x0;
vx=zeros(nT,1);
vx(1,:)=vx0;
fdx=zeros(nT,1);
fdx(1,:)=0.5*rho*vx0^2*Cd*pi*R^2;
ax=zeros(nT,1);
ax(1,:)=(F_N(1,3)+F_N(1,3)-fdx(1))/mx;
y=zeros(nT,1);
y(1,:)=y0;
vy=zeros(nT,1);
vy(1,:)=vy0;
ay=zeros(nT,1);
fdy=zeros(nT,1);
fdy(1,:)=0.5*rho*vy0^2*Cd2*2*pi*R*L;
ay(1,:)=(F_N(1,1)-F_N(1,1)-fdy(1))/my;

z=zeros(nT,1);
z(1,:)=z0;
vz=zeros(nT,1);
vz(1,:)=vz0;
fdz=zeros(nT,1);
fdz(1,:)=0.5*rho*vz0^2*Cd2*2*pi*R*L;
az=zeros(nT,1);
az(1,:)=(F_N(1,2)+F_N(1,2))/mz;




for n=2:(ntt+1)
    x(n,:)=x(n-1,:)+vx(n-1,:)*dt; % robot position
    vx(n,:)=vx(n-1,:)+ax(n-1,:)*dt;
    fdx(n,:)=0.5*rho*(vx(n-1,:))^2*Cd*pi*R^2;
    ax(n,:)=(F_N(n,3)+F_N(n,3)-fdx(n,:))/mx;
    
    y(n,:)=y(n-1,:)+vy(n-1,:)*dt; % robot position
    vy(n,:)=vy(n-1,:)+ay(n-1,:)*dt;
    fdy(n,:)=0.5*rho*(vy(n-1,:))^2*Cd2*2*pi*R*L;
    ay(n,:)=(F_N(n,1)-F_N(n,1)-fdy(n,:))/my;
    
    z(n,:)=z(n-1,:)+vz(n-1,:)*dt; % robot position
    vz(n,:)=vz(n-1,:)+az(n-1,:)*dt;
    fdz(n,:)=0.5*rho*(vz(n-1,:))^2*Cd2*2*pi*R*L;
    az(n,:)=(F_N(n,2)+F_N(n,2))/mz;
end

X=sqrt(x.^2+y.^2+z.^2);
V=sqrt(vx.^2+vy.^2);
A=sqrt(ax.^2+ay.^2+az.^2);

figure(1)
plot now
subplot(2,1,1)
plot(timeVec,q(:,2),'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,q(:,3),'Color',[0 0 0]+0.6);
plot(timeVec,q(:,4),':','Color',[0 0 0]);
hold off
legend('q_{1}','q_{2}','q_{3}')
ytixr = get(gca, 'YTick');
ytixlbl = regexp(sprintf('%.1f°\n',ytixr*180/pi), '\n', 'split');
set(gca, 'YTick', ytixr, 'YTickLabel',ytixlbl);
grid on 
ylabel('angle [\circ]')
title('(a) motion profile of the rays (for right and left fin)');
subplot(2,1,2)
plot(timeVec,FX,'linewidth',1.5,'Color',[0 0 0]);
hold on
plot(timeVec,FY,'-.','Color',[0 0 0]);
plot(timeVec,FZ,':','Color',[0 0 0]);
hold off
legend('F_{Surge}','F_{Heave}','F_{Sway}')
title('(b) hydrodynamic forces');
%ylim([-6 6])
%yticks([-10 -8 -4 -2 0 2 4 8 10])
%yticklabels({'-4','-2','0','2','4'})
%for figure1
%yticks([-5 -4 -2 0 2 4 5])
grid on
ylabel('Force [N]')
xlabel('time[s]')

figure(2)
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
%xlabel('time[s]')
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
title('(c) Robot accleartion in Surge (X), Sway (Y), Heave (Z)');
xlabel('time[s]')

