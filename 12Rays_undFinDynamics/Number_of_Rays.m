% this function is used to find the relation between rays phase-shift and
% number of the rays N and it effect on Sway component.
clear all
close all
%----------------------------------------------------------------------
global Cn rho hmin
% Traveling wave parameters(sinusoidal path)
N = 3;%number of rays 
A = deg2rad(30);%Amplitude
NMax = 6;
f = 2;%frequency
Phi = deg2rad(40);%the inter-ray phase shift 
X= deg2rad(0);%angular offset 
%----------------------------------------------------------------------
%time
tMax = 2;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%----------------------------------------------------------------------
%inter-ray phase shift
PhiMax = deg2rad(80);
dPhi = deg2rad(0.8);
PhiVec = [0:dPhi:PhiMax].';
nPhi = length(PhiVec);
%----------------------------------------------------------------------
%The membrane parameters
hmax = 0.175; %[m]
Hmax = 0.2;
dH = 0.05;
HVec = [dH:dH:Hmax].';
nH = length(HVec);
hmin = 0; %[m]
%----------------------------------------------------------------------
D = 0.10; %[m]
Dmax = 0.2;
dD = 0.05;
DVec = [dD:dD:Dmax].';
nD = length(DVec);
%----------------------------------------------------------------------
l = 0.08;
v = 0.5;
EK = 360; %[N/m]
alpha=1;
%----------------------------------------------------------------------

%hydrodynamic parameters
Cn=2.65;
rho=1022;%[kg/m3]
%------------------------------------------------------------------------

for indPhi = 1:nPhi
    Phi = PhiVec(indPhi);
for indT = 1:nT
    t = timeVec(indT);
    
    
    for i = 1:NMax
   
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*t);
    end
FN1=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2));
FN_1(indT, :)=FN1;
FN2=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3));
FN_2(indT, :)=FN2;
FN3=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(D,hmax,q_curr(3:4),dq_curr(3:4));
FN_3(indT, :)=(FN3);
FN4=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(D,hmax,q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(D,hmax,q_curr(4:5),dq_curr(4:5));
FN_4(indT, :)=FN4;
FN5=hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3))+hydrodynamicforces_Func(D,hmax,q_curr(3:4),dq_curr(3:4))+hydrodynamicforces_Func(D,hmax,q_curr(4:5),dq_curr(4:5))+hydrodynamicforces_Func(D,hmax,q_curr(5:6),dq_curr(5:6));
FN_5(indT, :)=FN5;
end
F1(indPhi,:)=max(FN_1);
F2(indPhi,:)=max(FN_2);
F3(indPhi,:)=max(FN_3);
F4(indPhi,:)=max(FN_4);
F5(indPhi,:)=max(FN_5);
F40(51,:)=[max(FN_1) max(FN_2) max(FN_3) max(FN_4) max(FN_5)];

Fheave=[F1(:,1) F2(:,1) F3(:,1) F4(:,1) F5(:,1)];
Fsway=[F1(:,2) F2(:,2) F3(:,2) F4(:,2) F5(:,2)];
Fsurge=[F1(:,3) F2(:,3) F3(:,3) F4(:,3) F5(:,3)];
end

figure(1)
plot(PhiVec,Fheave(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Fheave(:,2),':','Color',[0 0 0]);
plot(PhiVec,Fheave(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Fheave(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Fheave(:,5),'-o','Color',[0 0 0]);
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
plot(PhiVec,Fsway(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Fsway(:,2),':','Color',[0 0 0]);
plot(PhiVec,Fsway(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Fsway(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Fsway(:,5),'-o','Color',[0 0 0]);
legend('\it N=2','\it N=3','\it N=4','\it N=5','\it N=6');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{Sway}[N]')
xlabel('\it \phi [\circ]');

figure(3)
plot(PhiVec,Fsurge(:,1),'Color',[0 0 0]);
hold on
plot(PhiVec,Fsurge(:,2),':','Color',[0 0 0]);
plot(PhiVec,Fsurge(:,3),'-.','Color',[0 0 0]);
plot(PhiVec,Fsurge(:,4),'--','Color',[0 0 0]);
plot(PhiVec,Fsurge(:,5),'-o','Color',[0 0 0]);
legend('\it N=2','\it N=3','\it N=4','\it N=5','\it N=6');
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhiMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('F_{Surge}[N]')
xlabel('\it \phi [\circ]');

figure (4)
%phi=40
plot(NVec,Fheave(51,:),'Color',[0 0 0]);
hold on
plot(NVec,Fsway(51,:),':','Color',[0 0 0]);
plot(NVec,Fsurge(51,:),'-.','Color',[0 0 0]);
hold off
legend('F_{Sway}','F_{Heave}','F_{Surge}')
xticks([2 3 4 5 6])
xticklabels({'2','3','4','5','6'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it Number of Rays','FontName', 'Times' )

figure (5)
%phi=60
plot(NVec,Fheave(76,:),'Color',[0 0 0]);
hold on
plot(NVec,Fsway(76,:),':','Color',[0 0 0]);
plot(NVec,Fsurge(76,:),'-.','Color',[0 0 0]);
hold off
legend('F_{Sway}','F_{Heave}','F_{surge}')
xticks([2 3 4 5 6])
xticklabels({'2','3','4','5','6'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it Number of Rays','FontName', 'Times' )