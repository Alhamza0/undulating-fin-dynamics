%this code is for Implementing simulation results for a fixed two-ray fin
clear all;close all;

global D l v

g = 9.81; %gravitational accelaration [m/s^2]
N = 3;
f = 2;
t = 0;
phi0 = deg2rad(40);
%time
%tMax = 1;
%dt = 0.001;
%timeVec = [0:dt:tMax].';
%nT = length(timeVec);


%oscillation amplitude,
Amax = deg2rad(40);
dA = 0.001;
AmVec = [0:dA:Amax].';
nA = length(AmVec);

%inter-ray phase shift
%Phax = deg2rad(90);
%dPh = 0.001;
%PhVec = [0:dPh:PhMax].';
%nPh = length(PhVec);

%q = zeros(nT,N);

%%Geometric parameters
hmax = 0.12; %[m]
hmin = 0; %[m]
D = 0.08; %[m]
l = 0.08;
v = 0.5;


for (indA = 1:nA)
    
    A = AmVec(indA);
    
    for (i = 1:N)
        
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*phi0);
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(phi0*(1 - i) + 2*pi*f*t);
        %ddq_curr(i) = -4*A*f^2*pi^2*sin(phi0*(1 - i) + 2*pi*f*t);
    end
    
    q(indA,:) = (q_curr);
    dq(indA,:) = (dq_curr);
    FN = hydrodynamicForce(q(indA,[1:N-1]),dq(indA,[1:N-1]))+hydrodynamicForce(q(indA,[N-1:N]),dq(indA,[N-1:N]));
    FN=abs(FN);
    F_N(indA, :)=FN;
    
    
end

figure(1)
plot(AmVec,F_N(:,1),'-k',AmVec,F_N(:,2),'-b',AmVec,F_N(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it A [\circ]','FontName', 'Times' )
xtixr = get(gca, 'XTick');
xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

