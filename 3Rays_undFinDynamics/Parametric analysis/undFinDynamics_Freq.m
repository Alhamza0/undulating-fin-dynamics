%this code is for Implementing simulation results for a fixed two-ray fin
clear all;
close all;

global D l v

g = 9.81; %gravitational accelaration [m/s^2]
N = 3;
A = deg2rad(30);
t = 0;
phi0 = deg2rad(40);
%time
%tMax = 1;
%dt = 6;
%timeVec = [0:dt:tMax].';
%nT = length(timeVec);

%oscillation frequency,

fMax = 4;
df = 0.001;
freqVec = [0:df:fMax].';
nF = length(freqVec);

%oscillation amplitude,
%Aax = deg2rad(40);
%dA = 0.001;
%AmVec = [0:dA:AMax].';
%nA = length(AmVec);

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


for (indF = 1:nF)
    
    f = freqVec(indF);
    
    for (i = 1:N)
        
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*phi0);
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(phi0*(1 - i) + 2*pi*f*t);
        %ddq_curr(i) = -4*A*f^2*pi^2*sin(phi0*(1 - i) + 2*pi*f*t);
    end
    
    q(indF,:) = (q_curr);
    dq(indF,:) = (dq_curr);
    FN = hydrodynamicForce(q(indF,[1:N-1]),dq(indF,[1:N-1]))+hydrodynamicForce(q(indF,[N-1:N]),dq(indF,[N-1:N]));
    FN=abs(FN);
    F_N(indF, :)=FN;
    
end

figure(1)
plot(freqVec,F_N(:,1),'-k',freqVec,F_N(:,2),'-b',freqVec,F_N(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX','FontName', 'Times')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it f [Hz]','FontName', 'Times' )

