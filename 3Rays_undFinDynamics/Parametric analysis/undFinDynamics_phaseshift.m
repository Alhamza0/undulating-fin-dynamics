%this code is for Implementing simulation results for a fixed two-ray fin
clear all;close all;

global D l v

g = 9.81; %gravitational accelaration [m/s^2]
N = 3;
f = 2;
t = 0;
A = deg2rad(30);

%%Geometric parameters
hmax = 0.12; %[m]
hmin = 0; %[m]
D = 0.08; %[m]
l = 0.08;
v = 0.5;

%inter-ray phase shift
PhMax = deg2rad(80);
dPh = 0.001;
PhVec = [0:dPh:PhMax].';
nPh = length(PhVec);





for (indPh = 1:nPh)
    
    Ph = PhVec(indPh);
    
    for (i = 1:N)
        
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Ph);
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Ph*(1 - i) + 2*pi*f*t);
        %ddq_curr(i) = -4*A*f^2*pi^2*sin(Ph*(1 - i) + 2*pi*f*t);
    end
    
    q(indPh,:) = (q_curr);
    dq(indPh,:) = (dq_curr);
    FN = hydrodynamicForce(q(indPh,[1:N-1]),dq(indPh,[1:N-1]))+hydrodynamicForce(q(indPh,[N-1:N]),dq(indPh,[N-1:N]));
    FN=abs(FN);
    F_N(indPh, :)=FN;
    
    
end

figure(1)
plot(PhVec,F_N(:,1),'-k',PhVec,F_N(:,2),'-b',PhVec,F_N(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it A [\circ]','FontName', 'Times' )
xtixr = get(gca, 'XTick');
xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

