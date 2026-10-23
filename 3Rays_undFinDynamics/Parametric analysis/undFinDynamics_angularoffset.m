%this code is for Implementing simulation results for a fixed two-ray fin
clear all;close all;

global D l v

g = 9.81; %gravitational accelaration [m/s^2]
N = 3;
f = 2;
t = 1;
A = deg2rad(30);
Ph=deg2rad(40);

%%Geometric parameters
hmax = 0.12; %[m]
hmin = 0; %[m]
D = 0.08; %[m]
l = 0.08;
v = 0.5;

%inter-ray phase shift
XMax = deg2rad(90);
dX = 0.01;
XVec = [deg2rad(13):dX:XMax].';
nX = length(XVec);





for (indX = 1:nX)
    
    X = XVec(indX);
    
    for (i = 1:N)
        
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Ph)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Ph*(1 - i) + 2*pi*f*t);
        %ddq_curr(i) = -4*A*f^2*pi^2*sin(Ph*(1 - i) + 2*pi*f*t);
    end
    
    q(indX,:) = (q_curr);
    dq(indX,:) = (dq_curr);
    FN = hydrodynamicForce(q(indX,[1:N-1]),dq(indX,[1:N-1]))+hydrodynamicForce(q(indX,[N-1:N]),dq(indX,[N-1:N]));
    FN=abs(FN);
    F_N(indX, :)=FN;
    
    
end

figure(1)
plot(XVec,F_N(:,1),'-k',XVec,F_N(:,2),'-b',XVec,F_N(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it A [\circ]','FontName', 'Times' )
xtixr = get(gca, 'XTick');
xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

