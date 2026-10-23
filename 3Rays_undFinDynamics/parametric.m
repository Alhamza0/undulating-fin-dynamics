clear all;
close all;

global D l v N Cn rho hmin hmax Ph A t f X

N = 3;
t = 0;
A = deg2rad(30);
X=deg2rad(0);
f=2;
Ph=deg2rad(60);
%-----------------------------------
%time
tMax = 1;
dt = 0.01;
timeVec = [0:dt:tMax].';
nT = length(timeVec);
%-----------------------------------
%%Geometric parameters
hmax = 0.165; %[m]
hmin = 0; %[m]
D = 0.07; %[m]
l = 0.08;
v = 0.5;
Cn=2.8;
rho=1022;

%-----------------------------------
%oscillation frequency,
fMax = 4;
df = 0.01;
freqVec = [0:df:fMax].';
nF = length(freqVec);

%-----------------------------------
%oscillation amplitude,
Amax = deg2rad(40);
dA = 0.1;
AmVec = [0:dA:Amax].';
nA = length(AmVec);
%-----------------------------------
%inter-ray phase shift
PhMax = deg2rad(80);
dPh = 0.1;
PhVec = [0:dPh:PhMax].';
nPh = length(PhVec);
%-----------------------------------
%Angular offset
XMax = deg2rad(90);
XMin=deg2rad(0);
dX = 0.1;
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
PhFN=PhiFN(nPh,PhVec,timeVec,nT);
%-------------------------------------------------------------------
FN_freq = freqFN(nF,freqVec,timeVec,nT);
%-------------------------------------------------------------------
N_Amp = AmpFN(nA,AmVec,timeVec,nT);
%-------------------------------------------------------------------
XFN = XFN(nX,XVec,timeVec,nT);
%-------------------------------------------------------------------
FN_N = NFN_func(NVec,nN,NMax,nT,timeVec);
%-------------------------------------------------------------------
figure (1)
subplot(2,3,1)
plot(freqVec,FN_freq(:,1),'-k',freqVec,FN_freq(:,2),'-b',freqVec,FN_freq(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 fMax])
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it f [Hz] ','FontName', 'Times' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

subplot(2,3,2)
plot(AmVec,N_Amp(:,1),'-k',AmVec,N_Amp(:,2),'-b',AmVec,N_Amp(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 Amax])
xticks([0 deg2rad(10) deg2rad(20) deg2rad(30) deg2rad(40)])
xticklabels({'0','10','20','30','40'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it A [\circ]','FontName', 'Times' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)


subplot(2,3,3)
[maxy,idx] = max(PhFN(:,3));
FX=PhFN(:,3);
plot(PhVec,PhFN(:,1),'-k',PhVec,PhFN(:,2),'-b',PhVec,PhFN(:,3),'r',PhVec(idx),FX(idx),'pr','LineWidth',2.0);
legend('\it FZ','\it FY','\it FX', sprintf(' \\Phi= %0.1f',rad2deg(PhVec(idx))) );
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([0 PhMax])
xticks([0 deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it \Phi [\circ]','FontName', 'Times' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)

subplot(2,3,4)
plot(XVec,XFN(:,1),'-k',XVec,XFN(:,2),'-b',XVec,XFN(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
%for figure1
%yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
xlim([XMin XMax])
xticks([deg2rad(0) deg2rad(20) deg2rad(40) deg2rad(60) deg2rad(80)])
xticklabels({'0','20','40','60','80'})
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it X [\circ]','FontName', 'Times' )
%xtixr = get(gca, 'XTick');
%xtixlbl = regexp(sprintf('%.1f°\n',xtixr*180/pi), '\n', 'split');
%set(gca, 'XTick', xtixr, 'XTickLabel',xtixlbl)


subplot(2,3,5)
plot(NVec,FN_N(:,1),'-k',NVec,FN_N(:,2),'-b',NVec,FN_N(:,3),'r','LineWidth',2.0);
legend('FZ','FY','FX')
grid on
ylabel('Force [N]','FontName', 'Times')
xlabel('\it Number of Rays','FontName', 'Times' )



