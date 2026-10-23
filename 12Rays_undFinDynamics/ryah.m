clear all
close all
clc
syms y(t) x(t) df1 df2 Y a w V0 g lamda0 g alpha phi
%df1=0;
%df2=0;
Vy=diff(y);
Vx=diff(x);
dVx=diff(x,2);
dVy=diff(y,2);
V=sqrt(Vx^2+Vy^2);
Eq1=dVy==g-(1/320*a*alpha)*V^2*(pi*Vy*(df1)+Vx*(df2));
Eq2=dVx==-(1/320*a*alpha)*V^2*(pi*Vx*(df1)-Vy*(df2));
[VF,Subs] = odeToVectorField(Eq1, Eq2);
ftotal = matlabFunction(VF, 'Vars',{t,Y,a,w,g,alpha});

%F1_phi0=8*cos(0)^3-24*cos(0)+(6*a*w/V)*sin(0)^4;
F1_phi=8*cos(phi)^3-24*cos(phi)+(6*a*w/V)*sin(phi)^4;
df1=F1_phi;
%F2_phi0=32*sin(0)^3-+(3*a*w/V)*(4*0-sin(4*0));
F2_phi=32*sin(phi)^3-+(3*a*w/V)*(4*phi-sin(4*phi));
df2=F2_phi;
phi=atan(Vy/Vx)+acos(a-y/a);
df1 = @(t)df1;
df2= @(t)df2;
a=0.1;
alpha= 7.8;
g= 9.81;
w=200;
V0=100;
lamda0=10;
ic=[0;V0*sin(lamda0);0;V0*cos(lamda0)];
tspan = [0 5];
[T,Y] = ode45(@(t,Y) ftotal(t,Y,,a,w,g,alpha),tspan, ic);
plot(T, Y)
grid
%ySol = ode45( @(t,Y)f1total(t,y,1,1,1),interval,y0);
%odes = [ode1; ode2];
%cond1= x(0)==0;
%cond2= y(0)==0;
%cond3= Vx(0)==V0*cos(lamda0);
%cond4= Vy(0)==V0*sin(lamda0);
%conds = [cond2 cond4];
%tspan = [0 5];
%[t,y] = ode45(@(t,y) ode1, tspan,0,0);

%ySol(t) = dsolve(ode1,conds)
