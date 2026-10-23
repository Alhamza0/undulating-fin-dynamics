clear all
close all
clc
syms y(t) x(t) Y a w  g alpha 
Vy=diff(y);%==tan(lamda)*diff(x);
Vx=diff(x);%==Vy/tan(lamda);
dVx=diff(x,2);
dVy=diff(y,2);
V=sqrt(diff(x)^2+diff(y)^2);
Eq1=dVy==g - (a*alpha*((2^(1/2)*(32*sin(acos(a - y/a) + atan(Vy/Vx))^3 - (3*a*w*(4*acos(a - y/a) - sin(4*acos(a - y/a) + 4*atan(Vy/Vx)) + 4*atan(Vy/Vx)))/V))/2 + (2^(1/2)*pi*(8*cos(acos(a - y/a) + atan(Vy/Vx))^3 - 24*cos(acos(a - y/a) + atan(Vy/Vx)) + (6*a*w*sin(acos(a - y/a) + atan(Vy/Vx))^4)/V))/2)*(tan(lamda)^2 + 1))/320;
Eq2=dVx==(a*alpha*((2^(1/2)*(32*sin(acos(a - y/a) + atan(Vy/Vx))^3 - (3*a*w*(4*acos(a - y/a) - sin(4*acos(a - y/a) + 4*atan(Vy/Vx)) + 4*atan(Vy/Vx)))/V))/2 - (2^(1/2)*pi*(8*cos(acos(a - y/a) + atan(Vy/Vx))^3 - 24*cos(acos(a - y/a) + atan(Vy/Vx)) + (6*a*w*sin(acos(a - y/a) + atan(Vy/Vx))^4)/V))/2)*(tan(lamda)^2 + 1))/320;
[VF,Subs] = odeToVectorField(Eq1, Eq2);
%[V] = odeToVectorField(diff(y, 2) == (Eq1);
ftotal = matlabFunction(VF, 'Vars',{t,Y,a,w,g,alpha});
%F1_phi0=8*cos(0)^3-24*cos(0)+(6*a*w/V)*sin(0)^4;
%F1_phi=8*cos(atan(Vy/Vx)+acos(a-y/a))^3-24*cos(atan(Vy/Vx)+acos(a-y/a))+(6*a*w/V)*sin(atan(Vy/Vx)+acos(a-y/a))^4;
%df1=F1_phi;
%F2_phi0=32*sin(0)^3-+(3*a*w/V)*(4*0-sin(4*0));
%F2_phi=32*sin(atan(Vy/Vx)+acos(a-y/a))^3-+(3*a*w/V)*(4*atan(Vy/Vx)+acos(a-y/a)-sin(4*atan(Vy/Vx)+acos(a-y/a)));
%df2=F2_phi;
%phi=atan(Vy/Vx)+acos(a-y/a);
%df1 = @(t)df1;
%df2= @(t)df2;
a=0.1;
alpha= 7.8;
g= 9.81;
w=200;
V0=100;
lamda0=10;
%lamda=atan(diff(y)/diff(x));
ic=[0;V0*sin(lamda0);0;V0*cos(lamda0)];
tspan = [0 0.1];
[T,Y] = ode45(@(t,Y) ftotal(t,Y,a,w,g,alpha),tspan, ic);
plot(T, Y)
%fplot(@(x)deval(T,x,1), [0, 0.1])
%grid
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
