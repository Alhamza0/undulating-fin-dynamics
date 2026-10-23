clear all
close all
clc
syms y(t) x(t) Y a w  g alpha
Vy=diff(y);
Vx=diff(x);
dVx=diff(x,2);
dVy=diff(y,2);
V=sqrt(Vx^2+Vy^2);
lamda=atan(diff(y)/diff(x));
Eq1=dVy==g - (V^2*a*alpha*(cos(lamda)*(32*sin(acos(a - y/a) + atan(Vy/Vx))^3 - (3*a*w*(4*acos(a - y/a) - sin(4*acos(a - y/a) + 4*atan(Vy/Vx)) + 4*atan(Vy/Vx)))/V) + pi*sin(lamda)*(8*cos(acos(a - y/a) + atan(Vy/Vx))^3 - 24*cos(acos(a - y/a) + atan(Vy/Vx)) + (6*a*w*sin(acos(a - y/a) + atan(Vy/Vx))^4)/V)))/320;
Eq2=dVx==-(V^2*a*alpha*(sin(lamda)*(32*sin(acos(a - y/a) + atan(Vy/Vx))^3 - (3*a*w*(4*acos(a - y/a) - sin(4*acos(a - y/a) + 4*atan(Vy/Vx)) + 4*atan(Vy/Vx)))/V) - pi*cos(lamda)*(8*cos(acos(a - y/a) + atan(Vy/Vx))^3 - 24*cos(acos(a - y/a) + atan(Vy/Vx)) + (6*a*w*sin(acos(a - y/a) + atan(Vy/Vx))^4)/V)))/320;
[V1] = odeToVectorField(Eq1);
