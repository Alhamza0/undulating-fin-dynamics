clear all
close all
syms h w1 w2 hmin hmax Cn q1 q2 q3 dq1 dq2 dq3 rho D
c1=cos(q1);
c2=cos(q2);
c3=cos(q3);
s1=sin(q1);
s2=sin(q2);
s3=sin(q3);
cq12=q1-q2;
cq23=q2-q3;
wmax12=sqrt(2*h^2*(1-cos(q1-q2))+D^2);
wmax23=sqrt(2*h^2*(1-cos(q2-q3))+D^2);
w1=w1/wmax12;
w2=w2/wmax23;
s12= sqrt(2*D^2*w1*(w1-1)*(1-cos(cq12))+h^2*sin(cq12)^2+D^2);
s23= sqrt(2*D^2*w2*(w2-1)*(1-cos(cq23))+h^2*sin(cq23)^2+D^2);
en12 = (D/s12)*[s1*(w1-1)-s2*w1; -c1*(w1-1)+c2*w1; -sin(cq12)*h/D ];
en23 = (D/s23)*[s2*(w2-1)-s3*w2; -c2*(w2-1)+c3*w2; -sin(cq23)*h/D ];
vn12 = D*h/s12*(dq1+(1-cos(cq12))*(dq1+dq2)*w1^2+(cos(cq12)*(dq1+dq2)-2*dq1)*w1);
vn23 = D*h/s23*(dq2+(1-cos(cq23))*(dq2+dq3)*w2^2+(cos(cq23)*(dq2+dq3)-2*dq2)*w2);
J12= [0 0 1; -s1*h*(1-w1) c1*h*(1-w1) 0; -s2*h*w1 c2*h*w1 0 ];
J23= [0 0 1; -s2*h*(1-w2) c2*h*(1-w2) 0; -s3*h*w2 c3*h*w2 0 ];
fn12=-0.5*rho*Cn*en12*vn12^2*sign(vn12);
fn23=-0.5*rho*Cn*en23*vn23^2*sign(vn23);
phi=J12*(-fn12);
phi22=J23*(-fn23);
wmax1=@(h) sqrt(2.*h.^2*(1-cos(q1-q2))+D.^2)
wmax2=@(h) sqrt(2.*h.^2*(1-cos(q2-q3))+D.^2)

f1_12=matlabFunction(fn12(1))
f2_12=matlabFunction(fn12(2))
f3_12=matlabFunction(fn12(3))



f1_23=matlabFunction(fn23(1))
f2_23=matlabFunction(fn23(2))
f3_23=matlabFunction(fn23(3))

phi1=matlabFunction(phi(2))
phi2=matlabFunction(phi(3))
phi3=matlabFunction(phi22(3))