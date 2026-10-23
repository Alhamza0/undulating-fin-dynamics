function FN= NFN_func(NVec,nN,NMax,nT,timeVec,D,hmax)
global   A X f Phi 
for indN = 1:nN
    N=NVec(indN);
for (indT = 1:nT)
    t = timeVec(indT);
    
for (i = 1:NMax)  
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(1 - i) + 2*pi*f*t);
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

%FX=[FN1(1);FN2(1);FN3(1);FN4(1);FN5(1)];
%FY=[FN1(2);FN2(2);FN3(2);FN4(2);FN5(2)];
%FZ=[FN1(3);FN2(3);FN3(3);FN4(3);FN5(3)];
%FN_N=[FX, FY, FZ];
%FN_N=abs(FN_N);
%FN(indT, :)=FN_N;
end
F1(indN,:)=max(FN_1);
F2(indN,:)=max(FN_2);
F3(indN,:)=max(FN_3);
F4(indN,:)=max(FN_4);
F5(indN,:)=max(FN_5);
FN=[F1(1,:);F2(1,:);F3(1,:);F4(1,:);F5(1,:)];
end