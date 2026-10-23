function Fheave = Phi_N_heave(nPhi,PhiVec,timeVec,nT,NMax,D,hmax)
% this function is used to find the relation between rays phase-shift and
% number of the rays N and it effect on Sway component.
global X A f
for indPhi = 1:nPhi
    Phi = PhiVec(indPhi);
for indT = 1:nT
    t = timeVec(indT);
    
    
    for i = 1:NMax
   
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*t);
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
end
F1(indPhi,:)=max(FN_1);
F2(indPhi,:)=max(FN_2);
F3(indPhi,:)=max(FN_3);
F4(indPhi,:)=max(FN_4);
F5(indPhi,:)=max(FN_5);
Fheave=[F1(:,2) F2(:,2) F3(:,2) F4(:,2) F5(:,2)];
end