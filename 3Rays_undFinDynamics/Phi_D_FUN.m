function FA = Phi_D_FUN(nPhi,PhiVec,timeVec,nT,nD,DVec,hmax)
global N X A f
for indD = 1:nD
    D=DVec(indD);
for indPhi = 1:nPhi
    
    Phi = PhiVec(indPhi);
for indT = 1:nT
    t = timeVec(indT);
    
    
    for i = 1:N
   
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*t);
    end
    FN_Ph = hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3));
    %FN_Ph=abs(FN_Ph);
    PhiFN(indT, :)=FN_Ph;
end
F(indPhi,:)=max(PhiFN);
Fx=F(:,3);
end
FA(:,indD)=Fx;
end
end