function B = freqFN(nF,freqVec,timeVec,nT,D,hmax)
global N  A Phi X
for (indF = 1:nF)
   
    f = freqVec(indF);
    
   for (indT = 1:nT)
    t = timeVec(indT);
    for i = 1:N
   
        q_curr(i) = A*sin(2*pi*f*t - (i-1)*Phi)+X;
        %send current desired position (q(indT,i) for motor i
        
        dq_curr(i) = 2*A*f*pi*cos(Phi*(i - 1) - 2*pi*f*t);
    end
    FN = hydrodynamicforces_Func(D,hmax,q_curr(1:2),dq_curr(1:2))+hydrodynamicforces_Func(D,hmax,q_curr(2:3),dq_curr(2:3));
    %FN_Ph=abs(FN_Ph);
    F_N(indT, :)=FN;
end
B(indF,:)=max(F_N);
end