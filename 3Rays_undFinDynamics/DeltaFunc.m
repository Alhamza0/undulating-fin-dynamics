function Delta = DeltaFunc(D,x,q)

global l v;

alpha = D/l;
p = sqrt(2*(1-cos(q(1)-q(2))));
kappa = x/D*p;

Delta = (x^3*alpha/(3*D)) * (1 + v*(8-3*alpha*sqrt(kappa^2+1))/4)- ...
    (D^2/(2*p^3))*(sqrt(kappa^4+ kappa^2) - asinh(kappa)) * ...
    (1 + v*(4 + alpha^2)/4);
Delta(isnan(Delta))=0;
end 
