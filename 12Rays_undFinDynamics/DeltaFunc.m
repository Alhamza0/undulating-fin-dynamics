function Delta = DeltaFunc(D,x,q)
%DELTAFUNC Compute elastic membrane deflection between two rays.
%
%   Delta = DELTAFUNC(D, x, q) computes the deformation of the elastic
%   membrane connecting two adjacent fin rays, based on their angular
%   separation and the membrane geometry.
%
%   The deflection depends on:
%     - The inter-ray distance D
%     - The evaluation height x along the membrane
%     - The angular difference (q1 - q2) between adjacent rays
%
%   Inputs:
%     D - Inter-ray joint spacing [m]
%     x - Evaluation height along the membrane [m]
%     q - [1x2] Joint angles of adjacent rays [rad]
%
%   Outputs:
%     Delta - Membrane deflection [m^3] (used in elastic torque calculation)
%
%   Global variables required: l (characteristic length), v (Poisson's ratio)
%
%   Note: Returns 0 when rays are parallel (q1 = q2) to avoid NaN.
%
%   See also: hydrodynamicforces_Func, hydrodynamictorques_Func
%
global l v;

alpha = D/l;
p = sqrt(2*(1-cos(q(1)-q(2))));
kappa = x/D*p;

Delta = (x^3*alpha/(3*D)) * (1 + v*(8-3*alpha*sqrt(kappa^2+1))/4)- ...
    (D^2/(2*p^3))*(sqrt(kappa^4+ kappa^2) - asinh(kappa)) * ...
    (1 + v*(4 + alpha^2)/4);
Delta(isnan(Delta))=0;
end 
