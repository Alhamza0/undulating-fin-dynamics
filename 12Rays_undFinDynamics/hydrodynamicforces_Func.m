function F=hydrodynamicforces_Func(D,hmax,q,dq)
%HYDRODYNAMICFORCES_FUNC Compute hydrodynamic forces on a single membrane panel.
%
%   F = HYDRODYNAMICFORCES_FUNC(D, hmax, q, dq) calculates the 3D 
%   hydrodynamic force vector by integrating the normal force distribution
%   over the membrane surface between two adjacent fin rays.
%
%   Reference: Thesis Equation (13) — Double integration of f_n over
%   the fin surface area using resistive force theory.
%
%   Inputs:
%     D    - Inter-ray joint spacing [m]
%     hmax - Maximum membrane height [m]
%     q    - [1x2] Joint angles of adjacent rays [rad]
%     dq   - [1x2] Joint angular velocities of adjacent rays [rad/s]
%
%   Outputs:
%     F    - [3x1] Force vector [F_sway; F_heave; F_surge] [N]
%
%   Global variables required: hmin, Cn, rho
%     hmin - Minimum membrane height [m]
%     Cn   - Normal drag coefficient [-]
%     rho  - Fluid density [kg/m^3]
%
%   Integration:
%     h: hmin → hmax (membrane height direction)
%     w: 0 → wmax(h) = sqrt(2*h^2*(1-cos(q1-q2)) + D^2) (width)
%
%   Note: The integrand expressions were generated using MATLAB's
%   Symbolic Math Toolbox from the analytical derivation.
%
%   See also: hydrodynamictorques_Func, DeltaFunc, compute_total_forces
%
global hmin Cn rho
wmax1 = @(h)sqrt(2.*h.^2*(1-cos(q(1)-q(2)))+D.^2);
f1_12 = @(h,w)(Cn.*D.^3.*h.^2.*rho.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*(sin(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*sin(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0;
f2_12 = @(h,w)Cn.*D.^3.*h.^2.*rho.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*(cos(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*cos(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0).*(-1.0./2.0);
f3_12 = @(h,w)Cn.*D.^2.*h.^3.*rho.*sin(q(1)-q(2)).*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0).*(-1.0./2.0);
Fx=quad2d(f1_12,hmin,hmax,0,wmax1);
Fy=quad2d(f2_12,hmin,hmax,0,wmax1);
Fz=quad2d(f3_12,hmin,hmax,0,wmax1);
F=[Fx; Fy; Fz];
end