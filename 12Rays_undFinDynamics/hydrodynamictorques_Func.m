function T=hydrodynamictorques_Func(D,hmax,q,dq)
%HYDRODYNAMICTORQUES_FUNC Compute joint torques from hydrodynamic loads.
%
%   T = HYDRODYNAMICTORQUES_FUNC(D, hmax, q, dq) calculates the torques
%   on each fin ray joint due to hydrodynamic pressure loads on the
%   membrane surface. Uses double integration (similar to force calculation).
%
%   Reference: Thesis Equation (15) — Joint torques due to hydrodynamic loads
%
%   Inputs:
%     D    - Inter-ray joint spacing [m]
%     hmax - Maximum membrane height [m]
%     q    - [1x3] Joint angles of 3 adjacent rays [rad]
%     dq   - [1x3] Joint angular velocities [rad/s]
%
%   Outputs:
%     T    - [3x1] Torque vector [T1; T2; T3] [N*m]
%            T1: torque on joint 1 (from membrane 1-2)
%            T2: torque on joint 2 (from membrane 1-2)
%            T3: torque on joint 3 (from membrane 2-3)
%
%   Global variables required: hmin, Cn, rho
%
%   Note: Currently hardcoded for 3 rays. For N-ray generalization,
%   this function would need to be extended.
%
%   See also: hydrodynamicforces_Func, DeltaFunc
%
global hmin Cn rho
wmax1 = @(h)sqrt(2.*h.^2*(1-cos(q(1)-q(2)))+D.^2);
wmax2 = @(h)sqrt(2.*h.^2*(1-cos(q(2)-q(3)))+D.^2);

Tx=@(h,w)Cn.*D.^3.*h.^3.*rho.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*cos(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*cos(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0).*(-1.0./2.0)-(Cn.*D.^3.*h.^3.*rho.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*sin(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(sin(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*sin(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0;
Ty=@(h,w)(Cn.*D.^3.*h.^3.*rho.*w.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*cos(q(2)).*(cos(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*cos(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0+(Cn.*D.^3.*h.^3.*rho.*w.*sign(-D.*h.*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0)).*sin(q(2)).*(sin(q(1)).*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0)-w.*sin(q(2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*(dq(1)-w.*(dq(1).*2.0-cos(q(1)-q(2)).*(dq(1)+dq(2))).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)+(w.^2.*(dq(1)+dq(2)).*(cos(q(1)-q(2))-1.0))./(h.^2.*(cos(q(1)-q(2))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(1)-q(2)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2)-1.0).*(cos(q(1)-q(2))-1.0).*1.0./sqrt(h.^2.*(cos(q(1)-q(2))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0;
Tz=@(h,w)(Cn.*D.^3.*h.^3.*rho.*w.*sign(-D.*h.*(dq(2)-w.*(dq(2).*2.0-cos(q(2)-q(3)).*(dq(2)+dq(3))).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)+(w.^2.*(dq(2)+dq(3)).*(cos(q(2)-q(3))-1.0))./(h.^2.*(cos(q(2)-q(3))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(2)-q(3)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0).*(cos(q(2)-q(3))-1.0).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*2.0)).*cos(q(3)).*(cos(q(2)).*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0)-w.*cos(q(3)).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*(dq(2)-w.*(dq(2).*2.0-cos(q(2)-q(3)).*(dq(2)+dq(3))).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)+(w.^2.*(dq(2)+dq(3)).*(cos(q(2)-q(3))-1.0))./(h.^2.*(cos(q(2)-q(3))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(2)-q(3)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0).*(cos(q(2)-q(3))-1.0).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0+(Cn.*D.^3.*h.^3.*rho.*w.*sign(-D.*h.*(dq(2)-w.*(dq(2).*2.0-cos(q(2)-q(3)).*(dq(2)+dq(3))).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)+(w.^2.*(dq(2)+dq(3)).*(cos(q(2)-q(3))-1.0))./(h.^2.*(cos(q(2)-q(3))-1.0).*2.0-D.^2)).*1.0./sqrt(h.^2.*sin(q(2)-q(3)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0).*(cos(q(2)-q(3))-1.0).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*2.0)).*sin(q(3)).*(sin(q(2)).*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0)-w.*sin(q(3)).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*(dq(2)-w.*(dq(2).*2.0-cos(q(2)-q(3)).*(dq(2)+dq(3))).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)+(w.^2.*(dq(2)+dq(3)).*(cos(q(2)-q(3))-1.0))./(h.^2.*(cos(q(2)-q(3))-1.0).*2.0-D.^2)).^2.*1.0./(h.^2.*sin(q(2)-q(3)).^2+D.^2-D.^2.*w.*(w.*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2)-1.0).*(cos(q(2)-q(3))-1.0).*1.0./sqrt(h.^2.*(cos(q(2)-q(3))-1.0).*-2.0+D.^2).*2.0).^(3.0./2.0))./2.0;

T1=quad2d(Tx,hmin,hmax,0,wmax1);
T2=quad2d(Ty,hmin,hmax,0,wmax1);
T3=quad2d(Tz,hmin,hmax,0,wmax2);
T=[T1; T2; T3];
end
