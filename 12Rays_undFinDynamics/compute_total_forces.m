function F_total = compute_total_forces(params, q, dq)
%COMPUTE_TOTAL_FORCES Sum hydrodynamic forces over all adjacent ray pairs.
%
%   F_total = COMPUTE_TOTAL_FORCES(params, q, dq) computes the total
%   hydrodynamic force on the fin by summing contributions from each
%   membrane panel (between adjacent rays).
%
%   For N rays, there are (N-1) membrane panels. The total force is:
%     F_total = sum_{i=1}^{N-1} F_hydro(q_i, q_{i+1}, dq_i, dq_{i+1})
%
%   This replaces the hardcoded summation of 11 calls that was repeated
%   in PhiFN.m, AmpFN.m, freqFN.m, XFN.m, etc.
%
%   Inputs:
%     params - struct from fin_params() (needs D, hmax, hmin, Cn, rho)
%     q      - joint angles [1 x N] [rad] at current time step
%     dq     - joint velocities [1 x N] [rad/s] at current time step
%
%   Outputs:
%     F_total - [3 x 1] = [F_sway; F_heave; F_surge] [N]
%
%   See also: hydrodynamicforces_Func, fin_params, compute_kinematics

    N = length(q);
    F_total = zeros(3, 1);
    
    for i = 1:(N-1)
        F_panel = hydrodynamicforces_Func(params.D, params.hmax, ...
                                          q(i:i+1), dq(i:i+1));
        F_total = F_total + F_panel;
    end
end
