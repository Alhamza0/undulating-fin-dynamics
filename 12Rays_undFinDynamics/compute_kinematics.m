function [q, dq, ddq] = compute_kinematics(params, t)
%COMPUTE_KINEMATICS Compute fin ray angular position, velocity, acceleration.
%
%   [q, dq, ddq] = COMPUTE_KINEMATICS(params, t) computes the sinusoidal
%   kinematic profiles for all N fin rays at time instant t.
%
%   Reference: Thesis Eq. (1)
%     q_i(t)   = A * sin(2*pi*f*t - (i-1)*Phi) + X
%     dq_i(t)  = 2*pi*A*f * cos((i-1)*Phi - 2*pi*f*t)
%     ddq_i(t) = 4*pi^2*A*f^2 * sin((i-1)*Phi - 2*pi*f*t)
%
%   Inputs:
%     params - struct from fin_params() containing N, A, f, Phi, X
%     t      - time instant [s] (scalar or column vector)
%
%   Outputs:
%     q   - joint angles   [length(t) x N] [rad]
%     dq  - joint velocities [length(t) x N] [rad/s]
%     ddq - joint accelerations [length(t) x N] [rad/s^2]
%
%   See also: fin_params, compute_total_forces

    N   = params.N;
    A   = params.A;
    f   = params.f;
    Phi = params.Phi;
    X   = params.X;
    
    nT = length(t);
    q   = zeros(nT, N);
    dq  = zeros(nT, N);
    ddq = zeros(nT, N);
    
    for i = 1:N
        phase = (i-1) * Phi;
        q(:,i)   = A * sin(2*pi*f*t - phase) + X;
        dq(:,i)  = 2*pi*A*f * cos(phase - 2*pi*f*t);
        ddq(:,i) = 4*pi^2*A*f^2 * sin(phase - 2*pi*f*t);
    end
end
