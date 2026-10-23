%% MAIN_ANALYSIS - Undulating Fin Hydrodynamic Analysis (Refactored)
%
% This is the refactored version of undFinDynamics.m and twelve_Rays.m
% demonstrating the clean architecture with fin_params struct.
%
% Part A: Time-domain force and torque analysis
% Part B: Parametric studies (Phi, A, f, X, N)
%
% Original code by: Hamzah (Master's Thesis)
% Refactored: 2026
%
% See also: fin_params, compute_kinematics, compute_total_forces,
%           parametric_sweep

clear; close all; clc;

%% ========== CONFIGURATION ==========
% Change parameters here — no globals needed
p = fin_params( ...
    'N',   3, ...              % Number of rays
    'A',   deg2rad(30), ...    % Amplitude
    'f',   2, ...              % Frequency [Hz]
    'Phi', deg2rad(40), ...    % Phase shift
    'X',   deg2rad(0), ...     % Angular offset
    'D',   0.10, ...           % Joint spacing [m]
    'hmax', 0.175, ...         % Membrane height [m]
    'tMax', 2, ...             % Simulation time [s]
    'dt',  0.01 ...            % Time step [s]
);

%% ========== PART A: TIME-DOMAIN ANALYSIS ==========
fprintf('Part A: Computing time-domain forces for %d-ray fin...\n', p.N);

% Preallocate output arrays
F_N    = zeros(p.nT, 3);    % [Sway, Heave, Surge]
T_phi  = zeros(p.nT, 3);    % Hydrodynamic torques
q_all  = zeros(p.nT, p.N);  % Joint angles
dq_all = zeros(p.nT, p.N);  % Joint velocities

% Time loop
for indT = 1:p.nT
    t = p.timeVec(indT);
    
    % Kinematics
    [q_row, dq_row, ~] = compute_kinematics(p, t);
    q_all(indT, :)  = q_row;
    dq_all(indT, :) = dq_row;
    
    % Total hydrodynamic forces (sum over all membrane panels)
    F_N(indT, :) = compute_total_forces(p, q_row, dq_row).';
    
    % Joint torques (uses original 3-ray function)
    if p.N >= 3
        T_phi(indT, :) = hydrodynamictorques_Func(p.D, p.hmax, ...
                          q_row(1:3), dq_row(1:3)).';
    end
end

%% Plot Part A Results
figure('Name', 'Time-Domain Analysis', 'Position', [100 100 800 700]);

subplot(3,1,1)
hold on; grid on;
for i = 1:p.N
    plot(p.timeVec, rad2deg(q_all(:,i)), 'LineWidth', 1.5);
end
ylabel('Angle [°]');
title(sprintf('(a) Ray motion profiles (N=%d, A=%.0f°, f=%.0fHz, \\Phi=%.0f°)', ...
    p.N, rad2deg(p.A), p.f, rad2deg(p.Phi)));
legend(arrayfun(@(i) sprintf('q_{%d}', i), 1:p.N, 'UniformOutput', false));

subplot(3,1,2)
plot(p.timeVec, F_N(:,1), 'k-', 'LineWidth', 1.5); hold on; grid on;
plot(p.timeVec, F_N(:,2), 'k--', 'LineWidth', 1.5);
plot(p.timeVec, F_N(:,3), 'k:', 'LineWidth', 1.5);
ylabel('Force [N]');
title('(b) Hydrodynamic force components');
legend('F_{Sway}', 'F_{Heave}', 'F_{Surge}');

subplot(3,1,3)
if p.N >= 3
    plot(p.timeVec, T_phi(:,1), 'k-', 'LineWidth', 1.5); hold on; grid on;
    plot(p.timeVec, T_phi(:,2), 'k--', 'LineWidth', 1.5);
    plot(p.timeVec, T_phi(:,3), 'k:', 'LineWidth', 1.5);
    ylabel('Torque [Nm]');
    title('(c) Joint torques due to hydrodynamic loads');
    legend('\tau_{\Phi,1}', '\tau_{\Phi,2}', '\tau_{\Phi,3}');
end
xlabel('Time [s]');

%% ========== PART B: PARAMETRIC STUDIES ==========
fprintf('Part B: Running parametric sweeps...\n');

% Define sweep ranges
PhiVec  = linspace(0, deg2rad(80), 101).';
freqVec = linspace(0, 4, 101).';
AmVec   = linspace(0, deg2rad(40), 101).';
XVec    = linspace(0, deg2rad(90), 101).';

% 1D Sweeps (uncomment to run — computationally intensive)
fprintf('  Sweeping phase shift (Phi)...\n');
r_Phi = parametric_sweep(p, 'Phi', PhiVec);

fprintf('  Sweeping frequency (f)...\n');
r_freq = parametric_sweep(p, 'f', freqVec);

fprintf('  Sweeping amplitude (A)...\n');
r_Amp = parametric_sweep(p, 'A', AmVec);

fprintf('  Sweeping angular offset (X)...\n');
r_X = parametric_sweep(p, 'X', XVec);

%% Plot Parametric Results
figure('Name', 'Parametric Studies', 'Position', [100 100 900 500]);

subplot(2,2,1)
plot(rad2deg(PhiVec), r_Phi.F_max(:,1), 'k-', ...
     rad2deg(PhiVec), r_Phi.F_max(:,2), 'k--', ...
     rad2deg(PhiVec), r_Phi.F_max(:,3), 'k:', 'LineWidth', 1.5);
xlabel('\Phi [°]'); ylabel('F_{max} [N]'); grid on;
legend('Sway', 'Heave', 'Surge'); title('(a) Phase shift effect');

subplot(2,2,2)
plot(freqVec, r_freq.F_max(:,1), 'k-', ...
     freqVec, r_freq.F_max(:,2), 'k--', ...
     freqVec, r_freq.F_max(:,3), 'k:', 'LineWidth', 1.5);
xlabel('f [Hz]'); ylabel('F_{max} [N]'); grid on;
legend('Sway', 'Heave', 'Surge'); title('(b) Frequency effect');

subplot(2,2,3)
plot(rad2deg(AmVec), r_Amp.F_max(:,1), 'k-', ...
     rad2deg(AmVec), r_Amp.F_max(:,2), 'k--', ...
     rad2deg(AmVec), r_Amp.F_max(:,3), 'k:', 'LineWidth', 1.5);
xlabel('A [°]'); ylabel('F_{max} [N]'); grid on;
legend('Sway', 'Heave', 'Surge'); title('(c) Amplitude effect');

subplot(2,2,4)
plot(rad2deg(XVec), r_X.F_max(:,1), 'k-', ...
     rad2deg(XVec), r_X.F_max(:,2), 'k--', ...
     rad2deg(XVec), r_X.F_max(:,3), 'k:', 'LineWidth', 1.5);
xlabel('\chi [°]'); ylabel('F_{max} [N]'); grid on;
legend('Sway', 'Heave', 'Surge'); title('(d) Angular offset effect');

sgtitle(sprintf('Parametric Analysis: %d-ray fin', p.N));

fprintf('Done.\n');
