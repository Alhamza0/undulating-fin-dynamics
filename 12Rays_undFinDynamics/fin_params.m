function params = fin_params(varargin)
%FIN_PARAMS Create a parameter structure for undulating fin simulation.
%
%   params = FIN_PARAMS() returns a struct with default parameters.
%   params = FIN_PARAMS('Name', Value, ...) overrides specific defaults.
%
%   This replaces the use of global variables throughout the codebase.
%   All parameters are grouped logically with clear documentation.
%
%   Example:
%       p = fin_params('N', 6, 'A', deg2rad(20), 'f', 1.5);
%       F = compute_total_forces(p, q, dq);
%
%   Parameter Groups:
%     Wave:         N, A, f, Phi, X (fin ray oscillation)
%     Geometry:     D, hmax, hmin, l (fin dimensions)
%     Fluid:        Cn, rho (hydrodynamic coefficients)
%     Material:     v, EK, alpha (membrane properties)
%     Mechanical:   mi, li, Izi, bi, kci, g (ray mechanical properties)
%
%   References:
%     - Thesis Eq. (1): q_i(t) = A*sin(2*pi*f*t - (i-1)*Phi) + X
%     - Thesis Eq. (13): Hydrodynamic force integration
%     - Thesis Eq. (15): Joint torque calculation
%
%   See also: compute_total_forces, compute_joint_torques,
%             compute_membrane_deflection

    %% Wave parameters (sinusoidal fin ray oscillation)
    params.N   = 3;              % Number of fin rays [-]
    params.A   = deg2rad(30);    % Oscillation amplitude [rad]
    params.f   = 2;              % Oscillation frequency [Hz]
    params.Phi = deg2rad(40);    % Inter-ray phase shift [rad]
    params.X   = deg2rad(0);     % Angular offset (bias angle) [rad]
    
    %% Fin geometry
    params.D    = 0.10;          % Joint spacing (inter-ray distance) [m]
    params.hmax = 0.175;         % Maximum membrane height [m]
    params.hmin = 0;             % Minimum membrane height [m]
    params.l    = 0.08;          % Characteristic length for DeltaFunc [m]
    
    %% Hydrodynamic coefficients
    params.Cn   = 2.65;          % Normal drag coefficient [-]
    params.rho  = 1022;          % Fluid density [kg/m^3]
    
    %% Membrane material properties
    params.v     = 0.5;          % Poisson's ratio [-]
    params.EK    = 360;          % Membrane stiffness [N/m]
    params.alpha = 1;            % Geometric ratio [-]
    
    %% Fin ray mechanical properties
    params.mi  = 0.01;           % Ray mass [kg]
    params.li  = 0.12;           % Ray length [m]
    params.Izi = 72e-6;          % Ray moment of inertia [kg*m^2]
    params.bi  = 0.02;           % Viscous friction coefficient [N*m*s]
    params.kci = 0.001;          % Coulomb friction coefficient [N*m]
    params.g   = 9.81;           % Gravitational acceleration [m/s^2]
    
    %% Time parameters
    params.tMax = 2;             % Maximum simulation time [s]
    params.dt   = 0.01;          % Time step [s]
    
    %% Override defaults with user-supplied values
    for k = 1:2:length(varargin)
        if isfield(params, varargin{k})
            params.(varargin{k}) = varargin{k+1};
        else
            error('fin_params:unknownParam', ...
                'Unknown parameter: %s', varargin{k});
        end
    end
    
    %% Derived quantities (computed from base parameters)
    params.timeVec = (0 : params.dt : params.tMax).';
    params.nT      = length(params.timeVec);
end
