function result = parametric_sweep(params, sweep_param, sweep_values, varargin)
%PARAMETRIC_SWEEP Generic parametric study for undulating fin forces.
%
%   result = PARAMETRIC_SWEEP(params, sweep_param, sweep_values)
%   performs a parametric sweep over one variable, computing max
%   hydrodynamic forces at each value.
%
%   result = PARAMETRIC_SWEEP(params, sweep_param, sweep_values, ...
%                             sweep_param2, sweep_values2)
%   performs a 2D parametric sweep (e.g., phase shift × amplitude).
%
%   This is a unified replacement for the separate functions:
%     PhiFN.m      → parametric_sweep(p, 'Phi', PhiVec)
%     AmpFN.m      → parametric_sweep(p, 'A', AmVec)
%     freqFN.m     → parametric_sweep(p, 'f', freqVec)
%     XFN.m        → parametric_sweep(p, 'X', XVec)
%     Phi_A_FUN.m  → parametric_sweep(p, 'Phi', PhiVec, 'A', AaVec)
%     Phi_D_FUN.m  → parametric_sweep(p, 'Phi', PhiVec, 'D', DVec)
%     Phi_F_FUN.m  → parametric_sweep(p, 'Phi', PhiVec, 'f', FaVec)
%     Phi_H_FUN.m  → parametric_sweep(p, 'Phi', PhiVec, 'hmax', HVec)
%
%   Inputs:
%     params       - struct from fin_params()
%     sweep_param  - string: parameter name to sweep ('Phi','A','f','X','N','D','hmax')
%     sweep_values - vector of values for sweep_param
%     (optional) sweep_param2, sweep_values2 for 2D sweep
%
%   Outputs (1D sweep):
%     result.F_max       - [n x 3] max [Sway, Heave, Surge] for each value
%     result.sweep_vals  - sweep values used
%     result.sweep_name  - name of swept parameter
%
%   Outputs (2D sweep):
%     result.F_surge_2D  - [nPhi x nParam2] max surge for each combination
%     result.sweep_vals1 - first sweep values (outer parameter)
%     result.sweep_vals2 - second sweep values (inner parameter)
%
%   Example:
%     p = fin_params('N', 12);
%     % 1D: sweep phase shift
%     r = parametric_sweep(p, 'Phi', linspace(0, deg2rad(80), 100));
%     plot(r.sweep_vals, r.F_max);
%     
%     % 2D: phase shift x amplitude
%     r2 = parametric_sweep(p, 'Phi', PhiVec, 'A', AaVec);
%     plot(PhiVec, r2.F_surge_2D);
%
%   See also: fin_params, compute_total_forces, compute_kinematics

    is2D = (nargin >= 5);
    
    nSweep = length(sweep_values);
    timeVec = params.timeVec;
    nT = params.nT;
    
    if ~is2D
        %% 1D Parametric Sweep
        F_max = zeros(nSweep, 3);
        
        for idx = 1:nSweep
            p_local = params;
            p_local.(sweep_param) = sweep_values(idx);
            
            F_time = zeros(nT, 3);
            for indT = 1:nT
                [q, dq, ~] = compute_kinematics(p_local, timeVec(indT));
                F_time(indT, :) = compute_total_forces(p_local, q, dq).';
            end
            F_max(idx, :) = max(F_time);
        end
        
        result.F_max      = F_max;
        result.sweep_vals = sweep_values;
        result.sweep_name = sweep_param;
        
    else
        %% 2D Parametric Sweep
        sweep_param2  = varargin{1};
        sweep_values2 = varargin{2};
        nSweep2       = length(sweep_values2);
        
        F_surge_2D = zeros(nSweep, nSweep2);
        
        for idx2 = 1:nSweep2
            for idx1 = 1:nSweep
                p_local = params;
                p_local.(sweep_param)  = sweep_values(idx1);
                p_local.(sweep_param2) = sweep_values2(idx2);
                
                F_time = zeros(nT, 3);
                for indT = 1:nT
                    [q, dq, ~] = compute_kinematics(p_local, timeVec(indT));
                    F_time(indT, :) = compute_total_forces(p_local, q, dq).';
                end
                F_surge_2D(idx1, idx2) = max(F_time(:, 3)); % Surge
            end
        end
        
        result.F_surge_2D  = F_surge_2D;
        result.sweep_vals1 = sweep_values;
        result.sweep_vals2 = sweep_values2;
        result.sweep_name1 = sweep_param;
        result.sweep_name2 = sweep_param2;
    end
end
