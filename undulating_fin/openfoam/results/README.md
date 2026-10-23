# CFD force histories

One CSV per case: `time_s, Fx_heave_N, Fy_sway_N, Fz_surge_N`, the
total force on the fin patch (pressure + viscous) at every time step.

Only the total is kept; OpenFOAM's `force.dat` also carries the
pressure/viscous split, which is ten times the size and unused by the
comparison.

Load with `openfoam.compare.load_openfoam_forces`.

| case | configuration | notes |
|---|---|---|
| `N2_phi0` | N=2, Φ=0°, A=30°, f=1 Hz | 3 cycles, 191 742 cells; not fully settled (9.5 % cycle-to-cycle drift) |
