# Undulating Fin Dynamics

Analytical and CFD-validated hydrodynamic models of bio-inspired undulating
fins for underwater robots. This repository archives the M.Sc. thesis
*Kinematic and Dynamic Modeling of Biomimetic Undulating Fins for Underwater
Robots* (Alhamzah Ihsan Ali, Bahçeşehir University, 2021; supervisor
Assoc. Prof. Dr. Berke Gür), including the original MATLAB models and COMSOL
FSI study notes. It also contains a modern Python port (`undulating_fin/`) of
the analytical model (fin kinematics, resistive-force hydrodynamics, membrane
deflection, full ray dynamics and a dual-fin robot), together with a
generator for runnable OpenFOAM cases used to check the analytical forces
against independent CFD.

**Related publication:** B. Gür, A. I. Ali, B. Akbulut, "Propelling URSULA:
bioinspired fins for robotic squids", in *Maritime Autonomous Vehicles*, IET,
2025, pp. 105–126. DOI: [10.1049/SBRA565E_ch5](https://doi.org/10.1049/SBRA565E_ch5)

**Author:** Alhamzah Ihsan Ali ·
[ORCID 0000-0001-5211-5588](https://orcid.org/0000-0001-5211-5588) ·
[Google Scholar](https://scholar.google.com/citations?user=-L6WpSEAAAAJ)

License: [MIT](LICENSE) · How to cite: [CITATION.cff](CITATION.cff)

---

# Undulating Fin Hydrodynamics - Master's Thesis Project

## نظرة عامة | Overview

مشروع أطروحة الماجستير: **النمذجة الحركية والديناميكية للزعانف المتموجة المستوحاة بيولوجياً للروبوتات المائية**

Master's thesis project: **Kinematic and Dynamic Modeling of Biomimetic
Undulating Fins for Underwater Robots**

> Ali, A.I. (2021), Bahcesehir University, Dept. of Mechatronics Engineering.
> Supervisor: Assoc. Prof. Dr. Berke Gür

This project develops an analytical hydrodynamic model for a bio-inspired undulating fin mechanism with N fin rays connected by elastic membranes, validated against CFD (COMSOL) simulations.

---

## Project Structure

```
undulating-fin-dynamics/
├── README.md                          ← This file
│
├── undulating_fin/                    ← ★ PYTHON PACKAGE (modern port)
│   ├── pyproject.toml
│   ├── README.md                     ← Python-specific docs
│   ├── src/undulating_fin/           ← Source code
│   │   ├── params.py                 ← Parameter dataclasses
│   │   ├── kinematics.py             ← Traveling wave (Eq. 1-3)
│   │   ├── hydrodynamics.py          ← Force integration (Eq. 12-16)
│   │   ├── membrane.py               ← Elastic deflection (Eq. 8-10)
│   │   ├── dynamics.py               ← Full dynamics (Eq. 4-11)
│   │   ├── robot.py                  ← Dual-fin robot (Ch. 6)
│   │   ├── parametric.py             ← Parametric sweeps (Ch. 5)
│   │   └── plotting.py              ← Visualization
│   └── examples/                     ← Ready-to-run scripts
│
├── 12Rays_undFinDynamics/             ← MATLAB analytical models (original)
│   ├── MATLAB_README.md              ← MATLAB code guide
│   ├── hydrodynamicforces_Func.m     ← Core force function
│   ├── undFinDynamics.m              ← Main 3-ray analysis
│   └── ...
│
└── 12-Comsol/                         ← COMSOL FSI simulations
    ├── COMSOL_README.md
    └── ...
```

---

## Mathematical Model Summary

### Fin Kinematics
Each fin ray follows a sinusoidal oscillation:
$$q_i(t) = A \sin(2\pi f t - (i-1)\Phi) + \chi$$

where:
- $A$ = amplitude [rad]
- $f$ = frequency [Hz]  
- $\Phi$ = inter-ray phase shift [rad]
- $\chi$ = angular offset [rad]
- $i$ = ray index (1 to N)

### Hydrodynamic Forces
Normal force on membrane element using resistive force theory:
$$F_N = \int_0^{h_{max}} \int_0^{w_{max}} f_n(h,w) \, dw \, dh$$

Force components: **Sway** (lateral), **Heave** (vertical), **Surge** (forward thrust)

### Dynamic Equation
$$\tau = M\ddot{q} + B\dot{q} + K_c \text{sign}(\dot{q}) + G(q) + \tau_\Phi + \tau_\Gamma$$

---

## Key Parameters

| Parameter | Symbol | Default Value | Unit |
|-----------|--------|---------------|------|
| Number of rays | N | 3 (or 12) | - |
| Amplitude | A | 30 | deg |
| Frequency | f | 1-2 | Hz |
| Phase shift | Φ | 40 | deg |
| Membrane height | h_max | 0.175 | m |
| Joint spacing | D | 0.10 | m |
| Normal force coeff. | C_n | 2.65 | - |
| Fluid density | ρ | 1022 | kg/m³ |

---

## How to Run

### MATLAB (Analytical Model)
1. Open MATLAB and navigate to `12Rays_undFinDynamics/`
2. Run `undFinDynamics.m` for the 3-ray fin analysis
3. Run `twelve_Rays.m` for the 12-ray fin analysis
4. Run `positionandVelocity.m` for robot body dynamics

### COMSOL (CFD Validation)
1. Open COMSOL Multiphysics
2. Load model files from `12-Comsol/`
3. Import parameters from the corresponding `.txt` files
4. Run FSI transient study

---

## Key Results
- Optimal phase shift for maximum surge: ~40° (3-ray), ~60° (12-ray)
- Surge force scales with f² and A²
- Increasing N (ray count) increases total thrust but with diminishing returns
- Analytical model validated against CFD within acceptable error margins
