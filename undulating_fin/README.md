# Undulating Fin — Python Implementation

**Kinematic and Dynamic Modeling of Biomimetic Undulating Fins for Underwater Robots**

Python port of the MATLAB codebase from:
> Ali, A.I. (2021) *"Kinematic and Dynamic Modeling of Biomimetic Undulating Fins for Underwater Robots"*, M.Sc. Thesis, Bahcesehir University.

---

## نظرة عامة | Overview

هذه الحزمة هي النسخة الحديثة من كود MATLAB الأصلي، مكتوبة بلغة Python مع:
- مرجعية كاملة لمعادلات الأطروحة (Eq. 1–19)
- هيكل قابل للتوسعة لمشروع الدكتوراه
- بدون متغيرات عامة (globals → dataclasses)
- تجميع جميع الدوال البارامترية في دالة واحدة عامة
- كود نظيف مع type hints وتوثيق احترافي

---

## Installation

```bash
cd undulating_fin
pip install -e .
```

## 🚀 تشغيل التطبيق | Running the Dashboard

التطبيق يستخدم **Streamlit** لتوفير واجهة تفاعلية في المتصفح.

### الخطوات | Steps

```bash
# 1. افتح Terminal وانتقل لمجلد المشروع
cd "C:\Users\Hamzah\Documents\Projects\My Thesis\undulating_fin"

# 2. (مرة واحدة فقط) تثبيت الحزمة والمتطلبات
pip install -e .

# 3. تشغيل التطبيق
python -m streamlit run app.py
```

سيفتح المتصفح تلقائياً على **http://localhost:8501**

### إيقاف التطبيق | Stopping
اضغط `Ctrl+C` في نافذة Terminal لإيقاف الخادم.

### التبويبات المتاحة | Available Tabs

| Tab | الوصف |
|-----|-------|
| 📈 **Time Domain** | محاكاة ديناميكية كاملة مع عرض القوى والعزوم + تحميل CSV |
| 🔄 **Parametric Sweep** | مسح بارامتري (Φ, A, f, X, D, hmax) أو مقارنة عدد الأشعة |
| 🤖 **Robot Simulation** | محاكاة روبوت ثنائي الزعنفة (أمام، خلف، دوران، مخصص) |
| 🌊 **Kinematics** | رسم تفاعلي فوري للزوايا والسرعات والتسارعات |
| ℹ️ **About** | مرجع المعادلات وشرح النموذج |

جميع المتغيرات يمكن التحكم بها من **الشريط الجانبي** (Sidebar) على يسار الصفحة.

---

## Project Structure

```
undulating_fin/
├── pyproject.toml                    # Package config (PEP 621)
├── README.md                         # This file
├── app.py                            # 🖥️ Streamlit interactive dashboard
├── src/undulating_fin/
│   ├── __init__.py                   # Package entry
│   ├── params.py                     # FinParams, RobotParams (dataclasses)
│   ├── kinematics.py                 # Eq. (1)-(3): Traveling wave
│   ├── hydrodynamics.py              # Eq. (12)-(16): Force integration
│   ├── membrane.py                   # Eq. (8)-(10): Elastic deflection
│   ├── dynamics.py                   # Eq. (4)-(11): Full dynamics
│   ├── robot.py                      # Ch. 6: Dual-fin robot
│   ├── parametric.py                 # Ch. 5: Parametric sweeps
│   └── plotting.py                   # Visualization
└── examples/
    ├── run_analysis.py               # Time-domain (Fig. 2-5, 2-6)
    ├── run_parametric.py             # Parametric sweeps (Ch. 5)
    └── run_robot.py                  # Robot motion (Ch. 6, Fig. 6-4)
```

## Quick Start

### 1. Time-Domain Analysis (Thesis Ch. 2)
```python
from undulating_fin import FinParams, SimulationConfig, compute_kinematics
from undulating_fin.dynamics import run_time_domain
from undulating_fin.plotting import plot_time_domain
import numpy as np

params = FinParams(N=2, A=np.deg2rad(30), f=1.0, Phi=np.deg2rad(50))
sim = SimulationConfig(t_max=3.0, dt=0.01)
results = run_time_domain(params, sim.time_vector)
plot_time_domain(results)
```

### 2. Parametric Sweep (Thesis Ch. 5)
```python
from undulating_fin.parametric import sweep_1d

# Phase shift sweep — replaces MATLAB PhiFN.m
r = sweep_1d(params, sim, 'Phi', np.linspace(0, np.deg2rad(80), 50))
print(f"Optimal ϕ for max surge: {np.rad2deg(r['values'][np.argmax(r['F_max'][:,2])]):.1f}°")
```

### 3. Robot Motion (Thesis Ch. 6)
```python
from undulating_fin.robot import MotionConfig, simulate_robot
from undulating_fin.params import RobotParams

robot = RobotParams(R=0.1, L=0.22, m=2.5)
motion = MotionConfig.forward(A=np.deg2rad(30), f=1.0, Phi=np.deg2rad(60), N=3)
results = simulate_robot(robot, motion, sim)
```

## Mapping: MATLAB → Python

| MATLAB File | Python Equivalent | Notes |
|---|---|---|
| `hydrodynamicforces_Func.m` | `hydrodynamics.compute_panel_forces()` | Eq. 12-16, from fundamentals |
| `hydrodynamictorques_Func.m` | `hydrodynamics.compute_panel_torques()` | Now N-ray capable |
| `DeltaFunc.m` | `membrane.membrane_deflection()` | Eq. 10 |
| `undFinDynamics.m` | `dynamics.run_time_domain()` | Unified Part A + B |
| `PhiFN.m` | `parametric.sweep_1d(p, s, 'Phi', ...)` | One function for all |
| `AmpFN.m` | `parametric.sweep_1d(p, s, 'A', ...)` | " |
| `freqFN.m` | `parametric.sweep_1d(p, s, 'f', ...)` | " |
| `XFN.m` | `parametric.sweep_1d(p, s, 'X', ...)` | " |
| `Phi_A_FUN.m` | `parametric.sweep_2d(p, s, 'Phi', ..., 'A', ...)` | " |
| `Phi_D_FUN.m` | `parametric.sweep_2d(p, s, 'Phi', ..., 'D', ...)` | " |
| `Phi_F_FUN.m` | `parametric.sweep_2d(p, s, 'Phi', ..., 'f', ...)` | " |
| `Phi_H_FUN.m` | `parametric.sweep_2d(p, s, 'Phi', ..., 'hmax', ...)` | " |
| `Number_of_Rays.m` | `parametric.sweep_ray_count()` | " |
| `positionandVelocity.m` | `robot.simulate_robot()` | 3-DOF body dynamics |
| 11 global variables | `FinParams` dataclass | Type-safe, documented |

## Key Improvements Over MATLAB

1. **No globals**: All parameters in `FinParams` dataclass with defaults and docs
2. **No code duplication**: `sweep_1d()` replaces 8 separate parametric functions
3. **N-ray generalization**: Loop over panels instead of hardcoded 11 calls
4. **Modern integration**: `scipy.integrate.dblquad` replaces deprecated `quad2d`
5. **Type safety**: Full type hints for IDE support
6. **Equation references**: Every function documents its thesis equation
7. **Extensible**: Easy to add new parameters, physics, or control strategies

## Thesis Equations Referenced

| Equation | Description | Location in Code |
|---|---|---|
| Eq. (1) | `q_i(t) = A·sin(2πft − (i−1)ϕ) + X` | `kinematics.py` |
| Eq. (2) | `q̇_i(t) = 2πfA·cos(...)` | `kinematics.py` |
| Eq. (3) | `wn = \|ϕ\|·(N−1)/(2π)` | `params.py:wave_count` |
| Eq. (7) | `Mq̈ + Bq̇ + G(q) + Kc·sgn(q̇) + Γ + Φ = τ` | `dynamics.py` |
| Eq. (10) | Membrane deflection `Λ(x)` | `membrane.py` |
| Eq. (12) | `f_n = −½ρCd·en·vn²·sgn(vn)` | `hydrodynamics.py` |
| Eq. (13) | `F_N = ∬ f_n ds` | `hydrodynamics.py` |
| Eq. (14) | Normal velocity `vn(h, w̄)` | `hydrodynamics.py` |
| Eq. (15) | Scale factor `δ(h, w̄)` | `hydrodynamics.py` |
| Eq. (16) | Normal vector `en(h, w̄)` | `hydrodynamics.py` |
