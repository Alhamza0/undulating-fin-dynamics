"""
Example: Robot body motion simulation.

Reproduces Thesis Ch. 6 — forward motion of dual-fin robot.

Usage
-----
    python examples/run_robot.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import matplotlib.pyplot as plt
from undulating_fin.params import FinParams, RobotParams, SimulationConfig
from undulating_fin.robot import MotionConfig, simulate_robot
from undulating_fin.plotting import plot_robot_motion


def main():
    # Robot body (Thesis Sec. 6.2)
    robot = RobotParams(R=0.1, L=0.22, m=2.5, rho=1022.0)

    # Forward motion (Thesis Sec. 6.3.1, Fig. 6-2)
    motion = MotionConfig.forward(
        A=np.deg2rad(30),
        f=1.0,
        Phi=np.deg2rad(60),
        N=3,
        D=0.10,
        hmax=0.175,
        hmin=0.0,
        Cn=2.65,
        rho=1022.0,
    )

    sim = SimulationConfig(t_max=6.0, dt=0.01)

    print("Simulating forward motion (dual fins)...")
    print(f"  Right fin: A={np.rad2deg(motion.right_fin.A):.0f}°, "
          f"f={motion.right_fin.f}Hz, Φ={np.rad2deg(motion.right_fin.Phi):.0f}°")
    print(f"  Left fin:  A={np.rad2deg(motion.left_fin.A):.0f}°, "
          f"f={motion.left_fin.f}Hz, Φ={np.rad2deg(motion.left_fin.Phi):.0f}°")

    results = simulate_robot(robot, motion, sim)

    # Print summary
    t = results["time"]
    pos = results["pos"]
    vel = results["vel"]
    print(f"\nFinal position: x={pos[-1,0]:.3f}m, y={pos[-1,1]:.4f}m, z={pos[-1,2]:.4f}m")
    print(f"Final velocity: vx={vel[-1,0]:.3f}m/s, vy={vel[-1,1]:.4f}m/s, vz={vel[-1,2]:.4f}m/s")
    print(f"Max surge force: {np.max(results['F_total'][:, 0]):.3f} N")

    fig = plot_robot_motion(results)
    fig.suptitle("Forward motion — Dual-fin robot (Fig. 6-4)", fontsize=14)

    plt.show()


if __name__ == "__main__":
    main()
