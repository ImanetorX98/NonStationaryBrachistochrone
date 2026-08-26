#!/usr/bin/env python3
"""Proof of concept for a slowly modulated draining-vortex wave experiment.

The model is the nondispersive shallow-water draining-bathtub flow
u_r=-D/r, u_phi=C/r.  Its acoustic horizon, ergosurface and the two unstable
circular-ray radii are plotted together with ideal Hamiltonian ray traces.
Real data must replace this profile through PIV before comparison with theory.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from experimental_style_codex import configure_style, panel_label, save_figure
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


C0 = 1.20  # circulation parameter, m^2/s in the dimensional DBT model
D0 = 1.00  # drain parameter, m^2/s
WAVE_SPEED = 1.00
TANK_RADIUS = 4.25


def characteristic_radii(circulation: np.ndarray | float, drain: np.ndarray | float):
    """Return horizon, ergosurface, prograde and retrograde circular-ray radii."""
    circulation = np.asarray(circulation, dtype=float)
    drain = np.asarray(drain, dtype=float)
    norm = np.sqrt(circulation**2 + drain**2)
    horizon = drain / WAVE_SPEED
    ergosurface = norm / WAVE_SPEED
    prograde = np.sqrt(2.0 * norm * (norm - circulation)) / WAVE_SPEED
    retrograde = np.sqrt(2.0 * norm * (norm + circulation)) / WAVE_SPEED
    return horizon, ergosurface, prograde, retrograde


def velocity_and_jacobian(x: float, y: float, circulation: float, drain: float):
    """Ideal DBT velocity and Cartesian spatial Jacobian outside the drain core."""
    radius2 = x * x + y * y
    radius2 = max(radius2, 1.0e-10)
    numerator_x = -drain * x - circulation * y
    numerator_y = -drain * y + circulation * x
    velocity = np.array([numerator_x / radius2, numerator_y / radius2])
    radius4 = radius2**2
    jacobian = np.array(
        [
            [(-drain * radius2 - 2.0 * x * numerator_x) / radius4,
             (-circulation * radius2 - 2.0 * y * numerator_x) / radius4],
            [(circulation * radius2 - 2.0 * x * numerator_y) / radius4,
             (-drain * radius2 - 2.0 * y * numerator_y) / radius4],
        ]
    )
    return velocity, jacobian


def ray_rhs(_time: float, state: np.ndarray, circulation: float, drain: float):
    x, y, kx, ky = state
    wavevector = np.array([kx, ky])
    norm_k = max(float(np.linalg.norm(wavevector)), 1.0e-12)
    velocity, jacobian = velocity_and_jacobian(x, y, circulation, drain)
    group_velocity = velocity + WAVE_SPEED * wavevector / norm_k
    wavevector_rate = -jacobian.T @ wavevector
    return np.r_[group_velocity, wavevector_rate]


def circular_ray(sign: float, radius: float, circulation: float, drain: float, duration: float):
    """Integrate a ray initialized on one of the analytic unstable circular rays."""
    radial_component = drain / (WAVE_SPEED * radius)
    tangential_component = sign * np.sqrt(max(1.0 - radial_component**2, 0.0))
    state0 = np.array([radius, 0.0, radial_component, tangential_component])

    horizon, _, _, _ = characteristic_radii(circulation, drain)

    def hit_core(_time, state):
        return np.hypot(state[0], state[1]) - 1.015 * float(horizon)

    def leave_tank(_time, state):
        return TANK_RADIUS - np.hypot(state[0], state[1])

    hit_core.terminal = True
    hit_core.direction = -1
    leave_tank.terminal = True
    leave_tank.direction = -1
    solution = solve_ivp(
        lambda t, y: ray_rhs(t, y, circulation, drain),
        (0.0, duration),
        state0,
        rtol=2.0e-10,
        atol=1.0e-12,
        max_step=0.01,
        events=(hit_core, leave_tank),
    )
    return solution.y[0], solution.y[1]


def draw_apparatus(ax: plt.Axes) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-5.4, 6.0)
    ax.set_ylim(-5.2, 5.4)
    ax.axis("off")
    ax.add_patch(Circle((0.0, 0.0), 4.15, ec="#334155", fc="#cffafe", alpha=0.75, lw=1.5))
    ax.add_patch(Circle((0.0, 0.0), 0.48, ec="#111827", fc="#475569"))
    ax.text(0.0, 0.0, "DRAIN", color="white", ha="center", va="center", fontsize=7)
    for angle in np.linspace(0.2, 5.8, 7):
        start = 2.1 * np.array([np.cos(angle), np.sin(angle)])
        end = 2.1 * np.array([np.cos(angle + 0.38), np.sin(angle + 0.38)])
        ax.add_patch(FancyArrowPatch(start, end, connectionstyle="arc3,rad=0.25", arrowstyle="-|>", color="#0284c7", mutation_scale=8, lw=1.0))
    ax.add_patch(FancyBboxPatch((-1.25, 4.55), 2.5, 0.55, boxstyle="round,pad=0.08", fc="#e5e7eb", ec="#111827"))
    ax.text(0.0, 4.82, "overhead camera + PIV", ha="center", va="center", fontsize=8)
    ax.add_patch(FancyArrowPatch((0.0, 4.55), (0.0, 3.65), arrowstyle="-|>", color="#22c55e", mutation_scale=10))
    ax.add_patch(FancyBboxPatch((4.55, -2.0), 1.15, 1.05, boxstyle="round,pad=0.08", fc="#e2e8f0", ec="#111827"))
    ax.text(5.12, -1.48, "PUMP", ha="center", va="center", fontsize=8)
    ax.add_patch(FancyArrowPatch((4.55, -1.4), (3.65, -0.7), arrowstyle="-|>", color="#475569", mutation_scale=10))
    ax.add_patch(FancyBboxPatch((-5.05, 0.65), 1.45, 0.72, boxstyle="round,pad=0.08", fc="#fef3c7", ec="#111827"))
    ax.text(-4.32, 1.01, "wave maker", ha="center", va="center", fontsize=8)
    ax.add_patch(FancyArrowPatch((-3.6, 1.0), (-2.8, 0.75), arrowstyle="-|>", color="#f59e0b", mutation_scale=10))
    ax.text(0.0, -4.65, "shallow water; closed recirculation loop; one-parameter slow pump ramp", ha="center", fontsize=8)
    ax.set_title("Laboratory layout")


def draw_field_and_rays(ax: plt.Axes) -> None:
    points = np.linspace(-TANK_RADIUS, TANK_RADIUS, 150)
    xx, yy = np.meshgrid(points, points)
    radius2 = xx**2 + yy**2
    mask = (radius2 > 1.02**2) & (radius2 < TANK_RADIUS**2)
    uu = np.full_like(xx, np.nan)
    vv = np.full_like(yy, np.nan)
    uu[mask] = (-D0 * xx[mask] - C0 * yy[mask]) / radius2[mask]
    vv[mask] = (-D0 * yy[mask] + C0 * xx[mask]) / radius2[mask]
    speed = np.hypot(uu, vv)
    ax.streamplot(points, points, uu, vv, color=speed, cmap="Blues", density=1.1, linewidth=0.75, arrowsize=0.65)

    horizon, ergo, r_pro, r_retro = characteristic_radii(C0, D0)
    for radius, linestyle, color, label in (
        (horizon, "-", "#111827", r"horizon $r_h$"),
        (ergo, "--", "#64748b", r"ergosurface $r_e$"),
        (r_pro, ":", "#f97316", r"prograde light ring"),
        (r_retro, "-.", "#db2777", r"retrograde light ring"),
    ):
        ax.add_patch(Circle((0.0, 0.0), float(radius), fill=False, ls=linestyle, ec=color, lw=1.5, label=label))

    pro_x, pro_y = circular_ray(+1.0, float(r_pro) * 1.002, C0, D0, 13.0)
    retro_x, retro_y = circular_ray(-1.0, float(r_retro) * 1.002, C0, D0, 24.0)
    ax.plot(pro_x, pro_y, color="#f97316", lw=2.0)
    ax.plot(retro_x, retro_y, color="#db2777", lw=2.0)
    ax.add_patch(Circle((0.0, 0.0), TANK_RADIUS, fill=False, ec="#334155", lw=1.0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-4.45, 4.45)
    ax.set_ylim(-4.45, 4.45)
    ax.set_xlabel("x / characteristic length")
    ax.set_ylabel("y / characteristic length")
    ax.set_title("Calibrated DBT field and Hamiltonian rays")
    ax.legend(loc="upper right", framealpha=0.90, fontsize=7)


def draw_slow_drift(ax: plt.Axes) -> None:
    time = np.linspace(0.0, 1.0, 400)
    smooth = time**2 * (3.0 - 2.0 * time)
    circulation = C0 * (1.0 + 0.18 * smooth)
    drain = D0 * (1.0 + 0.08 * smooth)
    horizon, ergo, prograde, retrograde = characteristic_radii(circulation, drain)
    ax.plot(time, prograde, color="#f97316", label="prograde light ring")
    ax.plot(time, retrograde, color="#db2777", label="retrograde light ring")
    ax.plot(time, horizon, color="#111827", ls="--", label="horizon")
    ax.plot(time, ergo, color="#64748b", ls=":", label="ergosurface")
    ax.fill_between(time, horizon, ergo, color="#94a3b8", alpha=0.16)
    ax.set_xlabel(r"normalized ramp time $t/T_{\rm ramp}$")
    ax.set_ylabel("radius / characteristic length")
    ax.set_title("Observable motion under a one-parameter slow ramp")
    ax.grid(alpha=0.20)
    ax.legend(frameon=False)
    ax.text(
        0.04, 0.05, "ideal profile; replace C(t), D(t) by PIV inference", transform=ax.transAxes,
        fontsize=7.5, bbox={"fc": "white", "ec": "#cbd5e1", "pad": 2.0}
    )


def main() -> None:
    configure_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.65))
    draw_apparatus(axes[0])
    draw_field_and_rays(axes[1])
    draw_slow_drift(axes[2])
    for ax, label in zip(axes, ("(a)", "(b)", "(c)")):
        panel_label(ax, label)
    fig.suptitle("Experimental proof of concept II — slowly modulated draining vortex", fontsize=12)
    save_figure(fig, "draining_vortex_poc_codex")


if __name__ == "__main__":
    main()
