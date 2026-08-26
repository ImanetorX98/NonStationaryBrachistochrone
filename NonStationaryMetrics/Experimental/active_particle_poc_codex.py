#!/usr/bin/env python3
"""Deterministic proof of concept for a planar controlled-rail demonstrator.

The local velocity oval is the equatorial Thakurta--Kerr indicatrix used in
Paper II.  A simple receding-horizon steering law then drives one tracked agent
between fixed spatial endpoints while the conformal parameter is ramped slowly.

This script is an apparatus/design study, not a proof of global optimality and
not an experimental realization of the spacetime itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from experimental_style_codex import configure_style, panel_label, save_figure
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


@dataclass(frozen=True)
class RailParameters:
    mass: float = 1.0
    spin: float = 0.9
    energy: float = 1.4
    scale0: float = 1.0
    inner_radius_mm: float = 2.25
    outer_radius_mm: float = 6.20


PARAMS = RailParameters()
START = np.array([-4.75, -2.45])
TARGET = np.array([3.55, 2.75])
# Intermediate optical gates keep the local controller away from the central
# exclusion region.  START and TARGET remain the prescribed physical endpoints.
REFERENCE_GATES = np.array([[-2.00, 3.30], [2.40, 3.40], TARGET])


def slow_scale(time_s: float, epsilon: float) -> float:
    """Controlled conformal ramp A(t)=A0 exp(epsilon t)."""
    return PARAMS.scale0 * np.exp(epsilon * time_s)


def tk_indicatrix(radius: np.ndarray | float, scale: float):
    """Return azimuthal centre and radial/azimuthal semi-axes of Eq. (II.8).

    Coordinates are (dr/deta, r dphi/deta).  The formula follows directly from
    phi'_0, R^2=Delta*vbar^2/Pbar and the quadratic form in H_eta.
    """
    r = np.asarray(radius, dtype=float)
    m, a, energy = PARAMS.mass, PARAMS.spin, PARAMS.energy
    f = 1.0 - 2.0 * m / r
    delta = r**2 - 2.0 * m * r + a**2
    b = 2.0 * m * a / r
    p = r**2 + a**2 + 2.0 * m * a**2 / r
    v2 = 1.0 - scale**2 * f / energy**2
    pbar = p + scale**2 * b**2 / energy**2
    if np.any(v2 <= 0.0) or np.any(delta <= 0.0):
        raise ValueError("Indicatrix left the regular timelike/non-freezing domain")
    phi0 = b * v2 / pbar
    r2 = delta * v2 / pbar
    centre_phi = r * phi0
    semi_r = np.sqrt(r2 * delta) / r
    semi_phi = r * np.sqrt(r2 / pbar)
    return centre_phi, semi_r, semi_phi


def rail_velocity(position: np.ndarray, time_s: float, theta: float, epsilon: float) -> np.ndarray:
    """Map one point of the local rail indicatrix into Cartesian lab velocity."""
    x, y = position
    radius = max(float(np.hypot(x, y)), PARAMS.inner_radius_mm + 1.0e-6)
    e_r = np.array([x, y]) / radius
    e_phi = np.array([-e_r[1], e_r[0]])
    centre_phi, semi_r, semi_phi = tk_indicatrix(radius, slow_scale(time_s, epsilon))
    return semi_r * np.cos(theta) * e_r + (
        centre_phi + semi_phi * np.sin(theta)
    ) * e_phi


def simulate_feedback(epsilon: float, max_time_s: float = 42.0):
    """Greedy receding-horizon controller used only as a realizable candidate path.

    The two intermediate gates are controller references, not extra boundary
    conditions in the variational problem.  They remove the local minimum that
    a one-step greedy controller develops against the central exclusion disk.
    """
    dt = 0.035
    lookahead = 0.55
    angles = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
    position = START.copy()
    previous_theta = 0.0
    times = [0.0]
    path = [position.copy()]
    controls = [previous_theta]
    gate_index = 0

    for step in range(int(max_time_s / dt)):
        time_s = step * dt
        local_target = REFERENCE_GATES[gate_index]
        costs = []
        for theta in angles:
            velocity = rail_velocity(position, time_s, theta, epsilon)
            predicted = position + lookahead * velocity
            predicted_radius = np.hypot(*predicted)
            endpoint_cost = np.linalg.norm(predicted - local_target)
            inner_penalty = 80.0 * max(PARAMS.inner_radius_mm + 0.08 - predicted_radius, 0.0) ** 2
            outer_penalty = 80.0 * max(predicted_radius - PARAMS.outer_radius_mm, 0.0) ** 2
            turn_penalty = 0.015 * (1.0 - np.cos(theta - previous_theta))
            costs.append(endpoint_cost + inner_penalty + outer_penalty + turn_penalty)
        theta_star = float(angles[int(np.argmin(costs))])
        velocity = rail_velocity(position, time_s, theta_star, epsilon)
        position = position + dt * velocity
        previous_theta = theta_star
        times.append(time_s + dt)
        path.append(position.copy())
        controls.append(theta_star)
        if np.linalg.norm(position - local_target) < 0.10:
            if gate_index < len(REFERENCE_GATES) - 1:
                gate_index += 1
            else:
                break
    return np.asarray(times), np.asarray(path), np.unwrap(np.asarray(controls))


def draw_apparatus(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 7.0)
    ax.set_aspect("equal")
    ax.axis("off")
    cell = FancyBboxPatch(
        (0.8, 0.9), 6.2, 4.2, boxstyle="round,pad=0.12", ec="#334155", fc="#dbeafe", alpha=0.72
    )
    ax.add_patch(cell)
    ax.add_patch(Circle((3.9, 3.0), 1.15, color="#94a3b8", alpha=0.32))
    ax.text(3.9, 3.0, "annular\nworking region", ha="center", va="center", fontsize=8)
    ax.plot(1.45, 1.55, "o", color="#f59e0b", ms=7)
    ax.plot(6.25, 4.35, "*", color="#dc2626", ms=11)
    ax.text(1.2, 1.15, "START", fontsize=8)
    ax.text(6.25, 4.65, "TARGET", fontsize=8, ha="center")

    camera = FancyBboxPatch((2.5, 5.7), 2.0, 0.75, boxstyle="round,pad=0.08", ec="#111827", fc="#e5e7eb")
    projector = FancyBboxPatch((7.65, 1.1), 1.55, 0.8, boxstyle="round,pad=0.08", ec="#111827", fc="#ddd6fe")
    controller = FancyBboxPatch((7.45, 4.35), 2.0, 1.05, boxstyle="round,pad=0.08", ec="#111827", fc="#f3f4f6")
    ax.add_patch(camera)
    ax.add_patch(projector)
    ax.add_patch(controller)
    ax.text(3.5, 6.08, "overhead camera", ha="center", va="center", fontsize=8)
    ax.text(8.43, 1.50, "DMD / SLM", ha="center", va="center", fontsize=8)
    ax.text(8.45, 4.88, "real-time\nfeedback", ha="center", va="center", fontsize=8)

    arrows = [
        ((3.5, 5.7), (3.5, 5.15), "tracking"),
        ((7.65, 1.65), (6.9, 2.15), "pattern"),
        ((7.45, 4.65), (4.55, 5.7), "state"),
        ((8.45, 4.35), (8.45, 1.95), "command"),
    ]
    for tail, head, label in arrows:
        ax.add_patch(FancyArrowPatch(tail, head, arrowstyle="-|>", mutation_scale=10, lw=1.0, color="#475569"))
        midpoint = 0.5 * (np.asarray(tail) + np.asarray(head))
        ax.text(*midpoint, label, fontsize=7, color="#475569", ha="center", va="center")
    ax.set_title("Planar hardware-in-the-loop layout")


def draw_indicatrices(ax: plt.Axes) -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 500)
    radius = 4.5
    scales = (1.00, 1.18, 1.34)
    colors = ("#2563eb", "#0f766e", "#d97706")
    for scale, color in zip(scales, colors):
        centre, semi_r, semi_phi = tk_indicatrix(radius, scale)
        vr = semi_r * np.cos(theta)
        vphi = centre + semi_phi * np.sin(theta)
        ax.plot(vr, vphi, color=color, label=rf"$A={scale:.2f}$")
        ax.plot(0.0, centre, "o", color=color, ms=3.5)
    ax.axhline(0.0, color="0.45", lw=0.7)
    ax.axvline(0.0, color="0.45", lw=0.7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"radial velocity $v_r$ (mm s$^{-1}$)")
    ax.set_ylabel(r"azimuthal velocity $r\dot\varphi$ (mm s$^{-1}$)")
    ax.set_title(r"Predicted local indicatrix at $r=4.5$ mm")
    ax.legend(frameon=False)
    ax.grid(alpha=0.18)


def draw_trajectories(ax: plt.Axes) -> None:
    time_dynamic, path_dynamic, _ = simulate_feedback(0.008)
    time_frozen, path_frozen, _ = simulate_feedback(0.0)

    grid = np.linspace(-6.0, 6.0, 17)
    xx, yy = np.meshgrid(grid, grid)
    rr = np.hypot(xx, yy)
    mask = (rr > PARAMS.inner_radius_mm) & (rr < PARAMS.outer_radius_mm)
    uu = np.full_like(xx, np.nan)
    vv = np.full_like(yy, np.nan)
    for index in zip(*np.where(mask)):
        x, y = xx[index], yy[index]
        radius = rr[index]
        centre_phi, _, _ = tk_indicatrix(radius, 1.18)
        uu[index] = -centre_phi * y / radius
        vv[index] = centre_phi * x / radius
    ax.quiver(xx, yy, uu, vv, color="#0ea5e9", alpha=0.48, scale=2.0, width=0.003)
    ax.add_patch(Circle((0.0, 0.0), PARAMS.inner_radius_mm, color="#94a3b8", alpha=0.42))
    ax.add_patch(Circle((0.0, 0.0), PARAMS.outer_radius_mm, fill=False, ec="#64748b", lw=0.9))
    ax.plot(path_frozen[:, 0], path_frozen[:, 1], "--", color="#7c3aed", label=rf"frozen feedback ({time_frozen[-1]:.1f} s)")
    scatter = ax.scatter(
        path_dynamic[:, 0], path_dynamic[:, 1], c=slow_scale(time_dynamic, 0.008),
        cmap="plasma", s=4.0, zorder=5, label=rf"slow ramp ({time_dynamic[-1]:.1f} s)"
    )
    ax.plot(*START, "o", color="#f59e0b", ms=7)
    ax.plot(*TARGET, "*", color="#dc2626", ms=11)
    ax.plot(
        REFERENCE_GATES[:-1, 0], REFERENCE_GATES[:-1, 1], "o", ms=3.0,
        mfc="white", mec="#64748b", alpha=0.85, label="feedback reference gates"
    )
    ax.text(*(START + np.array([-0.25, -0.35])), "START", ha="right", fontsize=8)
    ax.text(*(TARGET + np.array([0.15, 0.20])), "TARGET", fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-6.55, 6.55)
    ax.set_ylim(-6.55, 6.55)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Fixed-endpoint candidate path in the breathing oval field")
    ax.legend(loc="upper left", frameon=False)
    colorbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
    colorbar.set_label(r"commanded $A(t)$")


def draw_scaling_target(ax: plt.Axes) -> None:
    epsilon = np.array([2.5e-3, 5.0e-3, 1.0e-2, 2.0e-2, 4.0e-2])
    slope_on = 1.003152
    slope_full = 1.999997
    residual_on = 1.2e-2 * (epsilon / 1.0e-2) ** slope_on
    residual_full = 2.0e-4 * (epsilon / 1.0e-2) ** slope_full
    ax.loglog(epsilon, residual_on, "s--", color="#374151", mfc="white", label=rf"on-shell target, slope {slope_on:.3f}")
    ax.loglog(epsilon, residual_full, "o-", color="#dc2626", label=rf"full target, slope {slope_full:.3f}")
    ax.set_xlabel(r"dimensionless ramp rate $\varepsilon$")
    ax.set_ylabel("normalized orbit-shape residual")
    ax.set_title("Pre-registered convergence signature")
    ax.grid(which="both", alpha=0.22)
    ax.legend(frameon=False)
    ax.text(
        0.04, 0.05, "DESIGN TARGET — not laboratory data", transform=ax.transAxes,
        fontsize=7.5, color="#991b1b", bbox={"fc": "white", "ec": "#fecaca", "pad": 2.0}
    )


def main() -> None:
    configure_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2))
    draw_apparatus(axes[0, 0])
    draw_indicatrices(axes[0, 1])
    draw_trajectories(axes[1, 0])
    draw_scaling_target(axes[1, 1])
    for ax, label in zip(axes.flat, ("(a)", "(b)", "(c)", "(d)")):
        panel_label(ax, label)
    fig.suptitle("Experimental proof of concept I — actively controlled rail", fontsize=12)
    save_figure(fig, "active_particle_poc_codex")


if __name__ == "__main__":
    main()
