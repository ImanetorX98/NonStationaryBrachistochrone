#!/usr/bin/env python3
"""Proof of concept for a slowly expanding/contracting ring-BEC analogue.

The numerical observable is a linear phonon field on a ring with prescribed
radius R(t).  Each Fourier amplitude obeys

    q_n'' + gamma q_n' + [n c_s / R(t)]^2 q_n = 0.

It demonstrates measurable non-stationarity and redshift, but it is a wave
analogue rather than a direct realization of the massive controlled rail.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from experimental_style_codex import configure_style, panel_label, save_figure
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


RADIUS0 = 1.0
EXPANSION = 0.32
RAMP_TIME = 8.0
TOTAL_TIME = 13.0
SOUND_SPEED = 1.0
DAMPING = 0.018
MODES = np.array([6.0, 9.0])


def smoothstep(value: np.ndarray | float) -> np.ndarray:
    value = np.clip(np.asarray(value, dtype=float), 0.0, 1.0)
    return value**2 * (3.0 - 2.0 * value)


def radius(time: np.ndarray | float) -> np.ndarray:
    return RADIUS0 * (1.0 + EXPANSION * smoothstep(np.asarray(time) / RAMP_TIME))


def mode_rhs(time: float, state: np.ndarray) -> np.ndarray:
    count = len(MODES)
    q = state[:count]
    qdot = state[count:]
    omega = MODES * SOUND_SPEED / radius(time)
    return np.r_[qdot, -DAMPING * qdot - omega**2 * q]


def simulate_modes():
    initial_q = np.array([0.78, 0.38])
    initial_qdot = np.array([0.0, 0.45])
    times = np.linspace(0.0, TOTAL_TIME, 2200)
    solution = solve_ivp(
        mode_rhs,
        (times[0], times[-1]),
        np.r_[initial_q, initial_qdot],
        t_eval=times,
        rtol=1.0e-10,
        atol=1.0e-12,
        max_step=0.01,
    )
    return times, solution.y[: len(MODES)], solution.y[len(MODES) :]


def draw_apparatus(ax: plt.Axes) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-5.6, 6.0)
    ax.set_ylim(-4.8, 5.0)
    ax.axis("off")
    ax.add_patch(Circle((0.0, 0.0), 2.55, ec="#334155", fc="#e0f2fe", alpha=0.62, lw=1.3))
    ax.add_patch(Circle((0.0, 0.0), 1.38, fill=False, ec="#06b6d4", lw=8.0, alpha=0.72))
    ax.add_patch(Circle((0.0, 0.0), 1.02, fill=False, ec="#a855f7", lw=1.0, ls="--"))
    ax.add_patch(Circle((0.0, 0.0), 1.78, fill=False, ec="#f97316", lw=1.0, ls="--"))
    ax.text(0.0, -0.05, "ring BEC", ha="center", va="center", fontsize=9)

    elements = [
        ((-5.15, 2.65), "SLM", "#ddd6fe"),
        ((3.85, 2.65), "AOM", "#fee2e2"),
        ((-1.2, 4.15), "imaging", "#e5e7eb"),
    ]
    for (x, y), text, color in elements:
        patch = FancyBboxPatch((x, y), 1.6, 0.72, boxstyle="round,pad=0.08", fc=color, ec="#111827")
        ax.add_patch(patch)
        ax.text(x + 0.8, y + 0.36, text, ha="center", va="center", fontsize=8)
    ax.add_patch(FancyArrowPatch((-3.55, 3.0), (-1.85, 1.65), arrowstyle="-|>", color="#7c3aed", mutation_scale=10, lw=1.5))
    ax.add_patch(FancyArrowPatch((3.85, 3.0), (1.75, 1.7), arrowstyle="-|>", color="#dc2626", mutation_scale=10, lw=1.5))
    ax.add_patch(FancyArrowPatch((-0.4, 4.15), (-0.15, 2.65), arrowstyle="-|>", color="#475569", mutation_scale=10, lw=1.2))
    ax.text(0.0, -3.55, "three programmed trap radii\n+ absorption imaging of density waves", ha="center", fontsize=8)
    ax.set_title("Cold-atom apparatus concept")


def draw_kymograph(ax: plt.Axes, times: np.ndarray, amplitudes: np.ndarray) -> None:
    theta = np.linspace(-np.pi, np.pi, 480)
    phase = np.array([0.0, 0.75])
    density = np.zeros((len(times), len(theta)))
    for index, mode in enumerate(MODES.astype(int)):
        density += amplitudes[index, :, None] * np.cos(mode * theta[None, :] + phase[index])
    image = ax.imshow(
        density,
        origin="lower",
        aspect="auto",
        extent=[-np.pi, np.pi, times[0], times[-1]],
        cmap="RdBu_r",
        vmin=-1.15,
        vmax=1.15,
        interpolation="bilinear",
    )
    ax.axhline(RAMP_TIME, color="white", lw=0.9, ls="--")
    ax.text(3.0, RAMP_TIME + 0.15, "ramp complete", color="white", fontsize=7, ha="right")
    ax.set_xlabel(r"azimuth $\theta$")
    ax.set_ylabel("laboratory time")
    ax.set_title("Predicted density-wave kymograph")
    colorbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    colorbar.set_label(r"$\delta n$ (arb. units)")


def draw_observables(ax: plt.Axes, times: np.ndarray, amplitudes: np.ndarray) -> None:
    radii = radius(times)
    omega = MODES[:, None] * SOUND_SPEED / radii[None, :]
    axis_frequency = ax
    axis_radius = ax.twinx()
    axis_frequency.plot(times, omega[0], color="#2563eb", label=rf"$\omega_{{{int(MODES[0])}}}(t)$")
    axis_frequency.plot(times, omega[1], color="#7c3aed", label=rf"$\omega_{{{int(MODES[1])}}}(t)$")
    axis_radius.plot(times, radii, color="#f97316", ls="--", label=r"$R(t)$")
    axis_frequency.set_xlabel("laboratory time")
    axis_frequency.set_ylabel("instantaneous phonon frequency")
    axis_radius.set_ylabel("ring radius", color="#c2410c")
    axis_frequency.grid(alpha=0.20)
    lines = axis_frequency.lines + axis_radius.lines
    labels = [line.get_label() for line in lines]
    axis_frequency.legend(lines, labels, frameon=False, loc="upper right")
    axis_frequency.set_title("Directly calibratable redshift under the slow ramp")
    inset = axis_frequency.inset_axes([0.08, 0.10, 0.40, 0.26])
    inset.plot(times, amplitudes[0], color="#2563eb", lw=0.8)
    inset.plot(times, amplitudes[1], color="#7c3aed", lw=0.8)
    inset.set_title("mode amplitudes", fontsize=7)
    inset.set_xticks([])
    inset.set_yticks([])


def main() -> None:
    configure_style()
    times, amplitudes, _ = simulate_modes()
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.65))
    draw_apparatus(axes[0])
    draw_kymograph(axes[1], times, amplitudes)
    draw_observables(axes[2], times, amplitudes)
    for ax, label in zip(axes, ("(a)", "(b)", "(c)")):
        panel_label(ax, label)
    fig.suptitle("Experimental proof of concept III — slowly breathing ring BEC", fontsize=12)
    save_figure(fig, "ring_bec_poc_codex")


if __name__ == "__main__":
    main()
