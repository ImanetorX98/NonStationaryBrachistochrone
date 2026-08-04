#!/usr/bin/env python3
"""Regenerate the Paper-I off-shell validation plot without changing its solver.

The numerical data are obtained by executing the frozen validation script
``vaidya_first_order_offshell.py``.  This wrapper only updates the presentation,
exports a vector PDF, and reports the log--log regression diagnostics used in
the manuscript caption.
"""

from pathlib import Path
import runpy

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__).with_name("vaidya_first_order_offshell.py")
OUTPUT_DIRS = (ROOT / "paper" / "Immagini", ROOT / "PaperFigures")


def log_fit(epsilon: np.ndarray, residual: np.ndarray) -> dict[str, float]:
    """Return slope and diagnostics for an ordinary fit in natural-log space."""
    x = np.log(epsilon)
    y = np.log(residual)
    slope, intercept = np.polyfit(x, y, 1)
    delta = y - (slope * x + intercept)
    ss_res = float(np.dot(delta, delta))
    ss_tot = float(np.dot(y - y.mean(), y - y.mean()))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "l2": float(np.linalg.norm(delta)),
        "linf": float(np.max(np.abs(delta))),
        "r2": 1.0 - ss_res / ss_tot,
    }


def main() -> None:
    data = runpy.run_path(str(SOURCE), run_name="__codex_validation_source__")
    epsilon = np.asarray(data["epss"], dtype=float)
    residual_on = np.asarray(data["res_on"], dtype=float)
    residual_full = np.asarray(data["res"], dtype=float)

    fit_on = log_fit(epsilon, residual_on)
    fit_full = log_fit(epsilon, residual_full)

    plt.rcParams.update(
        {
            "font.size": 9.0,
            "axes.labelsize": 9.0,
            "legend.fontsize": 8.0,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "lines.linewidth": 1.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(4.9, 3.45), constrained_layout=True)

    ax.loglog(
        epsilon,
        residual_on,
        color="0.22",
        marker="s",
        markersize=4.8,
        markerfacecolor="white",
        linestyle="--",
        label=rf"on-shell truncation: slope {fit_on['slope']:.3f}",
    )
    ax.loglog(
        epsilon,
        residual_full,
        color="0.05",
        marker="o",
        markersize=4.8,
        markerfacecolor="0.65",
        linestyle="-",
        label=rf"complete Eq. (22): slope {fit_full['slope']:.3f}",
    )

    # Unit-normalized guide lines make the convergence orders readable without
    # obscuring the measured curves.
    guide = np.array([epsilon.min(), epsilon.max()])
    pivot = epsilon[2]
    ax.loglog(
        guide,
        residual_on[2] * (guide / pivot),
        color="0.55",
        linewidth=0.9,
        linestyle=":",
        label=r"$\mathcal{O}(\epsilon)$ guide",
    )
    ax.loglog(
        guide,
        residual_full[2] * (guide / pivot) ** 2,
        color="0.55",
        linewidth=0.9,
        linestyle="-.",
        label=r"$\mathcal{O}(\epsilon^2)$ guide",
    )

    ax.set_xlabel(r"accretion rate $\epsilon$")
    ax.set_ylabel(r"maximum orbit-shape residual on $8\leq r/M\leq 11$")
    ax.grid(which="major", color="0.82", linewidth=0.55)
    ax.grid(which="minor", color="0.91", linewidth=0.35, linestyle=":")
    ax.legend(loc="best", frameon=True, framealpha=0.96)

    for output_dir in OUTPUT_DIRS:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / "fig_vaidya_offshell.pdf", bbox_inches="tight")
        fig.savefig(output_dir / "fig_vaidya_offshell.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"source: {SOURCE}")
    print(f"epsilon samples: {epsilon.tolist()}")
    for name, fit in (("on-shell", fit_on), ("complete", fit_full)):
        print(
            f"{name}: slope={fit['slope']:.10f}, "
            f"log-fit L2={fit['l2']:.6e}, "
            f"log-fit Linf={fit['linf']:.6e}, R2={fit['r2']:.12f}"
        )


if __name__ == "__main__":
    main()
