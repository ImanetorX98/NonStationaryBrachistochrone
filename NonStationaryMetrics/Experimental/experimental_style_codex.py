#!/usr/bin/env python3
"""Shared plotting style for deterministic experimental proof-of-concept figures."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

# Matplotlib tries to use ~/.matplotlib by default, which is not guaranteed to
# be writable in a clean or sandboxed reproduction environment.
_MPL_CONFIG = Path(tempfile.gettempdir()) / "nonstationarymetrics-mplconfig"
_MPL_CONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"


def configure_style() -> None:
    """Use a compact, paper-compatible style with embedded TrueType fonts."""
    plt.rcParams.update(
        {
            "figure.constrained_layout.use": True,
            "font.size": 9.0,
            "axes.titlesize": 10.0,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "legend.handlelength": 2.0,
            "lines.linewidth": 1.5,
            "axes.linewidth": 0.8,
            "mathtext.fontset": "dejavusans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "figure.dpi": 150,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    """Save a figure as PNG, vector PDF and SVG inside Experimental/output."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for suffix in ("png", "pdf", "svg"):
        target = OUTPUT_DIR / f"{stem}.{suffix}"
        kwargs = {"dpi": 300} if suffix == "png" else {}
        fig.savefig(target, bbox_inches="tight", **kwargs)
        paths.append(target)
        print(f"saved {target.relative_to(HERE.parent)}")
    plt.close(fig)
    return paths


def panel_label(ax: plt.Axes, label: str) -> None:
    """Add a consistent panel label just inside the upper-left corner."""
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.0,
        fontweight="bold",
        zorder=20,
    )
