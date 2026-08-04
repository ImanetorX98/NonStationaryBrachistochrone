# -*- coding: utf-8 -*-
"""Generate separate rail-indicatrix figures for Paper I and Paper II.

The historical three-panel asset mixed FLRW, Vaidya and Thakurta-Kerr.  This
generator preserves both scientific parts as independent vector figures:

* fig_indicatrici / fig_indicatrici_vaidya: FLRW + Vaidya, used by Paper I;
* fig_indicatrici_thakurta_kerr: rotating Thakurta-Kerr, retained for Paper II.

Every asset is exported as PDF and PNG to both PaperFigures/ and
paper/Immagini/.  The unsuffixed fig_indicatrici name deliberately remains the
Paper-I include target.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRS = (ROOT / "PaperFigures", ROOT / "paper" / "Immagini")

E_VAIDYA = 1.2
SPIN = 0.9

COLORS_FLRW = ("#3B4CC0", "#168B8C", "#43B66D", "#C7D719")
COLORS_VAIDYA = ("#4B0082", "#CC4778", "#F89540")
COLORS_TK = ("#3B5B92", "#1F9E89", "#73D055")
LINESTYLES_4 = ("-", "--", "-.", ":")
LINESTYLES_3 = ("-", "--", ":")


def configure_style():
    plt.rcParams.update(
        {
            "figure.constrained_layout.use": True,
            "font.size": 9.5,
            "axes.titlesize": 10.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "legend.handlelength": 2.2,
            "lines.linewidth": 1.8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "mathtext.fontset": "dejavusans",
            "savefig.dpi": 300,
            "figure.dpi": 160,
        }
    )


def save_figure(fig, names):
    for output_dir in OUTPUT_DIRS:
        output_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            for extension in ("pdf", "png"):
                target = output_dir / f"{name}.{extension}"
                fig.savefig(target, bbox_inches="tight")
                print(f"saved {target.relative_to(ROOT)}")
    plt.close(fig)


def finish_axis(ax, xlim, ylim):
    ax.set_xlabel(r"$r'$")
    ax.set_ylabel(r"$r\,\phi'$")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.16, linewidth=0.5)


def make_paper_i_figure():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.35))

    # FLRW: the centered circle shrinks as a approaches the rail energy.
    ax = axes[0]
    for scale, color, linestyle in zip(
        (1.0, 1.5, 2.0, 2.75), COLORS_FLRW, LINESTYLES_4
    ):
        speed = np.sqrt(max(1.0 - scale**2 / 3.0**2, 0.0))
        ax.add_patch(
            plt.Circle(
                (0.0, 0.0),
                speed,
                fill=False,
                color=color,
                linestyle=linestyle,
                label=rf"$a={scale}$",
            )
        )
    ax.plot(0.0, 0.0, "k+", markersize=8, markeredgewidth=1.2)
    ax.set_title(r"FLRW ($\hat E=3$): shrinking centered circles")
    ax.legend(loc="upper right", framealpha=0.9)
    finish_axis(ax, (-1.15, 1.15), (-1.15, 1.15))
    ax.text(
        0.03,
        0.04,
        r"no wind; freezing as $a\to\hat E$",
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
    )

    # Vaidya in ingoing EF coordinates: a radial wind offsets the ellipses.
    ax = axes[1]
    for radius, color, linestyle in zip(
        (10.0, 4.0, 2.5), COLORS_VAIDYA, LINESTYLES_3
    ):
        f_value = 1.0 - 2.0 / radius
        w_value = E_VAIDYA**2 - f_value
        center = f_value - E_VAIDYA**2
        ax.add_patch(
            Ellipse(
                (center, 0.0),
                2.0 * E_VAIDYA * np.sqrt(w_value),
                2.0 * np.sqrt(w_value),
                fill=False,
                color=color,
                linestyle=linestyle,
                label=rf"$r={radius}$",
            )
        )
        ax.plot(center, 0.0, "o", color=color, markersize=4.2)
    ax.plot(0.0, 0.0, "k+", markersize=8, markeredgewidth=1.2)
    ax.axvline(0.0, color="0.2", linewidth=0.7)
    ax.set_title(r"Vaidya (EF, $\hat E=1.2$): ingoing radial wind")
    ax.legend(loc="upper right", framealpha=0.9)
    finish_axis(ax, (-3.6, 1.2), (-2.0, 2.0))
    ax.text(
        0.03,
        0.04,
        r"$f-\hat E^2<0$; $m(v)$ moves the oval",
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
    )

    save_figure(fig, ("fig_indicatrici", "fig_indicatrici_vaidya"))


def make_thakurta_kerr_figure():
    fig, ax = plt.subplots(figsize=(4.7, 3.9))
    radius = 2.2
    f_value = 1.0 - 2.0 / radius
    delta = radius**2 - 2.0 * radius + SPIN**2
    p_value = radius**2 + SPIN**2 + 2.0 * SPIN**2 / radius

    for scale, color, linestyle in zip(
        (1.0, 2.0, 3.0), COLORS_TK, LINESTYLES_3
    ):
        speed_sq = 1.0 - scale**2 * f_value / E_VAIDYA**2
        p_bar = p_value + scale**2 * (2.0 * SPIN / radius) ** 2 / E_VAIDYA**2
        phi_center = (2.0 * SPIN / radius) * speed_sq / p_bar
        ellipse_scale_sq = delta * speed_sq / p_bar
        ax.add_patch(
            Ellipse(
                (0.0, phi_center * radius),
                2.0 * np.sqrt(ellipse_scale_sq * delta) / radius,
                2.0 * np.sqrt(ellipse_scale_sq / p_bar) * radius,
                fill=False,
                color=color,
                linestyle=linestyle,
                label=rf"$A={scale}$",
            )
        )
        ax.plot(0.0, phi_center * radius, "o", color=color, markersize=4.2)

    ax.plot(0.0, 0.0, "k+", markersize=8, markeredgewidth=1.2)
    ax.axhline(0.0, color="0.2", linewidth=0.7)
    ax.set_title(
        rf"Thakurta-Kerr ($r={radius}$, $s={SPIN}$): angular frame-dragging wind"
    )
    ax.legend(loc="upper right", framealpha=0.9)
    finish_axis(ax, (-0.45, 0.45), (-0.18, 0.72))
    ax.text(
        0.03,
        0.04,
        r"$A(\eta)$ shrinks the offset ellipse",
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
    )
    save_figure(fig, ("fig_indicatrici_thakurta_kerr",))


def main():
    configure_style()
    make_paper_i_figure()
    make_thakurta_kerr_figure()


if __name__ == "__main__":
    main()
