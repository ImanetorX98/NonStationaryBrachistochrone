#!/usr/bin/env python3
"""
fig_nlo_delta90_semicerchio.py

Figura publication-ready: famiglia NLO di brachistocrone per Delta=90 deg,
zoom sul semicerchio destro per massimizzare la visibilita' dell'effetto.

Output:
  figures/fig_nlo_delta90_semicerchio.png  (300 DPI)
  figures/fig_nlo_delta90_semicerchio.pdf  (vettoriale)
"""

import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ---------------------------------------------------------------------------
# Parametri
# ---------------------------------------------------------------------------
DELTA_DEG = 90.0
N_CURVES  = 10       # curve NLO (escludendo Newton)
MU_MIN    = 0.05
MU_MAX    = 0.85     # NLO usato illustrativamente; accurato per mu piccolo
N_PTS     = 4000     # punti per curva

OUTDIR  = "figures"
OUTBASE = "fig_nlo_delta90_semicerchio"
DPI     = 300

# ---------------------------------------------------------------------------
# Funzioni NLO (da genera_figure_paper.py)
# ---------------------------------------------------------------------------
def q_nlo(delta_rad: float, mu: float) -> float:
    q0 = 1.0 - delta_rad / math.pi
    return q0 + mu * q0 * (15.0 - 7.0 * q0**2) / 32.0


def curva_brachistocrona(delta_rad: float, q: float, n: int = 4000):
    """
    Curva brachistocrona newtoniana riancorata agli estremi Delta_target.
    """
    q = max(0.0, min(0.999999, float(q)))
    nh = n // 2
    u = np.linspace(0.0, math.pi / 2.0, nh)

    rho = np.sqrt(q**2 + (1.0 - q**2) * np.sin(u)**2)
    phi = np.arctan2(np.sin(u), q * np.cos(u)) - q * u

    rho_left = rho[::-1];   phi_left = -phi[::-1]
    rho_right = rho[1:];    phi_right = phi[1:]

    rho_full = np.concatenate([rho_left, rho_right])
    phi_full = np.concatenate([phi_left, phi_right])

    # Rinormalizzazione agli estremi target
    phi_end = abs(phi_full[0])
    if phi_end > 1e-15:
        phi_full = phi_full * (0.5 * delta_rad / phi_end)

    x = rho_full * np.cos(phi_full)
    y = rho_full * np.sin(phi_full)
    return x, y

# ---------------------------------------------------------------------------
# Calcolo traiettorie
# ---------------------------------------------------------------------------
delta_rad = math.radians(DELTA_DEG)
mu_values = np.linspace(MU_MIN, MU_MAX, N_CURVES)

print("=== NLO trajectory family, Delta=90 deg, semicircle zoom ===")

q0 = 1.0 - delta_rad / math.pi  # Newtonian
x_newt, y_newt = curva_brachistocrona(delta_rad, q0, n=N_PTS)
print(f"Newton:  q={q0:.4f}")

curves = []
q_list = []
for mu in mu_values:
    q = q_nlo(delta_rad, float(mu))
    x, y = curva_brachistocrona(delta_rad, q, n=N_PTS)
    curves.append((x, y))
    q_list.append(q)
    print(f"mu={mu:.3f}  q_NLO={q:.4f}")

# ---------------------------------------------------------------------------
# Limiti del plot (zoom sul semicerchio destro)
# Gli estremi della curva Delta=90 si trovano a x=cos(45)=0.707, y=+-0.707
# Il punto piu' a sinistra e' q_Newton = 0.5; aggiungiamo margine sx
# ---------------------------------------------------------------------------
x_left_margin = q0 - 0.08   # 0.42: ampio margine sinistro per vedere le curve
x_right_lim   = 1.07        # leggermente fuori dalla superficie
y_lim         = 0.82        # copre i poli +-0.707 con margine

# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "text.usetex":      False,
    "font.family":      "serif",
    "font.size":        11,
    "axes.labelsize":   12,
    "legend.fontsize":  10,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "figure.dpi":       DPI,
})

# Figura con aspetto tale da far sembrare il cerchio circolare nel crop scelto
x_range = x_right_lim - x_left_margin   # ~0.65
y_range = 2 * y_lim                      # 1.64
fig_w = 5.5
fig_h = fig_w * (y_range / x_range)     # mantieni proporzioni fisiche
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# Arco del cerchio unita' (solo parte destra, phi in [-90, 90] deg)
phi_arc = np.linspace(-math.pi / 2, math.pi / 2, 600)
ax.plot(np.cos(phi_arc), np.sin(phi_arc), "k-", linewidth=1.8, zorder=2, label=None)

# Curva newtoniana
ax.plot(x_newt, y_newt, color="#1f77b4", linewidth=2.8, zorder=3,
        label=r"Newtonian ($\mu = 0$,  $q_0 = 0.500$)")

# Colormap NLO
cmap = plt.get_cmap("plasma")
norm = Normalize(vmin=MU_MIN, vmax=MU_MAX)

for (x, y), mu, q in zip(curves, mu_values, q_list):
    col = cmap(norm(float(mu)))
    ax.plot(x, y, color=col, linewidth=2.0, alpha=0.95, zorder=3)

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(r"$\mu = r_s/R$  (NLO)", labelpad=8)
# Tick ogni 0.1
cbar_ticks = np.arange(0.1, MU_MAX + 0.01, 0.1)
cbar.set_ticks(cbar_ticks)
cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

# Punti estremi della traiettoria (su superficie)
ax.plot(math.cos(math.radians(45)),  math.sin(math.radians(45)),
        "ko", markersize=6, zorder=5)
ax.plot(math.cos(math.radians(-45)), math.sin(math.radians(-45)),
        "ko", markersize=6, zorder=5)

# Griglia e assi
ax.set_xlabel(r"$x/R$")
ax.set_ylabel(r"$y/R$")
ax.set_xlim(x_left_margin, x_right_lim)
ax.set_ylim(-y_lim, y_lim)
ax.set_aspect("equal", adjustable="box")
ax.grid(True, linestyle=":", alpha=0.25, zorder=1)
ax.legend(loc="lower right", framealpha=0.85)

fig.tight_layout()

# ---------------------------------------------------------------------------
# Salvataggio
# ---------------------------------------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)
out_png = os.path.join(OUTDIR, OUTBASE + ".png")
out_pdf = os.path.join(OUTDIR, OUTBASE + ".pdf")
fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
fig.savefig(out_pdf, bbox_inches="tight")
plt.close(fig)
print(f"\nSalvato: {out_png}")
print(f"Salvato: {out_pdf}")
