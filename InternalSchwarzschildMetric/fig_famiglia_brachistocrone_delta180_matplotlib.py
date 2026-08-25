#!/usr/bin/env python3
"""
fig_famiglia_brachistocrone_delta180_matplotlib.py

Figura publication-ready (matplotlib): famiglia di brachistocrone a tempo
coordinato t, Delta = pi (antipodi), mu da 0 a 8/9.

Output:
  figures/fig_famiglia_brachistocrone_delta180_t.png  (300 DPI)
  figures/fig_famiglia_brachistocrone_delta180_t.pdf  (vettoriale)
"""

import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from fig_brachistocrone_t_antipodi_numerico import (
    seleziona_q_brachistocrona_t,
    traiettoria_da_q_antipodi,
)
from simula_brachistocrona_schwarzschild_interna import U_MAX

# ---------------------------------------------------------------------------
# Parametri
# ---------------------------------------------------------------------------
N_CURVES = 8        # curve relativistiche (escludendo mu=0)
MU_MIN   = 0.790    # mu minimo della famiglia
N_TRAJ   = 2600     # punti per traiettoria

OUTDIR  = "figures"
OUTBASE = "fig_famiglia_brachistocrone_delta180_t"
DPI     = 300

# Spaziatura logaritmica di (8/9 - mu): da (8/9 - 0.7) a 1e-3
_d_max   = U_MAX - MU_MIN       # ~0.1889
_d_min   = 1e-3
distances = np.geomspace(_d_max, _d_min, N_CURVES)
mu_values = U_MAX - distances   # crescenti verso 8/9

# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "text.usetex":      False,   # True se LaTeX installato
    "font.family":      "serif",
    "font.size":        11,
    "axes.labelsize":   12,
    "axes.titlesize":   12,
    "legend.fontsize":  10,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "lines.linewidth":  2.0,
    "figure.dpi":       DPI,
})

# ---------------------------------------------------------------------------
# Calcolo traiettorie
# ---------------------------------------------------------------------------
print("=== Calcolo traiettorie ===")

# mu = 0: diametro verticale
x_newt = np.array([0.0, 0.0])
y_newt = np.array([1.0, -1.0])
print("mu = 0.000000  (Newtonian diameter)")

curves = []
q_list = []
for mu in mu_values:
    q, mode, _ = seleziona_q_brachistocrona_t(float(mu))
    x, y = traiettoria_da_q_antipodi(u=float(mu), q=q, n=N_TRAJ)
    curves.append((x, y))
    q_list.append(q)
    print(f"mu = {mu:.6f}  q = {q:.6f}  mode = {mode}")

# ---------------------------------------------------------------------------
# Figura
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.1, 6.5))

# Cerchio unita' (superficie sfera)
theta = np.linspace(0, 2 * math.pi, 800)
ax.plot(np.cos(theta), np.sin(theta), "k-", linewidth=1.5, zorder=2)

# Curva newtoniana
ax.plot(x_newt, y_newt,
        color="#1f77b4", linewidth=2.5, zorder=3,
        label=r"Newtonian ($\mu = 0$)")

# Colormap per le curve relativistiche
cmap = plt.get_cmap("plasma")
norm = Normalize(vmin=MU_MIN, vmax=float(U_MAX))

for (x, y), mu in zip(curves, mu_values):
    col = cmap(norm(float(mu)))
    ax.plot(x, y, color=col, linewidth=2.0, alpha=0.95, zorder=3)

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(r"$\mu = r_s/R$", labelpad=8)
# Mostra 5 tick selezionati (evita overlap nella zona mu -> 8/9)
tick_idx = [0, 1, 2, 3, 4]
cbar.set_ticks(mu_values[tick_idx])
cbar.ax.set_yticklabels([f"{mu_values[i]:.4f}" for i in tick_idx], fontsize=9)
# Label del limite di Buchdahl in cima
cbar.ax.text(0.5, 1.01, r"$\rightarrow 8/9$", transform=cbar.ax.transAxes,
             ha="center", va="bottom", fontsize=9)

# Poli N/S
pole_kw = dict(color="k", marker="o", markersize=5, zorder=5, linestyle="none")
ax.plot(0, 1,  **pole_kw)
ax.plot(0, -1, **pole_kw)
ax.annotate(r"$N$", xy=(0, 1),  xytext=(0.05, 1.03),  fontsize=12, ha="left", va="bottom")
ax.annotate(r"$S$", xy=(0, -1), xytext=(0.05, -1.07), fontsize=12, ha="left", va="top")

# Assi e griglia
ax.set_xlabel(r"$x/R$")
ax.set_ylabel(r"$y/R$")
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.18, 1.18)
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
