#!/usr/bin/env python3
"""
fig_nlo_delta60_primo_quadrante.py

Figura publication-ready (doppia colonna): famiglia NLO di brachistocrone
per Delta=60 deg, ruotata di 45 deg nel primo quadrante (simmetria y=x).

Output:
  figures/fig_nlo_delta60_primo_quadrante.png  (300 DPI)
  figures/fig_nlo_delta60_primo_quadrante.pdf  (vettoriale)
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
DELTA_DEG = 60.0
N_CURVES  = 10
MU_MIN    = 0.05
MU_MAX    = 0.80
N_PTS     = 4000

OUTDIR  = "figures"
OUTBASE = "fig_nlo_delta60_primo_quadrante"
DPI     = 300

ROT     = math.pi / 4.0   # rotazione 45 deg: asse y=x nel primo quadrante

# ---------------------------------------------------------------------------
# Funzioni
# ---------------------------------------------------------------------------
def q_nlo(delta_rad: float, mu: float) -> float:
    q0 = 1.0 - delta_rad / math.pi
    return q0 + mu * q0 * (15.0 - 7.0 * q0**2) / 32.0


def curva_brachistocrona(delta_rad: float, q: float, n: int = 4000):
    """Curva newtoniana riancorata, simmetrica attorno a phi=0."""
    q = max(1e-6, min(0.999999, float(q)))
    nh = n // 2
    u = np.linspace(0.0, math.pi / 2.0, nh)

    rho = np.sqrt(q**2 + (1.0 - q**2) * np.sin(u)**2)
    phi = np.arctan2(np.sin(u), q * np.cos(u)) - q * u

    rho_l = rho[::-1]; phi_l = -phi[::-1]
    rho_r = rho[1:];   phi_r =  phi[1:]

    rho_f = np.concatenate([rho_l, rho_r])
    phi_f = np.concatenate([phi_l, phi_r])

    phi_end = abs(phi_f[0])
    if phi_end > 1e-15:
        phi_f = phi_f * (0.5 * delta_rad / phi_end)

    x = rho_f * np.cos(phi_f)
    y = rho_f * np.sin(phi_f)
    return x, y


def ruota(x: np.ndarray, y: np.ndarray, angle: float):
    """Rotazione 2D di angolo `angle` (radianti)."""
    c, s = math.cos(angle), math.sin(angle)
    return c * x - s * y, s * x + c * y

# ---------------------------------------------------------------------------
# Calcolo traiettorie
# ---------------------------------------------------------------------------
delta_rad = math.radians(DELTA_DEG)
q0_newt   = 1.0 - delta_rad / math.pi
mu_values = np.linspace(MU_MIN, MU_MAX, N_CURVES)

print(f"=== NLO Delta={DELTA_DEG} deg, primo quadrante (rot 45 deg) ===")
print(f"Newton: q0={q0_newt:.4f}")

x0, y0 = curva_brachistocrona(delta_rad, q0_newt, n=N_PTS)
x0r, y0r = ruota(x0, y0, ROT)

curves_rot = []
q_list     = []
for mu in mu_values:
    q = q_nlo(delta_rad, float(mu))
    xc, yc = curva_brachistocrona(delta_rad, q, n=N_PTS)
    xr, yr = ruota(xc, yc, ROT)
    curves_rot.append((xr, yr))
    q_list.append(q)
    print(f"mu={mu:.3f}  q_NLO={q:.4f}  bisettrice=({q/math.sqrt(2):.3f}, {q/math.sqrt(2):.3f})")

# ---------------------------------------------------------------------------
# Limiti del plot (primo quadrante, zoom sulla zona delle curve)
# Estremi rotati della curva: (cos15, sin15) e (sin15, cos15) ~ (0.966, 0.259)
# Punto minimo sulla bisettrice: q/sqrt(2) ~ (0.471..0.566, idem)
# ---------------------------------------------------------------------------
PAD  = 0.04
# Estremi delle curve sulla superficie: phi=15 deg -> x=cos15=0.966, y=sin15=0.259
# e phi=75 deg -> x=0.259, y=0.966. Serve includere 0.259 con margine.
XMIN = 0.17
XMAX = 1.0 + PAD    # ~1.04
YMIN = XMIN
YMAX = XMAX

# ---------------------------------------------------------------------------
# Matplotlib style (doppia colonna: ~7.1 in larghezza)
# Plot quadrato piu' colorbar: fig 6.5 x 5.5 in
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

fig, ax = plt.subplots(figsize=(6.5, 5.5))

# Quarto di cerchio (primo quadrante: phi 0 -> 90 deg)
phi_arc = np.linspace(0.0, math.pi / 2.0, 500)
ax.plot(np.cos(phi_arc), np.sin(phi_arc), "k-", linewidth=1.8, zorder=2)

# Bisettrice y=x tratteggiata (asse di simmetria)
bmax = 1.02
ax.plot([XMIN, bmax], [XMIN, bmax],
        color="gray", linewidth=0.9, linestyle="--", alpha=0.55, zorder=1,
        label=r"$y = x$  (symmetry axis)")

# Curva newtoniana
ax.plot(x0r, y0r, color="#1f77b4", linewidth=2.6, zorder=4,
        label=rf"Newtonian ($\mu=0$,  $q_0={q0_newt:.3f}$)")

# Curve NLO
cmap = plt.get_cmap("plasma")
norm = Normalize(vmin=MU_MIN, vmax=MU_MAX)

for (xr, yr), mu in zip(curves_rot, mu_values):
    col = cmap(norm(float(mu)))
    ax.plot(xr, yr, color=col, linewidth=1.8, alpha=0.95, zorder=3)

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(r"$\mu = r_s/R$  (NLO)", labelpad=8)
cbar_ticks = np.arange(0.1, MU_MAX + 0.01, 0.1)
cbar.set_ticks(cbar_ticks)
cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

# Punti estremi (su superficie, rotati)
ep_phi = math.radians(DELTA_DEG / 2.0)            # 30 deg
for sign in [+1, -1]:
    xe, ye = math.cos(ep_phi), sign * math.sin(ep_phi)
    xer, yer = ruota(np.array([xe]), np.array([ye]), ROT)
    ax.plot(xer, yer, "ko", markersize=6, zorder=6)

# Assi e griglia
ax.set_xlabel(r"$x/R$")
ax.set_ylabel(r"$y/R$")
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(YMIN, YMAX)
ax.set_aspect("equal", adjustable="box")
ax.grid(True, linestyle=":", alpha=0.25, zorder=1)
ax.legend(loc="upper left", framealpha=0.88, fontsize=9.5)

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
