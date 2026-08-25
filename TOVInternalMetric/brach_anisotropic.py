#!/usr/bin/env python3
"""
brach_anisotropic.py
====================
Computes and plots coordinate-time and proper-time brachistochrones
for isotropic vs Bowers-Liang anisotropic polytropic stars at fixed
compactness mu and fixed opening angle Delta.

Produces two publication-quality figures:
  figures/fig_brach_aniso_curves.png/.pdf   -- brachistochrone curves
  figures/fig_qstar_vs_lambda.png/.pdf      -- turning radii vs lambda_BL
"""

import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# ── Output ─────────────────────────────────────────────────────────────────────
OUTDIR = "figures"
DPI    = 300

# ── Global parameters ─────────────────────────────────────────────────────────
G  = 1.0
c  = 1.0
MU    = 0.20
GAMMA = 2.0
DELTA = np.pi / 2

LAMBDA_VALUES = [-0.4, -0.2, 0.0, 0.2, 0.4]

# Colorblind-safe diverging palette (Wong/Okabe-Ito)
COLORS = ["#0072B2", "#56B4E9", "#333333", "#E69F00", "#D55E00"]

# ── Matplotlib publication style ───────────────────────────────────────────────
matplotlib.rcParams.update({
    "text.usetex":        False,
    "font.family":        "serif",
    "font.size":          9,
    "axes.labelsize":     10,
    "axes.titlesize":     9,
    "legend.fontsize":    7.5,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "lines.linewidth":    1.5,
    "figure.dpi":         DPI,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.major.size":   3.5,
    "ytick.major.size":   3.5,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
})


# ── EOS: Gamma=2 polytrope ─────────────────────────────────────────────────────
def eos_epsilon(p, K=1.0):
    if p <= 0:
        return 0.0
    rho = np.sqrt(p / K)
    return rho + p / (GAMMA - 1)


# ── TOV + Bowers-Liang integration ────────────────────────────────────────────
def tov_bl_rhs(r, y, lam, K=1.0):
    m, p, Phi = y
    if p <= 0:
        return [0.0, 0.0, 0.0]
    eps  = eos_epsilon(p, K)
    r2   = r * r
    fact = 1.0 - 2.0 * G * m / (r * c**2)
    delta_p = lam * (eps + p) * (eps + 3*p) / fact * r2
    dm_dr   = 4.0 * np.pi * r2 * eps / c**2
    dp_dr   = (
        -G * (eps + p) / (c**2 * r2)
        * (m + 4.0 * np.pi * r**3 * p / c**2)
        / fact
        + 2.0 * delta_p / r
    )
    dPhi_dr = G * (m + 4.0 * np.pi * r**3 * p / c**2) / (r2 * c**2 * fact)
    return [dm_dr, dp_dr, dPhi_dr]


def integrate_star(p_c, lam, K=1.0, Nr=4000):
    r0  = 1e-6
    m0  = 4.0 / 3.0 * np.pi * r0**3 * eos_epsilon(p_c, K) / c**2
    y0  = [m0, p_c, 0.0]

    def event_p_zero(r, y, lam, K):
        return y[1]
    event_p_zero.terminal  = True
    event_p_zero.direction = -1

    r_max = 50.0
    sol = solve_ivp(
        tov_bl_rhs, [r0, r_max], y0,
        args=(lam, K),
        events=event_p_zero,
        method='RK45',
        rtol=1e-9, atol=1e-11,
        dense_output=True,
        max_step=r_max / Nr,
    )

    r_arr   = sol.t
    m_arr   = sol.y[0]
    p_arr   = sol.y[1]
    Phi_arr = sol.y[2]

    R  = r_arr[-1]
    M  = m_arr[-1]
    mu = 2.0 * G * M / (R * c**2)
    A0 = 1.0 - mu
    Phi_shift = 0.5 * np.log(A0) - Phi_arr[-1]
    Phi_arr  += Phi_shift

    return r_arr, m_arr, p_arr, np.exp(2.0 * Phi_arr), R, M, A0


def find_pc_for_mu(lam, mu_target, K=1.0, p_lo=1e-4, p_hi=10.0):
    def residual(pc):
        _, _, _, _, R, M, _ = integrate_star(pc, lam, K)
        return 2.0 * G * M / (R * c**2) - mu_target
    if residual(p_lo) * residual(p_hi) > 0:
        p_hi = 50.0
    return brentq(residual, p_lo, p_hi, xtol=1e-10, rtol=1e-10)


# ── Optical index interpolants ─────────────────────────────────────────────────
def make_indices(r_arr, A_arr, B_arr, A0):
    nt   = 1.0 / np.sqrt(np.maximum(A_arr * (A0 - A_arr), 1e-30))
    ntau = np.sqrt(np.maximum(A_arr / np.maximum(A0 - A_arr, 1e-30), 0.0))
    nt[-1]   = nt[-2]
    ntau[-1] = ntau[-2]
    return CubicSpline(r_arr, nt), CubicSpline(r_arr, ntau)


# ── Half-opening angle integral ────────────────────────────────────────────────
def half_opening(q_star, R, nfun_cs, Bfun_cs, Nr=800):
    nq  = float(nfun_cs(q_star))
    Pq2 = (nq * q_star)**2

    def integrand(t):
        x    = q_star + (R - q_star) * t**2
        dxdt = 2.0 * (R - q_star) * t
        nx   = float(nfun_cs(x))
        Bx   = float(Bfun_cs(x))
        denom = (nx * x)**2 - Pq2
        if denom <= 0:
            return 0.0
        return nq * q_star * np.sqrt(Bx) / (x * np.sqrt(denom)) * dxdt

    ts   = np.linspace(0, 1, Nr + 1)
    vals = np.array([integrand(t) for t in ts])
    return np.trapezoid(vals, ts)


def find_qstar(Delta, R, nfun_cs, Bfun_cs, q_lo=None, q_hi=None):
    if q_lo is None:
        q_lo = R * 0.01
    if q_hi is None:
        q_hi = R * 0.9999
    target = Delta / 2.0
    f  = lambda q: half_opening(q, R, nfun_cs, Bfun_cs) - target
    qs = np.linspace(q_lo, q_hi, 40)
    fs = [f(q) for q in qs]
    for i in range(len(fs) - 1):
        if fs[i] * fs[i + 1] < 0:
            return brentq(f, qs[i], qs[i + 1], xtol=R * 1e-8, rtol=1e-8)
    raise ValueError(f"No turning point found for Delta={Delta:.3f}")


def reconstruct_curve(q_star, R, nfun_cs, Bfun_cs, Nr=600):
    nq  = float(nfun_cs(q_star))
    Pq2 = (nq * q_star)**2
    rs  = np.linspace(q_star * 1.0001, R, Nr)
    phi = np.zeros(Nr)
    for i, r in enumerate(rs):
        Ns   = 400
        ts   = np.linspace(0, 1, Ns + 1)
        def integrand(t, r=r):
            x    = q_star + (r - q_star) * t**2
            dxdt = 2.0 * (r - q_star) * t
            nx   = float(nfun_cs(x))
            Bx   = float(Bfun_cs(x))
            denom = (nx * x)**2 - Pq2
            if denom <= 0:
                return 0.0
            return nq * q_star * np.sqrt(Bx) / (x * np.sqrt(denom)) * dxdt
        vals   = np.array([integrand(t) for t in ts])
        phi[i] = np.trapezoid(vals, ts)
    return rs, phi


def polar_to_cartesian(rs, phi_from_mid, R):
    x_r =  rs * np.sin(phi_from_mid) / R
    y_r =  rs * np.cos(phi_from_mid) / R
    x_l = -rs * np.sin(phi_from_mid) / R
    y_l =  rs * np.cos(phi_from_mid) / R
    x_full = np.concatenate([x_l[::-1], x_r])
    y_full = np.concatenate([y_l[::-1], y_r])
    return x_full, y_full


# ── Main computation ──────────────────────────────────────────────────────────
def compute_all(lambda_values=LAMBDA_VALUES, mu=MU, Delta=DELTA, K=1.0):
    results = {}
    for lam in lambda_values:
        print(f"  lambda_BL = {lam:+.2f} ...", end=' ', flush=True)
        pc = find_pc_for_mu(lam, mu, K)
        r_arr, m_arr, p_arr, A_arr, R, M, A0 = integrate_star(pc, lam, K)
        B_arr = 1.0 / np.maximum(1.0 - 2.0 * G * m_arr / (r_arr * c**2), 1e-10)
        Bfun  = CubicSpline(r_arr, B_arr)
        nt_cs, ntau_cs = make_indices(r_arr, A_arr, B_arr, A0)
        q_t   = find_qstar(Delta, R, nt_cs,   Bfun)
        q_tau = find_qstar(Delta, R, ntau_cs,  Bfun)
        r_t,   phi_t   = reconstruct_curve(q_t,   R, nt_cs,   Bfun)
        r_tau, phi_tau = reconstruct_curve(q_tau,  R, ntau_cs,  Bfun)
        results[lam] = dict(
            R=R, M=M, A0=A0, mu=2 * G * M / (R * c**2),
            q_t=q_t / R, q_tau=q_tau / R,
            r_t=r_t, phi_t=phi_t,
            r_tau=r_tau, phi_tau=phi_tau,
        )
        print(f"q_t/R={q_t/R:.4f}  q_tau/R={q_tau/R:.4f}")
    return results


# ── Figure 1: brachistochrone curves ─────────────────────────────────────────
def plot_curves(results, lambda_values=LAMBDA_VALUES, colors=COLORS, Delta=DELTA):
    """
    Two-panel figure stacked vertically (single column, 3.375 in).
    Top panel   : coordinate-time brachistochrones.
    Bottom panel: proper-time brachistochrones.
    Each panel has a zoom inset on the turning-point cluster.
    """
    # Calcolo altezza: axes_w ≈ 2.70 in (3.375 meno margini);
    # con aspect uguale: axes_h = axes_w × (y_range / x_range) = 2.70 × (1.22/2.26) ≈ 1.46 in
    # 2 pannelli + inter-panel + margini title/legend ≈ 4.2 in totali
    fig, axes = plt.subplots(2, 1, figsize=(3.375, 4.2))

    panels = [
        ("$n_t$  (coordinate time)",    "q_t",   "r_t",   "phi_t"),
        (r"$n_\tau$  (proper time)",     "q_tau", "r_tau", "phi_tau"),
    ]

    phi_end = Delta / 2.0   # surface half-angle

    for ax, (panel_title, q_key, r_key, phi_key) in zip(axes, panels):
        # Unit circle (stellar surface)
        th = np.linspace(0, 2 * np.pi, 500)
        ax.plot(np.cos(th), np.sin(th), '-', color='#888888', lw=0.8,
                alpha=0.55, zorder=0)

        # Surface endpoints (zorder basso: l'inset li copre)
        for sgn in [-1, 1]:
            ax.plot(sgn * np.sin(phi_end), np.cos(phi_end),
                    'ko', ms=3.5, zorder=2)

        # Brachistochrone curves
        q_vals = []
        for lam, col in zip(lambda_values, colors):
            res = results[lam]
            R   = res['R']
            xc, yc = polar_to_cartesian(res[r_key], res[phi_key], R)
            lw  = 2.2 if lam == 0.0 else 1.3
            ax.plot(xc, yc, '-', color=col, lw=lw, zorder=3, solid_capstyle='round')
            ax.plot(0, res[q_key], 'o', color=col, ms=3.5, zorder=5)
            q_vals.append(res[q_key])

        # Axis cosmetics
        ax.set_xlim(-1.13, 1.13)
        ax.set_ylim(-0.05, 1.17)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(r'$x/R$')
        ax.set_ylabel(r'$y/R$')
        ax.set_title(panel_title, pad=4)
        ax.grid(True, ls=':', lw=0.5, alpha=0.3, zorder=1)
        ax.tick_params(which='both', top=True, right=True)

        # ── Zoom inset on turning-point cluster ──────────────────────────────
        q_min, q_max = min(q_vals), max(q_vals)
        dq    = q_max - q_min            # actual spread (no artificial floor)
        pad_y = dq * 0.45               # tight vertical padding
        pad_x = 0.025                   # narrow x window → curves nearly horizontal

        # Inset axes in axes-fraction coordinates (upper-right of panel)
        axins = ax.inset_axes([0.56, 0.52, 0.42, 0.44])
        axins.set_facecolor('white')   # schermo bianco: copre i marker del parent

        for lam, col in zip(lambda_values, colors):
            res = results[lam]
            R   = res['R']
            xc, yc = polar_to_cartesian(res[r_key], res[phi_key], R)
            lw  = 2.2 if lam == 0.0 else 1.3
            axins.plot(xc, yc, '-', color=col, lw=lw, solid_capstyle='round')
            axins.plot(0, res[q_key], 'o', color=col, ms=4.5, zorder=5)

        axins.set_xlim(-pad_x, pad_x)
        axins.set_ylim(q_min - pad_y, q_max + pad_y)
        axins.set_aspect('auto')
        axins.tick_params(labelsize=6, which='both',
                          top=True, right=True, direction='in')
        axins.xaxis.set_major_locator(ticker.MultipleLocator(0.02))
        axins.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
        axins.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune='both'))
        axins.grid(True, ls=':', lw=0.4, alpha=0.4)

        # Box + connectors from zoom region to inset
        mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.55", lw=0.7)

    # ── Shared colour legend ──────────────────────────────────────────────────
    legend_handles = []
    for lam, col in zip(lambda_values, colors):
        lw = 2.2 if lam == 0.0 else 1.3
        legend_handles.append(
            Line2D([0], [0], color=col, lw=lw,
                   label=rf'$\lambda_{{\rm BL}}={lam:+.1f}$')
        )

    # Legend: 3 colonne su 2 righe (larghezza singola colonna)
    fig.legend(
        handles=legend_handles,
        ncol=3,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.005),
        framealpha=0.95,
        edgecolor='0.7',
        fontsize=7,
    )

    mu_val = results[0.0]['mu']
    fig.suptitle(
        rf'Bowers-Liang star  ($\mu={mu_val:.2f}$, $\Delta=\pi/2$, $\Gamma={GAMMA:.0f}$)',
        fontsize=8.5,
    )

    fig.tight_layout(pad=0.4, h_pad=0.6, rect=[0, 0.13, 1, 0.95])

    os.makedirs(OUTDIR, exist_ok=True)
    outbase = os.path.join(OUTDIR, "fig_brach_aniso_curves")
    fig.savefig(outbase + ".png", dpi=DPI, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outbase + ".pdf",           bbox_inches="tight", pad_inches=0.08)
    print(f"Saved {outbase}.png / .pdf")
    plt.close(fig)


# ── Figure 2: q* vs lambda_BL ────────────────────────────────────────────────
def plot_qstar(results, lambda_values=LAMBDA_VALUES, colors=COLORS):
    """
    Single-column figure (3.375 in).
    Shows q^(t) and q^(tau) vs lambda_BL with shaded t-tau splitting.
    """
    lams   = np.array(lambda_values)
    q_ts   = np.array([results[lam]['q_t']   for lam in lambda_values])
    q_taus = np.array([results[lam]['q_tau']  for lam in lambda_values])

    fig, ax = plt.subplots(figsize=(3.375, 2.9))

    # Shaded splitting
    ax.fill_between(lams, q_ts, q_taus,
                    color='#BBBBBB', alpha=0.45, zorder=1,
                    label=r'$t$–$\tau$ splitting')

    # q_t points
    ax.plot(lams, q_ts, 'o-',
            color='#0072B2', lw=1.6, ms=5, mec='white', mew=0.6,
            zorder=3, label=r'$q^{(t)}$  (coord.-time)')

    # q_tau points
    ax.plot(lams, q_taus, 's--',
            color='#D55E00', lw=1.6, ms=5, mec='white', mew=0.6,
            zorder=3, label=r'$q^{(\tau)}$  (proper-time)')

    # Isotropic reference line
    ax.axvline(0, color='k', lw=0.7, ls=':', alpha=0.6)

    ax.set_xlabel(r'$\lambda_{\rm BL}$')
    ax.set_ylabel(r'$q^{(i)}/R$')
    ax.set_xlim(lams[0] - 0.05, lams[-1] + 0.05)

    y_all = np.concatenate([q_ts, q_taus])
    dy = y_all.max() - y_all.min()
    ax.set_ylim(y_all.min() - 0.2 * dy, y_all.max() + 0.3 * dy)

    ax.legend(fontsize=7, framealpha=0.9, edgecolor='0.7', loc='best')
    ax.grid(True, ls=':', lw=0.5, alpha=0.35)
    ax.tick_params(which='both', top=True, right=True)

    mu_val = results[0.0]['mu']
    ax.set_title(
        rf'Turning radii vs anisotropy  ($\mu={mu_val:.2f}$, $\Delta=\pi/2$)',
        fontsize=8.5,
    )

    fig.tight_layout()

    os.makedirs(OUTDIR, exist_ok=True)
    outbase = os.path.join(OUTDIR, "fig_qstar_vs_lambda")
    fig.savefig(outbase + ".png", dpi=DPI, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outbase + ".pdf",           bbox_inches="tight", pad_inches=0.08)
    print(f"Saved {outbase}.png / .pdf")
    plt.close(fig)


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print(f"Computing brachistochrones: mu={MU}, Delta=pi/{np.pi/DELTA:.0f}, "
          f"Gamma={GAMMA}")
    print(f"Lambda_BL values: {LAMBDA_VALUES}\n")

    results = compute_all()

    print("\nGenerating figures...")
    plot_curves(results)
    plot_qstar(results)

    print("\nDone.")
    print(f"  {OUTDIR}/fig_brach_aniso_curves.png/.pdf")
    print(f"  {OUTDIR}/fig_qstar_vs_lambda.png/.pdf")
