# -*- coding: utf-8 -*-
"""
FLRW brachistochrones at arbitrary spatial curvature k = +1, 0, -1.

The rail selector is W = d_eta with |W| = a(eta) for every k, so the local
speed v = sqrt(1 - a^2/Ehat^2) and the clock relation dtau/deta = a^2/Ehat are
k-independent: the brachistochrone is the SPATIAL GEODESIC of the constant-
curvature slice, and curvature enters only through the geodesic distance.

Panel (a)  the three geodesics drawn in one comoving map -- geodesic polar
           coordinates (chi, phi) about a reference point O, same closest
           approach chi_min.  Flat -> straight line; closed bends toward O;
           open bends away.
Panel (b)  maximal reachable comoving distance in de Sitter,
           Delta_max = int v deta  from a = 1 to the freezing surface a = Ehat,
           against the antipodal distance pi of a closed slice of unit comoving
           curvature radius.  Where Delta_max > pi the antipode is reachable and
           the minimiser stops being unique (a Maxwell set).

Output: paper/Immagini/fig_flrw_curvatura.{pdf,png}
"""

import os
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig            # noqa: E402
import matplotlib.pyplot as plt                            # noqa: E402

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, '..', 'paper', 'Immagini'))

# Okabe-Ito subset: colourblind-safe, validated; linestyle is the secondary
# encoding so the panel survives greyscale printing.
STYLE = {
    +1: ('#0072B2', '-',  r'closed  $k=+1$'),
     0: ('#E69F00', '--', r'flat  $k=0$'),
    -1: ('#009E73', '-.', r'open  $k=-1$'),
}


def S_k(chi, k):
    """Curvature function: sin, identity, sinh."""
    if k > 0:
        return np.sin(chi)
    if k < 0:
        return np.sinh(chi)
    return chi


def dS_k(chi, k):
    if k > 0:
        return np.cos(chi)
    if k < 0:
        return np.cosh(chi)
    return np.ones_like(chi)


def phi_of_chi(chi, chi_min, k):
    """Azimuth swept by the geodesic with closest approach chi_min.

    Clairaut: L = S_k^2 dphi/ds with L = S_k(chi_min), and
        dphi/dchi = L / ( S_k sqrt(S_k^2 - L**2) ).
    The integrand has an inverse-square-root singularity at chi_min, removed by
    the substitution chi = chi_min + t^2.
    """
    L = S_k(chi_min, k)

    def integrand(t):
        c = chi_min + t * t
        s = S_k(c, k)
        rad = s * s - L * L
        if rad <= 0.0:
            # limit t -> 0 : integrand -> 2 L / (S sqrt(2 S S' )) * 1/sqrt(1)
            return 2.0 * L / (L * np.sqrt(2.0 * S_k(chi_min, k)
                                          * dS_k(chi_min, k)))
        return 2.0 * t * L / (s * np.sqrt(rad))

    out = np.empty_like(chi)
    for i, c in enumerate(chi):
        t_up = np.sqrt(max(c - chi_min, 0.0))
        out[i] = quad(integrand, 0.0, t_up, limit=200)[0]
    return out


def delta_max(Ehat, H):
    """de Sitter: comoving distance reachable before freezing, a: 1 -> Ehat."""
    return (np.sqrt(Ehat ** 2 - 1.0) + np.arcsin(1.0 / Ehat)
            - np.pi / 2.0) / (Ehat * H)


fig, (axA, axB) = plt.subplots(2, 1, figsize=(COL, 4.6),
                               gridspec_kw={'height_ratios': [1.0, 1.25]})

# ------------------------------------------------------------------ panel (a)
# The geodesics are symmetric about the closest approach, so only the upper
# half is drawn: that buys horizontal room at column width.
chi_min = 0.60
chi_out = 1.55
chi = np.linspace(chi_min, chi_out, 400)

# drawn with the closest approach pointing "up", so that the flat geodesic is a
# horizontal line at height chi_min and curvature shows as bending toward O
# (closed) or away from it (open)
for k in (+1, 0, -1):
    colour, ls, label = STYLE[k]
    phi = phi_of_chi(chi, chi_min, k)
    x, y = chi * np.sin(phi), chi * np.cos(phi)
    axA.plot(x, y, color=colour, ls=ls, lw=1.4, zorder=3)
    axA.annotate(label.replace('  ', ' '), xy=(x[-1], y[-1]),
                 xytext=(4, 0), textcoords='offset points',
                 color=colour, fontsize=6.5, ha='left', va='center')

axA.plot([0], [0], marker='o', ms=4, color='0.25', zorder=4)
axA.annotate('$O$', xy=(0, 0), xytext=(4, -2), textcoords='offset points',
             fontsize=7, color='0.25')
axA.plot([0, 0], [0, chi_min], color='0.6', lw=0.6, ls=':', zorder=1)
axA.annotate(r'$\chi_{\min}$', xy=(0, chi_min / 2), xytext=(3, -2),
             textcoords='offset points', fontsize=6.5, color='0.45')
axA.set_aspect('equal')
axA.set_xlim(-0.06, 2.24)
axA.set_ylim(-0.05, 1.02)
axA.set_xlabel(r'comoving map $\chi\sin\varphi$', labelpad=1)
axA.set_ylabel(r'$\chi\cos\varphi$')
axA.set_title(r'(a) brachistochrone $=$ spatial geodesic, common $\chi_{\min}$')
axA.grid(True, lw=0.3, color='0.9')
axA.set_axisbelow(True)

# ------------------------------------------------------------------ panel (b)
E = np.linspace(1.02, 12.0, 500)
for H, ls, shade in ((0.15, '-', '#0072B2'), (0.30, '--', '#56A9D6'),
                     (0.60, ':', '#9CC9E4')):
    axB.plot(E, delta_max(E, H), color=shade, ls=ls, lw=1.3,
             label=rf'$H R_{{\rm c}}={H}$')
axB.axhline(np.pi, color='0.25', lw=0.9)
axB.annotate(r'antipode, $d_\gamma=\pi$', xy=(11.6, np.pi), xytext=(0, 3),
             textcoords='offset points', ha='right', fontsize=6.5, color='0.25')
axB.fill_between(E, np.pi, 30, color='0.92', zorder=0)
axB.annotate('cut locus reachable:\nminimiser not unique', xy=(6.4, 6.6),
             fontsize=6.5, color='0.35')
axB.set_xlabel(r'rail energy $\hat E$')
axB.set_ylabel(r'reachable $\Delta_{\max}$  (curvature radii)')
axB.set_title(r'(b) de Sitter: reach before freezing $a=\hat E$')
axB.set_ylim(0, 8)
axB.set_xlim(1, 12)
axB.legend(loc='lower right', frameon=False)
axB.grid(True, lw=0.3, color='0.9')
axB.set_axisbelow(True)

savefig(fig, OUT, 'fig_flrw_curvatura')

# ---------------------------------------------------------------- consistency
print('\nchecks (half-sweep: azimuth from the closest approach to the far end):')
for k in (+1, 0, -1):
    far = np.pi - chi_min - 1e-9 if k > 0 else 60.0
    half = phi_of_chi(np.array([far]), chi_min, k)[0]
    print(f'  k={k:+d}: half-sweep = {half:.6f} rad   ({half/np.pi:.4f} pi)')
print('  expected: flat -> pi/2 as chi -> infinity; closed -> pi at the')
print('  antipodal far point chi = pi - chi_min (so the full great circle is')
print('  2 pi); open -> strictly less than pi/2, geodesics diverge.')
