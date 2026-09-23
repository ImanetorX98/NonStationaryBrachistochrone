# -*- coding: utf-8 -*-
"""
Brachistocrone t e tau a ESTREMI FISSI (problema ai limiti, BVP).

Fissati due punti A e B nel piano equatoriale (Schwarzschild statico,
ottica Riemanniana), per ciascun ramo si cerca la geodetica del PROPRIO
indice ottico che connette A e B (fastest-time curve), e se ne estrae il
raggio minimo r_min. t e tau connettono gli STESSI due punti ma sono
curve diverse -> confronto pulito di profondita' (chi affonda di piu').

Indici ottici (rotaia con invariante E, f=1-2M/r):
  n_tau = sqrt(f/(E^2-f)) ,   n_t = E/sqrt(f(E^2-f))
metrica spaziale base dl^2 = dr^2/f + r^2 dphi^2.
Geodetica (Beltrami/Clairaut):  N r^2/sqrt(r^2+r'^2/f) = J  (r'=dr/dphi)
=> dphi/dr = J /( r sqrt(f) sqrt(N^2 r^2 - J^2) ) ,  r_min: N(r_min) r_min = J.

Casi:
  (a) SIMMETRICO:  A=(r0,-Phi), B=(r0,+Phi)  (stesso raggio)
  (b) ASIMMETRICO: A=(rA,0),    B=(rB,dphi)  (raggi diversi)
Per ciascun ramo si spara J per far combaciare Delta_phi = phi_B - phi_A.
"""

import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory (CQG-116884, major 10)
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, E = 1.0, 1.2

def f(r):
    return 1 - 2 * M / r

def n_tau(r):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def rmin_of_J(N, J, r_hi):
    """r_min = radice ESTERNA di N(r) r = J sotto r_hi (scan verso l'interno)."""
    rs = np.linspace(r_hi * (1 - 1e-6), 2.0001 * M, 5000)
    g = N(rs) * rs - J
    for i in range(len(rs) - 1):
        if g[i] * g[i + 1] < 0:
            return brentq(lambda r: N(r) * r - J, rs[i + 1], rs[i])
    return None

def _leg(N, J, rm, rE, kernel):
    """Integrate kernel(r, S) dr from the periastron rm out to rE, where
    S = sqrt(N(r)^2 r^2 - J^2) vanishes like sqrt(r-rm) at the lower limit.

    The singularity is removed exactly, not stepped over.  Factor
    N^2 r^2 - J^2 = (Nr - J)(Nr + J); the first factor has a simple zero at rm,
    so with r = rm + u^2,

        Nr - J = u^2 psi(u^2),   psi(w) = (N(rm+w)(rm+w) - J)/w,
        S = u sqrt( psi(u^2) (Nr + J) ),

    and the 2u du of the substitution cancels the u in S.  The transformed
    integrand is bounded, so plain quad converges without `points` and without
    the extrapolation whose roundoff QUADPACK used to warn about.  The previous
    version integrated from rm + 1e-9 (angle) and rm + 1e-13 (cost) and relied
    on QAGP extrapolating the omitted sliver back: it did, but only by an
    undocumented accident of the epsilon algorithm -- at a nudge of 1e-6 the
    same call is short by 8.1e-4, which is not small.  Agreement of this routine
    with a 50-digit mpmath reproduction: 1e-11 in the aperture, 3e-9 in r_min.
    """
    dNr = (N(rm + 1e-6) * (rm + 1e-6) - N(rm - 1e-6) * (rm - 1e-6)) / 2e-6
    def g(u):
        w = u * u
        r = rm + w
        psi = (N(r) * r - J) / w if w > 1e-10 else dNr
        S = u * np.sqrt(max(psi * (N(r) * r + J), 1e-300))
        return 2.0 * u * kernel(r, S)
    return quad(g, 0.0, np.sqrt(max(rE - rm, 0.0)), limit=400)[0]

def dphi_seg(N, J, rm, r_hi):
    """angolo dal periasse rm a r_hi lungo la geodetica."""
    return J * _leg(N, J, rm, r_hi,
                    lambda r, S: 1.0 / (r * np.sqrt(f(r)) * S))

def cost_optical(N, J, rm, rA, rB):
    """Elapsed clock along the orbit: the optical length int N*sqrt(h(dx,dx)).

    With h = dr^2/f + r^2 dphi^2 and dphi/dr read off the orbit, this is the
    functional whose extremals solve_bvp is shooting for.  It reproduces the
    published pair of solutions at aperture 2.0389 to the quoted digits.
    """
    def integ(r, S):
        dphidr = J / (r * np.sqrt(f(r)) * S)
        return N(r) * np.sqrt(1.0 / f(r) + r**2 * dphidr**2)
    return sum(_leg(N, J, rm, rE, integ) for rE in (rA, rB))


def solve_bvp(N, rA, rB, dphi_target, report=False, nscan=2000):
    """Shoot on J for Delta_phi(A->r_min->B) = dphi_target.

    Enumerates by MONOTONE BRANCHES of the angle map, not by sign changes of the
    defect, and returns the cheapest solution found.

    Two defects have been fixed here, in that order.

    (1) The first version returned the first sign change and never compared
    costs.  The problem has more than one solution: at r_A=r_B=6, M=1, E=1.2 and
    aperture 2.0389 there are two, r_min = 2.585956 and 2.611228, and the
    deeper one -- which the outward scan reached first -- costs 1.82e-6 MORE.

    (2) Comparing costs is not enough, because a sign-change scan does not find
    every solution.  Near an interior extremum of J -> Phi(J) two roots sit on
    opposite sides of the turning value, and when they are closer together than
    the sample spacing the defect has the SAME sign at both bracketing samples:
    the pair is invisible.  Counterexample, same endpoints, aperture
    2.038982616263322: the roots are J = 1.133661398037 and 1.134023989527,
    3.6e-4 apart against a grid step of 2.7e-3 -- 0.13 of one cell -- and the
    scan returned no solution at all.

    The partial cure is to split the J-interval at the interior extrema of Phi
    and bracket on each resulting subinterval, so that a pair straddling an
    extremum is approached from both sides.  That recovers the case above.

    It is NOT a completeness guarantee, and must not be described as one.  The
    node is the vertex of the parabola through three samples: an APPROXIMATE
    extremum, not a solved Phi'=0.  When the target lies between the approximate
    node's angle and the true maximum, both roots stay on the same side of the
    node and are missed again -- at aperture 2.038982626262098, same endpoints,
    the roots J=1.1338407 and 1.1338447 are not found.  What this routine returns
    is the least-cost candidate found by a sampled search, under the cuts
    Jlo=1e-4*Jmax and r>=2.0001*M, and nothing stronger.
    """
    r_hi = min(rA, rB)

    def angle(J):
        rm = rmin_of_J(N, J, r_hi)
        if rm is None:
            return np.nan, None
        return dphi_seg(N, J, rm, rA) + dphi_seg(N, J, rm, rB), rm

    def gap(J):
        return angle(J)[0] - dphi_target

    Jmax = N(r_hi) * r_hi * (1 - 1e-6)
    # Lower cut: kept small but nonzero because r_min -> 2M as J -> 0 and the
    # turning point is lost to the quadrature there.  It bounds the SEARCH, not
    # the class, and is stated rather than silent.
    Jlo = 1e-4 * Jmax
    Js = np.linspace(Jlo, Jmax, nscan)
    phis = np.array([angle(J)[0] for J in Js])
    ok = np.isfinite(phis)

    # Work run by run: the angle is undefined wherever no turning point exists
    # below r_hi, and for the t-branch that is most of the interval.  Inside each
    # maximal finite run, split at the interior extrema of Phi; a root pair
    # straddling an extremum then lands on opposite monotone pieces, where the
    # endpoint values alone reveal it.
    runs, i = [], 0
    while i < len(Js):
        if not ok[i]:
            i += 1; continue
        j = i
        while j + 1 < len(Js) and ok[j + 1]:
            j += 1
        if j > i:
            runs.append((i, j))
        i = j + 1

    sols, seen = [], []
    for i0, i1 in runs:
        knots = [Js[i0]]
        for k in range(i0 + 1, i1):
            d1, d2 = phis[k] - phis[k - 1], phis[k + 1] - phis[k]
            if d1 * d2 < 0:
                denom = d1 - d2
                shift = 0.5 * (d1 + d2) / denom if denom != 0 else 0.0
                shift = max(-0.9, min(0.9, shift))          # stay inside the run
                knots.append(Js[k] + shift * (Js[1] - Js[0]))
        knots.append(Js[i1])
        knots = sorted(set(knots))
        for a, b in zip(knots[:-1], knots[1:]):
            ga, gb = gap(a), gap(b)
            if not (np.isfinite(ga) and np.isfinite(gb)) or ga * gb > 0:
                continue
            try:
                J = brentq(gap, a, b, xtol=1e-14)
            except ValueError:
                continue
            rm = rmin_of_J(N, J, r_hi)
            if rm is None or any(abs(J - j) < 1e-9 for j in seen):
                continue
            seen.append(J)
            sols.append((cost_optical(N, J, rm, rA, rB), J, rm))
    if not sols:
        return None, None
    sols.sort()
    if report and len(sols) > 1:
        print(f"      {len(sols)} solutions with these endpoints; costs "
              + ", ".join(f"{c:.9f}" for c, _, _ in sols)
              + f"  -> taking r_min={sols[0][2]:.6f}")
    return sols[0][1], sols[0][2]

def curva(N, J, rA, rB, phiA):
    """(r,phi) della geodetica da A(rA,phiA) a B, via periasse r_min."""
    rm = rmin_of_J(N, J, min(rA, rB))
    r_in = np.linspace(rA, rm + 1e-7, 300)
    phi_in = phiA + np.array([dphi_seg(N, J, rm, rA)
                              - dphi_seg(N, J, rm, rr)
                              for rr in r_in])
    phi_peri = phiA + dphi_seg(N, J, rm, rA)
    r_out = np.linspace(rm + 1e-7, rB, 300)
    phi_out = phi_peri + np.array([dphi_seg(N, J, rm, rr)
                                   for rr in r_out])
    return (np.concatenate([r_in, r_out]),
            np.concatenate([phi_in, phi_out]), rm)

print("=" * 66)
print("Schwarzschild statico, M=1, E=1.2: brachistocrone t/tau a estremi fissi")
print("=" * 66)

# ---- (a) simmetrico ----
r0, Phi = 6.0, 0.9
print(f"\n(a) SIMMETRICO: A=({r0},-{Phi}), B=({r0},+{Phi})")
res_a = {}
for nome, N in (('tau', n_tau), ('t', n_t)):
    J, rm = solve_bvp(N, r0, r0, 2 * Phi, report=True)
    res_a[nome] = (J, rm)
    print(f"    ramo {nome}: J_opt={J:.4f}, r_min={rm:.4f}")
print(f"    => r_min^t - r_min^tau = "
      f"{res_a['t'][1]-res_a['tau'][1]:+.4f} (t piu' fondo se <0)")

# ---- (b) asimmetrico ----
rA, rB, dphi = 10.0, 6.0, 1.9
print(f"\n(b) ASIMMETRICO: A=({rA},0), B=({rB},{dphi})")
res_b = {}
for nome, N in (('tau', n_tau), ('t', n_t)):
    J, rm = solve_bvp(N, rA, rB, dphi, report=True)
    res_b[nome] = (J, rm)
    print(f"    ramo {nome}: J_opt={J:.4f}, r_min={rm:.4f}")
print(f"    => r_min^t - r_min^tau = "
      f"{res_b['t'][1]-res_b['tau'][1]:+.4f}")

# --------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 6.4))
th = np.linspace(0, 2 * np.pi, 200)
for ax in (a1, a2):
    ax.plot(2 * M * np.cos(th), 2 * M * np.sin(th), 'k-', lw=0.8)  # r=2M
    ax.set_aspect('equal')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

# (a) simmetrico
for nome, col in (('tau', 'C3'), ('t', 'C0')):
    J, rm = res_a[nome]
    r, ph, _ = curva(n_tau if nome == 'tau' else n_t, J, r0, r0, -Phi)
    a1.plot(r * np.cos(ph), r * np.sin(ph), col, lw=1.6,
            label=rf'${nome}$: $r_{{\min}}={rm:.3f}$')
a1.plot(r0 * np.cos(-Phi), r0 * np.sin(-Phi), 'ks', ms=6)
a1.plot(r0 * np.cos(Phi), r0 * np.sin(Phi), 'ks', ms=6, label='fixed A, B')
a1.set_title(f'(a) symmetric endpoints (same $r_0={r0}$, $\\pm\\Phi$):\n'
             '$t$ and $\\tau$ join the SAME A,B, differ in depth')
a1.legend(fontsize=6, loc='lower left')

# (b) asimmetrico
for nome, col in (('tau', 'C3'), ('t', 'C0')):
    J, rm = res_b[nome]
    r, ph, _ = curva(n_tau if nome == 'tau' else n_t, J, rA, rB, 0.0)
    a2.plot(r * np.cos(ph), r * np.sin(ph), col, lw=1.6,
            label=rf'${nome}$: $r_{{\min}}={rm:.3f}$')
a2.plot(rA, 0, 'ks', ms=6)
a2.plot(rB * np.cos(dphi), rB * np.sin(dphi), 'ks', ms=6, label='fixed A, B')
a2.set_title(f'(b) asymmetric endpoints ($r_A={rA}$, $r_B={rB}$):\n'
             'genuine two-point brachistochrones')
a2.legend(fontsize=6, loc='lower left')
savefig(fig, HERE, 'fig_bvp_estremi_fissi')
print("\nFATTO.")
