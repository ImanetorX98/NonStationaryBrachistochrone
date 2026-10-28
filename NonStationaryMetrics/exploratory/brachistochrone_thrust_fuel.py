# -*- coding: utf-8 -*-
"""
Proper thrust and fuel budget along a brachistochrone.

The rail force a^mu = Du^mu/dtau is exactly the PROPER acceleration (what an
on-board accelerometer reads; gravity is already in the covariant derivative,
so a geodesic has a=0). Hence it is directly the thrust per unit mass that a rocket
must supply to follow the brachistochrone:

    |a|  = sqrt(g_ab a^a a^b)            proper thrust per unit mass
    Dv   = int |a| dtau                  characteristic velocity (Delta-v)
    fuel fraction = 1 - exp(-Dv/v_e)     (relativistic Tsiolkovsky)

Verified properties: a.u = 0 (magnetic type: no work on the particle,
mass shell u.u=-1 preserved); in a static field (W Killing) a.W = 0, so the
thrust is pure STEERING; in the non-stationary case (Vaidya/Thakurta) there is in addition
a.W = -u^a u^b nabla_(a W_b) != 0, the power that maintains Ehat against the drift.

CAUTION (physics): the brachistochrone minimises TIME at fixed Ehat, NOT
fuel. The minimum-fuel path is the geodesic (free fall, a=0,
zero fuel). The brachistochrone trades fuel for speed. Towards
the horizon / the freezing surface |a| -> infinity, so the fuel
diverges: the deep stretch is the prohibitive one.

Example: Schwarzschild, tau branch. (The non-stationary case only adds the
component a.W to the same scheme.)
"""
import numpy as np
import sympy as sp

M, E, J, r0 = 1.0, 1.4, 4.0, 10.0        # geometrised G=c=1, M in mass units

r = sp.symbols('r', positive=True)
DE = (E**2 - 1) * r + 2 * M
S = r * (r - 2 * M) * DE * (r**2 * (r - 2 * M) - J**2 * DE)   # tau-branch sextic (a=0)
F = J * DE / sp.sqrt(S)                                        # dphi/dr (shape)
f = 1 - 2 * M / r
Fn = sp.lambdify(r, F, 'numpy'); fn = sp.lambdify(r, f, 'numpy'); Sn = sp.lambdify(r, S, 'numpy')


def uvec(rv):
    """Brachistochrone 4-velocity (u^t,u^r,u^phi) from Ehat (rail) + shape + normalisation."""
    ff, Ff = fn(rv), Fn(rv)
    ur = -np.sqrt((E**2 - ff) / (1 + ff * rv**2 * Ff**2))     # ingoing
    return np.array([E / ff, ur, Ff * ur])


def christ(rv):
    ff, fp = fn(rv), 2 * M / rv**2                            # f, df/dr
    G = np.zeros((3, 3, 3))
    G[0, 0, 1] = G[0, 1, 0] = fp / (2 * ff)
    G[1, 0, 0] = ff * fp / 2; G[1, 1, 1] = -fp / (2 * ff); G[1, 2, 2] = -rv * ff
    G[2, 1, 2] = G[2, 2, 1] = 1 / rv
    return G


def gmet(rv):
    return np.diag([-fn(rv), 1 / fn(rv), rv**2])


def accel(rv, h=1e-6):
    """a^mu = u^r du^mu/dr + Gamma^mu_ab u^a u^b ;  returns (|a|, a.u, u^r)."""
    u = uvec(rv); du = (uvec(rv + h) - uvec(rv - h)) / (2 * h)
    G = christ(rv)
    a = u[1] * du + np.array([sum(G[m, i, j] * u[i] * u[j] for i in range(3) for j in range(3))
                              for m in range(3)])
    g = gmet(rv)
    return np.sqrt(max(a @ g @ a, 0)), a @ g @ u, u[1]


rs = np.linspace(2.01, r0, 4000)
rmin = rs[(Sn(rs) > 0) & (E**2 - fn(rs) > 0)].min() + 0.05
rg = np.linspace(r0, rmin, 300)
amags = np.array([accel(rv)[0] for rv in rg])
urs = np.array([abs(uvec(rv)[1]) for rv in rg])
dv = np.trapezoid(amags / urs, r0 - rg)                       # Dv/c = int |a| dtau
au_max = max(abs(accel(rv)[1]) for rv in rg[::30])

print(f"Schwarzschild M={M}, tau brachistochrone: E={E}, J={J}, r: {r0} -> {rmin:.2f}")
print(f"  check a.u = {au_max:.1e}  (must be 0: magnetic-type force)")
print(f"  proper |a|: min={amags.min():.4f}  max={amags.max():.4f}   [in 1/M]")
print(f"  Delta v / c = int|a|dtau = {dv:.4f}")
print("  fuel fraction (Tsiolkovsky, 1-exp(-Dv/v_e)):")
for ve, nome in [(1.0, "photonic  v_e=c"), (0.03, "ion       v_e=0.03c"), (1e-5, "chemical  v_e~3km/s")]:
    print(f"     {nome:22s} -> {1 - np.exp(-dv / ve):.3g}")
print("  comparison: geodesic (free fall) a=0 -> Delta v=0 -> zero fuel.")
print(f"  physical scale: 1/M(M_sun) ~ 2e13 m/s^2 ~ 2e12 g; towards the horizon |a|->inf.")
