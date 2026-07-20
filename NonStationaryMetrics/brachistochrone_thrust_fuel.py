# -*- coding: utf-8 -*-
"""
Spinta propria e budget di carburante lungo una brachistocrona.

La forza di rotaia a^mu = Du^mu/dtau e' esattamente l'accelerazione PROPRIA (cio'
che legge un accelerometro a bordo; la gravita' e' gia' nella derivata covariante,
quindi una geodetica ha a=0). Percio' e' direttamente la spinta/massa che un razzo
deve fornire per seguire la brachistocrona:

    |a|  = sqrt(g_ab a^a a^b)            spinta propria per unita' di massa
    Dv   = int |a| dtau                  characteristic velocity (Delta-v)
    frazione carburante = 1 - exp(-Dv/v_e)   (Tsiolkovsky relativistico)

Proprieta' verificate: a.u = 0 (magnetic-type: nessun lavoro sulla particella,
mass-shell u.u=-1 preservata); in un campo statico (W Killing) a.W = 0, quindi la
spinta e' pura STERZATA; nel caso non-stazionario (Vaidya/Thakurta) c'e' in piu'
a.W = -u^a u^b nabla_(a W_b) != 0, la potenza che mantiene Ehat contro il drift.

ATTENZIONE (fisica): la brachistocrona minimizza il TEMPO a Ehat fisso, NON il
carburante. Il percorso a carburante minimo e' la geodetica (caduta libera, a=0,
zero carburante). La brachistocrona baratta carburante per velocita'. Verso
l'orizzonte / la superficie di freezing |a| -> infinito, quindi il carburante
diverge: il tratto profondo e' quello proibitivo.

Esempio: Schwarzschild, ramo tau. (Il caso non-stazionario aggiunge solo la
componente a.W al medesimo schema.)
"""
import numpy as np
import sympy as sp

M, E, J, r0 = 1.0, 1.4, 4.0, 10.0        # geometrizzato G=c=1, M in unita' di massa

r = sp.symbols('r', positive=True)
DE = (E**2 - 1) * r + 2 * M
S = r * (r - 2 * M) * DE * (r**2 * (r - 2 * M) - J**2 * DE)   # sestica ramo tau (a=0)
F = J * DE / sp.sqrt(S)                                        # dphi/dr (shape)
f = 1 - 2 * M / r
Fn = sp.lambdify(r, F, 'numpy'); fn = sp.lambdify(r, f, 'numpy'); Sn = sp.lambdify(r, S, 'numpy')


def uvec(rv):
    """4-velocita' brachistocrona (u^t,u^r,u^phi) da Ehat (rail) + shape + norma."""
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
    """a^mu = u^r du^mu/dr + Gamma^mu_ab u^a u^b ; ritorna (|a|, a.u, u^r)."""
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

print(f"Schwarzschild M={M}, brachistocrona tau: E={E}, J={J}, r: {r0} -> {rmin:.2f}")
print(f"  check a.u = {au_max:.1e}  (deve 0: forza magnetic-type)")
print(f"  |a| propria: min={amags.min():.4f}  max={amags.max():.4f}   [in 1/M]")
print(f"  Delta v / c = int|a|dtau = {dv:.4f}")
print("  frazione carburante (Tsiolkovsky, 1-exp(-Dv/v_e)):")
for ve, nome in [(1.0, "fotonico  v_e=c"), (0.03, "ionico    v_e=0.03c"), (1e-5, "chimico   v_e~3km/s")]:
    print(f"     {nome:22s} -> {1 - np.exp(-dv / ve):.3g}")
print("  confronto: geodetica (caduta libera) a=0 -> Delta v=0 -> zero carburante.")
print(f"  scala fisica: 1/M(M_sole) ~ 2e13 m/s^2 ~ 2e12 g; verso l'orizzonte |a|->inf.")
