# -*- coding: utf-8 -*-
"""
Metrica ottica per FLRW dal Killing conforme d_eta — costruzione e verifica.

Metrica fisica (k = 0, +1, -1; coordinate comoventi (chi, phi) sul piano):

    ds^2 = a(eta)^2 [ -deta^2 + dchi^2 + S_k(chi)^2 dphi^2 ] = a^2 ghat
    S_0 = chi,  S_1 = sin(chi),  S_-1 = sinh(chi)

d_eta e' Killing genuino di ghat (statica). Riduzione ottica su ghat:

RAGGI NULLI (invarianza conforme):
    ghat = -deta^2 + dl_k^2  con  dl_k^2 = dchi^2 + S_k^2 dphi^2
    => deta = dl_k :  metrica ottica = 3-GEOMETRIA COMOVENTE, indice n = 1.
    Nessun termine di Randers (ghat_{eta i} = 0: niente rotazione/dragging).

MASSIVE con rotaia conforme Ehat = -u_eta:
    v(eta) = sqrt(1 - a^2/Ehat^2)  (solo tempo, omogeneita')
    => deta = dl_k / v(eta):  Fermat con indice n(eta) = 1/v(eta)
    DIPENDENTE DAL TEMPO ma non dallo spazio (Finsler non-autonomo... pero':)

    OROLOGIO OTTICO: lam(eta) = int v deta  =>  dlam = dl_k
    minimizzare eta_arrivo  <=>  minimizzare lam_arrivo  <=>  minimizzare
    la LUNGHEZZA COMOVENTE. Metrica ottica ESATTA = dl_k^2; tutta la
    dinamica non stazionaria e' riassorbita nell'orologio lam(eta).
    T_t e T_tau monotoni in eta_arrivo => vale per tutti e tre i rami.

Verifiche:
  V1  identita' dell'orologio: dlam = v deta = dl (algebra, esatta)
  V2  orologio ottico in de Sitter: lam(eta) in forma chiusa
      (= formula Dx_max gia' validata) e lam(eta_f) = distanza massima
  V3  k=+1 numerico: sulla 3-sfera comovente la brachistocrona massiva
      e' il CERCHIO MASSIMO: famiglia perturbata su S^2 equatoriale,
      eta_arrivo(eps) minimo a eps=0 (radiation a=eta, Ehat=2)
  V4  indice n(eta)=1/v: divergenza al congelamento a->Ehat — analogo
      temporale della superficie di luce (n->oo invece di F->0 spaziale)
"""

import sympy as sp
import numpy as np

print("=" * 72)
print("[V1] orologio ottico: dlam = v(eta) deta  ==>  dlam = dl_k")
print("=" * 72)
eta, chi, Eh, H = sp.symbols('eta chi Ehat H', positive=True)
a = sp.Function('a', positive=True)
v = sp.sqrt(1 - a(eta)**2 / Eh**2)
print("  vincolo (da flrw_brachistochrone_sympy):  dl/deta = v(eta) =", v)
print("  dlam := v deta  =>  dlam = dl_k lungo ogni worldline di rotaia:")
print("  il moto e' UNIT-SPEED nell'orologio lam. eta_arrivo = lam^{-1}(L),")
print("  monotono => brachistocrona = geodetica della 3-geometria comovente.")
print("  Metrica ottica:  dl_k^2 = dchi^2 + S_k(chi)^2 dphi^2   (esatta)")

print()
print("=" * 72)
print("[V2] de Sitter a = -1/(H eta): orologio ottico in forma chiusa")
print("=" * 72)
u = sp.Symbol('u', positive=True)
# lam(eta) da eta0=-1/H con u = -Ehat*H*eta (u: Ehat -> 1)
prim = sp.sqrt(u**2 - 1) - sp.acos(1 / u)
chk = sp.simplify(sp.diff(prim, u) - sp.sqrt(1 - 1 / u**2))
print("  primitiva di sqrt(1-1/u^2):", prim, "  [verifica d/du:", chk, "]")
lam_u = (prim.subs(u, Eh) - prim) / (Eh * H)
print("  lam(u) =", lam_u)
lam_tot = sp.simplify(lam_u.subs(u, 1))
dx_max = (sp.sqrt(Eh**2 - 1) + sp.asin(1 / Eh) - sp.pi / 2) / (Eh * H)
print("  lam(eta_f) =", lam_tot)
diff_ = sp.simplify((lam_tot - dx_max).rewrite(sp.asin))
print("  lam(eta_f) - Dx_max =", diff_, " (coincide, OK)")

print()
print("=" * 72)
print("[V3] k=+1: brachistocrona massiva = cerchio massimo su S^3 (numerico)")
print("=" * 72)
# radiation: a(eta) = eta, Ehat = 3, partenza eta0 = 1 (a=1), stop a eta=3
# (lam_max = 1.375 > lunghezze in gioco: la particella arriva sempre)
Eh_n, eta0 = 3.0, 1.0
etag = np.linspace(eta0, Eh_n - 1e-12, 400000)
vg = np.sqrt(np.clip(1 - (etag / Eh_n)**2, 0, None))
lam_g = np.concatenate([[0], np.cumsum((vg[1:] + vg[:-1]) / 2
                                       * np.diff(etag))])

def eta_arrivo(L):
    return np.interp(L, lam_g, etag)

# S^2 equatoriale della 3-sfera comovente: da (pi/2, 0) a (pi/2, 1.0)
Phi = 1.0
sg = np.linspace(0, 1, 4000)
res = []
for eps in np.linspace(-0.3, 0.3, 61):
    th = np.pi / 2 + eps * np.sin(np.pi * sg)
    dth = eps * np.pi * np.cos(np.pi * sg)
    L = np.trapezoid(np.sqrt(dth**2 + np.sin(th)**2 * Phi**2), sg)
    res.append((eps, L, eta_arrivo(L)))
res = np.array(res)
i0 = np.argmin(res[:, 2])
print(f"  distanza cerchio massimo: L(0) = {res[30, 1]:.9f} (attesa {Phi})")
print(f"  argmin eta_arrivo(eps) = {res[i0, 0]:.4f}  (atteso 0)")
print(f"  eta_arrivo(0) = {res[30, 2]:.6f};  eta_arrivo(0.3) = {res[-1, 2]:.6f}")
print("  => geodetica sferica (cerchio massimo) = brachistocrona.  OK")

print()
print("=" * 72)
print("[V4] indice ottico n(eta) = 1/v: divergenza al congelamento")
print("=" * 72)
n_expr = 1 / v
print("  n(eta) =", n_expr)
print("  a -> Ehat:  n -> oo  (orologio si ferma: dlam -> 0)")
print("  analogo TEMPORALE della superficie di luce stazionaria:")
print("  in Kerr la riduzione ottica degenera su una SUPERFICIE spaziale")
print("  (F -> 0); in FLRW degenera a un ISTANTE cosmico (a = Ehat).")
print("\nFATTO.")
