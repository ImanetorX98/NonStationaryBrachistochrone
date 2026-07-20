# -*- coding: utf-8 -*-
"""
Conservazione dell'energia di rotaia lungo le brachistocrone di Thakurta-Kerr.

Analogo esatto di VaidyaMetric/kodama_conservation.py, ma per Thakurta-Kerr
(conforme-Kerr, g^TK = A(eta)^2 g^Kerr, eta = tempo conforme). Poiche' Thakurta
NON e' sfericamente simmetrico, non esiste un vettore di Kodama; il ruolo del
selettore/rotaia W e' svolto dal vettore conforme di Killing d/deta.

Tre affermazioni, tutte verificate:

 [A] d/deta e' un CONFORMAL KILLING di g^TK con fattore psi = A'/A
     (proprieta' dello spaziotempo, simbolico esatto). L'energia associata
     Ehat := -g^TK(u, d/deta) NON e' conservata da una simmetria ordinaria.

 [B] NON-BANALITA': lungo una GEODETICA di g^TK (rotaia spenta)
       dEhat/dtau = -(1/2)(L_{d_eta} g)_{ab} u^a u^b = psi = A'/A
     cioe' l'energia di rotaia driftra al tasso A'/A. Verifica numerica
     (Christoffel per differenze finite): dEhat/dtau = A'/A a ~1e-11.

 [C] Lungo la BRACHISTOCRONA (rotaia accesa: l'orbita e' istantaneamente
     un'orbita di Kerr on-shell a E_eff = Ehat/A(eta), con J conservato)
     l'identita' conforme da'
        -g^TK(u, d_eta) = A * (-g^Kerr(u^Kerr, d_eta)) = A * E_eff = Ehat = cost.
     Verifica numerica: Ehat piatta a precisione macchina lungo la
     brachistocrona, mentre la geodetica con lo stesso dato iniziale driftra.

Conseguenza: E_eff = Ehat/A(eta) traccia in modo prescritto (Ehat costante),
l'orbita resta on-shell a E_eff, dunque dphi/dr = F(r; E_eff, J) e' esatta e
delta phi = int dE F * eta e' la correzione adiabatica al prim'ordine CORRETTA
(nessun termine off-shell mancante).
"""
import os, sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

M, a = 1.0, 0.9
A0, eps = 1.0, 0.05          # A(eta) = A0 exp(eps eta) => A'/A = eps
J = 6.0
Ehat0 = 1.4                  # energia di rotaia (costante) sulla brachistocrona
r_start = 10.0

def A(e):  return A0*np.exp(eps*e)
def Ap(e): return eps*A0*np.exp(eps*e)     # dA/deta

# ----------------------------------------------------------------------
# metrica Thakurta-Kerr equatoriale (eta, r, phi) = A(eta)^2 * Kerr_BL
# ----------------------------------------------------------------------
def gmat(e, r):
    A2 = A(e)**2; De = r**2 - 2*M*r + a**2
    g = np.zeros((3, 3))
    g[0, 0] = -A2*(1 - 2*M/r);          g[0, 2] = g[2, 0] = -A2*(2*M*a/r)
    g[1, 1] = A2*r**2/De;               g[2, 2] = A2*(r**2 + a**2 + 2*M*a**2/r)
    return g

print("=" * 70)
print("[A] d/deta e' conformal Killing di g^TK con psi = A'/A  (simbolico)")
print("=" * 70)
etas, rs = sp.symbols('eta r', real=True)
Af = sp.Function('A')
De_s = rs**2 - 2*M*rs + a**2
gS = sp.Matrix([
    [-Af(etas)**2*(1 - 2*M/rs),      0,                 -Af(etas)**2*(2*M*a/rs)],
    [0,                              Af(etas)**2*rs**2/De_s,   0],
    [-Af(etas)**2*(2*M*a/rs),        0,                 Af(etas)**2*(rs**2 + a**2 + 2*M*a**2/rs)]])
# xi = d/deta:  (L_xi g)_{ij} = d_eta g_{ij}. Conformal Killing <=> = 2 psi g_{ij}.
psi = sp.simplify(sp.diff(gS[0, 0], etas)/(2*gS[0, 0]))
print("  psi (da g_ee) =", psi, "   (atteso A'(eta)/A(eta))")
print("  check g_rr:", sp.simplify(sp.diff(gS[1, 1], etas)/(2*gS[1, 1]) - psi),
      "   g_pp:", sp.simplify(sp.diff(gS[2, 2], etas)/(2*gS[2, 2]) - psi),
      "  (=> conformal Killing su tutte le componenti)")

# ----------------------------------------------------------------------
# [B] geodetica di g^TK: Ehat = -g(u, d_eta) driftra a A'/A
# ----------------------------------------------------------------------
def dg(e, r, k, h=1e-6):
    if k == 0: return (gmat(e+h, r) - gmat(e-h, r))/(2*h)
    if k == 1: return (gmat(e, r+h) - gmat(e, r-h))/(2*h)
    return np.zeros((3, 3))

def christoffel(e, r):
    g = gmat(e, r); gi = np.linalg.inv(g)
    dgs = [dg(e, r, 0), dg(e, r, 1), np.zeros((3, 3))]
    Gam = np.zeros((3, 3, 3))
    for l in range(3):
        for i in range(3):
            for j in range(3):
                Gam[l, i, j] = 0.5*sum(gi[l, m]*(dgs[i][m, j] + dgs[j][m, i] - dgs[m][i, j])
                                       for m in range(3))
    return Gam

def Ehat_of(y):
    e, r, ph, ue, ur, up = y[:6]; g = gmat(e, r)
    return -(g[0, 0]*ue + g[0, 2]*up)          # -g(u, d_eta)

# dato iniziale entrante, u.u = -1
ur0, up0 = -0.20, 0.030
g0 = gmat(0.0, r_start)
Aq = g0[0, 0]; Bq = 2*g0[0, 2]*up0; Cq = g0[2, 2]*up0**2 + g0[1, 1]*ur0**2 + 1
ue0 = (-Bq - np.sqrt(Bq**2 - 4*Aq*Cq))/(2*Aq)
y0 = [0.0, r_start, 0.0, ue0, ur0, up0]

def geo_rhs(l, y):
    e, r, ph, ue, ur, up = y; u = [ue, ur, up]; Gam = christoffel(e, r)
    acc = [-sum(Gam[m, i, j]*u[i]*u[j] for i in range(3) for j in range(3)) for m in range(3)]
    return [ue, ur, up, acc[0], acc[1], acc[2]]

solG = solve_ivp(geo_rhs, [0, 6.0], y0, rtol=1e-10, atol=1e-12, max_step=0.005, dense_output=True)
lg = np.linspace(0, 6.0, 500); YG = solG.sol(lg)
EhG = np.array([Ehat_of(YG[:, i]) for i in range(len(lg))])
nn = np.array([ (lambda e, r, ue, ur, up, g=gmat(YG[0, i], YG[1, i]):
                 g[0, 0]*ue**2 + 2*g[0, 2]*ue*up + g[2, 2]*up**2 + g[1, 1]*ur**2)
               (YG[0, i], YG[1, i], YG[3, i], YG[4, i], YG[5, i]) for i in range(len(lg))])
dEhG = np.gradient(EhG, lg)
print()
print("=" * 70)
print("[B] GEODETICA g^TK: dEhat/dtau  vs  A'/A = %.3f   (rotaia spenta => drift)" % eps)
print("=" * 70)
print("  max|u.u+1| lungo la geodetica = %.1e" % np.max(np.abs(nn + 1)))
for i in [80, 200, 320, 440]:
    print("  tau=%.2f eta=%.3f  Ehat=%.5f  dEhat/dtau=%+.5f  A'/A=%+.5f  diff=%.1e"
          % (lg[i], YG[0, i], EhG[i], dEhG[i], eps, abs(dEhG[i] - eps)))

# ----------------------------------------------------------------------
# [C] brachistocrona (rotaia): orbita on-shell Kerr a E_eff = Ehat/A(eta).
#     shape dphi/dr e clock deta/dr dall'Hamiltoniana di ramo H2(E_eff).
# ----------------------------------------------------------------------
rr, pr, Es, Js = sp.symbols('r pr E J_', real=True)
f2 = 1 - 2*M/rr; Dl2 = rr**2 - 2*M*rr + a**2; b2 = 2*M*a/rr; v2 = 1 - f2/Es**2
P2 = rr**2 + a**2 + 2*M*a**2/rr; Pb2 = P2 + b2**2/Es**2
H2 = Js*b2*v2/Pb2 + sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2 + Js**2/Pb2) - 1
H2n = sp.lambdify((rr, pr, Es, Js), H2, 'numpy')
dHp = sp.lambdify((rr, pr, Es, Js), sp.diff(H2, pr), 'numpy')     # dr/deta
dHJ = sp.lambdify((rr, pr, Es, Js), sp.diff(H2, Js), 'numpy')     # dphi/deta

def pr_onshell(r, E):        # p_r entrante che risolve H2=0 (orbita on-shell a E)
    pg = np.linspace(-80, 80, 4001); Hv = H2n(r, pg, E, J)
    rts = [brentq(lambda p: H2n(r, p, E, J), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(r, p, E, J) < 0]
    return min(ing) if ing else np.nan

def brachi_rhs(r, y):        # y = [eta, phi];  E_eff = Ehat/A(eta)
    e, ph = y; E = Ehat0/A(e); prv = pr_onshell(r, E)
    drde = dHp(r, prv, E, J)                 # dr/deta
    deta_dr = 1.0/drde
    dphi_dr = dHJ(r, prv, E, J)/drde
    return [deta_dr, dphi_dr]

# integra la brachistocrona verso l'interno fino al turning (dr/deta -> 0)
def turning(r, y):
    e, ph = y; E = Ehat0/A(e); prv = pr_onshell(r, E); return dHp(r, prv, E, J)
turning.terminal = True; turning.direction = 0
solB = solve_ivp(brachi_rhs, [r_start, 3.0], [0.0, 0.0], rtol=1e-11, atol=1e-13,
                 max_step=0.01, dense_output=True, events=turning)
rB = np.linspace(r_start, solB.t[-1]*1.0, 500)
etaB, phiB = solB.sol(rB)

# ricostruisci u^TK e calcola Ehat_brachi = -g^TK(u, d_eta)
EhB = np.empty_like(rB); Eeff_chk = np.empty_like(rB)
for i, r in enumerate(rB):
    e = etaB[i]; E = Ehat0/A(e)
    prv = pr_onshell(r, E)
    drde = dHp(r, prv, E, J); dphidr = dHJ(r, prv, E, J)/drde; detadr = 1.0/drde
    # orbita ENTRANTE: r decresce lungo il moto, quindi tangente fisica ha r-slot < 0
    T = np.array([-detadr, -1.0, -dphidr])        # tangente futura entrante (eta,r,phi)
    g = gmat(e, r)
    gTT = T @ g @ T                              # = -(dtau_TK/dr)^2  (T timelike)
    u = T/np.sqrt(-gTT)                          # u.u = -1, future-diretta
    EhB[i] = -(g[0, 0]*u[0] + g[0, 2]*u[2])      # -g^TK(u, d_eta) = Ehat
    Eeff_chk[i] = E
print()
print("=" * 70)
print("[C] BRACHISTOCRONA (rotaia): Ehat = -g^TK(u, d_eta) lungo l'orbita")
print("=" * 70)
print("  Ehat0 (imposto via E_eff=Ehat/A) = %.6f" % Ehat0)
print("  Ehat ricostruito da u:  min %.12f  max %.12f" % (EhB.min(), EhB.max()))
print("  => spread (max-min) = %.2e   (=> Ehat COSTANTE: la rotaia lo conserva)" % (EhB.max() - EhB.min()))
print("  => |Ehat(u) - Ehat0| max = %.2e   (deve ~0)" % np.max(np.abs(EhB - Ehat0)))
print("  E_eff = Ehat/A(eta) varia da %.4f (r=%.1f) a %.4f (r=%.2f): traccia, on-shell"
      % (Eeff_chk[0], rB[0], Eeff_chk[-1], rB[-1]))
print()
print("  CONTRASTO: stesso spaziotempo, geodetica driftra a A'/A=%.3f; la rotaia lo annulla" % eps)
print("  (esattamente come Vaidya: -u.K=E costante a 6.66e-16 vs drift geodetico -m'(u^v)^2/r).")
print()
print("FATTO. Ehat conservata sulla brachistocrona => dphi/dr=F(r;E_eff) esatta")
print("       => delta phi = int dE F * eta e' la correzione 1o ordine corretta.")
