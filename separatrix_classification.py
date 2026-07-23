# -*- coding: utf-8 -*-
"""
Classificazione delle degenerazioni della sestica: locus algebrico vs separatrice fisica.

Referee obiezioni #7, #8 (CRITICHE). Il paper chiamava "separatrice" un locus di
radice doppia della curva che, ai parametri dichiarati, cade a raggio NEGATIVO
(r_d<0): e' una degenerazione algebrica (genus 2->1), NON la separatrice fisica
(orbita circolare instabile, confine cattura/scattering), che richiede una radice
doppia ESTERNA r_d>r_+ dove due turning point reali si fondono.

Distinzione da mantenere nel paper:
  J_deg      : locus del discriminante Res_r(curva, d_r curva)=0 (genus-degeneration).
               r_d puo' essere <0, complesso o interno: il closed form (Weierstrass)
               resta valido, ma NON e' una separatrice fisica.
  J_cap(v0)  : soglia dinamica di penetrazione (Sec. 4.2), oggetto FISICO distinto.
  J_c^phys   : radice doppia ESTERNA (r_d>r_+) = separatrice fisica (merger di due
               turning point reali).

Uscita: per ogni ramo, TUTTE le radici doppie reali classificate; per il ramo t
di Thakurta-Kerr la separatrice fisica esterna e i suoi b_i (sorgente corretta
E d_E + J d_J, cfr. obiezione #5).
"""
import sympy as sp
import numpy as np
import mpmath as mp

r, J, Es = sp.symbols('r J Es', real=True)
M = 1
mp.mp.dps = 25


def classify(curve, aval, label):
    """Res_r(curve, d_r curve)=0 -> tutte le J con radice doppia; classifica r_d."""
    rp = 1 + np.sqrt(1 - aval**2) if aval < 1 else 2.0     # orizzonte esterno
    Res = sp.resultant(curve, sp.diff(curve, r), r)
    print(f"\n===== {label}  (r_+={rp:.4f}) =====")
    rows = []
    for js in sp.solve(sp.Eq(Res, 0), J):
        jc = complex(js)
        if abs(jc.imag) > 1e-8 or abs(jc.real) < 1e-6:
            continue
        jv = jc.real
        rts = np.roots([complex(c) for c in sp.Poly(curve.subs(J, jv), r).all_coeffs()])
        rd = None
        for i in range(len(rts)):
            for k in range(i + 1, len(rts)):
                if abs(rts[i] - rts[k]) < 1e-4 and abs(rts[i].imag) < 1e-6:
                    rd = (rts[i].real + rts[k].real) / 2
        if rd is None:
            continue
        # separatrice fisica = merger di due turning point reali esterni
        real_ext = sorted([z.real for z in rts if abs(z.imag) < 1e-6 and z.real > rp])
        phys = rd > rp and any(abs(x - rd) < 1e-3 for x in real_ext)
        region = "ESTERNO r>r_+" if rd > rp else ("interno" if rd > 0 else "NEGATIVO")
        rows.append((jv, rd, region, phys))
    for jv, rd, region, phys in sorted(rows):
        tag = "  <== SEPARATRICE FISICA" if phys else "  (degenerazione algebrica J_deg)"
        print(f"  J = {jv:+.5f}   r_d = {rd:+.5f}   [{region}]{tag}")
    return rows


DE = (Es**2 - 1) * r + 2 * M
Dl = lambda av: r**2 - 2 * M * r + av**2

# ---- #7: Schwarzschild tau (a=0), E=1.4 ----
S_schw = sp.expand(r * (r - 2 * M) * DE * (r * Dl(0) - J**2 * DE)).subs(Es, sp.Rational(7, 5))
classify(S_schw, 0.0, "#7 Schwarzschild tau, E=1.4")

# ---- #8: Thakurta-Kerr t-branch, a=0.9, E=1.2 ----
a = sp.Rational(9, 10)
Q2 = (2*Es**2*J**2*M*r - Es**2*J**2*r**2 - 4*Es**2*J*M*a*r + 2*Es**2*M*a**2*r + Es**2*a**2*r**2
      + Es**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6 = sp.expand(r * Q2 * DE).subs(Es, sp.Rational(6, 5))
classify(R6, 0.9, "#8 Thakurta-Kerr t-branch, a=0.9, E=1.2")

# ---- separatrice FISICA del ramo t: raffina (J_c, r_d) esterni e calcola b_i corretti ----
print("\n=== ramo t: SEPARATRICE FISICA esterna, b_i con sorgente (E d_E + J d_J) ===")
R6p = sp.lambdify((r, J), R6, 'mpmath'); R6r = sp.lambdify((r, J), sp.diff(R6, r), 'mpmath')
sol = mp.findroot(lambda rd, jc: [R6p(rd, jc), R6r(rd, jc)], (mp.mpf('3.514'), mp.mpf('-8.0535')))
rd_v, Jc_v = float(sol[0]), float(sol[1])
print(f"  raffinato: J_c^phys = {Jc_v:.8f}   r_d = {rd_v:.8f}  (esterno, r_+=1.4359)")

Dl9 = r**2 - 2 * M * r + sp.Rational(81, 100)
Q2s = (2*Es**2*J**2*M*r - Es**2*J**2*r**2 - 4*Es**2*J*M*a*r + 2*Es**2*M*a**2*r + Es**2*a**2*r**2
       + Es**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6s = sp.expand(r * Q2s * ((Es**2 - 1) * r + 2 * M))
Kts = r * ((Es**2 - 1) * r + 2 * M) * (J * (r - 2 * M) + 2 * M * a) / Dl9
Fs = Kts / sp.sqrt(R6s)
Ntot = sp.expand(sp.cancel((Es * sp.diff(Fs, Es) + J * sp.diff(Fs, J)) * R6s**sp.Rational(3, 2)).subs(Es, sp.Rational(6, 5)))
S6 = R6s.subs(Es, sp.Rational(6, 5))


def bcoeffs(S, N, rd):
    d = lambda ex, n: sp.diff(ex, r, n).subs(r, rd)
    Q4rd, Q4prd, Q4pprd = d(S, 2) / 2, d(S, 3) / 6, d(S, 4) / 12
    s = sp.sqrt(Q4rd); a1 = Q4prd / (4 * s); a2 = Q4pprd / 12
    Nr, Np, Npp = N.subs(r, rd), sp.diff(N, r).subs(r, rd), sp.diff(N, r, 2).subs(r, rd)
    h0 = Nr / Q4rd; Fp = (Np * Q4rd - Nr * Q4prd) / Q4rd**2
    Fpp = ((Npp * Q4rd - Nr * Q4pprd) * Q4rd - 2 * Q4prd * (Np * Q4rd - Nr * Q4prd)) / Q4rd**3
    h1 = Fp * s; h2 = sp.Rational(1, 2) * (Fpp * s**2 + Fp * (Q4prd / 2))
    return {'b3': h0 / s**3, 'b2': (h1 - 3 * a1 * h0) / s**3,
            'b1': (h2 - 3 * a1 * h1 + (6 * a1**2 - 3 * a2) * h0) / s**3}


bc = bcoeffs(S6.subs(J, Jc_v), Ntot.subs(J, Jc_v), rd_v)
for k in ['b1', 'b2', 'b3']:
    print(f"  {k} = {float(bc[k]):+.6f}")
print("\nEsito: la separatrice FISICA del ramo t e' J_c^phys=-8.054, r_d=+3.514 (esterna),")
print("       non J=+19.089/r_d=-6.62 (degenerazione algebrica). Schwarzschild tau (E=1.4)")
print("       NON ha separatrice fisica esterna: il suo r_d=-3.36 e' solo J_deg.")
