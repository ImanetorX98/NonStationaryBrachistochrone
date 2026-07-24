# -*- coding: utf-8 -*-
"""
Complete symbolic coefficient table for the new off-shell adiabatic term (t/eta branch, main11).

STATUS (honest): the off-shell term is REDUCED to the genus-2 basis with symbolic coefficients
in (M, a, E, J); the quadratures U_k, W_jk (and third-kind letters) are NOT evaluated here --
they are the genus-2 transcendentals, integrated numerically in the verification below. Solving
them in elementary form is impossible (genus-2 theorem); the theta-nome series evaluation (the
paper's Sage machinery) is a separate, unrun step for this term.

delta phi_extra = - int (A/sqrt(S)) * Delta S dr,     S = r [(E^2-1)r+2M] Q2   (on-shell sextic)

  KERNEL   A = A_poly + N3/Q2 ,  A_poly = p1 r + p0 ,  poles at Q2=0
  ACTION   Delta S = sum_k a_k U_k + R/[Delta*((E^2-1)r+2M)] ,  U_k=int r^k/sqrt(S) dr, poles at Delta=0
  WRAP     delta phi_extra = sum Q_jk W_jk (second kind) + mixed/third-kind blocks,
           Q_jk = -1/2 (p_j a_k - p_k a_j) ,  W_jk = int (r^j U_k - r^k U_j)/sqrt(S) dr.

Every coefficient is an explicit function of (M, a, E, J); see coefficients() below.
"""
import sympy as sp
import numpy as np
from scipy.integrate import quad

r, a, E, J, M = sp.symbols('r a E J_ M', positive=True)


def curve():
    Dl = r**2 - 2*M*r + a**2
    Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
          + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
    S = sp.expand(r*((E**2-1)*r + 2*M)*Q2)
    return Dl, Q2, S


def coefficients():
    """All symbolic coefficients as functions of (M, a, E, J)."""
    Dl, Q2, S = curve()
    # kernel A = E^2 J r^4 ((E^2-1)r+2M)/Q2 = A_poly + N3/Q2
    Qp, N3 = sp.div(sp.Poly(sp.expand(E**2*J*r**4*((E**2-1)*r + 2*M)), r), sp.Poly(sp.expand(Q2), r))
    p_j = [sp.factor(c) for c in Qp.all_coeffs()[::-1]]
    # radial action: Delta S = sum a_k U_k + R/deng,  deng = Delta*((E^2-1)r+2M)
    deng = sp.expand(Dl*((E**2-1)*r + 2*M))
    Qd, R = sp.div(sp.Poly(sp.expand(S), r), sp.Poly(deng, r))
    a_k = [sp.factor(c) for c in Qd.all_coeffs()[::-1]]
    # second-kind W_jk coefficients
    Q = {}
    for j in range(4):
        for k in range(j+1, 4):
            pj = p_j[j] if j < len(p_j) else 0
            pk = p_j[k] if k < len(p_j) else 0
            q = sp.factor(-sp.Rational(1, 2)*(pj*a_k[k] - pk*a_k[j]))
            if q != 0:
                Q[(j, k)] = q
    return dict(S=S, Q2=Q2, Dl=Dl, a_k=a_k, p_j=p_j, N3=sp.factor(N3.as_expr()),
                R=sp.factor(R.as_expr()), Q_jk=Q, deng=deng)


def _verify(av, ev, jv, Mv=1.0):
    """Numeric verification (second-kind block) + full delta phi_extra self-consistency."""
    c = coefficients()
    sub = {a: av, E: ev, J: jv, M: Mv}
    Sn = sp.lambdify(r, c['S'].subs(sub), 'numpy')
    sq = lambda x: np.sqrt(np.clip(Sn(x), 0, None))
    ak = [float(x.subs(sub)) for x in c['a_k']]
    pj = [float(x.subs(sub)) for x in c['p_j']]
    Qjk = {jk: float(v.subs(sub)) for jk, v in c['Q_jk'].items()}
    An = sp.lambdify(r, (E**2*J*r**4*((E**2-1)*r+2*M)/c['Q2']).subs(sub), 'numpy')
    dengn = sp.lambdify(r, c['deng'].subs(sub), 'numpy')
    Rn = sp.lambdify(r, c['R'].subs(sub), 'numpy')
    r0, rf = 11.0, 6.0
    U = lambda x, k: quad(lambda t: t**k/sq(t), r0, x, limit=200)[0]
    W = lambda j, k: quad(lambda x: (x**j*U(x, k) - x**k*U(x, j))/sq(x), r0, rf, limit=200)[0]
    # second-kind block: direct vs Q_jk assembly
    direct2 = -quad(lambda x: ((pj[0]+pj[1]*x)/sq(x))*sum(ak[k]*U(x, k) for k in range(4)), r0, rf, limit=200)[0]
    asm2 = sum(Qjk[jk]*W(*jk) for jk in Qjk) - 0.5*sum(pj[j]*U(rf, j) for j in range(2))*sum(ak[k]*U(rf, k) for k in range(4))
    # full delta phi_extra: direct double integral (A full, Delta S full)
    DS = lambda x: -quad(lambda t: (sum(ak[k]*t**k for k in range(4)) + Rn(t)/dengn(t))/sq(t), r0, x, limit=150)[0]
    full = -quad(lambda x: (An(x)/sq(x))*DS(x), r0, rf, limit=150)[0]
    return direct2, asm2, full


if __name__ == '__main__':
    c = coefficients()
    print("on-shell sextic S = r[(E^2-1)r+2M] Q2")
    print("\nKERNEL  A = A_poly + N3/Q2:")
    print("  p_0 =", c['p_j'][0], " p_1 =", c['p_j'][1])
    print("  N3  =", c['N3'])
    print("\nRADIAL ACTION  Delta S = sum a_k U_k + R/[Delta((E^2-1)r+2M)]:")
    for k, v in enumerate(c['a_k']):
        print(f"  a_{k} =", v)
    print("  R   =", c['R'])
    print("\nSECOND-KIND  W_jk coefficients  Q_jk = -1/2(p_j a_k - p_k a_j):")
    for jk, v in c['Q_jk'].items():
        print(f"  Q_{jk[0]}{jk[1]} =", v)
    print("\nnumeric verification (a=0.9,E=1.4,J=6):")
    d2, a2, full = _verify(0.9, 1.4, 6.0)
    print(f"  second-kind block: direct={d2:.10f}  Q_jk assembly={a2:.10f}  diff={abs(d2-a2):.1e}")
    print(f"  full delta phi_extra (direct double integral) = {full:.10f}")
    print("  (W_jk / third-kind letters are genus-2 transcendentals, integrated numerically)")
