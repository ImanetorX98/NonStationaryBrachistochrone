# -*- coding: utf-8 -*-
# The conformal-time (eta) brachistochrone in Thakurta-Kerr coincides, as a SPATIAL CURVE,
# with the coordinate-time (t) brachistochrone -- frozen and (by the monotone argument) adiabatic.
#
# Analytic basis:  dt = A(eta) d.eta  =>  F_t = A F_eta  =>  at frozen A the geodesics of F_eta and
# A*F_eta coincide (constant Finsler rescaling)  =>  same brachistochrone curve, with the conserved
# charge relabelled by the conformal factor  J_t = A J_eta.
# Under a drift, t_f = int_0^{eta_f} A d.eta is strictly increasing in eta_f (A>0), so argmin over
# fixed-endpoint paths coincides: eta- and t-brachistochrones are the same spatial family.
#
# This script verifies the frozen equivalence numerically: the eta-orbit at J_eta and the t-orbit at
# the matched charge J_t = A J_eta trace the SAME curve phi(r) to ODE precision, whereas the naive
# same-J comparison does not (different shells).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

M, a, Ehat, A = 1.0, 0.9, 1.4, 1.3     # frozen conformal factor A (const)
rr, pr, pp = sp.symbols('r pr pphi')
Dl = rr**2 - 2*M*rr + a**2; f = 1 - 2*M/rr; vb2 = 1 - A**2*f/Ehat**2
P = rr**2 + a**2 + 2*M*a**2/rr; Pb = P + A**2*(2*M*a/rr)**2/Ehat**2
phi0 = (2*M*a/rr)*vb2/Pb; R = sp.sqrt(Dl*vb2/Pb); kin = sp.sqrt((Dl/rr**2)*pr**2 + pp**2/Pb)
H_eta = pp*phi0 + R*kin - 1          # eq:Heta  (unit conformal-time cost)
H_t   = pp*phi0 + R*kin - A          # eq:Ht-tk (coordinate-time cost)

def mk(H):
    return (sp.lambdify((rr, pr, pp), sp.diff(H, pr), 'numpy'),
            sp.lambdify((rr, pr, pp), sp.diff(H, pp), 'numpy'),
            sp.lambdify((rr, pr, pp), H, 'numpy'))
HpE, HphE, HnE = mk(H_eta)
HpT, HphT, HnT = mk(H_t)
r0, r1 = 12.0, 7.0

def curve(Hp, Hph, Hn, J):
    def dphidr(r, y):
        pg = np.linspace(-80, 80, 2001); Hv = Hn(r, pg, J)
        rts = [brentq(lambda p: Hn(r, p, J), pg[i], pg[i+1]) for i in range(len(pg)-1)
               if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
        ing = [p for p in rts if Hp(r, p, J) < 0]
        return [Hph(r, min(ing), J)/Hp(r, min(ing), J)] if ing else [np.nan]
    rs = np.linspace(r0, r1, 300)
    so = solve_ivp(dphidr, [r0, r1], [0.0], t_eval=rs, rtol=1e-9, atol=1e-11)
    return so.t, so.y[0]

J_eta = 6.0
rE, phE = curve(HpE, HphE, HnE, J_eta)

def mismatch(Jt):
    rT, phT = curve(HpT, HphT, HnT, Jt); n = min(len(phE), len(phT))
    return np.nanmax(np.abs(phE[:n] - phT[:n]))

res = minimize_scalar(mismatch, bounds=(3.0, 9.0), method='bounded', options={'xatol': 1e-6})
print("eta-brachistochrone J_eta = %.4f" % J_eta)
print("best-matching t-brachistochrone J_t = %.5f   (A*J_eta = %.5f)" % (res.x, A*J_eta))
print("max|phi_eta(r) - phi_t(r)| at matched charge = %.2e  [same curve]" % res.fun)
print("naive same-J (J_t=J_eta) mismatch            = %.2e  [wrong comparison]" % mismatch(J_eta))
