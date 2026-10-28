# -*- coding: utf-8 -*-
"""
The single adiabatic WKB curve: A(eta) varies slowly ALONG the orbit.

Equation (leading WKB order): the same closed form, but with E_eff running
INSIDE the integral:
    phi_ad(r) = int_{r0}^r F(r'; Ehat/A(eta(r')), J/A(...)) dr'
The adiabatic curve does NOT stay on the initial frozen orbit: it DRIFTS across
the family, interpolating between the frozen orbit at the initial A and the one at the final A.
Deviation from the initial frozen orbit = O(A'/A).

t branch. A(eta(r)) modelled as monotone along the orbit (a proxy for the elapsed time).
"""
import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a = 1.0, 0.9
Ehat = 1.4
J = 6.0
r0 = 12.0


def Dl(r): return r**2 - 2*M*r + a**2
def R6(r, E):
    Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
          + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
          - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
    return r*Q2*((E**2-1)*r + 2*M)
def dphidr(r, E):
    K = r*((E**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/Dl(r)
    return K/np.sqrt(R6(r, E))
def rmin(E):
    rs = np.linspace(2.1, r0, 40000); v = R6(rs, E)
    idx = np.where(np.diff(np.sign(v)))[0]
    rts = [brentq(lambda r: R6(r, E), rs[i], rs[i+1]) for i in idx]
    return max([r for r in rts if r > 2.0])


def frozen(A):
    E = Ehat/A
    rmn = rmin(E)
    rr = np.linspace(r0, rmn+1e-3, 500)
    ph = np.array([quad(lambda x: dphidr(x, E), r, r0, limit=200)[0] for r in rr])
    rr = np.concatenate([rr, rr[::-1]]); ph = np.concatenate([ph, 2*ph[-1]-ph[::-1]])
    return rr, ph


def adiabatic(A_i, A_f):
    """A(eta(r)) grows from A_i to A_f towards the turning point (monotone proxy)."""
    rmn = rmin(Ehat/((A_i+A_f)/2))   # approximate turning point with the mean E
    def Aloc(r, going_in=True):
        s = max(0.0, min(1.0, (r0-r)/(r0-rmn)))
        return A_i + (A_f-A_i)*s
    rr = np.linspace(r0, rmn+1e-3, 500)
    ph = np.array([quad(lambda x: dphidr(x, Ehat/Aloc(x)), r, r0, limit=200)[0]
                   for r in rr])
    rr = np.concatenate([rr, rr[::-1]]); ph = np.concatenate([ph, 2*ph[-1]-ph[::-1]])
    return rr, ph


A_i, A_f = 1.0, 1.2
rr_i, ph_i = frozen(A_i)
rr_f, ph_f = frozen(A_f)
rr_a, ph_a = adiabatic(A_i, A_f)

fig, ax = plt.subplots(figsize=(COL, COL*0.95))
th = np.linspace(0, 2*np.pi, 200)
ax.plot(rr_i*np.cos(ph_i), rr_i*np.sin(ph_i), 'C0-', lw=1.6,
        label=rf'frozen $A={A_i}$ ($E_{{\rm eff}}={Ehat/A_i:.2f}$)')
ax.plot(rr_f*np.cos(ph_f), rr_f*np.sin(ph_f), 'C1-', lw=1.6,
        label=rf'frozen $A={A_f}$ ($E_{{\rm eff}}={Ehat/A_f:.2f}$)')
ax.plot(rr_a*np.cos(ph_a), rr_a*np.sin(ph_a), 'k--', lw=1.8,
        label=r'adiabatic $A:%.1f\to%.1f$ (running)' % (A_i, A_f))
ax.plot(2*M*np.cos(th), 2*M*np.sin(th), 'b:', lw=0.7, label='$r_e$')
ax.plot(r0, 0, 'ks', ms=3)
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title('adiabatic WKB curve: $A$ runs slowly along the orbit\n'
             r'$\phi_{\rm ad}(r)=\int F(r;\hat E/A(\eta(r)))\,dr$'
             ' threads the family', fontsize=6.4)
ax.legend(fontsize=5.4, loc='upper left', framealpha=0.9)
savefig(fig, HERE, 'fig_adiabatic_curve')
print(f"frozen A=1 turning r={rr_i.min():.3f}, A=1.2 turning r={rr_f.min():.3f}, "
      f"adiabatic turning r={rr_a.min():.3f}")
print("=> the adiabatic orbit starts like A=1 (outside) and turns like ~A=1.2 (inside): it threads the family.")
print('DONE.')
