# -*- coding: utf-8 -*-
"""Paper I, section 4.4: from the proper-time Hamiltonian H_tau of eq. (Htau-vaidya),
solving H_tau = 0 for p_r and forming

    phi' = H_J/H_pr,   dv/dr = 1/H_pr,   dtau/dr = (dtau/dv)(dv/dr),  dtau/dv = (f - r')/E,

returns exactly the three spectral expressions on the sextic S = r(r-2m) D_E [r^2(r-2m) - J^2 D_E]:

    phi' = J D_E/sqrt(S),   dv/dr = E r^3/sqrt(S) + r/(r-2m),   dtau/dr = r^2(r-2m)/sqrt(S),

on the orientation of increasing r that makes sqrt(S) positive.  Checked at 50 digits,
m = 1, E = 7/5, J = 5/2, at several radii.
"""
import sys
import mpmath as mp

mp.mp.dps = 50
m, E, J = mp.mpf(1), mp.mpf(7)/5, mp.mpf(5)/2
worst = mp.mpf(0)
for r in [mp.mpf(x) for x in ('5', '6.5', '8', '9', '12')]:
    f = 1 - 2*m/r; w = E**2 - f
    H = lambda p, j: p*(f - E**2) - E + mp.sqrt(w)*mp.sqrt((E*p + 1)**2 + j**2/r**2)
    Hp = lambda p: mp.diff(lambda q: H(q, J), p)
    HJ = lambda p: mp.diff(lambda j: H(p, j), J)
    DE = (E**2 - 1)*r + 2*m
    S = r*(r - 2*m)*DE*(r**2*(r - 2*m) - J**2*DE)
    best = None
    for guess in (mp.mpf(-3), mp.mpf(-1), mp.mpf(0), mp.mpf(1), mp.mpf(3)):
        try:
            p = mp.findroot(lambda q: H(q, J), guess)
        except (ValueError, ZeroDivisionError):
            continue
        if abs(H(p, J)) < mp.mpf(10)**-40 and Hp(p) > 0:
            best = p; break
    assert best is not None, "no outgoing root at r=%s" % r
    p = best
    rdot = Hp(p)                               # dr/dv
    got = [HJ(p)/Hp(p), 1/Hp(p), ((f - rdot)/E)/Hp(p)]
    want = [J*DE/mp.sqrt(S), E*r**3/mp.sqrt(S) + r/(r - 2*m), r**2*(r - 2*m)/mp.sqrt(S)]
    err = max(abs(g - t)/abs(t) for g, t in zip(got, want))
    worst = max(worst, err)
    print(f"  r = {mp.nstr(r, 4):>5}: max relative difference {mp.nstr(err, 3)}")
ok = worst < mp.mpf(10)**-25
print("ALL CHECKS PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
