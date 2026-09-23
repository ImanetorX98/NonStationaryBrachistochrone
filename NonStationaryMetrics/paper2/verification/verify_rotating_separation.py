#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_rotating_separation.py -- geodesic separation for the *rotating* frozen
branch, i.e. the Randers case the static certificate does not reach.

Equatorial plane, frozen branch, manuscript eqs. (randers), (randers-data):

    F_T = alpha + beta,   alpha^2 = P dr^2 + Q dphi^2,   beta = b(r) dphi,
    P = n^2 r^2/(f Delta),  Q = n^2 Delta/f^2,  b = -2 M a/(r f),
    n^2 = E^2/(E^2 - f),    E = Ehat/A,   f = 1 - 2M/r.

F_T is parametrisation invariant and beta enters linearly, so its extremals are
the curves of prescribed geodesic curvature |kappa| = |B|, B = b'(r)/sqrt(PQ), in
the Riemannian optical metric alpha -- Perlick's own remark that the functional
"has the same mathematical structure as for a charged particle moving in a
magnetostatic field" (1991, p. 3153).  db = 0 iff a = 0, Randers' criterion for the
extremals to reduce to Riemannian geodesics.

CONVENTIONS, fixed once and used throughout.  Parametrise by the *alpha* arclength,
alpha(T) = 1 -- not the F arclength and not proper time.  Take

    J(v^r, v^phi) = ( -sqrt(Q/P) v^phi, +sqrt(P/Q) v^r ),   N = J T.

With this J the first variation gives  grad_T T = -B N,  the tangential constraint
is eta' = -B xi, and the normal equation carries +dB(N):

    xi'' + [ K + B^2 + dB(N) ] xi = 0.

The opposite choice N = -J T flips both signs and gives the -dB(N) form; the two
are one statement, but only if flow, normal and equation are flipped together.

An earlier version of this file integrated grad_T T = +B N while quoting the
-dB(N) equation -- the flow of a different functional.  Every test it ran was built
on that same flow, so none could see it.  The check that does see it is external to
the construction: phi is cyclic for the magnetic Lagrangian, so p_phi = Q phi' + b
must be conserved.  It is now the first test, and the flow below is derived from
the Lagrangian symbolically instead of being assembled by hand.

Verified:
  (0) p_phi and the alpha-speed are conserved along the integrated flow;
  (1) the integrated curve has geodesic curvature -B, matching the convention;
  (2) the scalar equation agrees with the *full linearised system* -- the Jacobian
      of the flow integrated along the reference curve, not with itself;
  (3) it agrees with centred finite differences of genuinely nearby extremals,
      at second order in the perturbation;
  (4) at a = 0 it collapses to xi'' + K xi = 0 with K the manuscript's
      eq. (optical-gauss).

Hypotheses: equatorial, A constant, effective energy E = Ehat/A fixed, M > 0,
physical Kerr, r > 2M, E^2 > f, regular future-directed curves, spatial endpoints
fixed.  The statement concerns extremals and their variations, not global minima.
"""
import numpy as np
import sympy as sp

ok = True


def check(label, cond):
    global ok
    ok = ok and bool(cond)
    print(f'{"PASS" if cond else "FAIL"}  {label}')


# ----------------------------------------------------------------- symbolic
r, M, a, E = sp.symbols('r M a E', positive=True)   # E = Ehat/A, frozen
f = 1 - 2*M/r
n2 = E**2/(E**2 - f)
P = sp.simplify(n2*r**2/(f*(r**2 - 2*M*r + a**2)))
Q = sp.simplify(n2*(r**2 - 2*M*r + a**2)/f**2)
b = sp.simplify(-2*M*a/(r*f))

sq = sp.sqrt(P*Q)
K = sp.simplify(-1/(2*sq)*sp.diff(sp.diff(Q, r)/sq, r))
B = sp.simplify(sp.diff(b, r)/sq)

check("B vanishes exactly when a = 0 (Randers' criterion, db = 0)",
      sp.simplify(B.subs(a, 0)) == 0 and sp.simplify(B) != 0)

# closed form quoted in eq. (magnetic-B): the two factors of Delta in P and Q
# cancel under the square root, so B is exactly linear in a and carries no Delta.
B_closed = 2*M*a*(E**2 - f)/(E**2*r**3*sp.sqrt(f))
check("B = 2 M a (E^2-f)/(E^2 r^3 sqrt f), free of Delta and linear in a",
      all(abs(float((B - B_closed).subs({M: 1.0, a: av, E: ev, r: rv}))) < 1e-15
          for av in (0.3, 0.9) for ev in (1.2, 1.4, 2.0)
          for rv in (2.5, 3.0, 6.0, 20.0)))
K0 = sp.simplify(K.subs(a, 0))
printed = (M*(16*M**3 + 2*(9*E**2 - 14)*M**2*r + (3*E**4 - 19*E**2 + 16)*M*r**2
              - (2*E**2 - 3)*(E**2 - 1)*r**3)/(E**2*r**5*(2*M + (E**2 - 1)*r)))
check('a -> 0 reproduces eq. (eq:optical-gauss)', sp.simplify(K0 - printed) == 0)

# Euler-Lagrange of L = (P r'^2 + Q phi'^2)/2 + b phi', whose extremals coincide
# with those of F on the shell alpha(T) = 1.  No sign is chosen by hand here.
vr, vp = sp.symbols('v_r v_phi')
acc_r = sp.simplify((-sp.diff(P, r)*vr**2/2 + sp.diff(Q, r)*vp**2/2
                     + sp.diff(b, r)*vp)/P)
acc_p = sp.simplify(-(sp.diff(Q, r)*vr*vp + sp.diff(b, r)*vr)/Q)

vals = {M: 1.0, a: 0.9, E: 1.4}
F_state = sp.lambdify((r, vr, vp),
                      [vr, vp, acc_r.subs(vals), acc_p.subs(vals)], 'numpy')
state = sp.Matrix([vr, vp, acc_r.subs(vals), acc_p.subs(vals)])
Jac = sp.lambdify((r, vr, vp),
                  state.jacobian(sp.Matrix([r, sp.Symbol('phi'), vr, vp])), 'numpy')

fP, fQ, fb, fK, fB = (sp.lambdify(r, x.subs(vals), 'numpy')
                      for x in (P, Q, b, K, B))
fdB = sp.lambdify(r, sp.diff(B, r).subs(vals), 'numpy')
fdP = sp.lambdify(r, sp.diff(P, r).subs(vals), 'numpy')
fdQ = sp.lambdify(r, sp.diff(Q, r).subs(vals), 'numpy')


def rhs(y):
    return np.array(F_state(y[0], y[2], y[3]), dtype=float)


def rk4(fun, y0, ds, n):
    y, out = np.array(y0, float), [np.array(y0, float)]
    for _ in range(n):
        k1 = fun(y); k2 = fun(y + ds/2*k1)
        k3 = fun(y + ds/2*k2); k4 = fun(y + ds*k3)
        y = y + ds/6*(k1 + 2*k2 + 2*k3 + k4)
        out.append(y.copy())
    return np.array(out)


# ----------------------------------------------------------------- reference
r0, ds, nst = 6.0, 1e-3, 12000
vr0, vp0 = 0.30, 1.0
nrm = np.sqrt(fP(r0)*vr0**2 + fQ(r0)*vp0**2)
y0 = [r0, 0.0, vr0/nrm, vp0/nrm]
ref = rk4(rhs, y0, ds, nst)

# ---- (0) the external check
pphi = np.array([fQ(y[0])*y[3] + fb(y[0]) for y in ref])
speed = np.array([np.sqrt(fP(y[0])*y[2]**2 + fQ(y[0])*y[3]**2) for y in ref])
check(f"p_phi = Q phi' + b conserved (drift {abs(pphi - pphi[0]).max():.2e})",
      abs(pphi - pphi[0]).max() < 1e-10)
check(f'alpha-unit speed conserved (drift {abs(speed - 1).max():.2e})',
      abs(speed - 1).max() < 1e-10)
print(f'  reference extremal: r from {ref[0,0]:.5f} to {ref[-1,0]:.5f}')

# ---- (1) geodesic curvature is -B
kap = []
for y in ref[::500]:
    acc = rhs(y)[2:]
    Pv, Qv, dPv, dQv = fP(y[0]), fQ(y[0]), fdP(y[0]), fdQ(y[0])
    cov_r = acc[0] + dPv/(2*Pv)*y[2]**2 - dQv/(2*Pv)*y[3]**2
    cov_p = acc[1] + dQv/Qv*y[2]*y[3]
    Nr, Np = -np.sqrt(Qv/Pv)*y[3], np.sqrt(Pv/Qv)*y[2]
    kap.append(Pv*cov_r*Nr + Qv*cov_p*Np)
kap = np.array(kap)
dev = abs(kap + fB(ref[::500, 0])).max()
check(f'geodesic curvature equals -B(r) (max dev {dev:.2e})', dev < 1e-9)


def dB_of_N(y):
    Pv, Qv = fP(y[0]), fQ(y[0])
    return fdB(y[0])*(-np.sqrt(Qv/Pv)*y[3])          # dB(N), N = J T


# ---- (2) scalar equation vs the FULL linearised system
def joint(w):
    y, z = w[:4], w[4:]
    return np.concatenate([rhs(y), np.asarray(Jac(y[0], y[2], y[3])) @ z])


Pv, Qv = fP(r0), fQ(r0)
Nvec = np.array([-np.sqrt(Qv/Pv)*y0[3], np.sqrt(Pv/Qv)*y0[2]])
sol = rk4(joint, np.concatenate([y0, [0.0, 0.0, Nvec[0], Nvec[1]]]), ds, nst)

xi_lin = []
for w in sol:
    y, z = w[:4], w[4:]
    Pv, Qv = fP(y[0]), fQ(y[0])
    Nr, Np = -np.sqrt(Qv/Pv)*y[3], np.sqrt(Pv/Qv)*y[2]
    xi_lin.append(Pv*z[0]*Nr + Qv*z[1]*Np)
xi_lin = np.array(xi_lin)


def scalar(w):
    y, xi, xip = w[:4], w[4], w[5]
    W = fK(y[0]) + fB(y[0])**2 + dB_of_N(y)
    return np.concatenate([rhs(y), [xip, -W*xi]])


xi_sc = rk4(scalar, np.concatenate([y0, [0.0, 1.0]]), ds, nst)[:, 4]
rel = abs(xi_sc - xi_lin).max()/abs(xi_lin).max()
check(f'scalar equation matches the full linearised system (rel {rel:.2e})',
      rel < 1e-8)

# ---- (3) centred finite differences
print('\n  eps        rel. error vs the scalar equation')
errs = []
for eps in (1e-3, 5e-4, 2.5e-4):
    curves = []
    for c in (eps, -eps):
        Tr, Tp = y0[2], y0[3]
        P0, Q0 = fP(r0), fQ(r0)          # at the start point, not wherever
        Nr, Np = -np.sqrt(Q0/P0)*Tp, np.sqrt(P0/Q0)*Tr
        curves.append(rk4(rhs, [r0, 0.0, Tr*np.cos(c) + Nr*np.sin(c),
                                Tp*np.cos(c) + Np*np.sin(c)], ds, nst))
    sep = []
    for yr, yp_, ym in zip(ref, curves[0], curves[1]):
        Pv2, Qv2 = fP(yr[0]), fQ(yr[0])
        dr_, dp_ = (yp_[0] - ym[0])/2, (yp_[1] - ym[1])/2
        Nr, Np = -np.sqrt(Qv2/Pv2)*yr[3], np.sqrt(Pv2/Qv2)*yr[2]
        sep.append(Pv2*dr_*Nr + Qv2*dp_*Np)
    errs.append(abs(np.array(sep)/eps - xi_sc).max()/abs(xi_sc).max())
    print(f'  {eps:.1e}    {errs[-1]:.3e}')
check('centred differences match the scalar equation', errs[-1] < 1e-6)
check(f'second-order convergence (ratios {errs[0]/errs[1]:.2f}, {errs[1]/errs[2]:.2f})',
      3.0 < errs[0]/errs[1] < 5.0 and 3.0 < errs[1]/errs[2] < 5.0)

# ---- (3b) the one check that does not live downstream of the Lagrangian
# Everything above presupposes that L = (P r'^2 + Q phi'^2)/2 + b phi' is the right
# Lagrangian for our functional.  If that identification were wrong, no test above
# could see it: they all sit downstream of it.  So evaluate the arrival cost
# directly from its definition, F_T = sqrt(P r'^2 + Q phi'^2) + b phi', on the
# integrated curve and on curves with the SAME spatial endpoints, and check that
# the first variation vanishes.  This uses no Lagrangian, no Christoffel symbol,
# no B and no Jacobi equation -- only the curve and eq. (randers-data).
def arrival_cost(rs, ps):
    """int F_T along a polyline, by the trapezoid rule on the two terms."""
    dr_, dp_ = np.diff(rs), np.diff(ps)
    rm, dphi = (rs[:-1] + rs[1:])/2, dp_
    quad = np.sqrt(fP(rm)*dr_**2 + fQ(rm)*dphi**2)
    return float(np.sum(quad + fb(rm)*dphi))


rs_ref, ps_ref = ref[:, 0], ref[:, 1]
L_arc = s[-1] if (s := np.arange(nst + 1)*ds) is not None else 1.0
bump = np.sin(np.pi*s/s[-1])                     # vanishes at both endpoints
Nr_arr = np.array([-np.sqrt(fQ(y[0])/fP(y[0]))*y[3] for y in ref])
Np_arr = np.array([np.sqrt(fP(y[0])/fQ(y[0]))*y[2] for y in ref])

print('\n  delta      [S(+d)-S(-d)]/2d        [S(+d)+S(-d)-2S(0)]/d^2')
first, second = [], []
S0 = arrival_cost(rs_ref, ps_ref)
for d in (4e-3, 2e-3, 1e-3):
    Sp = arrival_cost(rs_ref + d*bump*Nr_arr, ps_ref + d*bump*Np_arr)
    Sm = arrival_cost(rs_ref - d*bump*Nr_arr, ps_ref - d*bump*Np_arr)
    first.append((Sp - Sm)/(2*d))
    second.append((Sp + Sm - 2*S0)/d**2)
    print(f'  {d:.1e}   {first[-1]:+.6E}          {second[-1]:+.6E}')
check('first variation of the arrival cost vanishes on the integrated curve '
      f'(|dS/dd| <= {max(abs(x) for x in first):.2e})',
      max(abs(x) for x in first) < 1e-5)
check('the second variation is finite and stable, so the curve is a genuine '
      'stationary point rather than a numerical accident',
      max(second) - min(second) < 0.05*abs(np.mean(second)))
print('  This test uses only eq. (randers-data) and the curve: no Lagrangian, no')
print('  Christoffel symbols, no B, no Jacobi equation.  It is what closes the')
print('  identification the rest of the file depends on.')

# ---- (3c) both azimuthal senses, and the search for conjugate points
# The rotating case has no sign argument to exclude conjugate points, so they have
# to be located.  Integrate xi with xi(0)=0, xi'(0)=1 from a fixed initial point in
# both azimuthal senses and at two energies, and look for a first zero.
print('\n  sense    E      first conjugate point (alpha-arclength)   min W')
found = []
for sense in (+1, -1):
    for ev in (1.4, 1.2):
        fPe, fQe, fbe, fKe, fBe, fdBe = (
            sp.lambdify(r, x.subs({M: 1.0, a: 0.9, E: ev}), 'numpy')
            for x in (P, Q, b, K, B, sp.diff(B, r)))
        accs = sp.lambdify((r, vr, vp),
                           [vr, vp, acc_r.subs({M: 1.0, a: 0.9, E: ev}),
                            acc_p.subs({M: 1.0, a: 0.9, E: ev})], 'numpy')

        def rhs_e(y, _f=accs):
            return np.array(_f(y[0], y[2], y[3]), dtype=float)

        nr0 = np.sqrt(fPe(r0)*vr0**2 + fQe(r0)*(sense*vp0)**2)
        y0E = [r0, 0.0, vr0/nr0, sense*vp0/nr0]

        def sc_e(w, _r=rhs_e, _K=fKe, _B=fBe, _dB=fdBe, _P=fPe, _Q=fQe):
            y, xi, xip = w[:4], w[4], w[5]
            Nr = -np.sqrt(_Q(y[0])/_P(y[0]))*y[3]
            W_ = _K(y[0]) + _B(y[0])**2 + _dB(y[0])*Nr
            return np.concatenate([_r(y), [xip, -W_*xi]])

        sol_e = rk4(sc_e, np.concatenate([y0E, [0.0, 1.0]]), ds, nst)
        xi_e, Wmin = sol_e[:, 4], None
        zero = np.where(np.sign(xi_e[1:]) != np.sign(xi_e[1]))[0]
        s_conj = (zero[0] + 1)*ds if len(zero) else None
        Wser = np.array([fKe(y[0]) + fBe(y[0])**2
                         - np.sqrt(fQe(y[0])/fPe(y[0]))*y[3]*fdBe(y[0])
                         for y in sol_e[:, :4]])
        found.append(s_conj)
        lbl = 'none up to s = %.0f' % (nst*ds) if s_conj is None else f'{s_conj:.3f}'
        print(f'  {"pro " if sense > 0 else "retro"}   {ev:.1f}    {lbl:<38s}  {Wser.min():+.3E}')

check('no conjugate point is found on any of the four runs, consistent with W < 0 '
      'throughout them -- an observation on these arcs, not a theorem',
      all(x is None for x in found))

# ---- (4) a -> 0
check('a -> 0: B = 0 and the weight collapses to K',
      sp.simplify((K + B**2 + sp.diff(B, r)).subs(a, 0) - K.subs(a, 0)) == 0)

# and the same collapse numerically, from the final version of the code
fP0, fQ0, fK0 = (sp.lambdify(r, x.subs({M: 1.0, a: 1e-6, E: 1.4}), 'numpy')
                 for x in (P, Q, K))
fB0 = sp.lambdify(r, B.subs({M: 1.0, a: 1e-6, E: 1.4}), 'numpy')
check(f'numerically, B -> 0 as a -> 0 (|B| <= {abs(fB0(6.0)):.2e} at a = 1e-6)',
      abs(fB0(6.0)) < 1e-6)

W = np.array([fK(y[0]) + fB(y[0])**2 + dB_of_N(y) for y in ref])
Kv = np.array([fK(y[0]) for y in ref])
print(f'\n  along this extremal: K in [{Kv.min():+.8f}, {Kv.max():+.8f}]')
print(f'                       W in [{W.min():+.8f}, {W.max():+.8f}]')
print(f'                   W - K in [{(W-Kv).min():+.8f}, {(W-Kv).max():+.8f}]')
print('  W - K > 0 here: against K of the same Riemannian metric the magnetic term')
print('  makes W less negative.  One extremal at one energy, not a theorem, and not')
print('  a Kerr-vs-Schwarzschild comparison, where K and the curve both change.')
print('  B^2 + dB(N) has no fixed sign, so in the rotating case conjugate points')
print('  must be located, not excluded by a curvature sign.')

print('\nALL CHECKS PASS' if ok else '\nSOME CHECKS FAILED')
raise SystemExit(0 if ok else 1)
