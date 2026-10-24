"""verify_optical_planes.py

Central planes and the optical (Randers) data of frozen Thakurta-Kerr, g = A^2 g_Kerr with constant A > 0
(introduction, `sec:intro` bullet on the equatorial plane; appendix `app:rail-surfaces`, remark after
`lem:rail-surfaces`).

At fixed charge Ehat = -g(d_eta, u) the arrival functional is d eta = n alpha_K + beta_K with
  N0^2 = -g_tt, gamma0_ij = g_ij - g_ti g_tj / g_tt (seed Kerr), alpha_K^2 = gamma0 / N0^2,
  beta_K = -(g_ti / g_tt) dx^i, n = Ehat / sqrt(Ehat^2 - A^2 N0^2),
which reduces on the equator to eq:randers and eq:randers-data. The script checks, symbolically:
  1. the solution of the two conditions g(u,u) = -1 and -g(d_eta,u) = Ehat;
  2. the equatorial reduction;
  3. n alpha_K is even in a and beta_K is odd;
  4. reflection symmetry: the theta-derivatives of the data vanish on the equator, so the equator is
     totally geodesic for the Riemannian part and the magnetic force of d beta_K is tangent to it;
  5. at a = 0 the Riemannian part is spherically symmetric, hence at first order in a every central plane is
     totally geodesic for it;
  6. the normal component of the magnetic force on a central plane inclined by iota, at first order in a:
     d beta_K(T, m) = 4 M a cos(s) sin(iota) / ((rho - 2M) rho^2): identically zero on the plane only in the equatorial
     case (M a != 0); on an inclined plane the displayed contraction also vanishes where cos(s) = 0 (T the unit in-plane
     tangent of the circle of radius rho at angle s from the line of steepest inclination, m the plane normal);
  7. the proper-time data (N^2/Ehat) n alpha_K carry no one-form and are even in a;
  8. ordinary (Riemannian) umbilicity: alpha_K^2 = Omega^2 hbar, hbar = dr^2 + Delta dtheta^2 + Delta^2 sin^2(theta) dphi^2/Xi,
     Xi = Sigma - 2Mr, so h_t and h_tau are conformal to hbar; hbar is diagonal, independent of phi and reflection symmetric,
     so equatorial and meridional planes are totally geodesic for h_t and h_tau exactly in a; on a central plane inclined
     by 0 < iota < pi/2 the trace-free second fundamental form of hbar at order a^2 has the minors
     -(r-M) cos(iota) sin(iota) sin(s)/B^(3/2) and -sin(2 iota) cos(s)/sqrt(B), B = r(r-2M), never zero together
     (scaled by Omega^3 for h_t and h_tau); sampled at 40 digits in all quadrants of s, and checked independently with
     the exact metric h_t in the Cartesian chart.
Then, numerically (independent route): extremals of int(n alpha_K + beta_K) integrated as magnetic geodesics,
H = (1/2) a^{ij}(p - beta)_i (p - beta)_j at H = 1/2, launched tangent to a central plane; the out-of-plane
displacement vanishes on the equator and for a = 0; on an inclined plane its odd part is linear in a and its
even part quadratic.
Coverage: 2-7 symbolic; 8 symbolic for the conformal form, sampled (40 digits) for the minors; 1 by exact substitution at one rational point evaluated to 40 digits (the general identity is
derived in the appendix); sampled for the numerical check (A = M = Ehat = 1, one plane, one arc).

Run: python3 verify_optical_planes.py   (exit status 1 if a check fails)
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

r, th, M, a, s, iota = sp.symbols('r theta M a s iota', real=True)
rho = sp.symbols('rho', positive=True)
A, Eh = sp.symbols('A Ehat', positive=True)
checks = []


def check(name, expr):
    ok = sp.simplify(expr) == 0
    checks.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name)


Sig = r**2 + a**2*sp.cos(th)**2
Del = r**2 - 2*M*r + a**2
gtt = -(1 - 2*M*r/Sig)
gtp = -2*M*a*r*sp.sin(th)**2/Sig
gpp = (r**2 + a**2 + 2*M*a**2*r*sp.sin(th)**2/Sig)*sp.sin(th)**2
grr = Sig/Del
gthth = Sig
N02 = -gtt
gam0 = sp.diag(grr, gthth, gpp - gtp**2/gtt)          # coordinates (r, theta, phi)
bphi = -gtp/gtt
n = Eh/sp.sqrt(Eh**2 - A**2*N02)

# 1. the arrival functional from the two conditions, for g = A^2 g_Kerr (A frozen)
dr, dth, dph, ds = sp.symbols('dr dtheta dphi dsigma', real=True)
dx = sp.Matrix([dr, dth, dph])
gam_len = sp.sqrt((dx.T*gam0*dx)[0])
deta = n*gam_len/sp.sqrt(N02) + bphi*dph                # n alpha_K + beta_K along dx
# u = (d eta, dx)/d tau with d tau = A^2 N0 gam_len / sqrt(Ehat^2 - A^2 N0^2) (proper time of the rail)
dtau = A**2*sp.sqrt(N02)*gam_len/sp.sqrt(Eh**2 - A**2*N02)
G = A**2*sp.Matrix([[gtt, 0, 0, gtp], [0, grr, 0, 0], [0, 0, gthth, 0], [gtp, 0, 0, gpp]])  # (eta, r, theta, phi)
U = sp.Matrix([deta, dr, dth, dph])/dtau
num = {r: sp.Rational(37, 3), th: sp.Rational(7, 10), M: 1, a: sp.Rational(3, 7), A: sp.Rational(5, 4), Eh: sp.Rational(13, 10),
       dr: sp.Rational(1, 3), dth: -sp.Rational(1, 5), dph: sp.Rational(2, 9)}
norm = (U.T*G*U)[0] + 1
charge = -(G*U)[0] - Eh
for name, e in [('normalisation g(u,u) = -1', norm), ('charge -g(d_eta,u) = Ehat', charge)]:
    ok = bool(abs(sp.N(e.subs(num), 40)) < 1e-30)
    checks.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name + ' (exact data, 40 digits)')

# 2. equatorial reduction to eq:randers-data
eq = {th: sp.pi/2}
f = 1 - 2*M/r
check('alpha_K^2 rr on the equator = r^2/(f Delta)', (gam0[0, 0]/N02).subs(eq) - r**2/(f*Del))
check('alpha_K^2 phiphi on the equator = Delta/f^2', ((gam0[2, 2])/N02).subs(eq) - Del/f**2)
check('beta_K on the equator = -2Ma/(r f)', bphi.subs(eq) + 2*M*a/(r*f))
check('n on the equator = Ehat/sqrt(Ehat^2 - A^2 f)', n.subs(eq) - Eh/sp.sqrt(Eh**2 - A**2*f))

# 3. parity in a
h = n**2*gam0/N02
check('n^2 alpha_K^2 even in a', sp.simplify(h - h.subs(a, -a)).norm())
check('beta_K odd in a', bphi + bphi.subs(a, -a))

# 4. reflection: theta-derivatives vanish on the equator; no r-theta or theta-phi terms
for i, j, lab in [(0, 0, 'rr'), (1, 1, 'thth'), (2, 2, 'phiphi')]:
    check(f'd_theta (n^2 alpha_K^2)_{lab} = 0 on the equator', sp.diff(h[i, j], th).subs(eq))
check('d_theta beta_phi = 0 on the equator', sp.diff(bphi, th).subs(eq))
check('off-diagonal r-theta, theta-phi terms vanish', h[0, 1]**2 + h[1, 2]**2 + h[0, 2]**2)

# 5. spherical symmetry of the Riemannian part at a = 0
h0 = h.subs(a, 0)
check('a = 0: h_thth = r^2 n^2 / f', h0[1, 1] - r**2*n.subs(a, 0)**2/f)
check('a = 0: h_phiphi = h_thth sin^2(theta)', h0[2, 2] - h0[1, 1]*sp.sin(th)**2)
check('a = 0: n depends on r only', sp.diff(n.subs(a, 0), th))

# 6. magnetic force normal to an inclined central plane, first order in a
b1 = sp.simplify(sp.series(bphi, a, 0, 2).removeO())
check('beta_K at first order = -2Ma sin^2(theta)/(r - 2M)', b1 + 2*M*a*sp.sin(th)**2/(r - 2*M))
X, Y, Z = sp.symbols('X Y Z', real=True)
R = sp.sqrt(X**2 + Y**2 + Z**2)
c = sp.simplify(b1/(r**2*sp.sin(th)**2)).subs(r, R)     # beta = c(R) (X dY - Y dX)
beta = [-c*Y, c*X, 0]
xs = [X, Y, Z]
Fm = sp.Matrix(3, 3, lambda i, j: sp.diff(beta[j], xs[i]) - sp.diff(beta[i], xs[j]))
e1 = sp.Matrix([sp.cos(iota), 0, -sp.sin(iota)]); e2 = sp.Matrix([0, 1, 0]); m = sp.Matrix([sp.sin(iota), 0, sp.cos(iota)])
pt = rho*(sp.cos(s)*e1 + sp.sin(s)*e2); T = -sp.sin(s)*e1 + sp.cos(s)*e2
force = sp.simplify((T.T*Fm.subs(dict(zip(xs, list(pt))))*m)[0])
print('normal magnetic force:', force)
check('normal magnetic force d beta(T, m) = 4Ma cos(s) sin(iota)/((rho-2M) rho^2)', force - 4*M*a*sp.cos(s)*sp.sin(iota)/((rho - 2*M)*rho**2))
check('normal magnetic force vanishes on the equator (iota = 0)', force.subs(iota, 0))
check('on an inclined plane the contraction also vanishes where cos(s) = 0', force.subs(s, sp.pi/2))

# 7. proper-time data: (N^2/Ehat) n alpha_K for g = A^2 g_Kerr is (A^2 N0^2/Ehat) n alpha_K, no one-form
htau = (A**2*N02/Eh)**2*h
check('proper-time data even in a', sp.simplify(htau - htau.subs(a, -a)).norm())

# numerical, independent route: magnetic geodesics of n alpha_K + beta_K (A = M = Ehat = 1)
subsN = {A: 1, M: 1, Eh: 1}
hN = h.subs(subsN); bN = bphi.subs(subsN)
hinv = sp.simplify(hN.inv())
pr, pt_, pp = sp.symbols('p_r p_theta p_phi', real=True)
P = sp.Matrix([pr, pt_, pp - bN])
Ham = sp.Rational(1, 2)*(P.T*hinv*P)[0]
q = [r, th, sp.Symbol('phi')]
rhs = [sp.diff(Ham, v) for v in (pr, pt_, pp)] + [-sp.diff(Ham, v) for v in (r, th)] + [0]
F = sp.lambdify((r, th, pr, pt_, pp, a), rhs, 'numpy')
Hf = sp.lambdify((r, th, pr, pt_, pp, a), Ham, 'numpy')
hf = sp.lambdify((r, th, a), hN, 'numpy'); bf = sp.lambdify((r, th, a), bN, 'numpy')


def launch(av, inc, r0=50.0):
    """Start at the point of the plane of inclination inc that lies on the x-axis, tangent to the plane."""
    x0 = np.array([r0, 0., 0.]); tvec = np.array([0., np.cos(inc), np.sin(inc)])   # plane spanned by x and tvec
    # spherical velocity: theta = pi/2 at x0, d theta = -t_z/r, d phi = t_y/r
    v = np.array([0., -tvec[2]/r0, tvec[1]/r0])
    hm = np.array(hf(r0, np.pi/2, av), dtype=float)
    v = v/np.sqrt(v @ hm @ v)                                   # unit alpha-speed
    p = hm @ v + np.array([0., 0., bf(r0, np.pi/2, av)])        # p = h v + beta
    return np.array([r0, np.pi/2, 0., *p])


def out_of_plane(av, inc, length=60.0):
    y0 = launch(av, inc)
    sol = solve_ivp(lambda t, y: np.array(F(y[0], y[1], y[3], y[4], y[5], av), dtype=float)[[0, 1, 2, 3, 4, 5]],
                    (0, length), y0, rtol=1e-12, atol=1e-14, method='DOP853')
    rr, tt, ph = sol.y[0, -1], sol.y[1, -1], sol.y[2, -1]
    x = rr*np.array([np.sin(tt)*np.cos(ph), np.sin(tt)*np.sin(ph), np.cos(tt)])
    mvec = np.array([0., -np.sin(inc), np.cos(inc)])            # normal of span{x, tvec}
    drift = Hf(sol.y[0, -1], sol.y[1, -1], sol.y[3, -1], sol.y[4, -1], sol.y[5, -1], av) - 0.5
    return float(x @ mvec), float(drift)


inc = np.pi/4
d0, e0 = out_of_plane(0.0, inc)
dp, ep = out_of_plane(1e-3, inc)
dm, em = out_of_plane(-1e-3, inc)
d2, _ = out_of_plane(2e-3, inc)
d2m, _ = out_of_plane(-2e-3, inc)
deq, _ = out_of_plane(1e-3, 0.0)
print(f'numerical: out-of-plane displacement after the arc (M units): a=0 {d0:.3e}, a=+1e-3 {dp:.6e}, a=-1e-3 {dm:.6e}, '
      f'a=2e-3 {d2:.6e}; equator with a=1e-3 {deq:.3e}; |H - 1/2| <= {max(abs(e0), abs(ep), abs(em)):.1e}')
for name, ok in [('numerical: a = 0 stays in the inclined plane', abs(d0) < 1e-9),
                 ('numerical: equator invariant for a != 0', abs(deq) < 1e-9),
                 ('numerical: odd part dominates (even part < 1% at a = 1e-3)', abs(dp + dm) < 1e-2*abs(dp)),
                 ('numerical: even part quadratic in a (ratio 4)', abs((d2 + d2m)/(dp + dm) - 4) < 0.05),
                 ('numerical: displacement linear in a (ratio 2)', abs(d2/dp - 2) < 1e-2),
                 ('numerical: displacement nonzero for a = 1e-3', abs(dp) > 1e3*max(abs(d0), 1e-12)),
                 ('numerical: Hamiltonian conserved', max(abs(e0), abs(ep), abs(em)) < 1e-9)]:
    checks.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name)


# 8. ordinary (Riemannian) umbilicity of central coordinate planes for h_t = n^2 alpha_K^2 and h_tau
#    (no notion of umbilicity is adopted for the Randers data themselves)
Xi = Sig - 2*M*r
hbar = sp.diag(1, Del, Del**2*sp.sin(th)**2/Xi)
Om2 = Sig**2/(Del*Xi)
check('alpha_K^2 = Omega^2 hbar with Omega^2 = Sigma^2/(Delta Xi), Xi = Sigma - 2Mr',
      sum(sp.simplify(e)**2 for e in (gam0/N02 - Om2*hbar)))
check('hbar invariant under theta -> pi - theta (and diagonal, independent of phi)',
      sum(sp.simplify(e)**2 for e in (hbar.subs(th, sp.pi - th) - hbar)))
ph_ = sp.Symbol('varphi', real=True)
Xc = [r, th, ph_]
h0b = hbar.subs(a, 0).applyfunc(sp.simplify)
h2b = hbar.applyfunc(lambda e: sp.simplify(sp.diff(e, a, 2).subs(a, 0)/2))
B = r*(r - 2*M)
check('hbar = dr^2 + B dOmega^2 + a^2 [dOmega^2 + (sin^2 theta dphi)^2] + O(a^4)',
      sp.simplify(h0b[1, 1] - B)**2 + sp.simplify(h0b[2, 2] - B*sp.sin(th)**2)**2 + sp.simplify(h2b[1, 1] - 1)**2
      + sp.simplify(h2b[2, 2] - sp.sin(th)**2 - sp.sin(th)**4)**2 + h2b[0, 0]**2)
h0i = h0b.inv()
G0 = [[[sum(h0i[k, l]*(sp.diff(h0b[l, j], Xc[i]) + sp.diff(h0b[l, i], Xc[j]) - sp.diff(h0b[i, j], Xc[l]))/2 for l in range(3))
        for j in range(3)] for i in range(3)] for k in range(3)]
dG = [[[sum(h0i[k, l]*(sp.diff(h2b[l, j], Xc[i]) + sp.diff(h2b[l, i], Xc[j]) - sp.diff(h2b[i, j], Xc[l]))/2 for l in range(3))
        - sum(h0i[k, l]*h2b[l, q]*G0[q][i][j] for l in range(3) for q in range(3))
        for j in range(3)] for i in range(3)] for k in range(3)]
Sf = sp.sin(iota)*sp.sin(th)*sp.cos(ph_) + sp.cos(iota)*sp.cos(th)
dS = [sp.diff(Sf, v) for v in Xc]
tg = [[1, 0, 0], [0, -sp.sin(iota)*sp.sin(ph_), sp.cos(iota)/sp.sin(th)**2]]   # d_rho and d_s on the plane
hess0 = [[sp.diff(Sf, Xc[i], Xc[j]) - sum(G0[k][i][j]*dS[k] for k in range(3)) for j in range(3)] for i in range(3)]
nrm0 = sp.sqrt(sum(h0i[i, i]*dS[i]**2 for i in range(3)))
II2 = [[-sum(dG[k][i][j]*dS[k]*tg[p][i]*tg[q][j] for k in range(3) for i in range(3) for j in range(3))/nrm0
        for q in range(2)] for p in range(2)]
H0 = [[sum(hess0[i][j]*tg[p][i]*tg[q][j] for i in range(3) for j in range(3)) for q in range(2)] for p in range(2)]
g0 = [[sum(h0b[i, j]*tg[p][i]*tg[q][j] for i in range(3) for j in range(3)) for q in range(2)] for p in range(2)]
m1 = II2[0][1]*g0[0][0] - II2[0][0]*g0[0][1]
m2 = II2[0][0]*g0[1][1] - II2[1][1]*g0[0][0]
m1t = -(r - M)*sp.cos(iota)*sp.sin(iota)*sp.sin(s)/B**sp.Rational(3, 2)
m2t = -sp.sin(2*iota)*sp.cos(s)/sp.sqrt(B)
fnum = sp.lambdify((r, th, ph_, iota, M), [H0[0][0], H0[0][1], H0[1][1], g0[0][1], m1, m2], 'mpmath')
ftgt = sp.lambdify((r, s, iota, M), [m1t, m2t], 'mpmath')
import mpmath as mp
mp.mp.dps = 40
worst_h0, worst_m = mp.mpf(0), mp.mpf(0)
for Mv in (1, 2):
    for rv in (mp.mpf(3)*Mv, mp.mpf(10)*Mv, mp.mpf(1000)*Mv):
        for io in (mp.pi/6, mp.pi/4, mp.pi/3):
            for sv in [mp.mpf(k)*mp.pi/6 + mp.mpf(1)/10 for k in range(12)]:
                thv = mp.acos(-mp.sin(io)*mp.cos(sv)); phv = mp.atan2(mp.sin(sv), mp.cos(io)*mp.cos(sv))
                v = fnum(rv, thv, phv, io, Mv); t = ftgt(rv, sv, io, Mv)
                worst_h0 = max(worst_h0, *[abs(x) for x in v[:4]])
                worst_m = max(worst_m, abs(v[4] - t[0])/abs(t[0]), abs(v[5] - t[1])/abs(t[1]))
print('umbilicity: largest tangential Hess_h0 S and g0_rs', mp.nstr(worst_h0, 3), '; largest relative deviation of the a^2 minors',
      mp.nstr(worst_m, 3), '(M = 1, 2; r/M = 3, 10, 1000; iota = 30, 45, 60 deg; s in all quadrants)')
checks.append(('a = 0: inclined central planes totally geodesic for hbar, induced metric diagonal (sampled, 40 digits)', worst_h0 < mp.mpf('1e-30')))
checks.append(('order a^2 minors of hbar: -(r-M) cos(i) sin(i) sin(s)/B^(3/2) and -sin(2i) cos(s)/sqrt(B) (sampled, 40 digits)', worst_m < mp.mpf('1e-30')))
fmut = sp.lambdify((r, s, iota, M), -(r - 2*M)*sp.cos(iota)*sp.sin(iota)*sp.sin(s)/B**sp.Rational(3, 2), 'mpmath')
thv = mp.acos(-mp.sin(mp.pi/4)*mp.cos(1)); phv = mp.atan2(mp.sin(1), mp.cos(mp.pi/4)*mp.cos(1))
checks.append(('mutation: factor (r - 2M) in place of (r - M) is rejected',
               abs(fnum(10, thv, phv, mp.pi/4, 1)[4]/fmut(10, 1, mp.pi/4, 1) - 1) > mp.mpf('1e-3')))
for name, ok in checks[-3:]:
    print(('PASS ' if ok else 'FAIL ') + name)


# independent route: exact metric h_t in the Cartesian chart built on Boyer-Lindquist, numerical Christoffel symbols,
# a^2 coefficient by Richardson; compare with Omega_t0^3 times the hbar minors (M = 1; Omega_t0 = Ehat r^(3/2)/((r-2) sqrt(D)),
# D = 2 + (Ehat^2 - 1) r). Ordinary umbilicity is conformally invariant: both minors scale by Omega^3.
def _metric_cart(x, av, Ev):
    X_, Y_, Z_ = x
    rr = mp.sqrt(X_*X_ + Y_*Y_ + Z_*Z_); rh = mp.sqrt(X_*X_ + Y_*Y_); tt = mp.acos(Z_/rr)
    Sg = rr**2 + av**2*mp.cos(tt)**2; Dl = rr**2 - 2*rr + av**2
    gtt_ = -(1 - 2*rr/Sg); gtp_ = -2*av*rr*mp.sin(tt)**2/Sg
    gpp_ = (rr**2 + av**2 + 2*av**2*rr*mp.sin(tt)**2/Sg)*mp.sin(tt)**2
    fac = Ev**2/(Ev**2 + gtt_)/(-gtt_)
    d = [fac*Sg/Dl, fac*Sg, fac*(gpp_ - gtp_**2/gtt_)]
    J = [[X_/rr, Y_/rr, Z_/rr], [X_*Z_/(rr*rr*rh), Y_*Z_/(rr*rr*rh), -rh/(rr*rr)], [-Y_/rh**2, X_/rh**2, 0]]
    return mp.matrix([[sum(d[k]*J[k][i]*J[k][j] for k in range(3)) for j in range(3)] for i in range(3)])


def _cart_minors(rv, io, sv, av, Ev):
    e1 = [mp.cos(io), 0, -mp.sin(io)]; e2 = [0, 1, 0]; mv = [mp.sin(io), 0, mp.cos(io)]
    x = [rv*(mp.cos(sv)*e1[i] + mp.sin(sv)*e2[i]) for i in range(3)]
    R = [mp.cos(sv)*e1[i] + mp.sin(sv)*e2[i] for i in range(3)]; T = [rv*(-mp.sin(sv)*e1[i] + mp.cos(sv)*e2[i]) for i in range(3)]
    g = _metric_cart(x, av, Ev); gi = g**-1; hh = mp.mpf('1e-10'); dg = []
    for k in range(3):
        def comp(t, k=k):
            y = list(x); y[k] += t
            return _metric_cart(y, av, Ev)
        dg.append((8*(comp(hh) - comp(-hh)) - (comp(2*hh) - comp(-2*hh)))/(12*hh))
    Gm = [[[sum(gi[k, l]*(dg[i][l, j] + dg[j][l, i] - dg[l][i, j])/2 for l in range(3)) for j in range(3)] for i in range(3)] for k in range(3)]
    mn = mp.sqrt(sum(gi[i, j]*mv[i]*mv[j] for i in range(3) for j in range(3)))
    II = lambda u, w: -sum(Gm[k][i][j]*mv[k]*u[i]*w[j] for i in range(3) for j in range(3) for k in range(3))/mn
    G = lambda u, w: sum(g[i, j]*u[i]*w[j] for i in range(3) for j in range(3))
    return [II(R, T)*G(R, R) - II(R, R)*G(R, T), II(R, R)*G(T, T) - II(T, T)*G(R, R)]


worst_c = mp.mpf(0)
for rv, Ev, io, sv in [(mp.mpf(10), mp.mpf(13)/10, mp.pi/4, mp.mpf(1)/3), (mp.mpf(4), mp.mpf(1), mp.pi/3, 7*mp.pi/6)]:
    c = []
    for av in (mp.mpf('1e-4'), mp.mpf('2e-4')):
        pp = _cart_minors(rv, io, sv, av, Ev); qq = _cart_minors(rv, io, sv, -av, Ev)
        c.append([(pp[i] + qq[i])/(2*av**2) for i in range(2)])
    c2 = [(4*c[0][i] - c[1][i])/3 for i in range(2)]
    Om = Ev*rv**mp.mpf(1.5)/((rv - 2)*mp.sqrt(2 + (Ev**2 - 1)*rv))
    t = ftgt(rv, sv, io, 1)
    worst_c = max(worst_c, *[abs(c2[i] - Om**3*t[i])/abs(Om**3*t[i]) for i in range(2)])
print('umbilicity, independent Cartesian route with the exact metric h_t: largest relative deviation', mp.nstr(worst_c, 3))
checks.append(('independent route: exact h_t minors = Omega_t0^3 x hbar minors (two points, Richardson)', worst_c < mp.mpf('1e-12')))
print(('PASS ' if checks[-1][1] else 'FAIL ') + checks[-1][0])

npass = sum(bool(ok) for _, ok in checks)
print(f'{npass}/{len(checks)} checks passed')
sys.exit(0 if npass == len(checks) else 1)
