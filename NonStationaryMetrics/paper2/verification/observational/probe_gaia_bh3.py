"""Clock thought experiment around Gaia BH3 (paragraph `par:obs-kerr-thought`).

Frozen Kerr, A = 1, G = c = M = 1 unless stated; M = 32.70 Msun, distance 590.6 pc and the companion orbit from the Gaia
Collaboration (2024): m* = 0.76 Msun, P = 4253.1 d, e = 0.7291, relative semi-major axis a = a1 (M + m*)/M with
a1 = 16.17 au (checked against Kepler's third law); illustrative spin a* = 0.5 (not measured).
1. Encounter clock effect. Two free probes on the same spherical orbit of radius r and inclination i, in opposite
   senses, compare clocks at each encounter. At first order in a, exact in r:
       Delta tau per encounter = 6 pi a cos(i) / sqrt(r (r - 3)).
   Symbolic derivation (E, L from R = R' = 0 at first order; dt/dlambda independent of theta; nodal periods and encounter
   interval pi r^(3/2)); independent check against the exact equatorial closed form.
2. Single-station circuit: the zero-spin Ehat = 1 arrival rail from A to B (turning radius r_min, endpoints at
   r0 = 2^(2/3) r_min), closed by the circular arc at r0, run by two probes in opposite senses. On a fixed curve
   dt = alpha + b dphi with alpha reversible and dtau = f alpha (Ehat = 1), so D = Delta t = 2 oint b exactly in a
   (equatorial plane), Delta tau = 0; the static station reads sqrt(f0) D. Benchmark r_min = 1.5e5 M; also 1e5, 2e5 M.
3. Companion: first-order shift of the arrival-time difference at the station caused by its Newtonian quadrupole tide,
   Phi = -(G m*/2R^3)(3 (x.n)^2 - x^2), in the black hole's freely falling frame, for two control laws:
   law B (no tangential thrust: the probe slides on the fixed circuit in the total field), with
       v dv = Phi(x(0),0) - Phi(x(t),t) + int_0^t d_t Phi dt',  dT_arr = -int dv/v dt
   (the static part cancels exactly between the two senses: checked); law A (ideal tangential feedforward compensation of
   the modelled tide, which keeps the charge with respect to the black hole exactly when the model is exact) leaves
   eps times the law-B shift if the modelled potential is (1 - eps) times the true one on the whole circuit, i.e. an error
   in the companion mass only, with its ephemeris, the circuit and the actuation held fixed; eps = 0.05/0.76 from the
   published mass 0.76 +- 0.05 Msun. This is a one-parameter illustration, not a joint uncertainty budget. Sampled over a
   deterministic grid of 48 orientations of the circuit plane, near periastron and at apastron; checked against a direct
   nonlinear integration of the sliding motion. Gaia BH1 and BH2 (El-Badry et al. 2023) are run for comparison.
4. Thrust budget of the circuit at a = 0: rail peak (2 g_s at r_min, relativistic rule of two) and arc (exact Christoffel
   symbols); I = int |b| dtau / c with onboard proper time and the two 120-degree turns at maintained local speed and
   charge (gamma beta theta each; the direct boost is printed as a lower bound); ideal photon-rocket ratio exp(I).
   Transfer from the Sun for scale only (two symmetric rapidity changes, total mass ratio 10 and 20, cruise at 590.6 pc).
Writes probe_gaia_bh3.json. Exit status 1 if a check fails.
"""
from pathlib import Path
import json
import math
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
C = 299792458.; GM_SUN = 1.32712440018e20; AU = 1.495978707e11; DAY = 86400.
MASS = 32.70; MSTAR = 0.76; PBIN = 4253.1*DAY; EBIN = 0.7291; ABIN = 16.17*AU*(MASS + MSTAR)/MASS; DIST_PC = 590.6
TG = MASS*GM_SUN/C**3; RG = MASS*GM_SUN/C**2
checks = []


def check(name, ok, value=None):
    checks.append(dict(name=name, status='PASS' if ok else 'FAIL', value=value))
    print(('PASS ' if ok else 'FAIL ') + name + ('' if value is None else f'  [{value}]'))


# ---------------------------------------------------------------- 1. encounter clock effect (symbolic)
r, a, ii, s_, th, lam = sp.symbols('r a i s theta lambda', real=True)
En, L, e1, l1, x = sp.symbols('E L e1 l1 x', real=True)
Rr = (En*(x**2 + a**2) - a*L)**2 - (x**2 - 2*x + a**2)*(x**2 + (L - a*En)**2 + L**2*sp.tan(ii)**2)
E0 = (1 - 2/r)/sp.sqrt(1 - 3/r); l0 = r/sp.sqrt(r - 3)
sub = {En: E0 + a*e1, L: s_*l0*sp.cos(ii) + a*l1}
eqs = [sp.series(e.subs(sub), a, 0, 2).removeO() for e in (Rr.subs(x, r), sp.diff(Rr, x).subs(x, r))]
assume = dict(r=sp.Q.positive(r - 3))
zero0 = [sp.simplify(e.subs(a, 0).subs(s_**2, 1)) for e in eqs]
check('spherical orbit at zeroth order: R = R\' = 0', all(sp.simplify(z.subs({r: sp.Rational(41, 3), ii: sp.Rational(3, 5), s_: 1})) == 0 and
                                                       sp.simplify(z.subs({r: sp.Rational(41, 3), ii: sp.Rational(3, 5), s_: -1})) == 0 for z in zero0))
sol1 = sp.solve([sp.expand(e).coeff(a, 1) for e in eqs], [e1, l1], dict=True)[0]
sol1 = {k: sp.simplify(v.subs(s_**2, 1)) for k, v in sol1.items()}
dtdl = (r**2 + a**2)/(r**2 - 2*r + a**2)*(En*(r**2 + a**2) - a*L) - a*(a*En*sp.sin(th)**2 - L)
G1 = sp.series(dtdl, a, 0, 2).removeO()
check('dt/dlambda independent of theta at first order', sp.simplify(sp.diff(G1, th)) == 0)


def first(expr, sg):
    return sp.series(expr.subs(sub).subs(sol1).subs(s_, sg), a, 0, 2).removeO()


rate = {sg: first(r**2/G1, sg) for sg in (1, -1)}                       # dtau/dt, Sigma = r^2 + O(a^2)
Lt = {sg: first(sg*L/sp.cos(ii), sg) for sg in (1, -1)}   # sqrt(L^2 + Q) = |L|/cos(i), Q = L^2 tan^2(i), sign(L) = s
Tn = {sg: sp.series(2*sp.pi/Lt[sg]*first(G1, sg), a, 0, 2).removeO() for sg in (1, -1)}
tmeet = sp.series(1/(1/Tn[1] + 1/Tn[-1]), a, 0, 2).removeO()
dtau = sp.simplify(tmeet*(rate[1] - rate[-1]))
pts = [{r: sp.Rational(41, 3), ii: sp.Rational(3, 5)}, {r: 200, ii: sp.Rational(1, 7)}, {r: sp.Rational(10**5), ii: 1}]
num_ok = all(abs(sp.N((Tn[1] - Tn[-1] - 12*sp.pi*a*sp.cos(ii)).subs(p).subs(a, sp.Rational(1, 10)), 30)) < 1e-20 for p in pts)
check('nodal periods T+ - T- = 12 pi a cos i (first order)', num_ok)
check('encounter interval pi r^(3/2) at first order', all(abs(sp.N((tmeet - sp.pi*r**sp.Rational(3, 2)).subs(p).subs(a, sp.Rational(1, 10)), 30)) < 1e-15 for p in pts))
target = 6*sp.pi*a*sp.cos(ii)/sp.sqrt(r*(r - 3))
check('Delta tau per encounter = 6 pi a cos i / sqrt(r(r-3))', all(abs(sp.N((dtau - target).subs(p).subs(a, sp.Rational(1, 10)), 30)) < 1e-20 for p in pts))
# exact equatorial closed form
Om = {sg: sg/(r**sp.Rational(3, 2) + sg*a) for sg in (1, -1)}
ut = {sg: (r**sp.Rational(3, 2) + sg*a)/(r**sp.Rational(3, 4)*sp.sqrt(r**sp.Rational(3, 2) - 3*sp.sqrt(r) + 2*sg*a)) for sg in (1, -1)}
tm_ex = 2*sp.pi/(Om[1] - Om[-1])
d_ex = tm_ex*(1/ut[1] - 1/ut[-1])
check('equatorial encounter interval exactly pi (r^3 - a^2)/r^(3/2)', sp.simplify(tm_ex - sp.pi*(r**3 - a**2)/r**sp.Rational(3, 2)) == 0)
check('exact equatorial Delta tau odd in a', abs(sp.N((d_ex + d_ex.subs(a, -a)).subs({r: sp.Rational(37, 3), a: sp.Rational(1, 7)}), 30)) < 1e-25)
lin = sp.series(d_ex, a, 0, 2).removeO().coeff(a, 1)
check('exact equatorial linear coefficient = 6 pi / sqrt(r(r-3))',
      all(abs(sp.N((lin - 6*sp.pi/sp.sqrt(r*(r - 3))).subs(r, rv), 40)) < 1e-35 for rv in (sp.Rational(37, 3), 200, 10**5)))
cub_ratio = {rv: float(sp.N((sp.series(d_ex, a, 0, 4).removeO().coeff(a, 3)/lin).subs(r, rv), 12)) for rv in (100, 10**5)}
enc = []
for rv in (1e4, 1e5, 1e6):
    d1 = 6*math.pi*0.5/math.sqrt(rv*(rv - 3))*TG
    taum = math.pi*rv**1.5*math.sqrt(1 - 3/rv)*TG
    enc.append(dict(r_M=rv, jump_ns=d1*1e9, fractional_rate=d1/taum, period_days=2*math.pi*rv**1.5*TG/DAY,
                    campaign_200_encounters_coordinate_days=200*math.pi*rv**1.5*TG/DAY))
check('encounter jump at 1e5 M = 4.46583 ns x (M/9.62 Msun) = 15.18 ns', abs(enc[1]['jump_ns'] - 4.46583*MASS/9.62) < 1e-4, round(enc[1]['jump_ns'], 5))
check('companion semi-major axis agrees with Kepler III', abs((GM_SUN*(MASS + MSTAR)*PBIN**2/(4*math.pi**2))**(1/3)/ABIN - 1) < 0.02, round(ABIN/AU, 3))

# ---------------------------------------------------------------- 2. single-station circuit
ra = sp.symbols('ra', positive=True)
P = lambda rv: rv**3/(2*(rv - 2)**2)           # Ehat = 1, a = 0
Q = lambda rv: rv**4/(2*(rv - 2))
rd, pd = sp.symbols('rd pd', real=True)
fr = 1 - 2/ra; Dl = ra**2 - 2*ra + a**2
alpha = sp.sqrt(ra**2/(fr*Dl*(2/ra))*rd**2 + Dl/(fr**2*(2/ra))*pd**2); bb = -2*a/(ra*fr)
check('fixed curve, opposite senses: dt(+) - dt(-) = 2 b dphi exactly in a',
      sp.simplify((alpha + bb*pd) - (alpha.subs({rd: -rd, pd: -pd}) - bb*pd) - 2*bb*pd) == 0)
check('fixed curve, opposite senses: dtau(+) = dtau(-)', sp.simplify(fr*alpha - fr*alpha.subs({rd: -rd, pd: -pd})) == 0)


def rail_integrals(rm, q=2**(2/3)):
    r0 = q*rm; J = math.sqrt(Q(rm)); S = math.sqrt(r0 - rm)

    def rp(rv):                                  # dphi/dr, singular as 1/sqrt(r - rm)
        return math.sqrt(P(rv))/math.sqrt(Q(rv)**2/J**2 - Q(rv))

    def integ(fun):
        return 2*quad(lambda u: fun(rm + u*u)*rp(rm + u*u)*2*u, 0, S, epsabs=0, epsrel=1e-12, limit=400)[0]
    Phi = integ(lambda rv: 1.0)
    Ig = integ(lambda rv: 1/(rv*(1 - 2/rv)))
    Tg = integ(lambda rv: Q(rv)/J)
    chi = math.atan((1/math.sqrt(1 - 2/r0))/(r0*rp(r0)))
    return dict(Phi=Phi, Ig=Ig, Ic=Phi/(r0*(1 - 2/r0)), Trail=Tg, Tarc=math.sqrt(Q(r0))*Phi, chi=chi, r0=r0, J=J, rp=rp, S=S)


circ = []
for rm in (1e5, 1.5e5, 2e5):
    R = rail_integrals(rm)
    k2 = -4*(R['Ig'] - R['Ic'])*TG                      # coordinate-time difference 2 oint beta_K per unit a*
    T = (R['Trail'] + R['Tarc'])*TG
    circ.append(dict(r_min_M=rm, k2_ps_per_unit_astar=k2*1e12, signal_ps_astar_half=0.5*k2*1e12,
                     station_proper_ps_astar_half=0.5*k2*1e12*math.sqrt(1 - 2/R['r0']), duration_s=T,
                     r_min_million_km=rm*RG/1e9, station_r0_million_km=R['r0']*RG/1e9,
                     requirement_y=0.1*abs(0.5*k2)/T, angle_tangents_deg=math.degrees(R['chi']),
                     velocity_turn_deg=180 - math.degrees(R['chi'])))
check('circuit k2 at 1e5 M per unit a* and per (M/9.62 Msun) = -655.27 ps (scale invariance)', abs(circ[0]['k2_ps_per_unit_astar']*9.62/MASS + 655.2713) < 1e-3, round(circ[0]['k2_ps_per_unit_astar'], 3))
check('benchmark 1.5e5 M: signal about -742 ps (a* = 0.5), duration about 9.4 h', abs(circ[1]['signal_ps_astar_half'] + 742) < 3 and 9.2 < circ[1]['duration_s']/3600 < 9.6,
      (round(circ[1]['signal_ps_astar_half'], 3), round(circ[1]['duration_s']/3600, 4)))
check('benchmark 1.5e5 M: turning point 7.24 and station 11.5 million km from the black hole',
      abs(circ[1]['r_min_million_km'] - 7.24) < 0.02 and abs(circ[1]['station_r0_million_km'] - 11.50) < 0.03,
      (round(circ[1]['r_min_million_km'], 4), round(circ[1]['station_r0_million_km'], 4)))
check('velocity turn at the junctions = 120 degrees', all(abs(z['velocity_turn_deg'] - 120) < 1e-3 for z in circ), round(circ[1]['velocity_turn_deg'], 6))


# ---------------------------------------------------------------- 3. companion: shift of the station arrival difference
from scipy.integrate import cumulative_trapezoid


def rotation(al, be, ga):
    def Rz(t): return np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1]])
    def Rx(t): return np.array([[1, 0, 0], [0, math.cos(t), -math.sin(t)], [0, math.sin(t), math.cos(t)]])
    return Rz(al) @ Rx(be) @ Rz(ga)


def make_system(M, ms, a, e, P):
    GMb = M*GM_SUN; GMs = ms*GM_SUN; rg = GMb/C**2

    def kep(t):
        Mm = 2*math.pi*t/P
        u = brentq(lambda z: z - e*math.sin(z) - Mm, Mm - 2, Mm + 2)
        return a*np.array([math.cos(u) - e, math.sqrt(1 - e**2)*math.sin(u), 0.])

    def circuit(rminM, n=20001):
        """Nominal Newtonian circuit at zero energy: rail r = r_min cos^(-2/3)(3 phi/2), |phi| <= 40 deg, arc at r0."""
        rmin = rminM*rg; r0 = 2**(2/3)*rmin
        ph = np.linspace(-math.radians(40), math.radians(40), n)
        r = rmin*np.cos(1.5*ph)**(-2/3); drdph = rmin*np.cos(1.5*ph)**(-5/3)*np.sin(1.5*ph)
        v = np.sqrt(2*GMb/r); t_rail = cumulative_trapezoid(np.sqrt(r**2 + drdph**2)/v, ph, initial=0.)
        pa = np.linspace(math.radians(40), -math.radians(40), n)[1:]; va = math.sqrt(2*GMb/r0)
        t_arc = t_rail[-1] + r0*(math.radians(40) - pa)/va
        PH = np.concatenate([ph, pa]); R = np.concatenate([r, np.full(pa.size, r0)])
        T = np.concatenate([t_rail, t_arc]); V = np.concatenate([v, np.full(pa.size, va)])
        X = R[:, None]*(np.cos(PH)[:, None]*np.array([0, 0, 1.]) + np.sin(PH)[:, None]*np.array([1., 0, 0]))
        return T, X, V

    def phi_v(X, Rc):
        d = np.linalg.norm(Rc, axis=1); nn = Rc/d[:, None]
        return -(GMs/(2*d**3))*(3*np.einsum('ij,ij->i', X, nn)**2 - np.einsum('ij,ij->i', X, X))

    def comps(tstart, times, static=False, h=60.):
        c = np.array([kep(tstart + (0. if static else t)) for t in times])
        if static:
            return c, c, c
        return c, np.array([kep(tstart + t + h) for t in times]), np.array([kep(tstart + t - h) for t in times])

    def delay(T, X, V, Rm, cc, h=60.):
        c, cp, cq = cc
        Phi = phi_v(X, c @ Rm.T); dtPhi = (phi_v(X, cp @ Rm.T) - phi_v(X, cq @ Rm.T))/(2*h)
        vdv = Phi[0] - Phi + cumulative_trapezoid(dtPhi, T, initial=0.)
        return -np.trapezoid(vdv/V**2, T)

    def shifts(rminM, tstart, grid, n=20001, static=False):
        T, X, V = circuit(rminM, n); Tm = T[-1] - T[::-1]
        cP = comps(tstart, T, static); cM = comps(tstart, Tm, static)
        return np.array([delay(T, X, V, rotation(*o), cP) - delay(Tm, X[::-1], V[::-1], rotation(*o), cM) for o in grid])
    return dict(shifts=shifts, kep=kep, rg=rg, GMb=GMb, GMs=GMs, P=P)


grid = [(al, be, 0.37) for al in np.linspace(0, 2*math.pi, 8, endpoint=False) for be in np.linspace(0.1, math.pi - 0.1, 6)]
RMB = 1.5e5                                              # benchmark circuit
sysB = make_system(MASS, MSTAR, ABIN, EBIN, PBIN)
st = sysB['shifts'](RMB, PBIN/2, grid[:6], 8001, static=True)
check('companion held fixed: the two senses cancel exactly', np.max(np.abs(st)) < 1e-20, float(np.max(np.abs(st))))
c1 = sysB['shifts'](RMB, PBIN/2, grid[:4], 8001); c2 = sysB['shifts'](RMB, PBIN/2, grid[:4], 20001)
check('quadrature converged (8001 vs 20001 points)', np.max(np.abs(c1 - c2)/np.abs(c2)) < 1e-5, float(np.max(np.abs(c1 - c2)/np.abs(c2))))
S_B = abs(circ[1]['signal_ps_astar_half'])*1e-12
EPS_MASS = 0.05/MSTAR                                  # relative uncertainty of the published companion mass
companion = {}
for label, ts in (('apastron', PBIN/2), ('periastron', -0.2*DAY)):
    d = np.abs(sysB['shifts'](RMB, ts, grid))
    companion[label] = dict(lawB_max_ps=d.max()*1e12, lawB_median_ps=float(np.median(d))*1e12,
                            lawA_eps_for_10pc_worst=0.1*S_B/d.max(), lawA_eps_for_10pc_median=0.1*S_B/float(np.median(d)),
                            lawA_eps_mass=EPS_MASS, lawA_residual_mass_only_max_ps=EPS_MASS*d.max()*1e12, lawA_residual_mass_only_median_ps=EPS_MASS*float(np.median(d))*1e12)
check('apastron, law B: shift of the station difference below the signal in every sampled orientation', companion['apastron']['lawB_max_ps'] < S_B*1e12, round(companion['apastron']['lawB_max_ps'], 2))
check('apastron, law A with a mass-only model error eps = 0.05/0.76: residual below 3% of the signal', companion['apastron']['lawA_residual_mass_only_max_ps'] < 0.03*S_B*1e12, round(companion['apastron']['lawA_residual_mass_only_max_ps'], 3))
check('periastron, law B: shift larger than the signal (campaign near apastron)', companion['periastron']['lawB_max_ps'] > S_B*1e12, round(companion['periastron']['lawB_max_ps'], 1))


# independent check: direct nonlinear integration of the sliding motion (one orientation near periastron, large shift)
def nonlinear_shift(sysd, rminM, Rm, tstart):
    rmin = rminM*sysd['rg']; r0 = 2**(2/3)*rmin; P40 = math.radians(40); GMb = sysd['GMb']; GMs = sysd['GMs']
    e1 = np.array([0, 0, 1.]); e2 = np.array([1., 0, 0])

    def gradphi(x, Rc):
        d = np.linalg.norm(Rc); n = Rc/d
        return -(GMs/d**3)*(3*(x @ n)*n - x)

    def arrival(sense, tide):
        def X(r, p): return r*(math.cos(p)*e1 + math.sin(p)*e2)
        def dX(r, dr, p): return dr*(math.cos(p)*e1 + math.sin(p)*e2) + r*(-math.sin(p)*e1 + math.cos(p)*e2)
        def rail(p): return rmin*math.cos(1.5*p)**(-2/3), rmin*math.cos(1.5*p)**(-5/3)*math.sin(1.5*p)
        def f_rail(t, y, sg):
            p, w = y; r, dr = rail(p); sp_ = math.hypot(r, dr)
            dPhi = GMb*dr/r**2 + (gradphi(X(r, p), Rm @ sysd['kep'](tstart + t)) @ dX(r, dr, p) if tide else 0.)
            return [sg*w/sp_, -sg*dPhi/sp_]
        def f_arc(t, y, sg):
            p, w = y
            dPhi = (gradphi(X(r0, p), Rm @ sysd['kep'](tstart + t)) @ dX(r0, 0., p) if tide else 0.)
            return [sg*w/r0, -sg*dPhi/r0]
        o = dict(rtol=1e-13, atol=1e-16, method='DOP853'); w0 = math.sqrt(2*GMb/r0)
        up = lambda t, y, s: y[0] - P40; up.terminal = True
        dn = lambda t, y, s: y[0] + P40; dn.terminal = True
        if sense > 0:
            a1 = solve_ivp(f_rail, (0, 1e8), [-P40, w0], args=(1,), events=up, **o)
            a2 = solve_ivp(f_arc, (a1.t_events[0][0], 1e8), [P40, a1.y_events[0][0][1]], args=(-1,), events=dn, **o)
        else:
            a1 = solve_ivp(f_arc, (0, 1e8), [-P40, w0], args=(1,), events=up, **o)
            a2 = solve_ivp(f_rail, (a1.t_events[0][0], 1e8), [P40, a1.y_events[0][0][1]], args=(-1,), events=dn, **o)
        return a2.t_events[0][0]
    return (arrival(1, True) - arrival(-1, True)) - (arrival(1, False) - arrival(-1, False))


dpk = sysB['shifts'](RMB, -0.2*DAY, grid)
kmax = int(np.argmax(np.abs(dpk)))
nl = nonlinear_shift(sysB, RMB, rotation(*grid[kmax]), -0.2*DAY)
check('nonlinear sliding motion reproduces the first-order shift (largest periastron case, within 2%)', abs(nl/dpk[kmax] - 1) < 0.02,
      (round(nl*1e12, 1), round(dpk[kmax]*1e12, 1)))
# comparison at the same r_min/M: Gaia BH1 and BH2 (El-Badry et al. 2023)
comparison = {}
for name, p in (('Gaia BH1', (9.62, 0.93, 1.40*AU, 0.451, 185.59*DAY)), ('Gaia BH2', (8.94, 1.07, 4.96*AU, 0.5176, 1276.7*DAY))):
    sy = make_system(*p)
    d = np.abs(sy['shifts'](RMB, p[4]/2, grid))
    comparison[name] = dict(signal_ps=abs(circ[1]['k2_ps_per_unit_astar'])*0.5*p[0]/MASS, apastron_lawB_max_ps=d.max()*1e12,
                            apastron_lawB_median_ps=float(np.median(d))*1e12)
check('Gaia BH1 at the same r_min/M: companion shift exceeds the signal by more than 100', comparison['Gaia BH1']['apastron_lawB_median_ps'] > 100*comparison['Gaia BH1']['signal_ps'],
      round(comparison['Gaia BH1']['apastron_lawB_median_ps']/comparison['Gaia BH1']['signal_ps'], 1))

# ---------------------------------------------------------------- 4. thrust budget of the circuit (ideal photon rocket)
# Proper acceleration of the controlled worldline at a = 0 (the spin changes it at O(a)), Ehat = 1, from the exact
# Christoffel symbols of equatorial Schwarzschild in (t, r, phi); integrated budget I = int |b| dtau / c over the rail and
# the arc, plus the two junction turns made at maintained local speed and charge, each costing gamma v theta
# (theta = turn angle); the direct boost psi between incoming and outgoing 4-velocities
# (cosh psi = gamma^2 (1 - v^2 cos theta)), the minimum of int |b| dtau for a pure turn, is printed as a lower bound.
rr, jj = sp.symbols('rr jj', positive=True)
fS = 1 - 2/rr
gS = sp.diag(-fS, 1/fS, rr**2)
giS = gS.inv()
xs3 = [sp.Symbol('t_'), rr, sp.Symbol('phi_')]
Gam = [[[sp.simplify(sum(giS[i, l]*(sp.diff(gS[l, k], xs3[h]) + sp.diff(gS[l, h], xs3[k]) - sp.diff(gS[k, h], xs3[l]))/2
                          for l in range(3))) for h in range(3)] for k in range(3)] for i in range(3)]
Pr, Qr = rr**3/(2*(rr - 2)**2), rr**4/(2*(rr - 2))
dphidr = sp.sqrt(Pr)/sp.sqrt(Qr**2/jj**2 - Qr)               # outgoing branch, jj = sqrt(Q(r_min))
ut_ = 1/fS                                                    # Ehat = 1
ur2 = (fS*ut_**2 - 1)/(1/fS + rr**2*dphidr**2)                # (u^r)^2 from the mass shell
ur_ = sp.sqrt(ur2); up_ = dphidr*ur_
U3 = [ut_, ur_, up_]
acc_rail = [ur_*sp.diff(U3[i], rr) + sum(Gam[i][k][h]*U3[k]*U3[h] for k in range(3) for h in range(3)) for i in range(3)]
b2_rail = sp.lambdify((rr, jj), sum(gS[i, i]*acc_rail[i]**2 for i in range(3)), 'mpmath')
dtaudr = sp.lambdify((rr, jj), 1/ur_, 'mpmath')
# circle at r0 with Ehat = 1: u^r = 0, u^phi^2 = (1/f - 1)/r^2
up_c = sp.sqrt((1/fS - 1)/rr**2)
acc_c = [sum(Gam[i][k][h]*[ut_, 0, up_c][k]*[ut_, 0, up_c][h] for k in range(3) for h in range(3)) for i in range(3)]
b_circ = sp.lambdify(rr, sp.sqrt(sum(gS[i, i]*acc_c[i]**2 for i in range(3))), 'mpmath')
# check: on a free circular geodesic, u^phi^2 = 1/(r^3 - 3 r^2)/..., the same formula with the geodesic velocity gives b = 0
ut_g = 1/sp.sqrt(1 - 3/rr); up_g = ut_g/rr**sp.Rational(3, 2)
acc_g = [sum(Gam[i][k][h]*[ut_g, 0, up_g][k]*[ut_g, 0, up_g][h] for k in range(3) for h in range(3)) for i in range(3)]
check('Christoffel check: circular geodesic has zero acceleration', all(sp.simplify(z) == 0 for z in acc_g))
import mpmath as mpm
mpm.mp.dps = 30
G0 = 9.80665
rows4 = []
for rm in (1e5, RMB, 2e5):
    R = rail_integrals(rm); r0 = R['r0']
    jv = mpm.sqrt(mpm.mpf(rm)**4/(2*(mpm.mpf(rm) - 2)))         # sqrt(Q(r_min)) at 30 digits
    # rail: integrate |b| dtau = |b| (dtau/dr) dr on both branches, r = rm + u^2
    rmm = mpm.mpf(rm); S0 = mpm.sqrt(mpm.mpf(r0) - rmm)
    Irail_c = 2*mpm.quad(lambda u: mpm.sqrt(b2_rail(rmm + u*u, jv))*dtaudr(rmm + u*u, jv)*2*u, [mpm.mpf('1e-12')*S0, S0])
    assert abs(mpm.im(Irail_c)) < 1e-20*abs(mpm.re(Irail_c)), 'complex rail budget'
    Irail = mpm.re(Irail_c)
    rsamp = [mpm.mpf(rm)*(1 + mpm.mpf('1e-20'))] + [mpm.mpf(rm) + (mpm.mpf(k)/400)**2*(r0 - rm) for k in range(1, 401)]
    bpk = max(float(mpm.sqrt(b2_rail(x, jv))) for x in rsamp)
    tau_arc = (1 - 2/r0)*R['Tarc']                             # onboard: dtau = f dt/Ehat on the arc, Ehat = 1
    Iarc = float(b_circ(r0))*tau_arc
    v = math.sqrt(2/r0); gam = 1/math.sqrt(1 - v*v); th_turn = math.pi - R['chi']
    psi = math.acosh(gam**2*(1 - v*v*math.cos(th_turn)))       # unconstrained direct boost (lower bound)
    turn = gam*v*th_turn                                        # rotation at maintained local speed and charge
    I = float(Irail) + Iarc + 2*turn
    to_si = C*C/(rm*0 + 1)/RG                                 # acceleration in M units times c^2/r_g -> m s^-2
    rows4.append(dict(r_min_M=rm, peak_rail_g=bpk*C*C/RG/G0, arc_g=float(b_circ(r0))*C*C/RG/G0, I_rail=float(Irail),
                      I_arc=Iarc, junction_turn_maintained_charge=turn, junction_direct_boost_lower_bound=psi,
                      I_total=I, photon_mass_ratio=math.exp(I), rule_of_two_peak_g=2/(rm**2*math.sqrt(1 - 2/rm))*C*C/RG/G0))
check('rail peak = 2 g_s at r_min (relativistic rule of two)', all(abs(z['peak_rail_g']/z['rule_of_two_peak_g'] - 1) < 1e-9 for z in rows4), [round(z['peak_rail_g'], 3) for z in rows4])
check('arc proper time: dtau/dt = f0 at Ehat = 1', abs(tau_arc/R['Tarc'] - (1 - 2/r0)) < 1e-15)
check('benchmark 1.5e5 M: rail peak about 17 g, arc about 3.4 g', 16.5 < rows4[1]['peak_rail_g'] < 17.5 and 3.2 < rows4[1]['arc_g'] < 3.6,
      (round(rows4[1]['peak_rail_g'], 3), round(rows4[1]['arc_g'], 3)))
# transfer, for scale only: total mass ratio R, two symmetric rapidity changes, D = 590.6 pc, cruise only
PC = 3.0856775814913673e16; YR = 365.25*DAY
transfer = []
for Rm in (10., 20.):
    beta = (Rm - 1)/(Rm + 1); gam = (math.sqrt(Rm) + 1/math.sqrt(Rm))/2
    transfer.append(dict(mass_ratio=Rm, beta=beta, gamma=gam, years_inertial=DIST_PC*PC/(beta*C)/YR, years_onboard=DIST_PC*PC/(gam*beta*C)/YR))
check('transfer R = 20: about 2129 years (907 on board)', abs(transfer[1]['years_inertial'] - 2129.0) < 1 and abs(transfer[1]['years_onboard'] - 906.8) < 1,
      (round(transfer[1]['years_inertial'], 1), round(transfer[1]['years_onboard'], 1)))

out = dict(model='frozen Kerr, A = 1; Gaia BH3 M = 32.70 Msun, companion from the Gaia Collaboration (2024); a* = 0.5 illustrative',
           encounter=dict(formula='Delta tau = 6 pi a cos(i)/sqrt(r(r-3)) (first order in a, exact in r)',
                          cubic_to_linear_ratio_over_a2=cub_ratio, rows=enc),
           circuit=circ, companion_station_shift_benchmark_sampled=companion, companion_comparison_same_rmin=comparison, thrust_budget=rows4, transfer_for_scale=transfer, checks=checks)
(HERE/'probe_gaia_bh3.json').write_text(json.dumps(out, indent=2, default=str))
for z in enc:
    print('encounters r=%.0eM: jump %.4f ns (a*=0.5, cos i=1), fractional %.3e, period %.4f d, 200 encounters %.3f d (coordinate time)' %
          (z['r_M'], z['jump_ns'], z['fractional_rate'], z['period_days'], z['campaign_200_encounters_coordinate_days']))
for z in circ:
    print('circuit r_min=%.3gM: k2 %.4f ps, signal %.4f ps (a*=0.5; station proper time %.4f ps), duration %.2f s, |y| < %.3e, turn %.4f deg' %
          (z['r_min_M'], z['k2_ps_per_unit_astar'], z['signal_ps_astar_half'], z['station_proper_ps_astar_half'], z['duration_s'], z['requirement_y'], z['velocity_turn_deg']))
for z in rows4:
    print('thrust r_min=%.3g M: rail peak %.3f g (rule of two %.3f g), arc %.3f g, I_rail %.5f, I_arc %.5f, junction turn %.5f each (direct-boost lower bound %.5f), I %.5f, exp(I) %.5f' % (z['r_min_M'], z['peak_rail_g'], z['rule_of_two_peak_g'], z['arc_g'], z['I_rail'], z['I_arc'], z['junction_turn_maintained_charge'], z['junction_direct_boost_lower_bound'], z['I_total'], z['photon_mass_ratio']))
for z in transfer:
    print('transfer R=%g: beta %.5f, %.1f years, %.1f on board' % (z['mass_ratio'], z['beta'], z['years_inertial'], z['years_onboard']))
for k, z in companion.items():
    print('companion, benchmark 1.5e5 M, %s (48 sampled orientations): law B shift max %.4g ps, median %.4g ps; law A with mass-only eps = %.4f: max %.3g ps (%.2f%% of the signal), median %.3g ps; eps for 10%%: %.2g (worst), %.2g (median)' % (k, z['lawB_max_ps'], z['lawB_median_ps'], z['lawA_eps_mass'], z['lawA_residual_mass_only_max_ps'], 100*z['lawA_residual_mass_only_max_ps']/(S_B*1e12), z['lawA_residual_mass_only_median_ps'], z['lawA_eps_for_10pc_worst'], z['lawA_eps_for_10pc_median']))
for k, z in comparison.items():
    print('comparison %s at r_min = 1.5e5 M, apastron: signal %.0f ps, law B shift max %.4g ps, median %.4g ps' % (k, z['signal_ps'], z['apastron_lawB_max_ps'], z['apastron_lawB_median_ps']))
npass = sum(c['status'] == 'PASS' for c in checks)
print(f'{npass}/{len(checks)} checks passed')
sys.exit(0 if npass == len(checks) else 1)
