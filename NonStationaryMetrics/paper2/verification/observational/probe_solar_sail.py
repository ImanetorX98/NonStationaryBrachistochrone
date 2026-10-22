"""Solar-sail clock comparison (section `sec:observational-context`).

Static Schwarzschild field of the Sun, G = c = M_sun = 1. A rail at fixed Ehat = sqrt(f) gamma has two
time-optimal problems with the same endpoints: the arrival (Killing-time) cost, optical metric
Lt^2 (dr^2/f + r^2 dphi^2) with Lt^2 = Ehat^2/(f(Ehat^2 - f)), and the onboard cost, Ltau^2 = f/(Ehat^2 - f).
Their ratio is f^2/Ehat^2, so the two rails coincide at Newtonian order and separate at order U.

At Ehat = 1 the rails are arcs symmetric about the pericentre r_min with endpoints at r0 = 2^(2/3) r_min,
where sin(alpha) = 1/2 at Newtonian order. Reports, for r_min = 0.1 and 0.05 au: the angle, the duration,
the separation of the two pericentres, the two inversions (each rail wins on its own clock), the thrust
2 g sin(alpha) and its Delta v, and the sail lightness beta = 2/sin(alpha) needed by an ideal sail whose
normal is the outward normal to the velocity (force beta g sin^2(alpha)); the critical areal density
L_sun/(2 pi c GM_sun) and the sail equilibrium temperature for absorptance 0.1 and total emissivity 1.
Daughters follow r(phi) = r_t(phi) + d sin^2(pi u), u = (phi + Phi/2)/Phi, with the same endpoints:
t(d) - t(0) is quadratic (the rail is stationary), tau(d) has its vertex at d* != 0, and X = t - tau is
linear and antisymmetric in d. The local energy sensitivity of X is -(Ehat - 1) times that of t
(weak-field series, symbolic). Writes probe_solar_sail.json.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
mp.mp.dps = 45
C = mp.mpf(299792458); GM = mp.mpf('1.32712440018e20'); AU = mp.mpf(149597870700)
L_SUN = mp.mpf('3.828e26'); SIGMA_SB = mp.mpf('5.670374419e-8'); R_SUN = mp.mpf(695700000)
LM = GM/C**2; TM = LM/C                                  # length and time units
checks = []


def check(name, value, atol):
    checks.append(dict(name=name, value=float(value), atol=atol, status='PASS' if abs(value) <= atol else 'FAIL'))


def symbolic():
    r, E, l, U, eps = sp.symbols('r E l U epsilon', positive=True)
    f = 1 - 2/r; gam = E/sp.sqrt(f); v2 = 1 - 1/gam**2
    lt2, ltau2 = 1/(v2*f), 1/(gam**2*v2)
    check('Lt^2 = E^2/(f(E^2-f))', sp.simplify(lt2 - E**2/(f*(E**2 - f))) != 0, 0)
    check('Ltau^2 = f/(E^2-f)', sp.simplify(ltau2 - f/(E**2 - f)) != 0, 0)
    check('Ltau^2/Lt^2 = f^2/E^2', sp.simplify(ltau2/lt2 - f**2/E**2) != 0, 0)
    x = (E - f)/(sp.sqrt(f)*sp.sqrt(E**2 - f))            # Lt - Ltau
    check('Lt - Ltau = (E-f)/(sqrt f sqrt(E^2-f))',
          abs(complex((sp.sqrt(lt2) - sp.sqrt(ltau2) - x).subs({r: sp.Rational(37, 3), E: sp.Rational(21, 20)}).evalf(40))), 1e-35)
    w = {r: 1/(l*U), E: 1 + l*eps}
    ratio = sp.series((sp.diff(x, E)/sp.diff(sp.sqrt(lt2), E)).subs(w), l, 0, 2).removeO()
    check('energy sensitivity of X / of t = -(E-1) + O(l^2)', sp.simplify(ratio + l*eps) != 0, 0)


class Rails:
    def __init__(self, rmin_au):
        self.rmin = rmin_au*AU/LM; self.r0 = mp.mpf(2)**(mp.mpf(2)/3)*self.rmin; self.e2 = mp.mpf(1)

    @staticmethod
    def f(r):
        return 1 - 2/r

    def l2(self, k, r):
        f = self.f(r)
        return self.e2/(f*(self.e2 - f)) if k == 't' else f/(self.e2 - f)

    def dphi(self, k, rm, r):
        j = mp.sqrt(self.l2(k, rm))*rm
        return j/(r*mp.sqrt(self.f(r))*mp.sqrt(self.l2(k, r)*r*r - j*j))

    def angle(self, k, rm):
        return 2*mp.quad(lambda s: self.dphi(k, rm, rm + s*s)*2*s, [0, mp.sqrt(self.r0 - rm)], method='gauss-legendre')

    def cost(self, a, b, rm):
        """cost in metric b along the rail of metric a"""
        def g(s):
            r = rm + s*s
            return mp.sqrt(self.l2(b, r))*mp.sqrt(1/self.f(r) + r*r*self.dphi(a, rm, r)**2)*2*s
        return 2*mp.quad(g, [0, mp.sqrt(self.r0 - rm)], method='gauss-legendre')

    def sin_alpha(self, r):
        p = r*self.dphi('t', self.rmin, r)
        return p/mp.sqrt(1/self.f(r) + p*p)

    def delta_v(self):
        def g(s):
            r = self.rmin + s*s
            return 2/r**2*self.sin_alpha(r)*mp.sqrt(self.l2('t', r))*mp.sqrt(1/self.f(r) + (r*self.dphi('t', self.rmin, r))**2)*2*s
        return 2*mp.quad(g, [0, mp.sqrt(self.r0 - self.rmin)], method='gauss-legendre')

    def energy_sensitivity(self, rel=mp.mpf(10)**-20):
        """k = (dX/dEhat^2)/(dt/dEhat^2) on the arrival rail with the endpoints fixed (symmetric differences)"""
        phi = self.angle('t', self.rmin); vals = []
        for e in (1 + rel, 1 - rel):
            self.e2 = e
            rm = mp.findroot(lambda x: self.angle('t', x) - phi, self.rmin, tol=mp.mpf(10)**-38)
            t = self.cost('t', 't', rm); vals.append((t, t - self.cost('t', 'tau', rm)))
        self.e2 = mp.mpf(1)
        return (vals[0][1] - vals[1][1])/(vals[0][0] - vals[1][0])

    def summary(self):
        phi = self.angle('t', self.rmin)
        rm_tau = mp.findroot(lambda x: self.angle('tau', x) - phi, self.rmin, tol=mp.mpf(10)**-38)
        tt, ttau = self.cost('t', 't', self.rmin), self.cost('tau', 't', rm_tau)
        taut, tautau = self.cost('t', 'tau', self.rmin), self.cost('tau', 'tau', rm_tau)
        self.phi = phi
        return dict(rmin_au=float(self.rmin*LM/AU), rmin_Rsun=float(self.rmin*LM/R_SUN), r0_au=float(self.r0*LM/AU),
                    phi_deg=float(phi*180/mp.pi), duration_h=float(tt*TM/3600), onboard_h=float(taut*TM/3600),
                    v_peri_kms=float(mp.sqrt(1 - self.f(self.rmin))*C/1000),
                    separation_km=float((rm_tau - self.rmin)*LM/1000),
                    inversion_t_ns=float((ttau - tt)*TM*1e9), inversion_tau_ns=float((taut - tautau)*TM*1e9),
                    g_peri=float(C**2/LM/self.rmin**2), delta_v_kms=float(self.delta_v()*C/1000),
                    sin_alpha_end=float(self.sin_alpha(self.r0)), beta_required=float(2/self.sin_alpha(self.r0)))

    # daughters: radial displacement d sin^2(pi u) at fixed polar angle, u = (phi + Phi/2)/Phi
    def daughters(self, d_km, degree=7):
        rm, half = self.rmin, self.phi/2
        S = mp.sqrt(self.r0 - rm)
        xs, ws = zip(*[(S*(x + 1)/2, S*w/2) for x, w in gauss_legendre(degree)])
        phis = [mp.quad(lambda s: self.dphi('t', rm, rm + s*s)*2*s, [0, x], method='gauss-legendre') for x in xs]
        dphids = [self.dphi('t', rm, rm + x*x)*2*x for x in xs]

        def costs(d):
            out = {'t': mp.mpf(0), 'tau': mp.mpf(0)}
            for x, w, ph, dp in zip(xs, ws, phis, dphids):
                for sign in (1, -1):                       # the two symmetric halves
                    u = (sign*ph + half)/(2*half)
                    rr = rm + x*x + d*mp.sin(mp.pi*u)**2
                    drdphi = sign*2*x/dp + d*mp.pi/(2*half)*mp.sin(2*mp.pi*u)
                    for k in out:
                        out[k] += w*dp*mp.sqrt(self.l2(k, rr))*mp.sqrt(drdphi**2/self.f(rr) + rr*rr)
            return out
        base = costs(0)
        rows = []
        for dk in d_km:
            c = costs(mp.mpf(dk)*1000/LM)
            dt, dtau = (c['t'] - base['t'])*TM, (c['tau'] - base['tau'])*TM
            rows.append(dict(d_km=dk, t_lag_s=float(dt), tau_lag_s=float(dtau), X_ns=float((dt - dtau)*1e9),
                             t_lag_over_d2=float(dt/dk**2)))
        return base, rows, costs


def gauss_legendre(degree):
    """3*2^(degree-1) Gauss-Legendre nodes and weights on [-1, 1] from mpmath"""
    return mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)


def main():
    symbolic()
    sigma_star = L_SUN/(2*mp.pi*C*GM)
    cases = []
    for q in ['0.1', '0.05']:
        R = Rails(mp.mpf(q)); row = R.summary()
        check(f'rmin={q} au: the two inversions agree (second order)', (row['inversion_t_ns'] - row['inversion_tau_ns'])/row['inversion_t_ns'], 1e-5)
        check(f'rmin={q} au: both inversions positive', 0 if row['inversion_t_ns'] > 0 and row['inversion_tau_ns'] > 0 else 1, 0)
        check(f'rmin={q} au: beta = 4 at the endpoints', row['beta_required'] - 4, 1e-6)
        S = mp.mpf(q)*AU
        row['sail_T_eq_K'] = float((mp.mpf('0.1')*L_SUN/(4*mp.pi*S*S)/SIGMA_SB)**mp.mpf(0.25))
        base, drows, costs = R.daughters([-1000, -100, -10, 10, 100, 1000])
        lag = [z['t_lag_over_d2'] for z in drows]
        check(f'rmin={q} au: t(d) quadratic (stationary rail)', (max(lag) - min(lag))/max(lag), 1e-3)
        D = {z['d_km']: z for z in drows}
        k = R.energy_sensitivity(); row['energy_k'] = float(k)
        check(f'rmin={q} au: X suppresses energy errors (|k| < 1e-7)', 0 if abs(k) < 1e-7 else 1, 0)
        # arrival shift allowed by a 10% measurement of X at d = 1000 km, and the clock stability it needs
        row['allowed_arrival_shift_s'] = float(abs(D[1000]['X_ns'])*1e-10/abs(k))
        row['clock_stability_10pct'] = float(abs(D[1000]['X_ns'])*1e-10/(row['duration_h']*3600))
        # along-track distance at the end between mother and the d = 1000 km daughter, ideal control
        row['end_separation_m'] = float(D[1000]['t_lag_s']*mp.sqrt(1 - R.f(R.r0))*C)
        check(f'rmin={q} au: X antisymmetric at 100 km', (D[-100]['X_ns'] + D[100]['X_ns'])/D[100]['X_ns'], 1e-3)
        slope = (D[-100]['X_ns'] - D[100]['X_ns'])/200
        # vertex of tau(d) from the quadratic through d = -10, 0, 10 km
        a2 = (D[-10]['tau_lag_s'] + D[10]['tau_lag_s'])/(2*10**2); a1 = (D[10]['tau_lag_s'] - D[-10]['tau_lag_s'])/20
        row.update(daughters=drows, X_slope_ns_per_km=slope, tau_vertex_km=-a1/(2*a2))
        cases.append(row)
    check('self-similarity at Ehat = 1: same separation', (cases[0]['separation_km'] - cases[1]['separation_km'])/cases[1]['separation_km'], 1e-6)
    check('self-similarity: same tau vertex', (cases[0]['tau_vertex_km'] - cases[1]['tau_vertex_km'])/cases[1]['tau_vertex_km'], 1e-4)
    check('tau vertex agrees with an independent minimisation (-0.6814747 km)', cases[1]['tau_vertex_km'] + 0.6814747, 1e-5)
    check('X slope scales as rmin^(-1/2)', cases[1]['X_slope_ns_per_km']/cases[0]['X_slope_ns_per_km'] - mp.sqrt(2), 1e-3)
    out = dict(model='static Schwarzschild field of the Sun, Ehat = 1, symmetric arcs r0 = 2^(2/3) rmin, ideal sail; '
                     'daughters r = r_t(phi) + d sin^2(pi u); a concept, not a mission design',
               sigma_critical_g_m2=float(sigma_star*1000), sigma_beta4_g_m2=float(sigma_star*250), cases=cases, checks=checks)
    (HERE/'probe_solar_sail.json').write_text(json.dumps(out, indent=2))
    print('checks', len(checks), 'PASS', sum(c['status'] == 'PASS' for c in checks))
    print('critical areal density %.4f g/m^2; beta = 4 -> %.4f g/m^2' % (sigma_star*1000, sigma_star*250))
    for z in cases:
        print('rmin %.3f au (%.2f Rsun): Phi %.5f deg, %.3f h, v_p %.2f km/s, separation %.4f km, inversions %.5f / %.5f ns, '
              'g_p %.3f m/s^2, dv %.1f km/s, beta %.4f, T_eq %.0f K' %
              (z['rmin_au'], z['rmin_Rsun'], z['phi_deg'], z['duration_h'], z['v_peri_kms'], z['separation_km'],
               z['inversion_t_ns'], z['inversion_tau_ns'], z['g_peri'], z['delta_v_kms'], z['beta_required'], z['sail_T_eq_K']))
        for y in z['daughters']:
            print('    d=%+5d km: t lag %.6e s, tau lag %.6e s, X %+.4f ns' % (y['d_km'], y['t_lag_s'], y['tau_lag_s'], y['X_ns']))
        print('    X slope %.5f ns/km, tau vertex %.5f km, energy k %.3e, allowed arrival shift %.2f s, clock %.2e, end separation %.0f m' %
              (z['X_slope_ns_per_km'], z['tau_vertex_km'], z['energy_k'], z['allowed_arrival_shift_s'], z['clock_stability_10pct'], z['end_separation_m']))
    for c in checks:
        if c["status"] != "PASS":
            print("FAIL", c)
    assert all(c["status"] == "PASS" for c in checks)


if __name__ == '__main__':
    main()
