"""Proposed Solar-System two-probe control comparison (section `sec:observational-context`).

Mission concept, not an engineering certification. Units: AU, Julian year.
  * Earth -> Saturn free transfer, perihelion 1 AU, aphelion 20 AU (patched conics).
  * Saturn flyby rotating the relative velocity by -85 deg (pericentre from the deflection law),
    propagated from pericentre to 1.01 Hill radii under Sun + Jupiter, Saturn, Uranus, Neptune.
  * Planets on prescribed planar circular orbits, J2000 mean longitudes from JPL
    "Approximate Positions of the Planets", Table 1; heliocentric indirect term included;
    planetary potential set to zero at the Sun.
  * Rail: maintained Newtonian energy e = v^2/2 + Phi(x,t), arrival-time extremal
    xdot = V e_v, e_v' = -(g - (e_v.g) e_v)/V, thrust -2 g_perp - (d_t Phi/V) e_v.
  * Mother probe on the rail from A (Hill exit) to B (first crossing of r_B = 4 or 3 AU).
    Daughter on x_R(u) + d sin^2(pi u) n: same A, B, initial tangent and energy, no release kick;
    d = 0.01, 0.03, 0.05, 0.10 AU. The lag grows as d^2 because the rail makes the arrival time stationary.
  * Phase: earliest circular-model window in 2030-2040 with Jupiter > 2.5 AU from the rail and
    two-body ring-plane crossings beyond 3 Saturn radii. Dates are model phases, not launch dates.
Writes probe_saturn_comparison.json.
"""
from pathlib import Path
from datetime import datetime, timedelta
import json
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
AU = 149597870700.; YEAR = 365.25*86400.
# Saturn GM: Jacobson et al., AJ 132, 2520 (2006), Table, 37931207.7 km^3/s^2; equatorial radius
# 60268 km and pole: Archinal et al., Celest. Mech. Dyn. Astron. 130, 22 (2018).
GM = 1.32712440018e20; GMS = 3.79312077e16; RS = 60268000.
MU = GM*YEAR**2/AU**3; MUS = GMS*YEAR**2/AU**3
VS = AU/YEAR; AC = AU/YEAR**2; SATR = 9.58
SPEC = [('Jupiter', 5.2, 34.39644051, 1.26686534e17), ('Saturn', SATR, 49.95424423, GMS),
        ('Uranus', 19.2, 313.23810451, 5.793939e15), ('Neptune', 30.1, -55.12002969, 6.836529e15)]
PLANETS = [(nm, R, math.radians(L), math.sqrt(MU/R**3), gm*YEAR**2/AU**3) for nm, R, L, gm in SPEC]
HILL = SATR*(MUS/(3*MU))**(1/3)
EARTH_L0 = math.radians(100.46457166)
checks = []


def check(name, value, tolerance):
    checks.append(dict(name=name, residual=float(value), tolerance=tolerance,
                       status='PASS' if abs(value) <= tolerance else 'FAIL'))


def kepler_state(Q, r):
    """(v_r, v_t), time from perihelion and true anomaly at radius r; perihelion 1 AU, aphelion Q."""
    a = (1 + Q)/2; eps = -MU/(2*a); h = math.sqrt(MU*2*Q/(1 + Q))
    vt = h/r; vr = math.sqrt(max(0., 2*eps + 2*MU/r - vt*vt))
    e = (Q - 1)/(Q + 1); E = math.acos(np.clip((1 - r/a)/e, -1, 1))
    tof = (E - e*math.sin(E))*math.sqrt(a**3/MU)
    phi = 2*math.atan2(math.sqrt(1 + e)*math.sin(E/2), math.sqrt(1 - e)*math.cos(E/2))
    return np.array([vr, vt]), tof, phi


def state(spec, t):
    name, R, L, n, gm = spec; phi = L + n*t
    x = R*np.array([math.cos(phi), math.sin(phi), 0.])
    return x, n*np.array([-x[1], x[0], 0.])


def potential(x, t, planets=True):
    R = np.linalg.norm(x); phi = -MU/R; g = -MU*x/R**3; pt = 0.
    if planets:
        for spec in PLANETS:
            name, r, L, n, gm = spec; xp, vp = state(spec, t); d = x - xp; dd = np.linalg.norm(d)
            phi += -gm/dd + gm/r + gm*np.dot(x, xp)/r**3
            g += -gm*d/dd**3 - gm*xp/r**3
            pt += -gm*np.dot(d, vp)/dd**3 + gm*np.dot(x, vp)/r**3
    return phi, g, pt


def windows():
    _, T, angle = kepler_state(20., SATR)
    ns = PLANETS[1][3]; ls = PLANETS[1][2]; ne = math.sqrt(MU)
    out = []
    for turn in range(-60, 61):
        launch = (angle + EARTH_L0 - ls - ns*T + 2*math.pi*turn)/(ns - ne)
        if 30 <= launch < 40:
            out.append((launch, launch + T))
    return sorted(out)


def pole(t):
    """Saturn pole (RA, Dec) of the IAU WGCCRE 2015 report (= JPL PCK00011), equatorial J2000 -> ecliptic J2000."""
    ra = math.radians(40.589 - .036*t/100); de = math.radians(83.537 - .004*t/100)
    eq = np.array([math.cos(de)*math.cos(ra), math.cos(de)*math.sin(ra), math.sin(de)])
    ob = math.radians(23.43928)
    return np.array([eq[0], math.cos(ob)*eq[1] + math.sin(ob)*eq[2], -math.sin(ob)*eq[1] + math.cos(ob)*eq[2]])


def flyby(epoch, theta=-85., tol=2e-11):
    vin, _, _ = kepler_state(20., SATR); vs = math.sqrt(MU/SATR)
    incoming = np.r_[vin - [0, vs], 0.]; th = math.radians(theta)
    rot = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1.]])
    outgoing = rot@incoming; vi = np.linalg.norm(incoming)
    rp = MUS/vi**2*(1/math.sin(abs(th)/2) - 1)   # turn angle of a hyperbolic flyby (Battin 1999, eqs. 9.1, 9.10)
    rpdir = (incoming - outgoing)/np.linalg.norm(incoming - outgoing)
    vpdir = (incoming + outgoing)/np.linalg.norm(incoming + outgoing)
    lam = PLANETS[1][2] + PLANETS[1][3]*epoch
    rotS = np.array([[math.cos(lam), -math.sin(lam), 0], [math.sin(lam), math.cos(lam), 0], [0, 0, 1.]])
    er = rotS@rpdir; et = rotS@vpdir
    xs, vsat = state(PLANETS[1], epoch)
    x = xs + rp*er; v = vsat + math.sqrt(vi*vi + 2*MUS/rp)*et

    def rhs(t, y):
        return np.r_[y[3:], potential(y[:3], epoch + t, True)[1]]

    def leave(t, y):
        return np.linalg.norm(y[:3] - state(PLANETS[1], epoch + t)[0]) - 1.01*HILL
    leave.direction = 1; leave.terminal = True
    sol = solve_ivp(rhs, [0, 1.], np.r_[x, v], events=leave, method='DOP853',
                    rtol=tol, atol=tol/30, max_step=.003)
    assert sol.success and len(sol.t_events[0]) == 1
    duration = float(sol.t_events[0][0]); y = sol.y_events[0][0]
    ecc = 1 + rp*vi*vi/MUS; p = rp*(1 + ecc); normal = pole(epoch)
    node0 = math.atan2(-np.dot(normal, er), np.dot(normal, et)); fmax = math.acos(-1/ecc)
    nodes = [p/(1 + ecc*math.cos(node0 + i*math.pi))*AU/RS for i in range(-2, 3)
             if -fmax < node0 + i*math.pi < fmax]
    v_inf_kms = vi*VS/1000
    return dict(theta_deg=theta, v_infinity_Saturn_km_s=v_inf_kms, rp_over_RS=rp*AU/RS,
                max_deflection_at_2p3RS_deg=math.degrees(2*math.asin(1/(1 + 2.3*RS/AU*vi*vi/MUS))),
                epoch_peri=epoch, exit_after_days=duration*365.25, exit_epoch=epoch + duration,
                exit_position_AU=y[:3].tolist(), exit_velocity_km_s=(y[3:]*VS/1000).tolist(),
                ring_plane_crossings_RS=nodes), y


def ray(A, vA, epoch, rB, planets=True, tol=3e-11, initial_direction=None, endtime=None):
    direction = vA/np.linalg.norm(vA) if initial_direction is None else initial_direction
    energy = np.dot(vA, vA)/2 + potential(A, epoch, planets)[0]

    def quantities(t, y):
        x = y[:3]; d = y[3:6]/np.linalg.norm(y[3:6]); phi, g, pt = potential(x, epoch + t, planets)
        vv = math.sqrt(2*(energy - phi)); gp = g - d*np.dot(d, g)
        return vv, d, gp, -2*gp - pt/vv*d, g

    def rhs(t, y):
        vv, d, gp, thr, g = quantities(t, y)
        return np.r_[vv*d, -gp/vv, np.linalg.norm(thr)]

    def stop(t, y):
        return np.linalg.norm(y[:3]) - rB
    stop.direction = -1; stop.terminal = True
    sol = solve_ivp(rhs, [0, 6. if endtime is None else endtime], np.r_[A, direction, 0.],
                    events=stop if endtime is None else None, dense_output=True,
                    rtol=tol, atol=tol/30, method='DOP853', max_step=.015)
    assert sol.success
    if endtime is None:
        T = float(sol.t_events[0][0]); yf = sol.y_events[0][0]
    else:
        T = endtime; yf = sol.y[:, -1]
    return sol, quantities, T, yf


def path_cost(position, derivative, second, A, vA, epoch, tol=1e-11):
    """Time, thrust budget and peak thrust along a prescribed path at the maintained energy."""
    energy = np.dot(vA, vA)/2 + potential(A, epoch)[0]

    def at(u, t):
        x = position(u); xp = derivative(u); xpp = second(u); ell = np.linalg.norm(xp)
        phi, g, pt = potential(x, epoch + t); V = math.sqrt(2*(energy - phi))
        tp = ell/V; dV = (-pt*tp + np.dot(g, xp))/V
        tpp = np.dot(xp, xpp)/(ell*V) - ell*dV/V**2
        return tp, (xpp - xp*tpp/tp)/tp**2 - g

    def rhs(u, y):
        tp, thrust = at(u, y[0]); return [tp, np.linalg.norm(thrust)*tp]
    sol = solve_ivp(rhs, [0, 1], [0, 0], dense_output=True, method='DOP853',
                    rtol=tol, atol=tol/30, max_step=.015)
    assert sol.success
    peak = max(np.linalg.norm(at(u, sol.sol(u)[0])[1]) for u in np.linspace(0, 1, 501))*AC
    jup = min(np.linalg.norm(position(u) - state(PLANETS[0], epoch + sol.sol(u)[0])[0])
              for u in np.linspace(0, 1, 501))
    return dict(T_years=float(sol.y[0, -1]), delta_v_km_s=float(sol.y[1, -1])*VS/1000,
                peak_thrust_m_s2=peak, min_Jupiter_distance_AU=jup)


def rail_summary(sol, quant, T, epoch):
    times = np.linspace(0, T, 501)
    peak = max(np.linalg.norm(quant(t, sol.sol(t))[3]) for t in times)*AC
    jup = [np.linalg.norm(sol.sol(t)[:3] - state(PLANETS[0], epoch + t)[0]) for t in times]
    j = int(np.argmin(jup))
    opt = minimize_scalar(lambda t: np.linalg.norm(sol.sol(t)[:3] - state(PLANETS[0], epoch + t)[0]),
                          bounds=(times[max(0, j - 1)], times[min(500, j + 1)]), method='bounded')
    return dict(T_years=T, delta_v_km_s=float(sol.sol(T)[6])*VS/1000, peak_thrust_m_s2=peak,
                min_Jupiter_distance_AU=min(min(jup), float(opt.fun)))


def experiment(fb, y, rB):
    A = y[:3]; vA = y[3:]; epoch = fb['exit_epoch']
    sol, quant, T, yf = ray(A, vA, epoch, rB)
    B = yf[:3]; disp = B - A
    out = rail_summary(sol, quant, T, epoch); out['B_AU'] = B.tolist()
    n = np.cross([0., 0., 1.], disp); n = n/np.linalg.norm(n)

    def base(u):
        t = T*u; yy = sol.sol(t); V, di, gp, thr, g = quant(t, yy)
        return yy[:3], T*V*di, T*T*(g + thr)
    daughters = []
    for d in [.01, .03, .05, .10]:
        row = path_cost(lambda u: base(u)[0] + d*math.sin(math.pi*u)**2*n,
                        lambda u: base(u)[1] + d*math.pi*math.sin(2*math.pi*u)*n,
                        lambda u: base(u)[2] + d*2*math.pi**2*math.cos(2*math.pi*u)*n, A, vA, epoch)
        row.update(amplitude_AU=d, lag_hours=(row['T_years'] - T)*YEAR/3600)
        row['lag_over_amplitude_squared_h_per_AU2'] = row['lag_hours']/d**2
        check(f'positive lag d={d} rB={rB}', min(0., row['lag_hours']), 0.)
        daughters.append(row)
    q = [z['lag_over_amplitude_squared_h_per_AU2'] for z in daughters]
    spread = (max(q) - min(q))/min(q)
    check(f'quadratic lag: spread of lag/d^2 over d = 0.01-0.10 AU below 0.2% (rB={rB})', spread, .002)
    out['lag_over_d2_relative_spread'] = spread
    out['daughters'] = daughters
    # Independent fixed-endpoint shooting: the rail from A to B is recovered.
    ini = vA/np.linalg.norm(vA); e1 = np.cross([0., 0., 1.], ini); e1 /= np.linalg.norm(e1); e2 = np.cross(ini, e1)

    def residual(z):
        di = ini + z[0]*e1 + z[1]*e2; di /= np.linalg.norm(di)
        return ray(A, vA, epoch, rB, initial_direction=di, endtime=T*math.exp(z[2]))[3][:3] - B
    rr = least_squares(residual, [.001, .0003, .0005], xtol=1e-12, ftol=1e-12, gtol=1e-12)
    check(f'fixed-endpoint shooting residual AU rB={rB}', np.linalg.norm(residual(rr.x)), 1e-8)
    return out


def main():
    options = []
    for launch, peri in windows():
        fb, y = flyby(peri)
        sol, quant, T, yf = ray(y[:3], y[3:], fb['exit_epoch'], 4.)
        jup = rail_summary(sol, quant, T, fb['exit_epoch'])['min_Jupiter_distance_AU']
        options.append((jup, launch, fb, y))
    eligible = [z for z in options if z[0] > 2.5 and min(z[2]['ring_plane_crossings_RS']) > 3.]
    jup, launch, fb, y = min(eligible, key=lambda z: z[1])
    check('rail starts outside Saturn Hill sphere', min(0., np.linalg.norm(
        np.array(fb['exit_position_AU']) - state(PLANETS[1], fb['exit_epoch'])[0]) - HILL), 0.)
    vec, tof, _ = kepler_state(20., SATR)
    launch_vinf = math.sqrt(MU*(2 - 1/10.5)) - math.sqrt(MU)
    results = {str(R): experiment(fb, y, R) for R in [4., 3.]}
    date = lambda t: (datetime(2000, 1, 1, 12) + timedelta(seconds=t*YEAR)).date().isoformat()
    out = dict(caveat='mission concept in a prescribed circular-planet model; not launch targeting, '
                      'ring/dust safety, power-system or navigation certification',
               transfer=dict(perihelion_AU=1, aphelion_AU=20, Earth_v_infinity_km_s=launch_vinf*VS/1000,
                             C3_km2_s2=(launch_vinf*VS/1000)**2, Earth_to_Saturn_years=tof),
               model_phase=dict(launch=date(launch), flyby=date(fb['epoch_peri']), rail_start=date(fb['exit_epoch'])),
               flyby=fb, results=results, checks=checks)
    (HERE/'probe_saturn_comparison.json').write_text(json.dumps(out, indent=2))
    print('checks', len(checks), 'PASS', sum(c['status'] == 'PASS' for c in checks))
    print('flyby rp/RS %.4f, ring crossings %s RS' % (fb['rp_over_RS'], [round(z, 3) for z in fb['ring_plane_crossings_RS']]))
    for key, row in results.items():
        print('B at %s AU: rail %.4f yr, dv %.4f km/s, peak %.3e m/s^2, Jupiter >= %.3f AU' %
              (key, row['T_years'], row['delta_v_km_s'], row['peak_thrust_m_s2'], row['min_Jupiter_distance_AU']))
        for z in row['daughters']:
            print('   d=%.2f AU: lag %.4f h, lag/d^2 %.2f h/AU^2, dv %.4f km/s, peak %.3e m/s^2' %
                  (z['amplitude_AU'], z['lag_hours'], z['lag_over_amplitude_squared_h_per_AU2'],
                   z['delta_v_km_s'], z['peak_thrust_m_s2']))
    assert all(c['status'] == 'PASS' for c in checks)


if __name__ == '__main__':
    main()
